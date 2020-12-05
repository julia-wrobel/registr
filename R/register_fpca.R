#' Register curves using constrained optimization and GFPCA
#'
#' Function combines constrained optimization and GFPCA to estimate warping functions for 
#' exponential family curves. See argument \code{family} for which families are
#' supported. Warping functions are calculated by the function \code{\link{registr}}.
#' The GFPCA step can be performed either using the variational EM-based GFPCA
#' approaches of Wrobel et al. (2019) (\code{fpca_type = "variationalEM"}, default)
#' or the mixed model-based two-step approach of Gertheiss et al. (2017)
#' (\code{fpca_type = "two-step"}). \cr \cr
#' Warping functions by default are forced to start and end on the diagonal to be
#' domain-preserving. This behavior can be changed by setting
#' \code{incompleteness} to some other value than NULL and a reasonable \code{lambda_inc} value.
#' For further details see the accompanying vignette. \cr \cr
#' By specifying \code{cores > 1} the registration call can be parallelized.
#' 
#' Requires input data \code{Y} to be a dataframe in long format with variables 
#' \code{id}, \code{index}, and \code{value} to indicate subject IDs, times, and observations, 
#' respectively.
#' 
#' One joint iteration consists of a GFPCA step and a registration step.
#' As preprocessing, one initial registration step is performed.
#' The template function for this registration step is defined by argument
#' \code{Y_template}.
#' After convergence or \code{max_iterations} is reached, one final GFPCA step
#' is performed.
#'
#' @param family One of \code{c("gaussian","binomial","gamma")}.
#' For \code{"gamma"}, the \code{fpca_type} is fixed to \code{"two-step"}.
#' Defaults to \code{"gaussian"}.
#' @param Y_template Optional dataframe with the same structure as \code{Y}.
#' Only used for the initial registration step. If NULL,
#' curves are registered to the overall mean of all curves in \code{Y} as template function.
#' If specified, the template function is taken as the mean
#' of all curves in \code{Y_template}. Defaults to NULL.
#' @param max_iterations Number of iterations for overall algorithm. Defaults to 10.
#' @param npc Number of principal components to calculate. Defaults to 1. 
#' @param fpca_type One of \code{c("variationalEM","two-step")}.
#' Defaults to \code{"variationalEM"}.
#' @param fpca_maxiter Only used if \code{fpca_type = "variationalEM"}. Number
#' to pass to the \code{maxiter} argument of `bfpca()` or `fpca_gauss()`. 
#' Defaults to 50.
#' @param fpca_seed Only used if \code{fpca_type = "variationalEM"}. Number to
#' pass to the \code{seed} argument of `bfpca()` or `fpca_gauss()`. Defaults to
#' 1988.
#' @param fpca_error_thresh Only used if \code{fpca_type = "variationalEM"}.
#' Number to pass to the \code{error_thresh} argument of `bfpca()` or
#' `fpca_gauss()`. Defaults to 0.0001.
#' @param fpca_index_relevantDigits Only used if \code{fpca_type = "two-step"}.
#' Positive integer \code{>= 2}, stating the number of relevant digits to which
#' the index grid should be rounded in the GFPCA step. Coarsening the index grid
#' is necessary since otherwise the covariance surface matrix explodes in size
#' in the presence of too many unique index values (which is the case after some
#' registration step). Defaults to 4. Set to \code{NULL} to prevent rounding.
#' @param ... Additional arguments passed to registr and to the gfpca functions
#' (if \code{fpca_type = "variationalEM"}).
#' @inheritParams registr
#'
#' @author Julia Wrobel \email{julia.wrobel@@cuanschutz.edu}
#' Jeff Goldsmith \email{ajg2202@@cumc.columbia.edu},
#' Alexander Bauer \email{alexander.bauer@@stat.uni-muenchen.de}
#' @export
#' 
#' @return An object of class \code{registration} containing:
#' \item{fpca_obj}{List of items from FPCA step.}
#' \item{Y}{The observed data plus variables \code{t_star} and \code{t_hat} which are the
#' unregistered grid and registered grid, respectively.}
#' \item{time_warps}{List of time values for each iteration of the algorithm. 
#' \code{time_warps[1]} is the original (observed) time
#' and \code{time_warps[n]} provides time values for the nth iteration.}
#' \item{loss}{Loss for each iteration of the algorithm, calculated in the registration step using an 
#' exponential family likelihood with natural parameter from the FPCA step.}
#' @return family \code{gaussian} or \code{binomial}.
#'  
#' @examples
#' library(ggplot2)
#' 
#' ### complete binomial curves
#' Y = simulate_unregistered_curves(I = 20, D = 200)
#' 
#' # estimation based on Wrobel et al. (2019)
#' reg = register_fpca(Y, npc = 2, family = "binomial",
#'                     fpca_type = "variationalEM", max_iterations = 5)
#' # estimation based on Gertheiss et al. (2017)
#' reg2 = register_fpca(Y, npc = 2, family = "binomial",
#'                      fpca_type = "two-step", max_iterations = 5,
#'                      fpca_index_relevantDigits = 4)
#' 
#' ggplot(reg$Y, aes(x = tstar, y = t_hat, group = id)) +
#'   geom_line(alpha = 0.2) + ggtitle("Estimated warping functions")
#' 
#' plot(reg$fpca_obj, response_function = function(x) { 1 / (1 + exp(-x)) })
#' 
#' 
#' \donttest{ 
#' # example using accelerometer data from nhanes 2003-2004 study
#' data(nhanes)
#' reg_nhanes = register_fpca(nhanes, npc = 2, family = "binomial", max_iterations = 5)
#' 
#' 
#' ### incomplete Gaussian curves
#' data(growth_incomplete)
#' 
#' # Force the warping functions to start and end on the diagonal
#' reg2a = register_fpca(growth_incomplete, npc = 2, family = "gaussian",
#'                       incompleteness = NULL, max_iterations = 5)
#' ggplot(reg2a$Y, aes(x = tstar, y = t_hat, group = id)) +
#'   geom_line(alpha = 0.2) +
#'   ggtitle("Estimated warping functions")
#' ggplot(reg2a$Y, aes(x = t_hat, y = value, group = id)) +
#'   geom_line(alpha = 0.2) +
#'   ggtitle("Registered curves")
#' 
#' # Allow the warping functions to not start / end on the diagonal.
#' # The higher lambda_inc, the more the starting points and endpoints are forced
#' # towards the diagonal.
#' reg2b = register_fpca(growth_incomplete, npc = 2, family = "gaussian",
#'                       incompleteness = "full", lambda_inc = 0.1,
#'                       max_iterations = 5)
#' ggplot(reg2b$Y, aes(x = tstar, y = t_hat, group = id)) +
#'   geom_line(alpha = 0.2) +
#'   ggtitle("Estimated warping functions")
#' ggplot(reg2b$Y, aes(x = t_hat, y = value, group = id)) +
#'   geom_line(alpha = 0.2) +
#'   ggtitle("Registered curves")
#' 
#' ### complete Gamma curves
#' Y             = simulate_unregistered_curves(I = 20, D = 100)
#' Y$value       = exp(Y$latent_mean)
#' registr_gamma = register_fpca(Y, npc = 2, family = "gamma", fpca_type = "two-step",
#'                               gradient = FALSE, max_iterations = 5)
#' }
#'
register_fpca = function(Y, Kt = 8, Kh = 4, family = "gaussian",
												 incompleteness = NULL, lambda_inc = NULL,
												 Y_template = NULL,
												 max_iterations = 10, npc = 1,
												 fpca_type = "variationalEM", fpca_maxiter = 50,
												 fpca_seed = 1988, fpca_error_thresh = 0.0001,
												 fpca_index_relevantDigits = 4L, cores = 1L, ...){
	
  if (!(family %in% c("gaussian","binomial","gamma"))) {
  	stop("Package currently handles only families 'gaussian', 'binomial' and 'gamma'.")
  }
	
	if (family == "gamma" && fpca_type == "variationalEM") {
		warning("fpca_type = 'variationalEM' is only available for families 'gaussian' and 'binomial'. Setting fpca_type = 'two-step'.")
		fpca_type = "two-step"
	}
		
  data = data_clean(Y)
  Y = data$Y
  rows = data$Y_rows
  
  time_warps = list(NA, max_iterations + 2)
  time_warps[[1]] = Y$index
  loss = rep(NA, max_iterations + 1)

  # first register values to the overall mean
  registr_step = registr(Y = Y, Kt = Kt, Kh = Kh, family = family,
  											 incompleteness = incompleteness,
  											 lambda_inc     = lambda_inc,
  											 Y_template     = Y_template,
  											 row_obj = rows, cores = cores, ...)
  time_warps[[2]] = registr_step$Y$index
  loss[1] = registr_step$loss
  
  iter = 1
  error = rep(NA, max_iterations)
  error[iter] = 100
  while( iter <= max_iterations && error[iter] > 0.01 ){
  	message("current iteration: ", iter)
  	
  	if (fpca_type == "variationalEM") { # GFPCA after Wrobel et al. (2019)
  		if (family == "binomial") {
  			fpca_step = bfpca(registr_step$Y, npc = npc, Kt = Kt, row_obj = rows, seed = fpca_seed, maxiter = fpca_maxiter, 
  												error_thresh = fpca_error_thresh, ...)
  		} else if (family == "gaussian") {
  			fpca_step = fpca_gauss(registr_step$Y, npc = npc, Kt = Kt, row_obj = rows, seed = fpca_seed, maxiter = fpca_maxiter,
  														 error_thresh = fpca_error_thresh, ...)
  		}
  		
  	} else if (fpca_type == "two-step") { # Two-step GFPCA after Gertheiss et al. (2017)
  		
  		# estimation details
  		estimation_accuracy = ifelse(iter == 1, "high", "low")
  		if (iter == 1) {
  			gamm4_startParams = NULL
  		} else {
  			gamm4_startParams = fpca_step$gamm4_theta # parameters from last step
  		}
  		
  		fpca_step = gfpca_twoStep(registr_step$Y, family = family, npc = npc,
  															Kt = Kt, row_obj = rows,
  															index_relevantDigits = fpca_index_relevantDigits,
  															estimation_accuracy  = estimation_accuracy,
  															start_params         = gamm4_startParams)
  	}
  	
  	registr_step = registr(obj = fpca_step, Kt = Kt, Kh = Kh, family = family, 
  												 incompleteness = incompleteness,
  												 lambda_inc     = lambda_inc,
  												 row_obj = rows, beta = registr_step$beta, cores = cores, ...)
  	
  	time_warps[[iter + 2]] = registr_step$Y$index
  	loss[iter + 1] = registr_step$loss
  	
  	## calculate error
  	error[iter + 1] = sum((time_warps[[iter + 2]]-time_warps[[iter + 1]])^2) 
  	iter = iter + 1
  }
  
  if(iter < max_iterations){
  	message("Registration converged.")
  } else{
  	warning("Convergence not reached. Try increasing max_iterations.")
  }

  # final fpca step
  if (fpca_type == "variationalEM") { # GFPCA after Wrobel et al. (2019)
  	if (family == "binomial") {
  		fpca_step = bfpca(registr_step$Y, npc = npc, Kt = Kt, row_obj = rows, seed = fpca_seed, maxiter = fpca_maxiter, 
  											error_thresh = fpca_error_thresh, ...)
  	} else if (family == "gaussian") {
  		fpca_step = fpca_gauss(registr_step$Y, npc = npc, Kt = Kt, row_obj = rows, seed = fpca_seed, maxiter = fpca_maxiter, 
  													 error_thresh = fpca_error_thresh, ...)
  	}
  	
  } else if (fpca_type == "two-step") { # Two-step GFPCA after Gertheiss et al. (2017)
  	fpca_step = gfpca_twoStep(registr_step$Y, family = family, npc = npc,
  														Kt = Kt, row_obj = rows,
  														index_relevantDigits = fpca_index_relevantDigits,
  														estimation_accuracy  = "high",
  														start_params         = fpca_step$gamm4_theta)
  }
  
  Y$tstar = time_warps[[1]]
  Y$t_hat = registr_step$Y$index
  
  beta = as.data.frame(t(registr_step$beta))
  beta$id = unique(Y$id)
  
	ret = list(fpca_obj = fpca_step, Y = Y, time_warps = time_warps[!is.na(time_warps)],
						 loss = loss[!is.na(loss)], family = family, beta = beta)
	class(ret) = "registration"
  return(ret)
}
