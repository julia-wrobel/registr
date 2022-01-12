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
#' The number of functional principal components (FPCs) can either be specified
#' directly (argument \code{npc}) or chosen based on the explained share of
#' variance in each iteration (argument \code{npc_criterion}). \cr \cr
#' By specifying \code{cores > 1} the registration call can be parallelized.
#' 
#' Requires input data \code{Y} to be a dataframe in long format with variables 
#' \code{id}, \code{index}, and \code{value} to indicate subject IDs, 
#' observation times on the domain, and observations, respectively.
#' 
#' One joint iteration consists of a GFPCA step and a registration step.
#' As preprocessing, one initial registration step is performed.
#' The template function for this registration step is defined by argument
#' \code{Y_template}.
#' After convergence or \code{max_iterations} is reached, one final GFPCA step
#' is performed.
#'
#' @param Kt Number of B-spline basis functions used to estimate mean functions
#' and functional principal components. Default is 8. If
#' \code{fpca_type = "variationalEM"} and \code{npc_criterion} is used,
#' \code{Kt} is set to 20.
#' @param family One of \code{c("gaussian","binomial","gamma","poisson")}.
#' Families \code{"gamma"} and \code{"poisson"} are only supported by
#' \code{fpca_type = "two-step"}. Defaults to \code{"gaussian"}.
#' @param Y_template Optional dataframe with the same structure as \code{Y}.
#' Only used for the initial registration step. If NULL,
#' curves are registered to the overall mean of all curves in \code{Y} as template function.
#' If specified, the template function is taken as the mean
#' of all curves in \code{Y_template}. Defaults to NULL.
#' @param max_iterations Number of iterations for overall algorithm. Defaults to 10.
#' @param npc,npc_criterion The number of functional principal components (FPCs)
#' has to be specified either directly as \code{npc} or based on their explained
#' share of variance. In the latter case, \code{npc_criterion} has to be set
#' to a number between 0 and 1. For \code{fpca_type = "two-step"}, it is also
#' possible to cut off potential tails of subordinate FPCs (see
#' \code{\link{gfpca_twoStep}} for details).
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
#' @param fpca_index_significantDigits Only used if \code{fpca_type = "two-step"}.
#' Positive integer \code{>= 2}, stating the number of significant digits to which
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
#' @importFrom magrittr %>%
#' 
#' @return An object of class \code{registration} containing:
#' \item{Y}{The observed data plus variables \code{t_star} and \code{t_hat} which are the
#' unregistered grid and registered grid, respectively.}
#' \item{fpca_obj}{List of items from FPCA step.}
#' \item{family}{Used exponential family.}
#' \item{index_warped}{List of the (warped) index values for each iteration.
#' Has \code{'convergence$iterations + 2'} elements since the first two elements
#' contain the original (observed) index and the warped index values from the
#' preprocessing registration step (see Details), respectively.}
#' \item{hinv_innerKnots}{List of inner knots for setting up the spline bases
#' for the inverse warping functions. Only contains \code{NULL} values for
#' \code{Kh <= 4}.}
#' \item{hinv_beta}{Matrix of B-spline basis coefficients used to construct the
#' subject-specific inverse warping functions. From the last performed
#' registration step. For details see \code{?registr}.}
#' \item{convergence}{List with information on the convergence of the joint
#' approach. Containing the following elements: \cr \cr
#' \emph{converged} \cr
#' Indicator if the joint algorithm converged or if not
#' (i.e., \code{max_iterations} was reached) \cr \cr
#' \emph{iterations} \cr
#' Number of joint iterations that were performed. \cr \cr
#' \emph{delta_index} \cr
#' Vector of mean squared differences between the (warped) index values
#' (scaled to [0,1] based on the size of the observed domain)
#' in the current and the previous iteration.
#' Convergence is reached if this measure drops below 0.0001. \cr \cr
#' \emph{registration_loss} \cr
#' Vector of the loss in each iteration of the algorithm.
#' Calculated in the registration step using the exponential family
#' likelihood with natural parameter from the FPCA step.
#' Has \code{'iterations + 1'} elements since the first element contains the
#' loss of the preprocessing registration step (see Details).
#' }
#'  
#' @examples
#' 
#' ### complete binomial curves
#' Y = simulate_unregistered_curves(I = 20, D = 200)
#' 
#' # estimation based on Wrobel et al. (2019)
#' reg = register_fpca(Y, npc = 2, family = "binomial",
#'                     fpca_type = "variationalEM", max_iterations = 5)

#' 
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   library(ggplot2)
#'   
#'   ggplot(reg$Y, aes(x = tstar, y = t_hat, group = id)) +
#'     geom_line(alpha = 0.2) + ggtitle("Estimated warping functions")
#'   
#'   plot(reg$fpca_obj, response_function = function(x) { 1 / (1 + exp(-x)) })
#' }
#' 
#' 
#' \donttest{
#' 
#' # estimation based on Gertheiss et al. (2017)
#' reg2 = register_fpca(Y, npc = 2, family = "binomial",
#'                      fpca_type = "two-step", max_iterations = 5,
#'                      fpca_index_significantDigits = 4)
#'                      
#' # example using accelerometer data from nhanes 2003-2004 study
#' data(nhanes)
#' nhanes_short = nhanes[nhanes$id %in% unique(nhanes$id)[1:5],]
#' reg_nhanes   = register_fpca(nhanes_short, npc = 2, family = "binomial", max_iterations = 5)
#' 
#' 
#' ### incomplete Gaussian curves
#' data(growth_incomplete)
#' 
#' # Force the warping functions to start and end on the diagonal
#' reg2a = register_fpca(growth_incomplete, npc = 2, family = "gaussian",
#'                       incompleteness = NULL, max_iterations = 5)
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   
#'   ggplot(reg2a$Y, aes(x = tstar, y = t_hat, group = id)) +
#'     geom_line(alpha = 0.2) +
#'     ggtitle("Estimated warping functions")
#'   ggplot(reg2a$Y, aes(x = t_hat, y = value, group = id)) +
#'     geom_line(alpha = 0.2) +
#'     ggtitle("Registered curves")
#' }
#' # Allow the warping functions to not start / end on the diagonal.
#' # The higher lambda_inc, the more the starting points and endpoints are forced
#' # towards the diagonal.
#' reg2b = register_fpca(growth_incomplete, npc = 2, family = "gaussian",
#'                       incompleteness = "full", lambda_inc = 0.1,
#'                       max_iterations = 5)
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   ggplot(reg2b$Y, aes(x = tstar, y = t_hat, group = id)) +
#'     geom_line(alpha = 0.2) +
#'     ggtitle("Estimated warping functions")
#'   ggplot(reg2b$Y, aes(x = t_hat, y = value, group = id)) +
#'     geom_line(alpha = 0.2) +
#'     ggtitle("Registered curves")
#' }
#' 
#' ### complete Gamma curves
#' Y             = simulate_unregistered_curves(I = 20, D = 100)
#' Y$value       = exp(Y$latent_mean)
#' registr_gamma = register_fpca(Y, npc = 2, family = "gamma", fpca_type = "two-step",
#'                               gradient = FALSE, max_iterations = 3)
#' }
#'
register_fpca = function(Y, Kt = 8, Kh = 4, family = "gaussian",
												 incompleteness = NULL, lambda_inc = NULL,
												 Y_template = NULL,
												 max_iterations = 10, npc = NULL, npc_criterion = NULL,
												 fpca_type = "variationalEM", fpca_maxiter = 50,
												 fpca_seed = 1988, fpca_error_thresh = 0.0001,
												 fpca_index_significantDigits = 4L, cores = 1L, 
												 verbose = 1, 
												 ...){
	
  index = NULL
  rm(list="index")
  if (!(family %in% c("gaussian","binomial","gamma","poisson"))) {
  	stop("Package currently handles only families 'gaussian', 'binomial', 'gamma' and 'poisson'.")
  }
	
	if (family %in% c("gamma","poisson") && fpca_type == "variationalEM") {
		warning("fpca_type = 'variationalEM' is only available for families 'gaussian' and 'binomial'. Calling variationalEM for 'gaussian' family.")
	}
  if (is.null(npc) & (is.null(npc_criterion) || length(npc_criterion) > 2)) {
    stop("Please either specify 'npc' or 'npc_criterion' appropriately.")
  }
  
  if (fpca_type == "variationalEM" & !is.null(npc_criterion)) {
    Kt = 20
  }
  
  if (verbose > 2) {
    message("Running data_clean")
  }
  data    = data_clean(Y)
  Y       = data$Y
  Y$tstar = Y$index
  rows    = data$Y_rows
  
  index_warped      = vector(mode = "list", length = max_iterations + 2)
  index_warped[[1]] = Y$index_scaled
  reg_loss          = rep(NA, max_iterations + 1)

  # first register values to the overall mean
  if (verbose > 0) {
    message("Running initial registration step")
  }
  registr_step = registr(Y = Y, Kt = Kt, Kh = Kh, family = family,
  											 incompleteness = incompleteness,
  											 lambda_inc     = lambda_inc,
  											 Y_template     = Y_template,
  											 row_obj = rows, cores = cores,
  											 verbose = verbose,
  											 ...)
  index_warped[[2]] = registr_step$Y$index_scaled
  reg_loss[1]       = registr_step$loss
  
  iter        = 0
  delta_index = rep(NA, max_iterations)
  convergence_threshold = 0.0001
  
  # main iterations
  while (iter == 0 ||
  			 (iter < max_iterations && delta_index[iter] > convergence_threshold)) {
  	
  	iter = iter + 1
  	if (verbose > 0) {
  	  message("current iteration: ", iter)
  	}
  	if (verbose > 1) {
  	  message("FPCA Step")
  	}  	
  	if (fpca_type == "variationalEM") { # GFPCA after Wrobel et al. (2019)
  	  if (family == "binomial") {
  	    fpca_step = bfpca(registr_step$Y, npc = npc, npc_varExplained = npc_criterion[1],
  	                      Kt = Kt, row_obj = rows, seed = fpca_seed, maxiter = fpca_maxiter, 
  	                      error_thresh = fpca_error_thresh, verbose = verbose, ...)
  	  } else {
  	    fpca_step = fpca_gauss(registr_step$Y, npc = npc, npc_varExplained = npc_criterion[1],
  	                           Kt = Kt, row_obj = rows, seed = fpca_seed, maxiter = fpca_maxiter,
  	                           error_thresh = fpca_error_thresh,  verbose = verbose, ...)
  	  }
  		
  	} else if (fpca_type == "two-step") { # Two-step GFPCA after Gertheiss et al. (2017)
  		
  		# estimation details
  		estimation_accuracy = ifelse(iter == 1, "high", "low")
  		if (iter == 1) {
  			gamm4_startParams = NULL
  		} else {
  			gamm4_startParams = fpca_step$gamm4_theta # parameters from last step
  		}
  		
  		fpca_step = gfpca_twoStep(registr_step$Y, family = family,
  		                          npc = npc, npc_criterion = npc_criterion,
  															Kt = Kt, row_obj = rows,
  															index_significantDigits = fpca_index_significantDigits,
  															estimation_accuracy     = estimation_accuracy,
  															start_params            = gamm4_startParams,
  															verbose                 = verbose)
  	}
  	if (verbose > 1) {
  	  message("Registration step")
  	}  	
  	registr_step = registr(obj = fpca_step, Kt = Kt, Kh = Kh, family = family, 
  												 incompleteness = incompleteness,
  												 lambda_inc     = lambda_inc,
  												 row_obj        = rows,
  												 beta           = registr_step$hinv_beta,
  												 cores          = cores,
  												 verbose        = verbose,
  												 ...)
  	
  	index_warped[[iter + 2]] = registr_step$Y$index_scaled
  	reg_loss[iter + 1]       = registr_step$loss
  	
  	# calculate how much the warping functions changed since the last iteration
  	delta_index[iter] = mean((index_warped[[iter + 2]] - index_warped[[iter + 1]])^2) 
  } # END while loop
  
  converged = (delta_index[iter] <= convergence_threshold)
  
  if(iter < max_iterations){
  	if (verbose > 1) {
  	  message("Registration converged.")
  	}
  } else{
  	warning("Convergence not reached. Try increasing max_iterations.")
  }

  if (verbose > 0) {
    message("Running final FPCA step")
  }  	
  # final FPCA step
  if (fpca_type == "variationalEM") { # GFPCA after Wrobel et al. (2019)
    if (family == "binomial") {
      fpca_step = bfpca(registr_step$Y, npc = npc, npc_varExplained = npc_criterion[1],
                        Kt = Kt, row_obj = rows, seed = fpca_seed, maxiter = fpca_maxiter, 
                        error_thresh = fpca_error_thresh, verbose = verbose, ...)
    } else {
      fpca_step = fpca_gauss(registr_step$Y, npc = npc, npc_varExplained = npc_criterion[1],
                             Kt = Kt, row_obj = rows, seed = fpca_seed, maxiter = fpca_maxiter,
                             error_thresh = fpca_error_thresh,  verbose = verbose, ...)
    }
  	
  } else if (fpca_type == "two-step") { # Two-step GFPCA after Gertheiss et al. (2017)
  	fpca_step = gfpca_twoStep(registr_step$Y, family = family,
  	                          npc = npc, npc_criterion = npc_criterion,
  														Kt = Kt, row_obj = rows,
  														index_significantDigits = fpca_index_significantDigits,
  														estimation_accuracy     = "high",
  														start_params            = fpca_step$gamm4_theta,
  														verbose                 = verbose)
  	
  	# restrict the Yhat element to the individual registered domains
  	t_max_warped_i = sapply(rows$id, function(x) { max(registr_step$Y$index[registr_step$Y$id == x], na.rm = TRUE) })
  	fpca_step$Yhat = fpca_step$Yhat %>%
  		dplyr::group_by(id) %>%
  	  dplyr::filter(index <= (t_max_warped_i[match(id[1], rows$id)] + 10^(-fpca_index_significantDigits))) %>%
  		ungroup()
  }
  
  Y$t_hat = registr_step$Y$index
  
  ret = list(Y               = Y,
  					 fpca_obj        = fpca_step,
  					 family          = family,
  					 index_warped    = index_warped[!is.na(index_warped)],
  					 hinv_innerKnots = registr_step$hinv_innerKnots,
  					 hinv_beta       = registr_step$hinv_beta,
  					 convergence     = list(converged        = converged,
  					 											 iterations        = iter,
  					 											 delta_index       = delta_index[!is.na(delta_index)],
  					 											 registration_loss = reg_loss[!is.na(reg_loss)]))
  class(ret) = "registration"
  return(ret)
}
