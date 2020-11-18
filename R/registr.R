#' Register (in)complete curves from exponential family
#' 
#' Function used in the registration step of an FPCA-based approach for 
#' registering exponential-family, potentially incomplete functional data,
#' called by \code{\link{register_fpca}}. 
#' This method uses constrained optimization to estimate spline 
#' coefficients for warping functions, where the objective function for optimization comes from 
#' maximizing the EF likelihood subject to monotonicity constraints on the warping functions. 
#' You have to either specify \code{obj}, which is a fpca 
#' object from an earlier step, or \code{Y}, a dataframe in long format with variables 
#' \code{id}, \code{index}, and \code{value} to indicate subject IDs, times, and observations, 
#' respectively. \cr \cr
#' Warping functions by default are forced to end on the diagonal to be
#' domain-preserving. This behavior can be changed by setting
#' \code{preserve_domain = FALSE} and a reasonable \code{lambda_endpoint} value.
#' For further details see the accompanying vignette. \cr \cr
#' By specifying \code{cores > 1} the registration call can be parallelized.
#' 
#' @param obj Current estimate of FPC object. 
#' Can be NULL only if Y argument is selected.
#' @param Y Dataframe. Should have values id, value, index.
#' @param Kt Number of B-spline basis functions used to estimate mean functions. Default is 8.
#' @param Kh Number of B-spline basis functions used to estimate warping functions \emph{h}. Default is 4.
#' @param family \code{gaussian} or \code{binomial}.
#' @param gradient if \code{TRUE}, uses analytic gradient to calculate derivative. 
#' If \code{FALSE}, calculates gradient numerically.
#' @param preserve_domain Indicator if the registration should preserve the
#' time domain, leading to warping functions that end on the diagonal.
#' Defaults to \code{TRUE}. Can only be set to \code{FALSE} when
#' \code{warping = "nonparametric"}.
#' @param lambda_endpoint Penalization parameter to control the amount of
#' deviation of the warping function endpoints from the diagonal. The higher
#' this lambda, the more the endpoints are forced towards the diagonal.
#' Only used if \code{preserve_domain = FALSE}.
#' @param beta Current estimates for beta for each subject. Default is NULL. 
#' @param t_min Minimum value to be evaluated on the time domain.
#' if `NULL`, taken to be minimum observed value.
#' @param t_max Maximum value to be evaluated on the time domain. 
#' if `NULL`, taken to be maximum observed value.
#' @param row_obj If NULL, the function cleans the data and calculates row indices. 
#' Keep this NULL if you are using standalone \code{registr} function.
#' @param periodic If \code{TRUE}, uses periodic b-spline basis functions. Default is \code{FALSE}.
#' @param warping If \code{nonparametric} (default), inverse warping functions are estimated nonparametrically. 
#' If \code{piecewise_linear2} they follow a piecewise linear function with 2 knots.
#' @param cores Number of cores to be used. If \code{cores > 1}, the registration
#' call is parallelized by using \code{parallel::mclapply} (for Unix-based
#' systems) or \code{parallel::parLapply} (for Windows). Defaults to 1,
#' no parallelized call.
#' @param ... additional arguments passed to or from other functions
#' 
#' @return An object of class \code{fpca} containing:
#' \item{Y}{The observed data. The variable index is the new estimated time domain.}
#' \item{loss}{Value of the loss function after registraton.}
#' \item{beta}{Matrix of B-spline basis coefficients used to construct subject-specific warping functions.}
#' 
#' @author Julia Wrobel \email{julia.wrobel@@cuanschutz.edu},
#' Erin McDonnell \email{eim2117@@cumc.columbia.edu},
#' Alexander Bauer \email{alexander.bauer@@stat.uni-muenchen.de}
#' @export
#' 
#' @importFrom stats glm coef quantile pbeta
#' @importFrom splines bs
#' @importFrom pbs pbs
#' @importFrom parallel mclapply makePSOCKcluster clusterExport clusterEvalQ parLapply stopCluster
#' 
#' @examples
#' ### complete binomial curves
#' Y = simulate_unregistered_curves()
#' register_step = registr(obj = NULL, Y = Y, Kt = 6, Kh = 4, family = "binomial", 
#'                         gradient = TRUE)
#' \donttest{
#' ### incomplete Gaussian curves
#' data(growth_incomplete)
#' library(ggplot2)
#' 
#' # Force the warping functions to end on the diagonal
#' register_step2a = registr(obj = NULL, Y = growth_incomplete, Kt = 6, Kh = 4,
#'                           family = "gaussian", gradient = TRUE,
#'                           preserve_domain = TRUE)
#' ggplot(register_step2a$Y, aes(x = tstar, y = index, group = id)) +
#'   geom_line(alpha = 0.2) +
#'   ggtitle("Estimated warping functions")
#' 	
#' # Allow the warping functions to not end on the diagonal.
#' # The higher lambda_endpoint, the more the endpoints are forced towards the diagonal.
#' register_step2b = registr(obj = NULL, Y = growth_incomplete, Kt = 6, Kh = 4,
#'                           family = "gaussian", gradient = TRUE,
#'                           preserve_domain = FALSE, lambda_endpoint = 1)
#' ggplot(register_step2b$Y, aes(x = tstar, y = index, group = id)) +
#'   geom_line(alpha = 0.2) +
#'   ggtitle("Estimated warping functions")
#' }
#'
registr = function(obj = NULL, Y = NULL, Kt = 8, Kh = 4, family = "binomial", gradient = TRUE,
									 preserve_domain = TRUE, lambda_endpoint = NULL,
									 beta = NULL, t_min = NULL, t_max = NULL, row_obj = NULL,
									 periodic = FALSE, warping = "nonparametric", cores = 1L, ...){
	
	if (!preserve_domain) {
		if (warping != "nonparametric") {
			stop("The functionality for 'preserve_domain = FALSE' is only available for 'warping = 'nonparametric''")
		}
		if ((is.null(lambda_endpoint) || (lambda_endpoint < 0))) {
			stop("For 'preserve_domain = FALSE' the penalization parameter 'lambda_endpoint' has to be set to some nonnegative value.")
		}
	}
	
  if(is.null(Y)) { 
  	Y = obj$Y
  }
	
	if(is.null(obj)) { 
		if(warping == "nonparametric"){
			Y$tstar = Y$index
			
		} else if(warping == "piecewise_linear2"){
			# scale time to 0, 1 for parametric warping
			if(is.null(Y$index_scaled)){
				stop("For warping = piecewise_linear2, need an index_scaled variable that ranges from 0 to 1.")
			}
			Y$tstar = Y$index_scaled
		}
	}
  
	if (is.null(row_obj)) {
		data = data_clean(Y)
		Y    = data$Y
		rows = data$Y_rows
		I    = data$I
	} else{
		rows = row_obj
		I    = nrow(rows)
	}
	
	if (Kh < 4) {
		stop("Kh must be greater than or equal to 4.")
	}
	if (Kt < 4) {
		stop("Kt must be greater than or equal to 4.")
	}
	
	if(!(warping %in% c("nonparametric", "piecewise_linear2"))){
		stop("warping argument can only take values of nonparametric or piecewise_linear2")
	}

	if (gradient & periodic){
		warning("gradient = TRUE is only available for periodic = FALSE. Setting gradient = FALSE.")
		gradient = FALSE
	}
	
	if (gradient & warping != "nonparametric"){
		warning("gradient = TRUE is only available for warping = nonparametric. Setting gradient = FALSE.")
		gradient = FALSE
	}
	
	tstar = Y$tstar
  if (is.null(t_min)) { t_min = min(tstar) }
  if (is.null(t_max)) { t_max = max(tstar) }

	if (!is.null(obj)) { # template function = GFPCA representation
		global_knots = obj$knots
		mean_coefs   = obj$subject_coefs
		
	} else { # template function = mean of all curves
  	
  	if (periodic) {
  		# if periodic, then we want more global knots, because the resulting object from pbs 
  		# only has (knots+intercept) columns.
  		global_knots = quantile(tstar, probs = seq(0, 1, length = Kt+1))[-c(1, Kt+1)]
  		mean_basis   = pbs(c(t_min, t_max, tstar), knots = global_knots, intercept = TRUE)[-(1:2),]
  		
  	} else {
  		# if not periodic, then we want fewer global knots, because the resulting object from bs
  		# has (knots+degree+intercept) columns, and degree is set to 3 by default.
  		global_knots = quantile(tstar, probs = seq(0, 1, length = Kt - 2))[-c(1, Kt - 2)]
  		mean_basis   =  bs(c(t_min, t_max, tstar), knots = global_knots, intercept = TRUE)[-(1:2),]
  	} 

    mean_coefs = coef(glm(Y$value ~ 0 + mean_basis, family = family))
    rm(mean_basis)
  }

  ### Calculate warping functions  
  arg_list <- list(obj             = obj,             Y               = Y,
  								 Kt              = Kt,              Kh              = Kh,
  								 family          = family,          gradient        = gradient,
  								 preserve_domain = preserve_domain, lambda_endpoint = lambda_endpoint,
  								 beta            = beta,            t_min           = t_min,
  								 t_max           = t_max,           rows            = rows,
  								 periodic        = periodic,        warping         = warping,
  								 global_knots    = global_knots,    mean_coefs      = mean_coefs)
  
  # main function call
  if (cores == 1) { # serial call
  	results_list <- lapply(1:I, registr_oneCurve, arg_list, ...)
  	
  } else if (.Platform$OS.type == "unix") { # parallelized call on Unix-based systems
  	results_list <- parallel::mclapply(1:I, registr_oneCurve, arg_list, ...,
  																		 mc.cores = cores)
  	
  } else { # parallelized call on Windows
  	local_cluster <- parallel::makePSOCKcluster(rep("localhost", cores)) # set up cluster
  	# export functions and packages to the cluster
  	parallel::clusterExport(cl = local_cluster, c("arg_list"), envir=environment())
  	parallel::clusterEvalQ(cl = local_cluster, c(library(registr), library(stats)))
  	
  	results_list <- parallel::parLapply(cl  = local_cluster, X = 1:I,
  																			fun = registr_oneCurve, arg_list, ...)
  	
  	parallel::stopCluster(cl = local_cluster) # close cluster
  }
  
  # gather the results
  beta_new      = sapply(results_list, function(x) x$beta_new)
  t_hat         = unlist(sapply(results_list, function(x) as.vector(x$t_hat), simplify = FALSE))
  loss_subjects = unlist(sapply(results_list, function(x) as.vector(x$loss),  simplify = FALSE))
  
  Y$index = t_hat
	Y$tstar = tstar
	
  return(list(Y    = Y,
  						loss = sum(loss_subjects),
  						beta = beta_new)) 
} 



#' Internal function to register one curve
#' 
#' This internal function is only to be used from within \code{registr}.
#' It performs the main optimization step with \code{constrOptim} for the
#' registration of one curve.
#' 
#' @param i Numeric index of the curve under focus.
#' @param arg_list Named list of all arguments necessary for the registration
#' step.
#' @param ... additional arguments passed to or from other functions
#' 
#' @return An list containing:
#' \item{beta_new}{Estimated parameter vector of the warping function.}
#' \item{t_hat}{Vector of registered time domain.}
#' \item{loss}{Loss of the optimal solution.}
#' 
#' @author Julia Wrobel \email{julia.wrobel@@cuanschutz.edu},
#' Erin McDonnell \email{eim2117@@cumc.columbia.edu},
#' Alexander Bauer \email{alexander.bauer@@stat.uni-muenchen.de}
#' 
#' @importFrom pbs pbs
#' @importFrom splines bs
#' @importFrom stats constrOptim
#' 
registr_oneCurve <- function(i, arg_list, ...) {
	
	t_max_i       = arg_list$Y$tstar[arg_list$rows$last_row[i]]
	Y_cropped     = arg_list$Y[arg_list$Y$tstar <= t_max_i,]
	tstar_cropped = Y_cropped$tstar
	rows_i        = which(Y_cropped$id == arg_list$rows$id[i])
	Y_i           = Y_cropped$value[rows_i]
	D_i           = length(Y_i)
	
	# spline basis on the curve-specific time interval.
	# Pre-calculate the knots similarly to splines::bs to make the final splines::bs call faster.
	tstar_bs_i      = c(tstar_cropped[rows_i], range(tstar_cropped))
	degree          = 3 # cubic splines
	n_innerKnots    = arg_list$Kh - (1 + degree)
	if (n_innerKnots <= 0) {
		inner_knots = NULL
	} else if (n_innerKnots > 0) {
		p_quantiles = seq.int(from = 0, to = 1, length.out = n_innerKnots + 2)[-c(1,n_innerKnots + 2)]
		inner_knots = quantile(tstar_cropped, probs = p_quantiles)
	}
	Theta_h_cropped = splines::bs(tstar_bs_i, df = NULL, knots = inner_knots, intercept = TRUE)
	Theta_h_i       = Theta_h_cropped[1:length(rows_i),]
	
	# start parameters
	if (is.null(arg_list$beta)) { # newly initialize start parameters
		t_vec_i     = Y_cropped$tstar[rows_i]
		beta_full_i = initial_params(warping = arg_list$warping,
																 K       = Theta_h_i,
																 t_vec   = t_vec_i)
		
		if (arg_list$warping == "nonparametric") {
			if (arg_list$preserve_domain) { # final param is fixed and shouldn't be estimated
				beta_i = beta_full_i[-c(1, length(beta_full_i))]
			} else { # final param should be estimated
				beta_i = beta_full_i[-1]
			}
		} else if (arg_list$warping == "piecewise_linear2") {
			beta_i <- beta_full_i
		}
		
	} else { # use the last estimates (from the joint approach) as start parameters
		beta_i = arg_list$beta[, i]
	}
	
	if (is.null(arg_list$obj)) { # template function = mean of all curves
		mean_coefs_i = arg_list$mean_coefs
		
	} else { # template function = current GFPCA representation
		
		if (all(!is.na(arg_list$mean_coefs))) { # GFPCA based on fpca_gauss or bfpca
			# In this case, the FPCA is explicitly based on the spline basis 'mean_basis'
			mean_coefs_i = arg_list$mean_coefs[, i]
		} else { # GFPCA based on gfpca_twoStep
			# In this case, the FPCA is not based on the spline basis 'mean_basis'.
			# Accordingly, smooth over the GFPCA representation of the i'th function
			# using 'mean_basis'.
			mean_dat_i = arg_list$obj$Yhat[arg_list$obj$Yhat$id == i,]

			if (arg_list$periodic) {
				mean_basis = pbs::pbs(c(arg_list$t_min, arg_list$t_max, mean_dat_i$index),
															knots = arg_list$global_knots, intercept = TRUE)[-(1:2),]
			} else {
				mean_basis = splines::bs(c(arg_list$t_min, arg_list$t_max, mean_dat_i$index),
																 knots = arg_list$global_knots, intercept = TRUE)[-(1:2),]
			} 
			
			mean_coefs_i = coef(glm(value ~ 0 + mean_basis, data = mean_dat_i))
		}
	}
	
	if (arg_list$gradient) { gradf = loss_h_gradient } else { gradf = NULL }
	
	# optimization constraints
	if (arg_list$preserve_domain) { # warping functions should all end on the diagonal
		if (arg_list$warping == "nonparametric") {
			constrs_i = constraints(arg_list$Kh, arg_list$t_min, t_max_i, warping = arg_list$warping)
			ui_i      = constrs_i$ui
			ui_i      = ui_i[-nrow(ui_i), -ncol(ui_i)]
			ci_i      = constrs_i$ci
			ci_i      = ci_i[-(length(ci_i) - 1)]
		} else if (arg_list$warping == "piecewise_linear2") {
			constrs_i = constraints(arg_list$Kh - 1, arg_list$t_min, t_max_i, warping = arg_list$warping)
			ui_i      = constrs_i$ui
			ci_i      = constrs_i$ci
		}
		
	} else { # warping functions not necessarily end on the diagonal
		constrs_i <- constraints(arg_list$Kh, arg_list$t_min, arg_list$t_max)
		ui_i      <- constrs_i$ui
		ci_i      <- constrs_i$ci
	}
	
	# workaround: sometimes constrOptim states that the starting values are not
	# inside the feasible region. This is only a numerical error, so let's simply
	# substract a minor value from ci
	# (source: https://stackoverflow.com/questions/50472525/constroptim-in-r-init-val-is-not-in-the-interior-of-the-feasible-region-error)
	ci_i <- ci_i - 1e-6
	
	# when an analytic gradient is used, constrOptim sometimes leads to beta
	# values slightly outside the possible domain, or to slightly nonmonotone beta
	# values that don't fulfill the constraints.
	# Correct these slight inconsistencies to ensure proper beta values.
	if (arg_list$warping != "piecewise_linear2") {
		beta_i <- ensure_proper_beta(beta  = beta_i,
																 t_min = arg_list$t_min,
																 t_max = ifelse(arg_list$preserve_domain, t_max_i, arg_list$t_max))
	}
	
	# main registration step	
	beta_optim = constrOptim(theta           = beta_i,
													 f               = loss_h,
													 grad            = gradf,
													 ui              = ui_i,
													 ci              = ci_i,
													 Y               = Y_i, 
													 Theta_h         = Theta_h_i,
													 mean_coefs      = mean_coefs_i, 
													 knots           = arg_list$global_knots, 
													 family          = arg_list$family,
													 preserve_domain = arg_list$preserve_domain,
													 lambda_endpoint = arg_list$lambda_endpoint,
													 t_min           = arg_list$t_min,
													 t_max           = arg_list$t_max,
													 t_max_curve     = t_max_i,
													 periodic        = arg_list$periodic,
													 Kt              = arg_list$Kt,
													 warping         = arg_list$warping,
													 ...)
	
	beta_new = beta_optim$par
	
	if (arg_list$warping == "nonparametric") {
		if (arg_list$preserve_domain) { # final param is fixed and wasn't estimated
			beta_full_i = c(arg_list$t_min, beta_new, t_max_i)
		} else { # final param was estimated
			beta_full_i = c(arg_list$t_min, beta_new)
		}
		#t_hat = as.vector(cbind(1, Theta_h_i) %*% beta_full_i)
		t_hat = as.vector(Theta_h_i %*% beta_full_i)
		
	} else if (arg_list$warping == "piecewise_linear2") {
		t_hat = piecewise_linear2_hinv(grid = seq(0, t_max_i, length.out = D_i),
																	 knot_locations = beta_new)
	}
	
	return(list(beta_new = beta_new,
							t_hat    = t_hat,
							loss     = beta_optim$value))
	
}
