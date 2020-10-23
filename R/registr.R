#' Register curves from exponential family
#' 
#' Function used in the registration step of an FPCA-based approach for 
#' registering exponential-family functional data, called by \code{\link{register_fpca}}. 
#' This method uses constrained optimization to estimate spline 
#' coefficients for warping functions, where the objective function for optimization comes from 
#' maximizing the EF likelihood subject to monotonicity constraints on the warping functions. 
#' You have to either specify \code{obj}, which is a fpca 
#' object from an earlier step, or \code{Y}, a dataframe in long format with variables 
#' \code{id}, \code{index}, and \code{value} to indicate subject IDs, times, and observations, 
#' respectively.
#' By specifying \code{cores > 1}, the registration call can be parallelized.
#' 
#' @param obj Current estimate of FPC object. 
#' Can be NULL only if Y argument is selected.
#' @param Y Dataframe. Should have values id, value, index.
#' @param Kt Number of B-spline basis functions used to estimate mean functions. Default is 8.
#' @param Kh Number of B-spline basis functions used to estimate warping functions \emph{h}. Default is 4.
#' @param family \code{gaussian} or \code{binomial}.
#' @param gradient if \code{TRUE}, uses analytic gradient to calculate derivative. 
#' If \code{FALSE}, calculates gradient numerically.
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
#' Y = simulate_unregistered_curves()
#' register_step = registr(obj = NULL, Y = Y, Kt = 6, Kh = 4, family = "binomial", 
#'    gradient = TRUE)
#' testthat::expect_error({
#' registr(obj = list(Y = Y), Kt = 6, Kh = 4, family = "binomial", 
#'    gradient = TRUE)
#' })
#' \donttest{
#' Y = simulate_unregistered_curves()
#' register_step = registr(obj = NULL, Y = Y, Kt = 6, Kh = 4, family = "binomial", 
#'    gradient = TRUE)
#' }
#'
registr = function(obj = NULL, Y = NULL, Kt = 8, Kh = 4, family = "binomial", gradient = TRUE,
									 beta = NULL, t_min = NULL, t_max = NULL, row_obj = NULL,
									 periodic = FALSE, warping = "nonparametric", cores = 1L, ...){
	
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
		Y = data$Y
		rows = data$Y_rows
		I = data$I
	} else{
		rows = row_obj
		I = nrow(rows)
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
  if (is.null(t_min)) {t_min = min(tstar)}
  if (is.null(t_max)) {t_max = max(tstar)}
  Theta_h = bs(tstar, df = Kh, intercept = TRUE) 

  if (is.null(obj)) {
    # define population mean
  	if(periodic){
  		# if periodic, then we want more global knots, because the resulting object from pbs 
  		# only has (knots+intercept) columns.
  		global_knots = quantile(tstar, probs = seq(0, 1, length = Kt+1))[-c(1, Kt+1)]
  		basis = pbs(c(t_min, t_max, tstar), knots = global_knots, intercept = TRUE)[-(1:2),]
  	}else{
  		# if not periodic, then we want fewer global knots, because the resulting object from bs
  		# has (knots+degree+intercept) columns, and degree is set to 3 by default.
  		global_knots = quantile(tstar, probs = seq(0, 1, length = Kt - 2))[-c(1, Kt - 2)]
  		basis =  bs(c(t_min, t_max, tstar), knots = global_knots, intercept = TRUE)[-(1:2),]
  	} 

    mean_coefs = coef(glm(Y$value ~ 0 + basis, family = family))
    rm(basis)
  } else {
  	# stopifnot(!all(c("knots", "subject_coefs") %in% names(obj)))
    global_knots = obj$knots
    mean_coefs = obj$subject_coefs
  }

  #### Define optimization constraints
  constrs = constraints(Kh, t_min, t_max, warping = warping)
  ui = constrs$ui
  ci = constrs$ci

  ### Calculate warping functions  
  t_hat = rep(NA, dim(Y)[1])
  loss_subjects = rep(NA, I)
  
  # initialize beta and create an empty matrix to store new values
  initial_beta = initial_params(warping = warping, Kh, t_min, t_max, I)
	
  arg_list <- list(obj          = obj,          Y            = Y,
  								 Kt           = Kt,           family       = family,
  								 gradient     = gradient,     beta         = beta,
  								 t_min        = t_min,        t_max        = t_max,
  								 rows         = rows,         periodic     = periodic,
  								 warping      = warping,      Theta_h      = Theta_h,
  								 initial_beta = initial_beta, global_knots = global_knots,
  								 mean_coefs   = mean_coefs,   ui           = ui,
  								 ci           = ci)
  
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
#' @author Alexander Bauer \email{alexander.bauer@@stat.uni-muenchen.de}
#' 
#' @importFrom stats constrOptim
#' 
registr_oneCurve <- function(i, arg_list, ...) {
	
	subject_rows = arg_list$rows$first_row[i]:arg_list$rows$last_row[i]
	
	Yi = arg_list$Y$value[subject_rows]
	Di = length(Yi)
	
	Theta_h_i = arg_list$Theta_h[subject_rows ,]
	
	if (is.null(arg_list$beta)) {beta_i = arg_list$initial_beta$beta_0} else {beta_i = arg_list$beta[, i]}
	if (is.null(arg_list$obj)) {mean_coefs_i = arg_list$mean_coefs} else {mean_coefs_i = arg_list$mean_coefs[, i]}
	if (arg_list$gradient) {gradf = loss_h_gradient} else {gradf = NULL}
	
	beta_optim = constrOptim(theta      = beta_i,
													 f          = loss_h,
													 grad       = gradf,
													 ui         = arg_list$ui,
													 ci         = arg_list$ci,
													 Y          = Yi, 
													 Theta_h    = Theta_h_i,
													 mean_coefs = mean_coefs_i, 
													 knots      = arg_list$global_knots, 
													 family     = arg_list$family,
													 t_min      = arg_list$t_min,
													 t_max      = arg_list$t_max, 
													 periodic   = arg_list$periodic,
													 Kt         = arg_list$Kt,
													 warping    = arg_list$warping,
													 ...)
	
	beta_new = beta_optim$par
	
	if (arg_list$warping == "nonparametric") {
		beta_full_i = c(arg_list$t_min, 	beta_new, arg_list$t_max)
		#t_hat = as.vector(cbind(1, Theta_h_i) %*% beta_full_i)
		t_hat = as.vector(Theta_h_i %*% beta_full_i)
		
	} else if (arg_list$warping == "piecewise_linear2") {
		t_hat = piecewise_linear2_hinv(grid = seq(0, arg_list$t_max, length.out = Di),
																	 knot_locations = beta_new)
	}
	
	return(list(beta_new = beta_new,
							t_hat    = t_hat,
							loss     = beta_optim$value))
	
}
