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
#' @param ... additional arguments passed to or from other functions
#' 
#' @return An object of class \code{fpca} containing:
#' \item{Y}{The observed data. The variable index is the new estimated time domain.}
#' \item{loss}{Value of the loss function after registraton.}
#' \item{beta}{Matrix of B-spline basis coefficients used to construct subject-specific warping functions.}
#' 
#' @author Julia Wrobel \email{jw3134@@cumc.columbia.edu},
#' Erin McDonnell \email{eim2117@@cumc.columbia.edu}
#' @export
#' 
#' @importFrom stats glm coef constrOptim quantile optim pbeta
#' @importFrom splines bs
#' @importFrom pbs pbs
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
									 periodic = FALSE, warping = "nonparametric", ...){
	
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
  
  if(warping == "nonparametric"){
	  beta_new = matrix(NA, Kh - 2, I) 
	  beta_0 = seq(t_min, t_max, length.out = Kh)[-c(1, Kh)] 
  } else if(warping == "piecewise_linear2"){
  	beta_new = matrix(NA, 4, I)
  	beta_0 = c(0.25, 0.3,  0.75, 0.8)
  	rownames(beta_new) = c("knot1_x", "knot1_y", "knot2_x", "knot2_y")
  }
	
  for (i in 1:I) {

    subject_rows = rows$first_row[i]:rows$last_row[i]
    
    Yi = Y$value[subject_rows]
    Di = length(Yi)
    
    Theta_h_i = Theta_h[subject_rows ,]
    
    if (is.null(beta)) {beta_i = beta_0} else {beta_i = beta[, i]}
    if (is.null(obj)) {mean_coefs_i = mean_coefs} else {mean_coefs_i = mean_coefs[, i]}
    if (gradient) {gradf = loss_h_gradient} else {gradf = NULL}

  	beta_optim = constrOptim(beta_i, loss_h, grad = gradf, ui = ui, ci = ci, Y = Yi, 
  													 Theta_h = Theta_h_i, mean_coefs = mean_coefs_i, 
  													 knots = global_knots, 
  													 family = family, t_min = t_min, t_max = t_max, 
  													 periodic = periodic, Kt = Kt, warping = warping, ...)
  	
  	beta_new[,i] = beta_optim$par

    if(warping == "nonparametric"){
    	beta_full_i = c(t_min, 	beta_new[,i], t_max)
    	#t_hat[subject_rows] = cbind(1, Theta_h_i) %*% beta_full_i
    	t_hat[subject_rows] = Theta_h_i %*% beta_full_i
    } else if(warping == "piecewise_linear2"){
    	t_hat[subject_rows] = piecewise_linear2_hinv(seq(0, t_max, length.out = Di),
    																							 beta_new[1, i], beta_new[2, i],
    																							 beta_new[3, i], beta_new[4, i])
    }
    
    loss_subjects[i] = beta_optim$value
  }
  Y$index = t_hat
	Y$tstar = tstar
  return(list(Y = Y, loss = sum(loss_subjects), beta = beta_new)) 
} 
