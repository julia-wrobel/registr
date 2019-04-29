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
#' @param parametric_warps If FALSE (default), inverse warping functions are 
#' estimated nonparametrically. If 'beta_cdf', they are assumed to have the form of a 
#' Beta(a,b) CDF. If 'piecewise' they follow a piecewise parameterized function. 
#' If 'piecewise_linear1' they follow a piecewise linear function with 1 knot.
#' If 'piecewise_linear2' they follow a piecewise linear function with 2 knots.
#' @param periodic If TRUE, uses periodic b-spline basis functions. Default is FALSE.
#' @param prior_1_x For parametric_warps = "piecewise_linear1" or "piecewise_linear2" only. 
#' If TRUE, will incorporate a prior Normal distribution for the first knot's x location into the loss function.
#' @param prior_1_x_mean Mean of the Normal distribution prior for the first knot's x location. 
#' @param prior_1_x_sd Standard deviation of the Normal distribution prior for the first knot's x location. 
#' @param prior_1_y For parametric_warps = "piecewise_linear1" or "piecewise_linear2" only. 
#' If TRUE, will incorporate a prior Normal distribution for the first knot's y location into the loss function.
#' @param prior_1_y_mean Mean of the Normal distribution prior for the first knot's y location. 
#' @param prior_1_y_sd Standard deviation of the Normal distribution prior for the first knot's y location. 
#' @param prior_2_x For parametric_warps = "piecewise_linear2" only. 
#' If TRUE, will incorporate a prior Normal distribution for the second knot's x location into the loss function.
#' @param prior_2_x_mean Mean of the Normal distribution prior for the second knot's x location. 
#' @param prior_2_x_sd Standard deviation of the Normal distribution prior for the second knot's x location. 
#' @param prior_2_y For parametric_warps = "piecewise_linear2" only. 
#' If TRUE, will incorporate a prior Normal distribution for the second knot's y location into the loss function.
#' @param prior_2_y_mean Mean of the Normal distribution prior for the second knot's y location. 
#' @param prior_2_y_sd Standard deviation of the Normal distribution prior for the second knot's y location. 
#' @param ... additional arguments passed to or from other functions
#' 
#' @return An object of class \code{fpca} containing:
#' \item{Y}{The observed data. The variable index is the new estimated time domain.}
#' \item{loss}{Value of the loss function after registraton.}
#' \item{beta}{Matrix of B-spline basis coefficients used to construct subject-specific warping functions.}
#' 
#' @author Julia Wrobel \email{jw3134@@cumc.columbia.edu}
#' @export
#' 
#' @importFrom stats glm coef constrOptim quantile optim pbeta
#' @importFrom splines bs
#' @importFrom pbs pbs
#' 
#' @examples
#' 
#' \dontrun{
#' Y = simulate_unregistered_curves()
#' register_step = registr(obj = NULL, Y = Y, Kt = 6, Kh = 3, family = "binomial", 
#'    gradient = TRUE)
#' }
#'
registr = function(obj = NULL, Y = NULL, Kt = 8, Kh = 4, family = "binomial", gradient = TRUE,
									 beta = NULL, t_min = NULL, t_max = NULL, row_obj = NULL,
									 parametric_warps = FALSE, periodic = FALSE,
									 prior_1_x = FALSE, prior_1_x_mean = 0.5, prior_1_x_sd = 1,
									 prior_1_y = FALSE, prior_1_y_mean = 0.5, prior_1_y_sd = 1,
									 prior_2_x = FALSE, prior_2_x_mean = 0.5, prior_2_x_sd = 1,
									 prior_2_y = FALSE, prior_2_y_mean = 0.5, prior_2_y_sd = 1, ...){
	
	print(paste("registr() Periodic = ", periodic))
  
  if(is.null(Y)) { Y = obj$Y}
	if(is.null(obj)) { Y$tstar = Y$index
		# scale time to 0, 1 for parametric warping
		if(!(parametric_warps == FALSE) ){
			Y$tstar = Y$index_scaled
		}
	}
  
  if(is.null(row_obj)){
    data = data_clean(Y)
    Y = data$Y
    rows = data$Y_rows
    I = data$I
  }else{
    rows = row_obj
    I = dim(rows)[1]
  }
	
	if(Kh < 3){
		stop("Kh must be greater than or equal to 3.")
	}
	if(Kt < 3){
		stop("Kt must be greater than or equal to 3.")
	}
  
	tstar = Y$tstar
  if (is.null(t_min)) {t_min = min(tstar)}
  if (is.null(t_max)) {t_max = max(tstar)}
  Theta_h = bs(tstar, df = Kh, intercept = FALSE)
  
  if(is.null(obj)){
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
  }else{
    global_knots = obj$knots
    mean_coefs = obj$subject_coefs
  }

  #### Define optimization constraints
  constrs = constraints(Kh, t_min, t_max, parametric_warps = parametric_warps)
  ui = constrs$ui
  ci = constrs$ci

  ### Calculate warping functions  
  t_hat = rep(NA, dim(Y)[1])
  loss_subjects = rep(NA, I)
  
  beta_new = matrix(NA, Kh - 1, I)
  beta_0 = seq(t_min, t_max, length.out = Kh + 1)[-c(1, Kh + 1)]
  if(!(parametric_warps == FALSE) ){
  	beta_new = matrix(NA, 2, I)
  	rownames(beta_new) = c("a", "b")
  	beta_0 = c(1, 1)
  	if(parametric_warps == "piecewise"){
  		beta_0 = c(.1, 0.5)
  		rownames(beta_new) = c("beta", "midpoint_percentile")
  	}else if(parametric_warps == "piecewise_linear1"){
  		beta_0 = c(0.5, 0.5)
  		rownames(beta_new) = c("knot_x", "knot_y")
  	}else if(parametric_warps == "piecewise_linear2"){
  		beta_new = matrix(NA, 4, I)
  		beta_0 = c(0.25, 0.3,  0.75, 0.8)
  		rownames(beta_new) = c("knot1_x", "knot1_y", "knot2_x", "knot2_y")
  	}
  }
	
  for (i in 1:I) {

    subject_rows = rows$first_row[i]:rows$last_row[i]
    
    Yi = Y$value[subject_rows]
    Di = length(Yi)
    
    Theta_h_i = Theta_h[subject_rows ,]
    
    if (is.null(beta)) {beta_i = beta_0} else {beta_i = beta[, i]}
    if (is.null(obj)) {mean_coefs_i = mean_coefs} else {mean_coefs_i = mean_coefs[, i]}
    if (gradient & !periodic) {gradf = loss_h_gradient} else {gradf = NULL}

    if(parametric_warps == "beta_cdf"){
    	beta_optim = optim(beta_i, loss_h, Y = Yi, Theta_h = Theta_h_i,
    										 mean_coefs = mean_coefs_i,knots = global_knots,
    										 family = family, t_min = t_min, t_max = t_max,
    										 parametric_warps = parametric_warps, periodic = periodic, 
    										 Kt = Kt,
    										 lower = 0.001, 
    										 upper = 5,
    										 method = "L-BFGS-B")
    	
    	beta_new[,i] = beta_optim$par
    	t_hat[subject_rows] = pbeta(seq(t_min, t_max, length.out = Di), 
    															beta_new[1, i], beta_new[2, i])
    
    }else if(parametric_warps == "piecewise"){
    	## these are not ideal endpoints.
    	beta_optim = constrOptim(beta_i, loss_h, grad = gradf, ui = ui, ci = ci, Y = Yi, 
    													 Theta_h = Theta_h_i, mean_coefs = mean_coefs_i, 
    													 knots = global_knots,
    													 parametric_warps = parametric_warps, 
    													 periodic = periodic, Kt = Kt,
    													 family = family, t_min = t_min, t_max = t_max)
    	
    	# beta_optim = optim(beta_i, loss_h, Y = Yi, Theta_h = Theta_h_i,
    	# 									 mean_coefs = mean_coefs_i,knots = global_knots,
    	# 									 family = family, t_min = t_min, t_max = t_max,
    	# 									 parametric_warps = parametric_warps, 
    	#										 periodic = periodic, Kt = Kt,
    	# 									 lower = c(0.01, 0.11), 
    	# 									 upper = c(10, 0.99),
    	# 									 method = "L-BFGS-B")
    	
    	beta_new[,i] = beta_optim$par
    	t_hat[subject_rows] = piecewise_parametric_hinv(seq(0, t_max, length.out = Di),
    																									beta_new[1, i], beta_new[2, i])
    	
    }else if(parametric_warps == "piecewise_linear1"){
    	beta_optim = constrOptim(beta_i, loss_h, grad = gradf, ui = ui, ci = ci, Y = Yi, 
    													 Theta_h = Theta_h_i, mean_coefs = mean_coefs_i, 
    													 knots = global_knots,
    													 parametric_warps = parametric_warps, 
    													 periodic = periodic, Kt = Kt,
    													 family = family, t_min = t_min, t_max = t_max,
    													 prior_1_x = prior_1_x, prior_1_x_mean = prior_1_x_mean, prior_1_x_sd = prior_1_x_sd,
    													 prior_1_y = prior_1_y, prior_1_y_mean = prior_1_y_mean, prior_1_y_sd = prior_1_y_sd)
    	
    	beta_new[,i] = beta_optim$par
    	t_hat[subject_rows] = piecewise_linear1_hinv(seq(0, t_max, length.out = Di),
    																							 beta_new[1, i], beta_new[2, i])
    	
    }else if(parametric_warps == "piecewise_linear2"){
    	beta_optim = constrOptim(beta_i, loss_h, grad = gradf, ui = ui, ci = ci, Y = Yi, 
    													 Theta_h = Theta_h_i, mean_coefs = mean_coefs_i, 
    													 knots = global_knots,
    													 parametric_warps = parametric_warps, 
    													 periodic = periodic, Kt = Kt,
    													 family = family, t_min = t_min, t_max = t_max,
    													 prior_1_x = prior_1_x, prior_1_x_mean = prior_1_x_mean, prior_1_x_sd = prior_1_x_sd,
    													 prior_1_y = prior_1_y, prior_1_y_mean = prior_1_y_mean, prior_1_y_sd = prior_1_y_sd,
    													 prior_2_x = prior_2_x, prior_2_x_mean = prior_2_x_mean, prior_2_x_sd = prior_2_x_sd,
    													 prior_2_y = prior_2_y, prior_2_y_mean = prior_2_y_mean, prior_2_y_sd = prior_2_y_sd)
    	
    	beta_new[,i] = beta_optim$par
    	t_hat[subject_rows] = piecewise_linear2_hinv(seq(0, t_max, length.out = Di),
    																							 beta_new[1, i], beta_new[2, i],
    																							 beta_new[3, i], beta_new[4, i])
    	
    }else{
    	beta_optim = constrOptim(beta_i, loss_h, grad = gradf, ui = ui, ci = ci, 
    													 Y = Yi, 
    													 Theta_h = Theta_h_i, mean_coefs = mean_coefs_i, 
    													 knots = global_knots, 
    													 periodic = periodic, Kt = Kt,
    													 family = family, t_min = t_min, t_max = t_max)
    	
    	
    	beta_new[,i] = beta_optim$par
    	
    	beta_full_i = c(t_min, 	beta_new[,i], t_max)
    	t_hat[subject_rows] = cbind(1, Theta_h_i) %*% beta_full_i
    }
    
    
    loss_subjects[i] = beta_optim$value
  }
  Y$index = t_hat
	Y$tstar = tstar
  return(list(Y = Y, loss = sum(loss_subjects), beta = beta_new)) 
} 
