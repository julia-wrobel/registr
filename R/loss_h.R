#' Loss function for registration step optimization
#'
#' @param Y vector of observed points.
#' @param Theta_h B-spline basis for inverse warping functions.
#' @param mean_coefs spline coefficient vector for mean curve.
#' @param knots knot locations for B-spline basis used to estimate mean and FPC basis function.
#' @param beta.inner spline coefficient vector to be estimated for warping function h.
#' @param family \code{gaussian} or \code{binomial}.
#' @param t_min minimum value to be evaluated on the time domain. 
#' @param t_max maximum value to be evaluated on the time domain. 
#' @param parametric_warps If FALSE (default), inverse warping functions are 
#' estimated nonparametrically. If 'beta_cdf', they are assumed to have the form of a 
#' Beta(a,b) CDF. If 'piecewise' they follow a piecewise parameterized function.
#' 
#' @importFrom stats plogis
#' @export
#'

loss_h = function(Y, Theta_h, mean_coefs, knots, beta.inner, family, t_min, t_max, 
									parametric_warps = FALSE, 
									prior_1_x = FALSE, prior_1_x_mean = 0.5, prior_1_x_sd = 1,
									prior_1_y = FALSE, prior_1_y_mean = 0.5, prior_1_y_sd = 1,
									prior_2_x = FALSE, prior_2_x_mean = 0.5, prior_2_x_sd = 1,
									prior_2_y = FALSE, prior_2_y_mean = 0.5, prior_2_y_sd = 1){
  
  if(parametric_warps == "beta_cdf"){
  	tstar = seq(0, 1, length.out = length(Y))
  	hinv_tstar = pbeta(tstar, beta.inner[1], beta.inner[2])
  	
  }else if(parametric_warps == "piecewise"){
  	# does not currently allow minimum values different from zero
  	tstar = seq(0, t_max, length.out = length(Y))
  	if(abs(beta.inner[1]) < 0.0001){beta.inner[1] = sign(beta.inner[1])*0.0001}
  	
  	hinv_tstar = piecewise_parametric_hinv(tstar, beta.inner[1], beta.inner[2])
  	
  }else if(parametric_warps == "piecewise_linear1"){
  	# does not currently allow minimum values different from zero
  	tstar = seq(0, t_max, length.out = length(Y))
  	hinv_tstar = piecewise_linear1_hinv(tstar, beta.inner[1], beta.inner[2])
  	
  }else if(parametric_warps == "piecewise_linear2"){
  	# does not currently allow minimum values different from zero
  	tstar = seq(0, t_max, length.out = length(Y))
  	hinv_tstar = piecewise_linear2_hinv(tstar, beta.inner[1], beta.inner[2], beta.inner[3], beta.inner[4])
  	
	}else{
  	beta = c(t_min, beta.inner, t_max)
  	hinv_tstar = cbind(1, Theta_h) %*% beta
  }
  
  Theta_phi = bs(hinv_tstar, knots = knots, intercept = TRUE)
  g_mu_t = Theta_phi %*% mean_coefs
  
  if (family == "gaussian") {
    return(sum((Y - g_mu_t) ^ 2))
  } else if (family == "binomial") {
    pi_h = plogis(g_mu_t)
    # Allows for a prior distribution on the knot locations
    loss = -1 * sum(Y * log(pi_h) + (1 - Y) * log(1 - pi_h))
  	if((parametric_warps == "piecewise_linear1" | parametric_warps == "piecewise_linear2") & prior_1_x == TRUE){
  		loss = loss + dnorm(x = beta.inner[1], mean = prior_1_x_mean, sd = prior_1_x_sd)
  	}
    if((parametric_warps == "piecewise_linear1" | parametric_warps == "piecewise_linear2") & prior_1_y == TRUE){
    	loss = loss + dnorm(x = beta.inner[2], mean = prior_1_y_mean, sd = prior_1_y_sd)
    }
    if(parametric_warps == "piecewise_linear2" & prior_2_x == TRUE){
    	loss = loss + dnorm(x = beta.inner[1], mean = prior_2_x_mean, sd = prior_2_x_sd)
    }
    if(parametric_warps == "piecewise_linear2" & prior_2_y == TRUE){
    	loss = loss + dnorm(x = beta.inner[2], mean = prior_2_y_mean, sd = prior_2_y_sd)
    }
	  return(loss)
  }
}
