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
#' @param warping If \code{nonparametric} (default), inverse warping functions are estimated nonparametrically. 
#' If \code{piecewise_linear2} they follow a piecewise linear function with 2 knots.
#' @param periodic If \code{TRUE} uses periodic b-spline basis functions. Default is \code{FALSE}.
#' @param Kt Number of B-spline basis functions used to estimate mean functions. Default is 8.
#' @param prior_1_x For \code{warping = "piecewise_linear2"} only. If \code{TRUE}, 
#' will incorporate a prior Normal distribution for the first knot's x location into the loss function.
#' @param prior_1_x_mean Mean of the Normal distribution prior for the first knot's x location. 
#' @param prior_1_x_sd Standard deviation of the Normal distribution prior for the first knot's x location. 
#' @param prior_1_y For \codePwarping = "piecewise_linear2"} only. If \code{TRUE}, 
#' will incorporate a prior Normal distribution for the first knot's y location into the loss function.
#' @param prior_1_y_mean Mean of the Normal distribution prior for the first knot's y location. 
#' @param prior_1_y_sd Standard deviation of the Normal distribution prior for the first knot's y location. 
#' @param prior_2_x For \code{warping = "piecewise_linear2"} only. If \code{TRUE}, 
#' will incorporate a prior Normal distribution for the second knot's x location into the loss function.
#' @param prior_2_x_mean Mean of the Normal distribution prior for the second knot's x location. 
#' @param prior_2_x_sd Standard deviation of the Normal distribution prior for the second knot's x location. 
#' @param prior_2_y For \code{warping = "piecewise_linear2"} only. If TRUE, 
#' will incorporate a prior Normal distribution for the second knot's y location into the loss function.
#' @param prior_2_y_mean Mean of the Normal distribution prior for the second knot's y location. 
#' @param prior_2_y_sd Standard deviation of the Normal distribution prior for the second knot's y location. 
#'  
#' 
#' @importFrom stats plogis
#' @importFrom splines bs
#' @importFrom pbs pbs
#' @export
#'

loss_h = function(Y, Theta_h, mean_coefs, knots, beta.inner, family, t_min, t_max, 
									warping = "nonparametric", periodic = FALSE, Kt = 8,
									prior_1_x = FALSE, prior_1_x_mean = 0.5, prior_1_x_sd = 1,
									prior_1_y = FALSE, prior_1_y_mean = 0.5, prior_1_y_sd = 1,
									prior_2_x = FALSE, prior_2_x_mean = 0.5, prior_2_x_sd = 1,
									prior_2_y = FALSE, prior_2_y_mean = 0.5, prior_2_y_sd = 1){
	
	if(warping == "nonparametric"){
  	beta = c(t_min, beta.inner, t_max)
  	hinv_tstar = cbind(1, Theta_h) %*% beta
	} else if(warping == "piecewise_linear2"){
  	# does not currently allow minimum values different from zero
  	tstar = seq(0, t_max, length.out = length(Y))
  	hinv_tstar = piecewise_linear2_hinv(tstar, beta.inner[1], beta.inner[2], beta.inner[3], beta.inner[4])
	}

  if(periodic){
  	Theta_phi = pbs(hinv_tstar, knots = knots, intercept = TRUE)
  }else{
  	Theta_phi =  bs(hinv_tstar, knots = knots, intercept = TRUE)
  }
	
  g_mu_t = Theta_phi %*% mean_coefs
  
  if (family == "gaussian") {
    return(sum((Y - g_mu_t) ^ 2))
  } else if (family == "binomial") {
    pi_h = plogis(g_mu_t)
    # Allows for a prior distribution on the knot locations
    loss = -1 * sum(Y * log(pi_h) + (1 - Y) * log(1 - pi_h))

  	if(warping == "piecewise_linear2"){
  		if(prior_1_x == TRUE){
  			loss = loss - log(dnorm(x = beta.inner[1], mean = prior_1_x_mean, sd = prior_1_x_sd))
  		}
    	if(prior_1_y == TRUE){
    		loss = loss - log(dnorm(x = beta.inner[2], mean = prior_1_y_mean, sd = prior_1_y_sd))
    	}
    	if(prior_2_x == TRUE){
    		loss = loss - log(dnorm(x = beta.inner[3], mean = prior_2_x_mean, sd = prior_2_x_sd))
    	}
    	if(prior_2_y == TRUE){
    		loss = loss - log(dnorm(x = beta.inner[4], mean = prior_2_y_mean, sd = prior_2_y_sd))
    	}
  	}
	  return(loss)
  }
}

