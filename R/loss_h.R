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
#' @param prior_1_x For \code{warping = "piecewise_linear2"} only. Optional mean for a Normal distribution prior on the first knot's x location. 
#' Default prior means are on the identity line. To use these priors, users must specify a value for `prior_sd`. 
#' @param prior_1_y For \code{warping = "piecewise_linear2"} only. Optional mean for a Normal distribution prior on the first knot's y location. 
#' Default prior means are on the identity line. To use these priors, users must specify a value for `prior_sd`. 
#' @param prior_2_x For \code{warping = "piecewise_linear2"} only. Optional mean for a Normal distribution prior on the second knot's x location. 
#' Default prior means are on the identity line. To use these priors, users must specify a value for `prior_sd`. 
#' @param prior_2_y For \code{warping = "piecewise_linear2"} only. Optional mean for a Normal distribution prior on the second knot's y location. 
#' Default prior means are on the identity line. To use these priors, users must specify a value for `prior_sd`. 
#' @param prior_sd For \code{warping = "piecewise_linear2"} only. Optional sd to apply to all 4 Normal priors. Default is `prior_sd = Inf` (no priors).
#' 
#' @author Julia Wrobel \email{jw3134@@cumc.columbia.edu},
#' Erin McDonnell \email{eim2117@@cumc.columbia.edu}
#' 
#' @return The scalar value taken by the loss function.
#' 
#' @importFrom stats plogis
#' @importFrom splines bs
#' @importFrom pbs pbs
#' @importFrom stats dnorm
#' @export
#'

loss_h = function(Y, Theta_h, mean_coefs, knots, beta.inner, family, t_min, t_max, 
									periodic = FALSE, Kt = 8, warping = "nonparametric", 
									prior_1_x = 0.25, prior_1_y = 0.25, prior_2_x = 0.75, prior_2_y = 0.75,
									prior_sd = NULL){
	
	if(warping != "piecewise_linear2" & !is.null(prior_sd)){
		stop("priors are only available for warping = piecewise_linear2")
	}
	if(family != "binomial" & !is.null(prior_sd)){
		stop("priors are only available for family = binomial")
	}

	if(warping == "nonparametric"){
  	beta = c(t_min, beta.inner, t_max)
  	#hinv_tstar = cbind(1, Theta_h) %*% beta
  	hinv_tstar = Theta_h %*% beta ## changed
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
  
  # Calculate the negative log likelihood
  # For gaussian, drop unnecessary constant terms and assume variance = 1
  if (family == "gaussian") {
  	loss = 1/2 * sum((Y - g_mu_t) ^ 2)
  } else if (family == "binomial") {
    pi_h = plogis(g_mu_t)
    loss = -1 * sum(Y * log(pi_h) + (1 - Y) * log(1 - pi_h))
  }

  # Add prior distributions on the knot locations as needed
	if(warping == "piecewise_linear2" & !is.null(prior_sd)){
		if(!is.null(prior_1_x)){ loss = loss - log(dnorm(x = beta.inner[1], mean = prior_1_x, sd = prior_sd)) }
		if(!is.null(prior_1_y)){ loss = loss - log(dnorm(x = beta.inner[2], mean = prior_1_y, sd = prior_sd)) }
		if(!is.null(prior_2_x)){ loss = loss - log(dnorm(x = beta.inner[3], mean = prior_2_x, sd = prior_sd)) }
		if(!is.null(prior_2_y)){ loss = loss - log(dnorm(x = beta.inner[4], mean = prior_2_y, sd = prior_sd)) }
	}
  return(loss)
}

