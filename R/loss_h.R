#' Loss function for registration step optimization
#'
#' @param Y vector of observed points.
#' @param Theta_h B-spline basis for inverse warping functions.
#' @param mean_coefs spline coefficient vector for mean curve.
#' @param knots knot locations for B-spline basis used to estimate mean and FPC basis function.
#' @param beta.inner spline coefficient vector to be estimated for warping function h.
#' @param family One of \code{c("gaussian","binomial","gamma","poisson")}.
#' For internal purposes, can also be set to \code{"gamma-varEM"} and
#' \code{"poisson-varEM"} if the preceding FPCA step in \code{register_fpca} was
#' performed with \code{fpca_type = "variationalEM"} which uses Gaussian family.
#' @param t_min,t_max minimum and maximum value to be evaluated on the time domain. 
#' @param t_min_curve,t_max_curve minimum and maximum value of the observed time domain of the
#' (potentially incomplete) curve.
#' @param warping If \code{nonparametric} (default), inverse warping functions are estimated nonparametrically. 
#' If \code{piecewise_linear2} they follow a piecewise linear function with 2 knots.
#' @param periodic If \code{TRUE} uses periodic b-spline basis functions. Default is \code{FALSE}.
#' @param Kt Number of B-spline basis functions used to estimate mean functions. Default is 8.
#' @param priors For \code{warping = "piecewise_linear2"} only. Logical indicator of whether to add Normal priors to pull the knots toward the identity line. 
#' @param prior_sd For \code{warping = "piecewise_linear2"} with \code{priors = TRUE} only. User-specified standard deviation for the Normal priors 
#' (single value applied to all 4 knot priors).
#' @inheritParams registr
#' 
#' @author Julia Wrobel \email{julia.wrobel@@cuanschutz.edu},
#' Erin McDonnell \email{eim2117@@cumc.columbia.edu},
#' Alexander Bauer \email{alexander.bauer@@stat.uni-muenchen.de}
#' 
#' @return The scalar value taken by the loss function.
#' 
#' @importFrom stats plogis
#' @importFrom splines bs
#' @importFrom pbs pbs
#' @importFrom stats dnorm
#' @importFrom utils tail
#' @export
#'
loss_h = function(Y, Theta_h, mean_coefs, knots, beta.inner, family, t_min, t_max, 
									t_min_curve, t_max_curve, incompleteness = NULL, lambda_inc = NULL,
									periodic = FALSE, Kt = 8, warping = "nonparametric", 
									priors = FALSE, prior_sd = NULL){
	
	if(warping != "piecewise_linear2" & priors == TRUE){
		stop("priors are only available for warping = piecewise_linear2")
	}
	if(family != "binomial" & priors == TRUE){
		stop("priors are only available for family = binomial")
	}
	if(priors == TRUE & is.null(prior_sd)){
		stop("priors = TRUE but no prior_sd supplied.")
	}
	if(priors == FALSE & !is.null(prior_sd)){
		warning("prior_sd supplied but priors = FALSE. No priors included.")
	}
	
	if (family %in% c("gamma","gamma-varEM")) { # the last element of beta.inner is the scale parameter
		scale      = tail(beta.inner, 1)
		beta.inner = beta.inner[1:(length(beta.inner) - 1)]
		if (scale <= 0) # ensure a positive scale value
			return(-1 * scale * 10^10)
	}
	
	
	# get the registered t values
	if (warping == "nonparametric") {
		if (is.null(incompleteness)) { # initial and final parameters are fixed
			beta = c(t_min_curve, beta.inner, t_max_curve)
		} else if (incompleteness == "leading") { # final parameter is fixed
			beta = c(beta.inner, t_max_curve)
		} else if (incompleteness == "trailing") { # initial parameter is fixed
			beta = c(t_min_curve, beta.inner)
		} else if (incompleteness == "full") { # no parameter is fixed
			beta = beta.inner
		}
  	#hinv_tstar = cbind(1, Theta_h) %*% beta
  	hinv_tstar = c(Theta_h %*% beta) # c() as a slightly faster version of as.vector()
  	
	} else if(warping == "piecewise_linear2"){
  	# does not currently allow minimum values different from zero
  	tstar = seq(0, t_max, length.out = length(Y))
  	#hinv_tstar = piecewise_linear2_hinv(tstar, beta.inner[1], beta.inner[2], beta.inner[3], beta.inner[4])
  	hinv_tstar = piecewise_linear2_hinv(tstar, beta.inner)
	}

	# evaluate the template curve at the registered t values
  if (periodic) {
  	Theta_phi = pbs::pbs(c(t_min, t_max, hinv_tstar),
  											 knots = knots, intercept = TRUE)[-(1:2),]
  } else {
  	Theta_phi = splines::bs(c(t_min, t_max, hinv_tstar),
  													knots = knots, intercept = TRUE)[-(1:2),]
  }
  g_mu_t = Theta_phi %*% mean_coefs

  # penalization term for the endpoint of the warping function
  pen_term = 0
  if (!is.null(incompleteness) && (lambda_inc != 0)) { # penalization
  	if (incompleteness == "leading") { # penalize the starting point only
  		pen_term_raw = (hinv_tstar[1] - t_min_curve)^2
  	} else if (incompleteness == "trailing") { # penalize the endpoint only
  		pen_term_raw = (utils::tail(hinv_tstar, 1) - t_max_curve)^2
  	} else if (incompleteness == "full") { # penalize overall dilation
  		pen_term_raw = ((utils::tail(hinv_tstar, 1) - hinv_tstar[1]) -
  											(t_max_curve - t_min_curve))^2
  	}
  	pen_term = lambda_inc * pen_term_raw
  }
  
  # Calculate the negative log likelihood
  if (family == "gaussian") { # drop unnecessary constant terms, assume variance = 1
  	loss = 1/2 * sum((Y - g_mu_t) ^ 2)
  	
  } else if (family == "binomial") {
    pi_h = plogis(g_mu_t)
    loss = -1 * sum(Y * log(pi_h) + (1 - Y) * log(1 - pi_h))
    
  } else if (family %in% c("gamma","gamma-varEM")) {
    mu_t <- if (family == "gamma") { as.vector(exp(g_mu_t)) } else { as.vector(g_mu_t) }
    if (family == "gamma-varEM") { # ensure positive values
      if (any(mu_t < 1e-8))
        mu_t[mu_t < 1e-8] <- 1e-8
    }
  	n     = length(Y)
  	alpha = mu_t * scale
  	loss  = -1 * (
  		sum((alpha - 1) * log(Y)) -
  			scale * sum(Y) +
  			log(scale) * scale * sum(mu_t) -
  			sum(log(gamma(alpha)))
  	)
  } else if (family %in% c("poisson","poisson-varEM")) {
    mu_t <- if (family == "poisson") { as.vector(exp(g_mu_t)) } else { as.vector(g_mu_t) }
    if (family == "poisson-varEM") { # ensure positive values
      if (any(mu_t < 1e-8))
        mu_t[mu_t < 1e-8] <- 1e-8
    }
    
    n    = length(Y)
  	loss = -1 * sum(Y * log(mu_t) - log(factorial(Y)) - mu_t)
  }

  # Add prior distributions on the knot locations as needed
	if(warping == "piecewise_linear2" & priors == TRUE){
		loss = loss - 
			log(dnorm(x = beta.inner[1], mean = 0.25, sd = prior_sd)) - 
			log(dnorm(x = beta.inner[2], mean = 0.25, sd = prior_sd)) - 
			log(dnorm(x = beta.inner[3], mean = 0.75, sd = prior_sd)) -
			log(dnorm(x = beta.inner[4], mean = 0.75, sd = prior_sd))
	}
  
  # compute the penalized log-likelihood
  loss_pen = loss + length(Y) * pen_term
  
  return(loss_pen)
}

