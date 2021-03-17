#' Gradient of loss function for registration step
#'
#' @param family One of \code{c("gaussian","binomial")}. Defaults to \code{"gaussian"}.
#' @param periodic If \code{TRUE}, uses periodic b-spline basis functions. Default is \code{FALSE}. 
#' \code{loss_h_gradient()} is currently only available for \code{periodic = FALSE}.
#' @param warping If \code{nonparametric} (default), inverse warping functions are estimated nonparametrically. 
#' If \code{piecewise_linear2} they follow a piecewise linear function with 2 knots.
#' \code{loss_h_gradient()} is currently only available for \code{warping = "nonparametric"}.
#' @inheritParams loss_h
#' 
#' @author Julia Wrobel \email{julia.wrobel@@cuanschutz.edu},
#' Alexander Bauer \email{alexander.bauer@@stat.uni-muenchen.de}
#' 
#' @importFrom stats plogis
#' @importFrom splines bs
#' 
#' @return A numeric vector of spline coefficients for the gradient of the loss function.
#' 
#' @export
#'
loss_h_gradient = function(Y, Theta_h, mean_coefs, knots, beta.inner, family = "gaussian",
                           incompleteness = NULL, lambda_inc = NULL,
                           t_min, t_max, t_min_curve, t_max_curve, Kt = 8, periodic = FALSE,
                           warping = "nonparametric"){
  
  if(periodic){
    stop("loss_h_gradient() only available for periodic = FALSE")
  }
  if(warping != "nonparametric"){
    stop("loss_h_gradient() only available for warping = nonparametric")
  }
  
	D_i        = length(Y)
	Kh         = dim(Theta_h)[2]
  mean_coefs = matrix(mean_coefs, ncol = 1)
  
	# get the registered t values
  if (is.null(incompleteness)) { # initial and final parameters are fixed
  	beta = c(t_min_curve, beta.inner, t_max_curve)
  } else if (incompleteness == "leading") { # final parameter is fixed
  	beta = c(beta.inner, t_max_curve)
  } else if (incompleteness == "trailing") { # initial parameter is fixed
  	beta = c(t_min_curve, beta.inner)
  } else if (incompleteness == "full") { # no parameter is fixed
  	beta = beta.inner
  }
	hinv_tstar = c(Theta_h %*% beta) # c() as a slightly faster version of as.vector()
  
	# evaluate the gradient at the registered t values
	Theta_phi       = splines::bs(c(t_min, t_max, hinv_tstar), 
	                              knots = knots, intercept = TRUE)[-(1:2),]
	boundary_knots  = c(t_min, t_max)
  Theta_phi_deriv = bs_deriv(hinv_tstar, knots, Boundary.knots = boundary_knots)
  if (family == "gaussian") {
  	varphi = 1
  	b_g_deriv = Theta_phi %*% mean_coefs
  } else if (family == "binomial") {
  	varphi = 1
  	b_g_deriv = plogis(Theta_phi %*% mean_coefs)
  } else {
  	stop("The gradient is currently only available for families 'gaussian' and 'binomial'.")
  }
  
  # gradient_mat = matrix(NA, Kh, D_i)
  # for(j in 1:D_i){
  # 	gradient_mat[, j] = (Y - b_g_deriv)[j] * (Theta_phi_deriv %*% mean_coefs)[j] * Theta_h[j,]
  # }
  gradient_mat = c((Y - b_g_deriv) * (Theta_phi_deriv %*% mean_coefs)) * Theta_h
  grad = colSums(gradient_mat)
  
  
  # derivatives of the penalization term
  pen_term = 0
  if (!is.null(incompleteness) && (lambda_inc != 0)) { # penalization
  	
    if (incompleteness == "leading") { # penalize the starting point only
    	theta_h_leading  = Theta_h[1,]
    	pen_term_raw     = 2 * (hinv_tstar[1] - t_min_curve) * theta_h_leading
    } else if (incompleteness == "trailing") { # penalize the endpoint only
    	theta_h_trailing  = Theta_h[D_i,]
    	pen_term_raw      = 2 * (hinv_tstar[D_i] - t_max_curve) * theta_h_trailing
    } else if (incompleteness == "full") { # penalize overall dilation
    	theta_h_leading  = Theta_h[1,]
    	theta_h_trailing = Theta_h[D_i,]
    	pen_term_raw     = 2 * ((hinv_tstar[D_i] - hinv_tstar[1]) - (t_max_curve - t_min_curve)) *
    		(theta_h_trailing - theta_h_leading)
    }
  	pen_term = lambda_inc * pen_term_raw
  }
  
  grad = 1/varphi * grad - length(Y) * pen_term
  
  if (is.null(incompleteness)) { # initial and final parameters are fixed
    grad.inner = grad[-c(1, length(grad))]
  } else if (incompleteness == "leading") { # final parameter is fixed
  	grad.inner = grad[-length(grad)]
  } else if (incompleteness == "trailing") { # initial parameter is fixed
    grad.inner = grad[-1]
  } else if (incompleteness == "full") { # no parameter is fixed
  	grad.inner = grad
  }
  
  return(-1 * grad.inner)
}
