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
                           preserve_domain = TRUE, lambda_endpoint = NULL,
                           t_min, t_max, t_max_curve, Kt = 8, periodic = FALSE,
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
	if (preserve_domain) { # the warping function should end on the diagonal
	  beta = c(t_min, beta.inner, t_max)
	} else { # the warping function not necessarily ends on the diagonal
	  beta = c(t_min, beta.inner)
	}
	hinv_tstar = c(Theta_h %*% beta) # c() as a slightly faster version of as.vector()
  
	# evaluate the gradient at the registered t values
	Theta_phi       = splines::bs(c(t_min, t_max, hinv_tstar), 
	                              knots = knots, intercept = TRUE)[-(1:2),]
	boundary_knots  = c(t_min, t_max)
  Theta_phi_deriv = bs_deriv(hinv_tstar, knots)
  if (family == "binomial") {
  	varphi = 1
  	b_g_deriv = plogis(Theta_phi %*% mean_coefs)
  } else if (family == "gaussian") {
  	varphi = 1
  	b_g_deriv = Theta_phi %*% mean_coefs
  } else {
  	stop("Package currently handles only 'binomial' or 'gaussian' families.")
  }
  
  gradient_mat = matrix(NA, Kh, D_i)
  for(j in 1:D_i){
  	gradient_mat[, j] = (Y - b_g_deriv)[j] * (Theta_phi_deriv %*% mean_coefs)[j] * Theta_h[j,]
  }
  
  # derivatives of the penalization term
  if (preserve_domain || (lambda_endpoint == 0)) { # no penalization
    pen_term = 0
  } else { # penalize the deviation of the endpoint from the diagonal
    theta_h  = Theta_h[D_i,]
    pen_term = 2 * (hinv_tstar[D_i] - t_max_curve) * theta_h
    pen_term = lambda_endpoint * pen_term
  }
  
  grad = 1/varphi * rowSums(gradient_mat) - pen_term
  
  if (preserve_domain) { # the warping function should end on the diagonal
    grad.inner = grad[-c(1, length(grad))]
  } else { # the warping function not necessarily ends on the diagonal
    grad.inner = grad[-1]
  }
  
  return(-1 * grad.inner)
}
