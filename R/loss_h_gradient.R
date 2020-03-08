#' Gradient of loss function for registration step
#'
#' @param Y vector of observed points.
#' @param Theta_h B-spline basis for inverse warping functions.
#' @param mean_coefs spline coefficient vector for mean curve.
#' @param knots knot locations for B-spline basis used to estimate mean and FPC basis function.
#' @param beta.inner spline coefficient vector to be estimated for warping function h.
#' @param family \code{gaussian} or \code{binomial}.
#' @param t_min minimum value to be evaluated on the time domain. 
#' @param t_max maximum value to be evaluated on the time domain. 
#' 
#' @author Julia Wrobel \email{jw3134@@cumc.columbia.edu}
#' 
#' @importFrom stats plogis
#' 
#' @return A numeric vector of spline coefficients for the gradient of the loss function.
#' 
#' @export
#'
loss_h_gradient = function(Y, Theta_h, mean_coefs, knots, beta.inner, family = "gaussian",
                           t_min, t_max){
  
	Di = length(Y)
	Kh = dim(Theta_h)[2]
  beta = c(t_min, beta.inner, t_max)

  Theta_h = cbind(1, Theta_h)
  
  hinv_tstar = Theta_h %*% beta
  mean_coefs = matrix(mean_coefs, ncol = 1)
  
  Theta_phi = bs(hinv_tstar, knots = knots, intercept = TRUE)
  Theta_phi_deriv = bs_deriv(hinv_tstar, knots)
  
  if(family == "binomial"){
  	varphi = 1
  	b_g_deriv = plogis(Theta_phi %*% mean_coefs)
  }else if (family == "gaussian"){
  	varphi = 1
  	b_g_deriv = Theta_phi %*% mean_coefs
  }else{
  	stop("Package currently handles only 'binomial' or 'gaussian' families.")
  }
  
  gradient_mat = matrix(NA, Kh + 1, Di)
  for(j in 1:Di){
  	gradient_mat[, j] = (Y - b_g_deriv)[j] * (Theta_phi_deriv %*% mean_coefs)[j] * Theta_h[j,]
  }
  
  grad = 1/varphi * rowSums(gradient_mat)
  
  grad.last = length(grad)
  grad.inner = grad[-c(1, grad.last)]
  return(-1 * grad.inner)
}