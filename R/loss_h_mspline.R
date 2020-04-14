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
#' @return The scalar value taken by the loss function.
#' 
#' @importFrom stats plogis
#' @export
#'

loss_h_mspline = function(Y, Theta_h, mean_coefs, knots, beta.inner, family, t_min, 
                          t_max){
  
  beta = c(t_min, beta.inner, t_max)
  hinv_tstar = cbind(1, Theta_h) %*% beta
  
  Theta_phi = bs(hinv_tstar, knots = knots, intercept = TRUE)
  g_mu_t = Theta_phi %*% mean_coefs
  
  if (family == "gaussian") {
    return(sum((Y - g_mu_t) ^ 2))
  } else if (family == "binomial") {
    pi_h = plogis(g_mu_t)
    return(-1 * sum(Y * log(pi_h) + (1 - Y) * log(1 - pi_h) ))
  }
}