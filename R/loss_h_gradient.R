#' Gradient of loss function for registration step
#'
#' @param Y vector of observed points.
#' @param Theta_phi B-spline basis for vector Y.
#' @param mean.coefs spline coefficient vector for mean curve.
#' @param knots knot locations for B-spline basis used to estimate mean and FPC basis function.
#' @param beta.inner spline coefficient vector to be estimated for warping function h.
#' @param family \code{gaussian} or \code{binomial}.
#' @param t_min minimum value to be evaluated on the time domain (useful if data are sparse and / or irregular). 
#' @param t_max maximum value to be evaluated on the time domain (useful if data are sparse and / or irregular). 
#' 
#' @author Julia Wrobel \email{jw3134@@cumc.columbia.edu}
#' @export
#'
loss_h_gradient = function(Y, Theta_phi, mean.coefs, knots, beta.inner, family = "gaussian",
                           t_min, t_max){
  
  beta = c(t_min, beta.inner, t_max)
  B.tstar = cbind(1, Theta_phi)
  htstar = B.tstar %*% beta
  mean.coefs = matrix(mean.coefs, ncol = 1)
  
  Theta_h = bs(htstar, knots = knots, intercept = TRUE)
  Theta_h_deriv = bs_deriv(htstar, knots)
  
  Theta_h_quad = crossprod(Theta_h, Theta_h_deriv)
  
  mu.t = Theta_h %*% mean.coefs
  mu.t.deriv = Theta_h_deriv %*% mean.coefs
  
  if(family == "binomial"){
    grad = crossprod(B.tstar, (Y * mu.t.deriv)) - 
      t(B.tstar) %*% (1/(1 + exp(mu.t)) * exp(mu.t) * mu.t.deriv )
  }else{
    D = length(Y)
    Y.diag = diag(Y)
    #grad = -2 * t(B.tstar) %*% (Y.diag %*% mu.t.deriv) + 2 * rowSums(t(B.tstar) %*% tcrossprod(mu.t.deriv))
    grad_A =  -2 * t(B.tstar) %*% (Y.diag %*% mu.t.deriv)
    grad_B = list(NA, D)
    for(t in 1:D){
      mu.t.deriv.temp = Theta_h_deriv[t,] %*% mean.coefs
      grad_B[[t]] = crossprod(mu.t.deriv.temp) * B.tstar[t,] 
    }
    grad_B_sum = Reduce("+", grad_B)
    grad = grad_A = grad_B_sum
  }
  
  grad.last = length(grad)
  grad.inner = grad[-c(1, grad.last)]
  return(-1 * grad.inner)
}