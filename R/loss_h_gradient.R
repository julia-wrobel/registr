#' Gradient of loss function for registration step
#'
#' @param Y vector of observed points.
#' @param basis.tstar B-spline basis for vector Y.
#' @param mean.coefs spline coefficient vector for mean curve.
#' @param knots knot locations for B-spline basis used to estimate mean and FPC basis function.
#' @param beta.inner spline coefficient vector to be estimated for warping function h.
#' @param family \code{gaussian} or \code{binomial}.
#' @param t.min minimum value to be evaluated on the time domain (useful if data are sparse and / or irregular). 
#' @param t.max maximum value to be evaluated on the time domain (useful if data are sparse and / or irregular). 
#' 
#' @author Julia Wrobel \email{jw3134@@cumc.columbia.edu}
#' @export
#'
loss_h_gradient = function(Y, basis.tstar, mean.coefs, knots, beta.inner, family = "gaussian",
                           t.min, t.max){
  
  Beta = c(t.min, beta.inner, t.max)
  B.tstar = cbind(1, basis.tstar)
  htstar = B.tstar %*% Beta
  mean.coefs = matrix(mean.coefs, ncol = 1)
  
  basis.h = bs(htstar, knots = knots, intercept = TRUE)
  basis.h.deriv = bs_deriv(htstar, knots)
  
  basis.h.quad = crossprod(basis.h, basis.h.deriv)
  
  mu.t = basis.h %*% mean.coefs
  mu.t.deriv = basis.h.deriv %*% mean.coefs
  
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
      mu.t.deriv.temp = basis.h.deriv[t,] %*% mean.coefs
      grad_B[[t]] = crossprod(mu.t.deriv.temp) * B.tstar[t,] 
    }
    grad_B_sum = Reduce("+", grad_B)
    grad = grad_A = grad_B_sum
  }
  
  grad.last = length(grad)
  grad.inner = grad[-c(1, grad.last)]
  return(-1 * grad.inner)
}