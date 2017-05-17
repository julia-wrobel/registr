#' Loss function for registration step optimization
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
#' @importFrom boot inv.logit
#' @export
#'

loss_h = function(Y, basis.tstar, mean.coefs, knots, beta.inner, family, t.min, t.max){
  
  Beta = c(t.min, beta.inner, t.max)
  
  htstar = cbind(1, basis.tstar) %*% Beta
  basis.h = bs(htstar, knots = knots, intercept = TRUE)
  mu.h = basis.h %*% mean.coefs
  
  
  if (family == "gaussian") {
    return(sum((Y - mu.h) ^ 2))
  } else if (family == "binomial") {
    pi.h = inv.logit(mu.h)
    return(-1 * sum(Y * log(pi.h) + (1 - Y) * log(1 - pi.h) ))
  }
}