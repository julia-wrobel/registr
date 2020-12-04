#' Nth derivative of spline basis 
#' 
#' This function gets derivative of a spline basis. Adapted from \code{bs()} function in \code{splines} package.
#'
#' @param x a numeric vector of values at which to evaluate the B-spline functions or derivatives.
#' @param knots the internal breakpoints that define the spline.
#' @param degree degree of the piecewise polynomial—default is 3 for cubic splines.
#' @param Boundary.knots boundary points at which to anchor the B-spline basis. 
#' Set to [0,1] if you want this to be your domain.
#' @param derivative a positive integer value that specifies which derivative to take. Defaults to 1 for 1st derivative.
#' Value of 0 returns the original set of b-spline basis functions. 
#' @param intercept if TRUE, an intercept is included in the basis; default is TRUE
#' 
#' @author Julia Wrobel \email{julia.wrobel@@cuanschutz.edu}
#' @importFrom splines splineDesign
#' 
#' @return A matrix containing:
#' \item{basis}{A B-spline basis that can be used to approximate the derivative of a function.}
#' 
#' @export
#'

bs_deriv = function(x, knots, degree = 3L, Boundary.knots = range(x), derivative = 1, intercept = TRUE){
  order = 1L + as.integer(degree)
  D = length(x)
  
  derivatives = rep(derivative, D)
  Aknots = sort(c(rep(Boundary.knots, order), knots))
  
  basis = splineDesign(Aknots, x, order, derivs = derivatives, outer.ok = TRUE)
  if (!intercept) 
    basis = basis[, -1L, drop = FALSE]
  
  return(basis)
}
