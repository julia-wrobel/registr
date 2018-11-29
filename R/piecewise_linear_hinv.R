#' Create two-parameter piecewise linear (inverse) warping functions
#'
#' This function uses a parametric model to calculate inverse warping functions for 
#' registration. The parameter \code{beta1} controls the slope from time 0 to some knot,
#' and the parameter \code{knot} controls where the knot takes place.
#' The designation (inverse) is intended to communicate that these functions take 
#' data from the unregistered space to the registered space, consistent with 
#' functional data literature on registration.
#' 
#' @param grid grid of values over which to evaluate the function.
#' @param beta1 parameter that controls the slope from time 0 to the knot. 
#' Result approaches identity warp as beta approaches 1.
#' @param knot controls where the result crosses the identity warp. 
#' Default is 0.5.
#'  
#' @author Julia Wrobel \email{jw3134@@cumc.columbia.edu}
#' @importFrom stats quantile
#' @export

piecewise_linear_hinv = function(grid, beta1 = 0.5, knot = 0.5){
	
	# knot*beta1 + (1 - knot)*beta2 = 1
	# (1 - knot)*beta2 = 1 - knot*beta1
	beta2 = (1 - knot*beta1) / (1 - knot)
	
	h_inv1 = grid[grid <  knot]*beta1
	h_inv2 = knot*beta1 + (grid[grid >= knot] - knot)*beta2
	
	print(c(h_inv1, h_inv2))

	return(c(h_inv1, h_inv2))
}
