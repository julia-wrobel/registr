#' Create two-parameter piecewise linear (inverse) warping functions
#'
#' This function uses a piecewise linear model to calculate inverse warping 
#' functions for #' registration. The parameter \code{knot_x} controls the 
#' x-location of the knot, #' and the parameter \code{knot_y} the y-location 
#' of the knot. The designation (inverse) is intended to communicate that these 
#' functions take data from the unregistered space to the registered space, 
#' consistent with #' functional data literature on registration.
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

piecewise_linear1_hinv = function(grid, knot_x = 0.5, knot_y = 0.5){
	
	beta1 = knot_y / knot_x
	beta2 = (1 - knot_y) / (1 - knot_x)
	
	h_inv1 = grid[grid <  knot_x]*beta1
	h_inv2 = knot_x*beta1 + (grid[grid >= knot_x] - knot_x)*beta2

	return(c(h_inv1, h_inv2))
}
