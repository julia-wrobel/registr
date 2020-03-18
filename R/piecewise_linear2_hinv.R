#' Create two-parameter piecewise linear (inverse) warping functions
#'
#' This function uses a 2-knot piecewise linear model to calculate inverse warping 
#' functions for registration. The parameters \code{knot1_x} and \code{knot1_y}
#' control the x and y locations of the first knot, and the parameters
#' \code{knot1_x} and \code{knot1_y} control the x and y locations of the second
#' knot. The designation (inverse) is intended to communicate that these 
#' functions take data from the unregistered space to the registered space, 
#' consistent with #' functional data literature on registration.
#' 
#' @param grid grid of values over which to evaluate the function.
#' @param knot1_x controls the x-location of the first knot. Defaults is \code{0.25}.
#' @param knot1_y controls the y-location of the first knot. Defaults is \code{0.3}. 
#' @param knot2_x controls the x-location of the second knot. Defaults is \code{0.75}. 
#' @param knot2_y controls the y-location of the second knot. Defaults is \code{0.9}. 
#'  
#' @author Erin McDonnell \email{eim2117@@cumc.columbia.edu}
#' @author Julia Wrobel \email{jw3134@@cumc.columbia.edu}
#' @importFrom stats quantile
#' @export

piecewise_linear2_hinv = function(grid, knot1_x = 0.25, knot1_y = 0.3, knot2_x = 0.75, knot2_y = 0.9){
	
	beta1 = knot1_y / knot1_x
	beta2 = (knot2_y - knot1_y) / (knot2_x - knot1_x)
	beta3 = (1 - knot2_y) / (1 - knot2_x)
	
	h_inv1 = grid[grid <  knot1_x]*beta1
	h_inv2 = knot1_x*beta1 + (grid[grid >= knot1_x & grid < knot2_x] - knot1_x)*beta2
	h_inv3 = knot1_x*beta1 + (knot2_x - knot1_x)*beta2 + (grid[grid >= knot2_x] - knot2_x)*beta3
	
	return(c(h_inv1, h_inv2, h_inv3))
}
