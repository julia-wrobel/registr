initial_params = function(warping, Kh, t_min, t_max, I){
	if(warping == "nonparametric"){
		beta_new = matrix(NA, Kh - 2, I) 
		beta_0 = seq(t_min, t_max, length.out = Kh)[-c(1, Kh)] 
	} else if(warping == "piecewise_linear2"){
		beta_new = matrix(NA, 4, I)
		beta_0 = c(0.25, 0.3,  0.75, 0.8)
		rownames(beta_new) = c("knot1_x", "knot1_y", "knot2_x", "knot2_y")
	}
	
	return(list(beta_new = beta_new, beta_0 = beta_0))
}


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
#' @param knot_locations controls the x and y locations of the two knots.
#'  
#' @author Erin McDonnell \email{eim2117@@cumc.columbia.edu}
#' @importFrom stats quantile
piecewise_linear2_hinv = function(grid, knot_locations = c(0.25, 0.3, 0.75, 0.9)){
	knot1_x = knot_locations[1]
	knot1_y = knot_locations[2]
	knot2_x = knot_locations[3]
	knot2_y = knot_locations[4]
	
	beta1 = knot1_y / knot1_x
	beta2 = (knot2_y - knot1_y) / (knot2_x - knot1_x)
	beta3 = (1 - knot2_y) / (1 - knot2_x)
	
	h_inv1 = grid[grid <  knot1_x]*beta1
	h_inv2 = knot1_x*beta1 + (grid[grid >= knot1_x & grid < knot2_x] - knot1_x)*beta2
	h_inv3 = knot1_x*beta1 + (knot2_x - knot1_x)*beta2 + (grid[grid >= knot2_x] - knot2_x)*beta3
	
	return(c(h_inv1, h_inv2, h_inv3))
}