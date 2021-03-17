#' Create initial parameters for (inverse) warping functions
#' 
#' Dependent on the specific type of warping functions, this function creates
#' a vector of initial parameters. For \code{"nonparametric"} warpings that
#' are based on a given spline basis matrix, the initial parameters are defined
#' s.t. the resulting (inverse) warping function equals a diagonal line.
#' For \code{"piecewise_linear2"} warpings a fixed parameter vector is returned.
#' 
#' @param K Spline basis matrix defined over the interval \code{c(t_min, t_max)}.
#' @param t_vec Vector of the observed and potentially irregular time grid.
#' @inheritParams registr
#' 
#' @importFrom MASS ginv
#' 
initial_params = function(warping = "nonparametric", K, t_vec){

	if (warping == "nonparametric") {
		# generalized inverse of non-square matrix
		K_inv  = MASS::ginv(K)
		beta_0 = as.vector(K_inv %*% t_vec)
		
	} else if (warping == "piecewise_linear2") {
		beta_0 = c(0.25, 0.3,  0.75, 0.8)
	}
	
	return(beta_0)
}



#' Create two-parameter piecewise linear (inverse) warping functions
#'
#' This function uses a 2-knot piecewise linear model to calculate inverse warping 
#' functions for registration. The parameters \code{knot1_x} and \code{knot1_y}
#' control the x and y locations of the first knot, and the parameters
#' \code{knot1_x} and \code{knot1_y} control the x and y locations of the second
#' knot. The designation (inverse) is intended to communicate that these 
#' functions take data from the unregistered space to the registered space, 
#' consistent with functional data literature on registration.
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