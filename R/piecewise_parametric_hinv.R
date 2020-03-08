#' One parameter parametric warping on (0, T)
#'
#' @param grid grid of values over which to evaluate the function.
#' @param t_max maximum value to be evaluated on the time domain. 
#' @param beta parameter that controls shape of warping. Result 
#' approaches identity warp as beta approaches zero.
#' 
#' @return A numeric vector containing values for a single warping function. 
#' 
#' @examples 
#' x = runif(100)
#' plot(x, type = 'l')
#' lines(registr:::h_inv_parametric(grid = x, beta = 0.5), col = "red")
# reference the Marron paper
h_inv_parametric = function(grid, t_max = 1, beta = 0.01){
	1/beta * log( (grid * (exp(beta * t_max)-1) + t_max)/t_max  )
}


#' Create two-parameter piecewise (inverse) warping functions
#'
#' This function uses a parametric model to calculate inverse warping functions for 
#' registration. The parameter \code{beta} controls the shape of warping,
#' and the parameter \code{midpoint_percentile} control where the warping function
#' crosses the identity line. The designation (inverse) is intended to 
#' communicate that these functions take data from the unregistered space to the 
#' registered space, consistent with functional data literature on registration.
#' 
#' @param grid grid of values over which to evaluate the function.
#' @param beta parameter that controls shape of warping. Result 
#' approaches identity warp as beta approaches zero.
#' @param midpoint_percentile controls where the result crosses the identity warp. 
#' Default is 0.5, which forces result to cross identity line at median of grid.
#'  
#' @author Julia Wrobel \email{jw3134@@cumc.columbia.edu}
#' @importFrom stats quantile
#' 
#' @return A numeric vector containing values for a single warping function. 
#' 
#' @export
#' @examples 
#' x = runif(100)
#' plot(x, type = 'l')
#' lines(piecewise_parametric_hinv(grid = x, beta = 0.5), col = "red")
piecewise_parametric_hinv = function(grid, beta = 0.01, midpoint_percentile = 0.5){
	
	midpoint = quantile(grid, midpoint_percentile, type = 1)[[1]] 
	
	h_inv1 = h_inv_parametric(grid[grid <= midpoint], t_max = midpoint, beta = beta)
	h_inv2 = h_inv_parametric(grid[grid >= midpoint] - midpoint, 
														t_max = (max(grid) - midpoint), beta = -beta)
	
	if(midpoint_percentile == 0){
		return(h_inv2 + midpoint)
	}else{
		return(c(h_inv1, (h_inv2 + midpoint)[-1]))
	}
}