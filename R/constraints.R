#' Define constraints for optimization of warping functions
#' 
#' Constraints ensure monotonicity of spline coefficients for warping functions 
#' for use with \code{constrOptim()} function.  
#'
#' @param Kh Number of B-spline basis functions used to estimate warping functions \emph{h}.
#' @param t_min Minimum value to be evaluated on the time domain. 
#' @param t_max Maximum value to be evaluated on the time domain. 
#' @param warping If \code{nonparametric} (default), inverse warping functions are estimated nonparametrically. 
#' If \code{piecewise_linear2} they follow a piecewise linear function with 2 knots.
#' 
#' @author Julia Wrobel \email{julia.wrobel@@cuanschutz.edu},
#' Erin McDonnell \email{eim2117@@cumc.columbia.edu}
#' 
#' @return An list containing:
#' \item{ui}{A constraint matrix.}
#' \item{ci}{A constraint vector.}
#' 
#' @export
constraints = function(Kh, t_min = 0, t_max = 1, warping = "nonparametric"){
	if(warping == "nonparametric"){
		ui = matrix(0, nrow = Kh, ncol = Kh - 1) 
		ui[1,1] = 1; 
		ui[Kh, Kh - 1] = -1
		if(Kh > 3){
			for(i in 2:(Kh-1)){ 
				ui[i, (i-1):i] = c(-1,1)  
			}
		}
		ci = rep(0, Kh)
		
		ci[1]  = t_min  
		ci[Kh] = -t_max 
		
	} else if(warping == "piecewise_linear2"){
		ui = matrix(c(1, -1,  0,  0,  0,  0,
									0,  0,  0,  1, -1,  0,
									0,  1, -1,  0,  0,  0,
									0,  0,  0,  0,  1, -1), 6, 4) 
		ci = c(0.1, 0.1, -0.9, 0.1, 0.1, -0.9)
		
		# Additional constraints ensuring that the slopes can be no bigger than 5
		#		ui = matrix(c(1, -1,  0,  0,  0,  0,  5, -5,  0,
		#									0,  0,  0,  1, -1,  0, -1,  1,  0,
		#									0,  1, -1,  0,  0,  0,  0,  5, -5,
		#									0,  0,  0,  0,  1, -1,  0, -1,  1), 9, 4) 
		#		ci = c(0.1, 0.1, -0.9, 0.1, 0.1, -0.9, 0, 0, -4)
	}
  return(list(ui = ui, ci = ci))
}
