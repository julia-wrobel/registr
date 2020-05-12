#' Define constraints for optimization of warping functions
#' 
#' Constraints ensure monotonicity of spline coefficients for warping functions 
#' for use with \code{constrOptim()} function.  
#'
#' @param Kh Number of B-spline basis functions used to estimate warping functions \emph{h}.
#' @param t_min Minimum value to be evaluated on the time domain. 
#' @param t_max Maximum value to be evaluated on the time domain. 
#' @param parametric_warps If FALSE (default), inverse warping functions are 
#' estimated nonparametrically. If 'beta_cdf', they are assumed to have the form of a 
#' Beta(a,b) CDF. If 'piecewise' they follow a piecewise parameterized function.
#' 
#' @return An list containing:
#' \item{ui}{A constraint matrix.}
#' \item{ci}{A constraint vector.}
#' 
#' @export
constraints = function(Kh, t_min = 0, t_max = 1, parametric_warps = FALSE){
	if(parametric_warps == "piecewise"){
		ui = matrix(c(1,-1,0,0,0,0,1,-1), 4, 2)
		ci = c(-25, -25, 0, -0.99)
	}else{
		ui = matrix(0, nrow = Kh -1 , ncol = Kh - 2) 
		ui[1,1] = 1; 
		ui[Kh-1, Kh - 2] = -1
		if(Kh > 3){
			for(i in 2:(Kh-2)){ 
				ui[i, (i-1):i] = c(-1,1)  
			}
		}
		ci = rep(0, Kh-1)
		
		ci[1] = t_min  
		ci[Kh-1] = -t_max 
	}
  return(list(ui = ui, ci = ci))
}
