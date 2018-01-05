#' Define constraints for optimization of warping functions
#' 
#' Constraints ensure monotonicity of spline coefficients for warping functions 
#' for use with \code{constrOptim()} function.  
#'
#' @param Kh Number of B-spline basis functions used to estimate warping functions \emph{h}.
#' @param t_min Minimum value to be evaluated on the time domain. 
#' @param t_max Maximum value to be evaluated on the time domain. 
#' @export
constraints = function(Kh, t_min = 0, t_max = 1){
  ui = matrix(0, nrow = Kh , ncol = Kh - 1) 
  ui[1,1] = 1; ui[Kh, Kh - 1] = -1
  for(i in 2:(Kh-1)){ ui[i, (i-1):i] = c(-1,1)  }
  
  ci = rep(0, Kh)
  
  ci[1] = t_min  
  ci[Kh] = -t_max  
  
  return(list(ui = ui, ci = ci))
}
