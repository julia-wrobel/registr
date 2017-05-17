#' Define constraints for optimization of warping functions
#' 
#' Constraints ensure monotonicty of spline coefficients for warping functions for use with \code{constrOptim()} function.  
#'
#' @param Kh number of B-spline basis functions used to estimate warping functions \emph{h}.
#' @export
constraints = function(Kh){
  ui = matrix(0, nrow = Kh , ncol = Kh - 1) 
  ui[1,1] = 1; ui[Kh, Kh - 1] = -1
  for(i in 2:(Kh-1)){ ui[i, (i-1):i] = c(-1,1)  }
  
  ci = rep(0, Kh); ci[Kh] = -1  
  
  return(list(ui = ui, ci = ci))
}
