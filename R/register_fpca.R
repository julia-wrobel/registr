#' Register curves from exponential family using constrained optimization and generalized FPCA
#'
#' Function combines constrained optimization and FPCA to estimate warping functions for exponential family data.
#'
#' @param Y dataframe. Should have values id, value, index.
#' @param Kt number of B-spline basis functions used to estimate mean functions. 
#' Defaults to 10.
#' @param Kh number of B-spline basis functions used to estimate warping functions \emph{h}. Defaults to 5.
#' @param family \code{gaussian} or \code{binomial}.
#' @param iterations number of iterations between fpca step and registration step.
#' @param npc defaults to 1. Number of principal components to calculate.
#' @param ... additional arguments passed to registr and fpca functions
#'
#' @author Julia Wrobel \email{jw3134@@cumc.columbia.edu}
#' @export
#' 
#' @return fpca_obj list of items from FPCA step
#' @return reg_object some registration stuff, should be cleaned up
#' @return time_warps list of time values for each iteration of the algorithm. time_warps[1] returns original (observed) time 
#' and time_warps[n] provides time values for the final iteration
#' @return loss Loss for each iteration of the algorithm. Loss is calculated in the registration step using an 
#' exponential family likelihood with natural parameter calculated in the FPCA step.
#' @return family \code{gaussian} or \code{binomial}.
#'  
#' @examples
#'
#' \dontrun{
#'  simulate_unregistered_curves()
#'  reg_sim = register_fpca(Y_sim, Kt = 8, Kh = 4, family = "binomial", iterations = 10, npc = 1)
#' }
#'
register_fpca <- function(Y, Kt = 10, Kh = 4, family = "binomial", iterations = 10, npc = 1, ...){
  # ... argument should take care of anything that has a default value, but I also should be change it if I want to
      # for example I should be able to put maxiter= 50 as an argument, if I want. Test this out.

  # should include all the returns for this function too
  # should have a plot that handles just registration, just fpca, and both

  ## clean data
  if( !(family %in% c("binomial", "gaussian")) ){
  	stop("Package currently handles only 'binomial' or 'gaussian' families.")
  }
	
  # save original tstar values and all other t values calculated
  time_warps = list(NA, iterations + 2)
  time_warps[[1]] = Y$index
  loss = rep(NA, iterations + 1)

  data = data_clean(Y)
  Y = data$Y
  rows = data$Y_rows
  

  # first register values to the overall mean
  registr_step = registr(Y = Y, Kt = Kt, Kh = Kh, family = family, row_obj = rows, ...)
  time_warps[[2]] = registr_step$Y$index
  loss[1] = registr_step$loss
  
 
  # iteratively do fpca and register to newly calculated subject-specific means.
  for(iter in 1:iterations){
    message("current iteration: ", iter)
  	
		if(family == "binomial"){
			fpca_step = bfpca(registr_step$Y, index = NULL, id = NULL, npc = npc, Kt = Kt, 
												row_obj = rows, seed = 1988 + iter, ...)
		}else if(family == "gaussian"){
			stop("'gaussian' family not yet implemented for fpca step")
		} 
  	
  	registr_step = registr(obj = fpca_step, Kt = Kt, Kh = Kh, family = family, 
    												 row_obj = rows, beta = registr_step$beta, ...)

    time_warps[[iter + 2]] = registr_step$Y$index
    loss[iter + 1] = registr_step$loss

  }
  
	ret = list(fpca_obj = fpca_step, reg_object = registr_step, time_warps = time_warps,
						 loss = loss, family = family)
	class(ret) <- "registration"
  return(ret)
} # end function
