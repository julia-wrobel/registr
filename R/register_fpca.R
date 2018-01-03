#' Register curves from exponential family using constrained optimization and generalized FPCA
#'
#' Function combines constrained optimization and FPCA to estimate warping functions for exponential family data.
#'
#' @param Y dataframe. Should have values id, value, index.
#' @param Kt number of B-spline basis functions used to estimate mean functions. 
#' Defaults to 10.
#' @param Kh number of B-spline basis functions used to estimate warping functions \emph{h}. Defaults to 5.
#' @param family \code{gaussian} or \code{binomial}.
#' @param max_iterations number of iterations between fpca step and registration step.
#' @param npc defaults to 1. Number of principal components to calculate.
#' @param ... additional arguments passed to registr and fpca functions
#'
#' @author Julia Wrobel \email{jw3134@@cumc.columbia.edu}
#' @export
#' 
#' @return fpca_obj list of items from FPCA step
#' @return Y dataframe of data plus unregistered grid t_star and registered grid t_hat
#' @return time_warps list of time values for each iteration of the algorithm. time_warps[1] returns original (observed) time 
#' and time_warps[n] provides time values for the final iteration
#' @return loss Loss for each iteration of the algorithm. Loss is calculated in the registration step using an 
#' exponential family likelihood with natural parameter calculated in the FPCA step.
#' @return family \code{gaussian} or \code{binomial}.
#'  
#' @examples
#'
#' \dontrun{
#'  Y = simulate_unregistered_curves(I = 50, D = 200)
#'  registr_object = register_fpca(Y, family = "binomial", max_iterations = 25)
#' }
#'
register_fpca <- function(Y, Kt = 10, Kh = 4, family = "binomial", max_iterations = 20, npc = 1, ...){
  ## clean data
  if( !(family %in% c("binomial", "gaussian")) ){
  	stop("Package currently handles only 'binomial' or 'gaussian' families.")
  }
	
  # save original tstar values and all other t values calculated
  time_warps = list(NA, max_iterations + 2)
  time_warps[[1]] = Y$index
  loss = rep(NA, max_iterations + 1)

  data = data_clean(Y)
  Y = data$Y
  rows = data$Y_rows

  # first register values to the overall mean
  registr_step = registr(Y = Y, Kt = Kt, Kh = Kh, family = family, row_obj = rows, ...)
  time_warps[[2]] = registr_step$Y$index
  loss[1] = registr_step$loss
  
  iter = 1
  error = rep(NA, max_iterations)
  error[iter] = 100
  while( iter < max_iterations && error[iter] > 0.01 ){
  	message("current iteration: ", iter)
 
  	if(family == "binomial"){
  		fpca_step = bfpca(registr_step$Y, npc = npc, Kt = Kt, row_obj = rows, seed = 1988 + iter, ...)
  	}else if(family == "gaussian"){
  		stop("'gaussian' family not yet implemented for fpca step")
  	}
  	
  	registr_step = registr(obj = fpca_step, Kt = Kt, Kh = Kh, family = family, 
  												 row_obj = rows, beta = registr_step$beta, ...)
  	
  	time_warps[[iter + 2]] = registr_step$Y$index
  	loss[iter + 1] = registr_step$loss
  	
  	## calculate error
  	error[iter + 1] = sum((time_warps[[iter + 2]]-time_warps[[iter + 1]])^2) 
  	iter = iter + 1
  	
  }

  # final fpca step
  if(family == "binomial"){
  	fpca_step = bfpca(registr_step$Y,npc = npc, Kt = Kt, row_obj = rows)
  }else if(family == "gaussian"){
  	stop("'gaussian' family not yet implemented for fpca step")
  }
  
  Y$tstar = time_warps[[1]]
  Y$t_hat = registr_step$Y$index
  
	ret = list(fpca_obj = fpca_step, Y = Y, time_warps = time_warps[!is.na(time_warps)],
						 loss = loss[!is.na(loss)], family = family)
	class(ret) <- "registration"
  return(ret)
} # end function
