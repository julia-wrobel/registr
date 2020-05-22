#' Register curves using constrained optimization and GFPCA
#'
#' Function combines constrained optimization and FPCA to estimate warping functions for 
#' exponential family curves. The FPCA step is performed through the function 
#' \code{\link{bfpca}} if \code{family = "binomial"} or the function 
#' \code{\link{fpca_gauss}} if \code{family = "gaussian"}. Warping functions are calculated 
#' by the function \code{\link{registr}}.
#' 
#' Requires input data \code{Y} to be a dataframe in long format with variables 
#' \code{id}, \code{index}, and \code{value} to indicate subject IDs, times, and observations, 
#' respectively. The code calls two 
#'
#' @param Y Dataframe. Should have values id, value, index.
#' @param Kt Number of B-spline basis functions used to estimate mean functions. Defaults to 8.
#' @param Kh Number of B-spline basis functions used to estimate warping functions \emph{h}. Defaults to 4.
#' @param family \code{gaussian} or \code{binomial}.
#' @param max_iterations Number of iterations for overall algorithm. Defaults to 10.
#' @param npc Number of principal components to calculate. Defaults to 1. 
#' @param fpca_maxiter Number to pass to the \code{maxiter} argument of `bfpca()` or `fpca_gauss()`. Default is 50.
#' @param ... Additional arguments passed to registr and fpca functions.
#'
#' @author Julia Wrobel \email{jw3134@@cumc.columbia.edu}
#' Jeff Goldsmith \email{ajg2202@@cumc.columbia.edu}
#' @export
#' 
#' @return An object of class \code{registration} containing:
#' \item{fpca_obj}{List of items from FPCA step.}
#' \item{Y}{The observed data plus variables \code{t_star} and \code{t_hat} which are the
#' unregistered grid and registered grid, respectively.}
#' \item{time_warps}{List of time values for each iteration of the algorithm. 
#' \code{time_warps[1]} is the original (observed) time
#' and \code{time_warps[n]} provides time values for the nth iteration.}
#' \item{loss}{Loss for each iteration of the algorithm, calculated in the registration step using an 
#' exponential family likelihood with natural parameter from the FPCA step.}
#' @return family \code{gaussian} or \code{binomial}.
#'  
#' @examples
#'
#' 
#'  Y = simulate_unregistered_curves(I = 20, D = 200)
#'  registr_object = register_fpca(Y, family = "binomial", max_iterations = 5)
#' 
#' \donttest{ 
#'  # example using accelerometer data from nhanes 2003-2004 study
#'  data(nhanes)
#'  register_nhanes = register_fpca(nhanes, npc = 2, family = "binomial", max_iterations = 5)
#' }
#'
register_fpca <- function(Y, Kt = 8, Kh = 4, family = "binomial", max_iterations = 10, 
													npc = 1, fpca_maxiter = 50, ...){
	
  if( !(family %in% c("binomial", "gaussian")) ){
  	stop("Package currently handles only 'binomial' or 'gaussian' families.")
  }
		
  data = data_clean(Y)
  Y = data$Y
  rows = data$Y_rows
  
  time_warps = list(NA, max_iterations + 2)
  time_warps[[1]] = Y$index
  loss = rep(NA, max_iterations + 1)

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
  		fpca_step = bfpca(registr_step$Y, npc = npc, Kt = Kt, row_obj = rows, seed = 1988 + iter, maxiter = fpca_maxiter, ...)
  	}else if(family == "gaussian"){
  		fpca_step = fpca_gauss(registr_step$Y, npc = npc, Kt = Kt, row_obj = rows, seed = 1988 + iter, maxiter = fpca_maxiter, ...)
  	}
  	
  	registr_step = registr(obj = fpca_step, Kt = Kt, Kh = Kh, family = family, 
  												 row_obj = rows, beta = registr_step$beta, ...)
  	
  	time_warps[[iter + 2]] = registr_step$Y$index
  	loss[iter + 1] = registr_step$loss
  	
  	## calculate error
  	error[iter + 1] = sum((time_warps[[iter + 2]]-time_warps[[iter + 1]])^2) 
  	iter = iter + 1
  }
  
  if(iter < max_iterations){
  	message("Registration converged.")
  } else{
  	warning("Convergence not reached. Try increasing max_iterations.")
  }

  # final fpca step
  if(family == "binomial"){
  	fpca_step = bfpca(registr_step$Y,npc = npc, Kt = Kt, row_obj = rows, maxiter = fpca_maxiter, ...)
  }else if(family == "gaussian"){
  	fpca_step = fpca_gauss(registr_step$Y,npc = npc, Kt = Kt, row_obj = rows, maxiter = fpca_maxiter, ...)
  }
  
  Y$tstar = time_warps[[1]]
  Y$t_hat = registr_step$Y$index
  
  beta = as.data.frame(t(registr_step$beta))
  beta$id = unique(Y$id)
  
	ret = list(fpca_obj = fpca_step, Y = Y, time_warps = time_warps[!is.na(time_warps)],
						 loss = loss[!is.na(loss)], family = family, beta = beta)
	class(ret) <- "registration"
  return(ret)
}
