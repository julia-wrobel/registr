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
#' library(tidyverse)
#' 
#' library(gridExtra)
#' library(refund.shiny)
#'
#' library(mvtnorm)
#' library(splines)
#' library(boot)
#'
#' ## design elements for simulated data
#' I = 50                            ## number of subjects
#' D = 100                           ## size of grid for observations
#' tstar = seq(0, 1, length = D)     ## shared time grid across subjects
#' Kt = 10                           ## dimension of bspline basis for mean and fpc functions
#' Kt_h = 5                          ## dimension of spline basis for warping functions
#' lambda.true = c(20)               ## variance of fpca scores
#'
#' ## functions to create mean and fpca functions on specified grid
#' mean.curve = function(grid) {
#'   3*(0 - sin(2*grid*pi) - cos(2*grid*pi) )
#' }
#' amp.curve = function(grid) {
#'   (0 - sin(2*grid*pi) - cos(2*grid*pi) )  / sqrt(322)
#' }
#'
#' ## function to create subject-specific time grid
#' grid.subj.create = function(coefs) {
#'   BS = bs(1:D, df = 3, intercept = FALSE, degree = 3)
#'   coefs = cumsum(coefs) / sum(coefs)
#'   (BS %*% coefs)
#' }
#'
#'
#' ###############################################################
#' ## registration for binary curves using logistic models for
#' ## fpca and registration
#' ###############################################################
#'
#' ## generate features for curves
#' c.true = rmvnorm(I, mean = rep(0, 1), sigma = diag(lambda.true, 1, 1))
#'
#' ## generate curves
#' Yi.obs = Yi.latent = pi.true = t.subj = matrix(NA, I, D)
#' Yi.regis.true = matrix(NA, I, D)
#' for (i in 1:I) {
#'   t.subj[i,] = runif(3, 0, 1) %>% grid.subj.create %>% as.vector
#'   Yi.latent[i,] = mean.curve(grid = t.subj[i,]) + c.true[i] * amp.curve(grid = t.subj[i,])
#'   pi.true[i,] = inv.logit(Yi.latent[i,])
#'   Yi.regis.true[i,] = mean.curve(grid = tstar) + c.true[i] * amp.curve(grid = tstar)
#'   for (j in 1:D) {
#'     Yi.obs[i,j] = rbinom(1, 1, pi.true[i,j])
#'   }
#' }
#'
#'  Y_sim = as_refundObj(Yi.obs)
#'  reg_sim = register_fpca(Y_sim, Kt = 8, Kh = 4, family = "binomial", iterations = 10, npc = 1)
#' }
#'
register_fpca <- function(Y, Kt = 10, Kh = 4, family = "binomial", iterations = 10, npc = 1, ...){
  # ... argument should take care of anything that has a default value, but I also should be change it if I want to
      # for example I should be able to put maxiter= 50 as an argument, if I want. Test this out.

  # should include all the returns for this function too
  # should have a plot that handles just registration, just fpca, and both

  ## clean data
  
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
  	
  	fpca_step = bfpca(registr_step$Y, index = NULL, id = NULL, npc = npc, Kt = Kt, 
  										row_obj = rows, seed = 1988 + iter, ...)
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
