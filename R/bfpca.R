#' Binary functional principal components analysis
#' 
#' Function used in the FPCA step for registering binary functional data,
#' called by \code{\link{register_fpca}} when \code{family = "binomial"}. 
#' This method uses a variational EM algorithm to estimate scores and principal components for 
#' binary functional data.
#'
#' @param Y Dataframe. Should have variables id, value, index. 
#' @param npc Default is 1. Number of principal components to calculate.
#' @param Kt Number of B-spline basis functions used to estimate mean functions. Default is 8.
#' @param maxiter Maximum number of iterations to perform for EM algorithm. Default is 50.
#' @param t_min Minimum value to be evaluated on the time domain. 
#' @param t_max Maximum value to be evaluated on the time domain.
#' @param print.iter Prints current error and iteration
#' @param row_obj If NULL, the function cleans the data and calculates row indices. 
#' Keep this NULL if you are using standalone \code{register} function.
#' @param seed Set seed for reproducibility. Default is 1988.
#' @param ... Additional arguments passed to or from other functions
#' 
#' @author Julia Wrobel \email{jw3134@@cumc.columbia.edu},
#' Jeff Goldsmith \email{ajg2202@@cumc.columbia.edu}
#' @importFrom splines bs
#' @importFrom stats rnorm quantile
#' 
#' @return An object of class \code{fpca} containing:
#' \item{knots}{Cutpoints for B-spline basis used to rebuild \code{alpha}.}
#' \item{efunctions}{\eqn{D \times npc} matrix of estimated FPC basis functions.}
#' \item{evalues}{Estimated variance of the FPC scores.}
#' \item{npc}{number of FPCs.}
#' \item{scores}{\eqn{I \times npc} matrix of estimated FPC scores.}
#' \item{alpha}{Estimated population-level mean.}
#' \item{mu}{Estimated population-level mean. Same value as \code{alpha} but included for compatibility
#' with \code{refund.shiny} package.}
#' \item{subject_coefs}{B-spline basis coefficients used to construct subject-specific means. 
#' For use in \code{registr()} function.}
#' \item{Yhat}{FPC approximation of subject-specific means.}
#' \item{Y}{The observed data.}
#' \item{family}{\code{binomial}, for compatibility with \code{refund.shiny} package.}
#' \item{error}{vector containing error for each iteration of the algorithm.}
#' @export
#' @references Jaakkola, T. S. and Jordan, M. I. (1997).
##' A variational approach to Bayesian logistic regression models and their extensions. 
##' \emph{Proceedings of the Sixth International Workshop on Artificial Intelligence 
##' and Statistics}.
#' 
#' Tipping, M. E. (1999). Probabilistic Visualisation of High-dimensional binary data.
#' \emph{Advances in neural information processing systems}, 592--598.
#' 
#' @examples
##' Y = simulate_functional_data()$Y
##' bfpca_object = bfpca(Y, npc = 2, print.iter = TRUE)
#'
bfpca <- function(Y, npc = 1, Kt = 8, maxiter = 50, t_min = NULL, t_max = NULL, 
                  print.iter = FALSE, row_obj= NULL,
									seed = 1988, ...){
   
  curr_iter = 1
  error = rep(NA, maxiter)
  error[1] = 100.0
  
  ## clean data
  if(is.null(row_obj)){
    data = data_clean(Y)
    Y = data$Y
    rows = data$Y_rows
    I = data$I
  }else{
    rows = row_obj
    I = dim(rows)[1]
  }
  
  if(Kt < 3){
  	stop("Kt must be greater than or equal to 3.")
  }
  
  time = Y$index
  
  ## check that data is binary
  if(any( !(Y$value %in% c(0, 1)))){
  	stop("'binomial' family requires data with binary values of 0 or 1")
  }
  
  ## construct theta matrix
  if (is.null(t_min)) {t_min = min(time)}
  if (is.null(t_max)) {t_max = max(time)}
  
  knots = quantile(time, probs = seq(0, 1, length = Kt - 2))[-c(1, Kt - 2)]
  Theta_phi = bs(c(t_min, t_max, time), knots = knots, intercept = TRUE)[-(1:2),] 
  
  ## initialize all your vectors
  set.seed(seed)
  xi = matrix(rnorm(dim(Y)[1]), ncol = 1) * 0.5
  alpha_coefs = matrix(coef(glm(Y$value ~ 0 + Theta_phi, family = "binomial")), Kt, 1)
  psi_coefs = matrix(rnorm(Kt * npc), Kt, npc) * 0.5
  
  temp_alpha_coefs = alpha_coefs
  temp_psi_coefs = psi_coefs
  
  phi_a = list(NA, I)
  phi_b = matrix(0, nrow = Kt * (npc+1), ncol = I)
  scores = matrix(NA, I, npc)
  
  while(curr_iter < maxiter && error[curr_iter] > 0.001){

    if(print.iter){
      message("current iteration: ", curr_iter)
      message("current error: ", error[curr_iter])
    }
    
    for(i in 1:I){
    
      subject_rows = rows$first_row[i]:rows$last_row[i]
      
      Yi = Y$value[subject_rows]
      Di = length(Yi)
      Theta_i = Theta_phi[subject_rows, ] 
      xi_i = xi[subject_rows,]
      Theta_i_quad = squareTheta(xi_i, Theta_i)
       
      ##### posterior scores
      mlist = expectedScores(Yi, temp_alpha_coefs, temp_psi_coefs, Theta_i, Theta_i_quad)
      
      Ci = mlist$Ci
      mi = mlist$mi
      mm = Ci + tcrossprod(mi)
      
      ##### variational parameter xi
      xi[ subject_rows, 1] = expectedXi(Theta_i, temp_alpha_coefs, mi, temp_psi_coefs, Ci)
      xi_i = xi[subject_rows, ]
      
      Theta_i_quad = squareTheta(xi_i, Theta_i)
      
      #***************** estimate phi by maximum likelihood
      si = rbind(mi, 1)
      ss = cbind(rbind(mm, t(mi)), si)
      
      phi_a[[i]] = 2.0 * kronecker(Theta_i_quad, ss)  
      phi_b[,i] = t((Yi - 0.5) %*% kronecker(Theta_i, t(si)) )
      
      scores[i,] = mi
      
    } # end loop over I
    
    phi_a_sum = Reduce("+", phi_a)
    
    phi_vec = -solve(phi_a_sum) %*% rowSums(phi_b)
    phi_mat = matrix(phi_vec, nrow = Kt, ncol = npc + 1, byrow = TRUE)
    
    alpha_coefs = phi_mat[, npc+1]
    psi_coefs = phi_mat[, 1:npc]
    
    if(npc == 1){ psi_coefs = matrix(psi_coefs, ncol = 1)}
    
    ## calculate error
    curr_iter = curr_iter + 1
    error[curr_iter] = sum((psi_coefs-temp_psi_coefs)^2) + sum((alpha_coefs-temp_alpha_coefs)^2)

    temp_psi_coefs = psi_coefs
    temp_alpha_coefs = alpha_coefs
    
  } ## end while loop

  fits = rep(NA, dim(Y)[1])
  subject_coef = alpha_coefs + tcrossprod(psi_coefs, scores)
 
   for(i in 1:I){
    subject_rows = rows$first_row[i]:rows$last_row[i]
    fits[subject_rows] = Theta_phi[subject_rows, ] %*% subject_coef[,i]
   }
  
  fittedVals = data.frame(id = Y$id, index = Y$index, value = fits)
  
  ## mean and eigenfunctions will have same number of grid points as last subject
  Theta_phi_mean = bs(seq(t_min, t_max, length.out = Di), knots = knots, intercept = TRUE) 
  
  # orthogonalize eigenvectors and extract eigenvalues
  psi_svd = svd(Theta_phi_mean %*% psi_coefs)
  efunctions = psi_svd$u
  evalues = ( psi_svd$d ) ^ 2
  scores = scores %*% psi_svd$v
  
  ret = list(
    "knots" = knots, 
    "alpha" = Theta_phi_mean %*% alpha_coefs,
    "mu" = Theta_phi_mean %*% alpha_coefs, # return this to be consistent with refund.shiny
    "efunctions" = efunctions, 
    "evalues" =  evalues,
    "npc" = npc,
    "scores" = scores,
    "subject_coefs" = subject_coef,
    "Yhat" = fittedVals, 
    "Y" = Y, 
    "family" = "binomial",
    "error" = error[!is.na(error)]
  )

  class(ret) = "fpca" 
  return(ret)
} 
