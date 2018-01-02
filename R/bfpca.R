#' Binary functional principal components analysis
#' 
#' Function used in the FPCA step for registering binary functional data. 
#' This method uses a variational EM algorithm to estimates scores and principal components for 
#' binary functional data. The returned values subject_coef are subject specific means.
#'
#' @param Y input data. Can be long format dataframe or wide format matrix. If you provide a matrix
#' @param index defaults to NULL. Indicates which column is the x-axis if Y input is a dataframe.
#' @param id defaults to NULL. Indicates which column gives the subject ids if Y input is a dataframe.
#' @param npc defaults to 1. Number of principal components to calculate.
#' @param Kt number of B-spline basis functions used to estimate mean functions. Defaults to 10.
#' @param maxiter maximum number of iterations to perform.
#' @param t_min minimum value to be evaluated on the time domain (useful if data are sparse and / or irregular). 
#' if `NULL`, taken to be minimum observed value.
#' @param t_max maximum value to be evaluated on the time domain (useful if data are sparse and / or irregular). 
#' if `NULL`, taken to be maximum observed value.
#' @param print.iter prints current error and iteration
#' @param row_obj if NULL, the function cleans the data and calculates row indices. Keep this NULL if you are using 
#' standalone \code{register} function.
#' @param seed set seed for reproducibility. Seed value defaults to 1988.
#' @param ... additional arguments passed to or from other functions
#' 
#' @author Julia Wrobel \email{jw3134@@cumc.columbia.edu}
#' @importFrom splines bs
#' @importFrom stats rnorm quantile
#' 
#' @export
#'

bfpca <- function(Y,index = NULL, id = NULL, npc = 1, Kt = 10, maxiter = 50, t_min = NULL, t_max = NULL, 
                  print.iter = FALSE, row_obj= NULL, seed = 1988, ...){
   
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
  psi_coefs = matrix(rnorm(Kt * npc), Kt, npc) * 0.5
  alpha_coefs = matrix(coef(glm(Y$value ~ 0 + Theta_phi, family = "binomial")), Kt, 1)
  xi = matrix(rnorm(dim(Y)[1]), ncol = 1) * 0.5
  
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
    curr_iter = curr_iter + 1;
    error[curr_iter] = sum((psi_coefs-temp_psi_coefs)^2) + sum((alpha_coefs-temp_alpha_coefs)^2);

    temp_psi_coefs = psi_coefs
    temp_alpha_coefs = alpha_coefs
    
  } ## end while loop

  fits = rep(NA, dim(Y)[1])
  subject_coef = alpha_coefs + tcrossprod(psi_coefs, scores)
 
  ## vectorize
   for(i in 1:I){
    subject_rows = rows$first_row[i]:rows$last_row[i]
   
    fits[subject_rows] = Theta_phi[subject_rows, ] %*% subject_coef[,i]
   }
  
  fittedVals = data.frame(id = Y$id, index = Y$index, value = fits)
  
  ## mean and eigenfuntions will have same number of grid points as last subject
  Theta_phi_mean = bs(seq(t_min, t_max, length.out = Di), knots = knots, intercept = TRUE) 
  
  # orthogonalize eigenvectors and extract eigenvalues
  psi_svd = svd(Theta_phi_mean %*% psi_coefs)
  efunctions = psi_svd$u
  evalues = ( psi_svd$d ) ^ 2
  scores = scores %*% psi_svd$v
  
  ret = list(
    "knots" = knots, 
    "alpha" = Theta_phi_mean %*% alpha_coefs,#
    "mu" = Theta_phi_mean %*% alpha_coefs, # return this to be consistent with refund.shiny
    "efunctions" = efunctions, 
    "evalues" =  evalues,
    "npc" = npc,
    "scores" = scores,
    "error" = error,
    "subject_coefs" = subject_coef,
    "Yhat" = fittedVals, 
    "Y" = Y, #
    "family" = "binomial"
  )

  class(ret) = "fpca" 
  return(ret)
} 
