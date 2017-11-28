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
#' @param t.min minimum value to be evaluated on the time domain (useful if data are sparse and / or irregular). 
#' if `NULL`, taken to be minimum observed value.
#' @param t.max maximum value to be evaluated on the time domain (useful if data are sparse and / or irregular). 
#' if `NULL`, taken to be maximum observed value.
#' @param print.iter prints current error and iteration
#' @param row_obj if NULL, the function cleans the data and calculates row indices. Keep this NULL if you are using 
#' standalone \code{register} function.
#' 
#' @author Julia Wrobel \email{jw3134@@cumc.columbia.edu}
#' @importFrom splines bs
#' @importFrom stats rnorm quantile
#' 
#' @export
#'

bfpca <- function(Y,index = NULL, id = NULL, npc = 1, Kt = 10, maxiter = 20, t.min = NULL, t.max = NULL, 
                  print.iter = FALSE, row_obj= NULL){
   
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
  
  ## construct theta matrix
  if (is.null(t.min)) {t.min = min(time)}
  if (is.null(t.max)) {t.max = max(time)}
  
  
  knots = quantile(time, probs = seq(0, 1, length = Kt - 2))[-c(1, Kt - 2)]
  Theta = bs(c(t.min, t.max, time), knots = knots, intercept = TRUE)[-(1:2),] 

  ## initialize all your vectors
  mu_coef = matrix(rnorm(Kt), Kt, 1)
  psi_coef = matrix(rnorm(Kt * npc), Kt, npc) * 0.5
  xi = matrix(rnorm(dim(Y)[1]), ncol = 1) * 0.5
  
  temp_mu_coef = mu_coef
  temp_psi_coef = psi_coef
  
  phi_a = list(NA, I)
  phi_b = matrix(0, nrow = Kt * (npc+1), ncol = I)
  scores = matrix(NA, I, npc)
  
  while( curr_iter < maxiter){

    if(print.iter){
      message("current iteration: ", curr_iter)
      message("current error: ", error[curr_iter])
    }
    
    for(i in 1:I){
    
      subject_rows = rows$first_row[i]:rows$last_row[i]
      
      Yi = Y$value[subject_rows]
      Di = length(Yi)
      Theta_i = Theta[subject_rows, ] 
      xi_i = xi[subject_rows,]
      Theta_i_quad = squareTheta(xi_i, Theta_i)
       
      ##### posterior scores
      mlist = expectedScores(Yi, temp_mu_coef, temp_psi_coef, Theta_i, Theta_i_quad)
      
      Ci = mlist$Ci
      mi = mlist$mi
      mm = Ci + tcrossprod(mi)
      
      ##### variational parameter xi
      xi[ subject_rows, 1] = expectedXi(Theta_i, temp_mu_coef, mi, temp_psi_coef, Ci)
      xi_i = xi[subject_rows, ]
      
      Theta_i_quad = squareTheta(xi_i, Theta_i)
      
      #***************** estimate phi by maximum likelihood
      si = rbind(mi, 1)
      ss = cbind(rbind(mm, t(mi)), si)
      
      phi_a[[i]] = 2.0 * kronecker(Theta_i_quad, ss)  
      phi_b[,i] = t((Yi - 0.5) %*% kronecker(Theta_i, t(si)) )
      
      if(curr_iter == (maxiter - 1) ){
        scores[i,] = mi
      }
      
    } # end loop over I
    
    phi_a_sum = Reduce("+", phi_a)
    
    phi_vec = -solve(phi_a_sum) %*% rowSums(phi_b)
    phi_mat = matrix(phi_vec, nrow = Kt, ncol = npc + 1, byrow = TRUE)
    
    mu_coef = phi_mat[, npc+1]
    psi_coef = phi_mat[, 1:npc]
    
    if(npc == 1){ psi_coef = matrix(psi_coef, ncol = 1)}
    
    ## calculate error
    curr_iter = curr_iter + 1;
    error[curr_iter] = sum((psi_coef-temp_psi_coef)^2) + sum((mu_coef-temp_mu_coef)^2);

    temp_psi_coef = psi_coef
    temp_mu_coef = mu_coef
    
  } ## end while loop

  fits = rep(NA, dim(Y)[1])
  subject_coef = mu_coef + tcrossprod(psi_coef, scores)
 
  ## vectorize
   for(i in 1:I){
    subject_rows = rows$first_row[i]:rows$last_row[i]
    fits[subject_rows] = Theta[subject_rows, ] %*% subject_coef[,i]
   }
  
  fittedVals = data.frame(id = Y$id, index = Y$index, value = fits)
  
  ### what do we do if unevenly spaced grids?
  Theta2 = bs(seq(t.min, t.max, length.out = Di), knots = knots, intercept = TRUE) 
  
  ret = list(
    "knots" = knots, 
    "mu" = Theta2 %*% mu_coef,#
    "efunctions" = Theta2 %*% psi_coef, #
    "evalues" =  rep(1, npc),#
    "npc" = npc,#
    "scores" = scores,#
    "error" = error,
    "subject_coefs" = subject_coef,
    "Yhat" = fittedVals, #
    "Y" = Y, #
    "family" = "binomial"
  )

  class(ret) = "fpca" 
  return(ret)
} 
