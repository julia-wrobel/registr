#' Functional principal components analysis via variational EM
#' 
#' Function used in the FPCA step for registering functional data,
#' called by \code{\link{register_fpca}} when \code{family = "gaussian"}. 
#' Parameters estimated based on probabilistic PCA framework originally 
#' introduced by Tipping and Bishop in 1999. \cr \cr
#' The number of functional principal components (FPCs) can either be specified
#' directly (argument \code{npc}) or chosen based on the explained share of
#' variance (\code{npc_varExplained}). For the latter, the solution is based
#' on running the FPCA with \code{npc = 20} (and correspondingly \code{Kt = 20})
#' before reducing the solution to the relevant number of FPCs. Doing so,
#' we approximate the overall variance in the data \code{Y} with the variance
#' represented by the FPC basis with 20 FPCs.
#'
#' @param Y Dataframe. Should have variables id, value, index. 
#' @param npc,npc_varExplained The number of functional principal components (FPCs)
#' has to be specified either directly as \code{npc} or based on their explained
#' share of variance. In the latter case, \code{npc_criterion} has to be set
#' to a share between 0 and 1.
#' @param Kt Number of B-spline basis functions used to estimate mean functions
#' and functional principal components. Default is 8. If \code{npc_varExplained}
#' is used, \code{Kt} is set to 20.
#' @param maxiter Maximum number of iterations to perform for EM algorithm. Default is 50.
#' @param t_min Minimum value to be evaluated on the time domain. 
#' @param t_max Maximum value to be evaluated on the time domain.
#' @param print.iter Prints current error and iteration
#' @param row_obj If NULL, the function cleans the data and calculates row indices. 
#' Keep this NULL if you are using standalone \code{register} function.
#' @param seed Set seed for reproducibility. Defaults to 1988.
#' @param periodic If TRUE, uses periodic b-spline basis functions. Default is FALSE.
#' @param error_thresh Error threshold to end iterations. Defaults to 0.0001.
#' @param ... Additional arguments passed to or from other functions
#' @param verbose Can be set to integers between 0 and 4 to control the level of
#' detail of the printed diagnostic messages. Higher numbers lead to more detailed
#' messages. Defaults to 1.
#' @param subsample if the number of rows of the data is greater than 
#' 10 million rows, the `id` values are subsampled to get the mean coefficients.
#' 
#' @author Julia Wrobel \email{julia.wrobel@@cuanschutz.edu},
#' Jeff Goldsmith \email{ajg2202@@cumc.columbia.edu},
#' Alexander Bauer \email{alexander.bauer@@stat.uni-muenchen.de}
#' @importFrom splines bs
#' @importFrom pbs pbs
#' @importFrom stats quantile
#' 
#' @return An object of class \code{fpca} containing:
#' \item{fpca_type}{Information that FPCA was performed with the 'variationEM' approach,
#' in contrast to registr::gfpca_twoStep.}
#' \item{t_vec}{Time vector over which the mean \code{mu} was evaluated.}
#' \item{knots}{Cutpoints for B-spline basis used to rebuild \code{alpha}.}
#' \item{efunctions}{\eqn{D \times npc} matrix of estimated FPC basis functions.}
#' \item{evalues}{Estimated variance of the FPC scores.}
#' \item{evalues_sum}{Approximation of the overall variance in \code{Y}, based
#' on an initial run of the FPCA with \code{npc = 20}. Is \code{NULL} if
#' \code{npc_varExplained} was not specified.}
#' \item{npc}{number of FPCs.}
#' \item{scores}{\eqn{I \times npc} matrix of estimated FPC scores.}
#' \item{alpha}{Estimated population-level mean.}
#' \item{mu}{Estimated population-level mean. Same value as \code{alpha} but included for compatibility
#' with \code{refund.shiny} package.}
#' \item{subject_coefs}{B-spline basis coefficients used to construct subject-specific means. 
#' For use in \code{registr()} function.}
#' \item{Yhat}{FPC approximation of subject-specific means.}
#' \item{Y}{The observed data.}
#' \item{family}{\code{gaussian}, for compatibility with \code{refund.shiny} package.}
#' \item{sigma2}{Estimated error variance}
#' @export
#' 
#' @references Tipping, M. E. and Bishop, C (1999). Probabilistic Principal Component Analysis.
#' \emph{Journal of the Royal Statistical Society Series B,}, 592--598.
#' 
#' @examples
#' data(growth_incomplete)
#' 
#' fpca_obj = fpca_gauss(Y = growth_incomplete, npc = 2)
#' plot(fpca_obj)
#' 
fpca_gauss = function(Y, npc = NULL, npc_varExplained = NULL, Kt = 8, maxiter = 20,
                      t_min = NULL, t_max = NULL, 
											print.iter = FALSE, row_obj= NULL, seed = 1988, periodic = FALSE, 
											error_thresh = 0.0001, subsample = TRUE, verbose = 1, 
											...){
	
  if (is.null(npc) & is.null(npc_varExplained))
    stop("Please either specify 'npc' or 'npc_varExplained'.")
  if (!is.null(npc_varExplained) && ((npc_varExplained < 0) | (npc_varExplained > 1)))
    stop("'npc_varExplained' must be a number between 0 and 1.")
  if (!is.null(npc) & !is.null(npc_varExplained))
    message("Ignoring argument 'npc' since 'npc_varExplained' is specified.")
  
  if (!is.null(npc_varExplained)) {
    npc = 20
    Kt  = 20
  }
  
	## clean data
	if (is.null(row_obj)) {
		data = data_clean(Y)
		Y    = data$Y
		rows = data$Y_rows
		I    = data$I
	} else {
		rows = row_obj
		I    = dim(rows)[1]
	}
	
	if (Kt < 3) {
		stop("Kt must be greater than or equal to 3.")
	}
	
	time = Y$index
	
	## construct theta matrix
	if (is.null(t_min)) { t_min = min(time) }
	if (is.null(t_max)) { t_max = max(time) }
	
	if (periodic) {
		# if periodic, then we want more global knots, because the resulting object from pbs 
		# only has (knots+intercept) columns.
		knots = quantile(time, probs = seq(0, 1, length = Kt + 1))[-c(1, Kt + 1)]
		Theta_phi = pbs(c(t_min, t_max, time), knots = knots, intercept = TRUE)[-(1:2),]
	} else {
		# if not periodic, then we want fewer global knots, because the resulting object from bs
		# has (knots+degree+intercept) columns, and degree is set to 3 by default.
		knots = quantile(time, probs = seq(0, 1, length = Kt - 2))[-c(1, Kt - 2)]
		Theta_phi =  bs(c(t_min, t_max, time), knots = knots, intercept = TRUE)[-(1:2),]
	}
	
	## initialize all parameters
	set.seed(seed)
	nrows_basis = nrow(Theta_phi)
	if (subsample && nrows_basis > 1000000) {
	  if (verbose > 2) {
	    message("fpca_gauss: Running Sub-sampling")
	  }       
	  uids = unique(Y$id)
	  nids = length(uids)
	  avg_rows_per_id = nrows_basis / nids
	  size = round(1000000 / avg_rows_per_id)
	  if(nids > size){
	    ids = sample(uids, size = size, replace = FALSE)
	  } else {
	    ids = uids  
	  } 
	  subsampling_index = which(Y$id %in% ids)      
	  rm(uids,nids)
	  rm(ids)
	} else {
	  subsampling_index = 1:nrows_basis
	}
	if (verbose > 2) {
	  message("fpca_gauss: running GLM")
	}    
	if (requireNamespace("fastglm", quietly = TRUE)) {
	  glm_obj = fastglm::fastglm(y = Y$value[subsampling_index], x = Theta_phi[subsampling_index,], 
	                             family = "gaussian", method=2)
	} else {
	  glm_obj = glm(Y$value[subsampling_index] ~ 0 + Theta_phi[subsampling_index,], family = "gaussian",
	                control = list(trace = verbose > 3))
	}
	rm(subsampling_index)
	if (verbose > 2) {
	  message("fpca_gauss: GLM finished")
	}    
	alpha_coefs = coef(glm_obj)
	alpha_coefs = matrix(alpha_coefs, Kt, 1)
	
	arg_list <- list(npc              = npc,
	                 npc_varExplained = npc_varExplained,
	                 Kt               = Kt,
	                 maxiter          = maxiter,
	                 print.iter       = print.iter,
	                 seed             = seed,
	                 periodic         = periodic,
	                 error_thresh     = error_thresh,
	                 verbose          = verbose,
	                 Y                = Y,
	                 rows             = rows,
	                 I                = I,
	                 knots            = knots,
	                 Theta_phi        = Theta_phi,
	                 alpha_coefs      = alpha_coefs)
	
	# main optimization step
	fpca_resList = do.call(fpca_gauss_optimization, args = arg_list)
  
  ret = list(
    fpca_type     = "variationalEM",
  	t_vec         = fpca_resList$t_vec,
  	knots         = knots, 
  	alpha         = fpca_resList$Theta_phi_mean %*% alpha_coefs,
  	mu            = fpca_resList$Theta_phi_mean %*% alpha_coefs, # return this to be consistent with refund.shiny
  	efunctions    = fpca_resList$efunctions, 
  	evalues       = fpca_resList$evalues,
    evalues_sum   = fpca_resList$evalues_sum,
  	npc           = fpca_resList$npc,
  	scores        = fpca_resList$scores,
  	subject_coefs = fpca_resList$subject_coef,
  	Yhat          = fpca_resList$fittedVals, 
  	Y             = Y, 
  	family        = "gaussian",
  	sigma2        = fpca_resList$sigma2
  )
  
  class(ret) = "fpca" 
  return(ret)
  
} # end function



#' Internal main optimization for fpca_gauss
#' 
#' @inheritParams fpca_gauss
#' @param Y,rows,I,knots,Theta_phi,alpha_coefs Internal objects created in
#' \code{fpca_gauss}.
#' 
#' @importFrom splines bs
#' @importFrom pbs pbs
#' @importFrom stats rnorm
#' 
#' @return list with elements \code{t_vec}, \code{Theta_phi_mean},
#' \code{efunctions}, \code{evalues}, \code{evalues_sum}, \code{scores},
#' \code{subject_coef}, \code{fittedVals}, \code{sigma2}. See documentation of
#' \code{\link{fpca_gauss}} for details.
fpca_gauss_optimization <- function(npc, npc_varExplained, Kt, maxiter, print.iter,
                                    seed, periodic, error_thresh, verbose,
                                    Y, rows, I, knots, Theta_phi, alpha_coefs) {
  
  curr_iter = 1
  error     = rep(NA, maxiter)
  error[1]  = 100.0
  
  set.seed(seed)
  psi_coefs   = matrix(rnorm(Kt * npc), Kt, npc) * 0.5
  sigma2      = 1
  
  temp_alpha_coefs = alpha_coefs
  temp_psi_coefs   = psi_coefs
  temp_sigma2      = sigma2
  
  phi_a     = list(NA, I)
  phi_b     = matrix(0, nrow = Kt * (npc+1), ncol = I)
  scores    = matrix(NA, I, npc)
  sigma_vec = rep(NA, I)
  ##### update parameters
  while (curr_iter < maxiter && error[curr_iter] > error_thresh) {
    
    if (print.iter) {
      message("current iteration: ", curr_iter)
      message("current error: ", error[curr_iter])
    }
    
    for (i in 1:I) {
      subject_rows = rows$first_row[i]:rows$last_row[i]
      
      Yi           = Y$value[subject_rows]
      Theta_i      = Theta_phi[subject_rows, ] 
      Theta_i_quad = crossprod(Theta_i)
      
      # posterior scores
      Ci = solve(1/sigma2 * crossprod(psi_coefs, Theta_i_quad) %*% psi_coefs + diag(npc))
      
      mi_inner = 1/sigma2 * (crossprod(Yi, Theta_i)  - crossprod(alpha_coefs, Theta_i_quad)) %*% psi_coefs
      mi = tcrossprod(Ci, mi_inner)
      mm = Ci + tcrossprod(mi)
      
      # estimate sigma2 by maximum likelihood
      sigma_vec[i] =  -2 * t(mi) %*% t(psi_coefs) %*% t(Theta_i) %*% (Yi - Theta_i %*% alpha_coefs) +
        crossprod((Yi - Theta_i %*% alpha_coefs)) +
        sum(diag( crossprod(psi_coefs, Theta_i_quad) %*% psi_coefs %*% Ci)) + 
        t(mi) %*% crossprod(psi_coefs, Theta_i_quad) %*% psi_coefs %*% mi
      
      # estimate phi by maximum likelihood
      si = rbind(mi, 1)
      ss = cbind(rbind(mm, t(mi)), si)
      
      phi_a[[i]] = kronecker(Theta_i_quad, ss)  
      phi_b[,i]  =  t(Yi) %*% kronecker(Theta_i, t(si))
      
      scores[i,] = mi
      
    } # end loop over subjects
    
    sigma2    = 1/length(Y$value) * sum(sigma_vec)
    phi_a_sum = Reduce("+", phi_a)
    
    phi_vec = solve(phi_a_sum) %*% rowSums(phi_b)
    phi_mat = matrix(phi_vec, nrow = Kt, ncol = npc + 1, byrow = TRUE)
    
    alpha_coefs = phi_mat[, npc+1]
    psi_coefs   = phi_mat[, 1:npc, drop = TRUE]
    
    ## calculate error
    curr_iter = curr_iter + 1
    error[curr_iter] = sum((psi_coefs - temp_psi_coefs)^2) +
      sum((alpha_coefs - temp_alpha_coefs)^2) +
      (sigma2 - temp_sigma2)^2
    
    temp_psi_coefs   = psi_coefs
    temp_alpha_coefs = alpha_coefs
    temp_sigma2      = sigma2
    
  } # end while loop
  
  ## mean and eigenfunctions will have the same grid as the subject with most measurements
  id_mostObserved = names(sort(table(Y$id), decreasing = TRUE))[1]
  t_vec           = sort(Y$index[Y$id == id_mostObserved])
  if (periodic) {
    Theta_phi_mean = pbs(t_vec, knots = knots, intercept = TRUE)
  } else {
    Theta_phi_mean =  bs(t_vec, knots = knots, intercept = TRUE)
  }
  
  # orthogonalize eigenvectors and extract eigenvalues
  scores_unorthogonal = scores
  psi_svd    = svd(Theta_phi_mean %*% psi_coefs)
  efunctions = psi_svd$u
  evalues    = ( psi_svd$d ) ^ 2
  
  # choose npc adaptively, if 'npc_varExplained' is specified
  if (is.null(npc_varExplained)) { # no adaptive choice of npc
    evalues_sum = NULL
  } else { # adaptive choice of npc
    evalues_sum = sum(evalues)
    npc         = which(cumsum(evalues) / evalues_sum > npc_varExplained)[1]
    if (verbose > 0) {
      message(paste0("Using the first ",npc," FPCs which explain ",
                     round(sum(evalues[1:npc]) / evalues_sum * 100, 1),"% of the (approximated) total variance."))
    }
    # restrict some elements to the first npc dimensions
    efunctions = efunctions[, 1:npc, drop=FALSE]
    evalues    = evalues[1:npc]
  }
  
  # orthogonalize scores
  d_diag     = if (length(psi_svd$d) == 1) { matrix(psi_svd$d) } else { diag(psi_svd$d) }
  scores     = scores_unorthogonal %*% psi_svd$v %*% d_diag
  
  if (!is.null(npc_varExplained)) { # restrict further elements to npc dimensions
    scores_unorthogonal = scores_unorthogonal[, 1:npc, drop=FALSE]
    scores              = scores[, 1:npc, drop=FALSE]
    psi_coefs           = psi_coefs[, 1:npc, drop=FALSE]
  }
  
  # calculate the fitted values on Theta_phi for compatibility with registr()
  fits         = rep(NA, dim(Y)[1])
  subject_coef = alpha_coefs + tcrossprod(psi_coefs, scores_unorthogonal)
  for (i in 1:I) {
    subject_rows       = rows$first_row[i]:rows$last_row[i]
    fits[subject_rows] = Theta_phi[subject_rows, ] %*% subject_coef[,i]
  }
  
  fittedVals = data.frame(id = Y$id, index = Y$index, value = fits)
  
  
  return(list(t_vec          = t_vec,
              Theta_phi_mean = Theta_phi_mean,
              npc            = npc,
              efunctions     = efunctions,
              evalues        = evalues,
              evalues_sum    = evalues_sum,
              scores         = scores,
              subject_coef   = subject_coef,
              fittedVals     = fittedVals,
              sigma2         = sigma2))
}
