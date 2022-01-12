#' Binary functional principal components analysis
#' 
#' Function used in the FPCA step for registering binary functional data,
#' called by \code{\link{register_fpca}} when \code{family = "binomial"}. 
#' This method uses a variational EM algorithm to estimate scores and principal components for 
#' binary functional data. \cr \cr
#' The number of functional principal components (FPCs) can either be specified
#' directly (argument \code{npc}) or chosen based on the explained share of
#' variance (\code{npc_varExplained}). In the latter case, the explained share of
#' variance and accordingly the number of FPCs is estimated before the main
#' estimation step by once running the FPCA with \code{npc = 20} (and
#' correspondingly \code{Kt = 20}). Doing so, we approximate the overall
#' variance in the data \code{Y} with the variance represented by the FPC basis
#' with 20 FPCs.
#'
#' @inheritParams fpca_gauss
#' 
#' @author Julia Wrobel \email{julia.wrobel@@cuanschutz.edu},
#' Jeff Goldsmith \email{ajg2202@@cumc.columbia.edu},
#' Alexander Bauer \email{alexander.bauer@@stat.uni-muenchen.de}
#' @importFrom stats quantile
#' 
#' @return An object of class \code{fpca} containing:
#' \item{fpca_type}{Information that FPCA was performed with the 'variationEM' approach,
#' in contrast to registr::gfpca_twoStep.}
#' \item{t_vec}{Time vector over which the mean \code{mu} and the functional principal
#' components \code{efunctions} were evaluated.}
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
#' \item{Yhat}{FPC approximation of subject-specific means, before applying the
#' response function.}
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
#' Y = simulate_functional_data()$Y
#' 
#' # estimate 2 FPCs
#' bfpca_obj = bfpca(Y, npc = 2, print.iter = TRUE, maxiter = 25)
#' 
#' 
#' \donttest{
#' plot(bfpca_obj)
#' 
#' # estimate npc adaptively, to explain 90% of the overall variation
#' bfpca_obj2 = bfpca(Y, npc_varExplained = 0.9, print.iter = TRUE, maxiter = 30)
#' plot(bfpca_obj2)
#'}
bfpca = function(Y, npc = NULL, npc_varExplained = NULL, Kt = 8, maxiter = 50,
                 t_min = NULL, t_max = NULL, print.iter = FALSE, row_obj= NULL,
                 seed = 1988, periodic = FALSE, error_thresh = 0.0001,
                 verbose = 1, subsample=TRUE,
                 ...){
  
  if (is.null(npc) & is.null(npc_varExplained))
    stop("Please either specify 'npc' or 'npc_varExplained'.")
  if (!is.null(npc_varExplained) && ((npc_varExplained < 0) | (npc_varExplained > 1)))
    stop("'npc_varExplained' must be a number between 0 and 1.")
  if (!is.null(npc) & !is.null(npc_varExplained))
    message("Ignoring argument 'npc' since 'npc_varExplained' is specified.")
  
  if (!is.null(npc_varExplained)) {
    Kt = 20
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
  
  ## check that data is binary
  if (any( !(Y$value %in% c(0, 1)))) {
  	stop("'binomial' family requires data with binary values of 0 or 1")
  }
  
  if (is.null(t_min)) { t_min = min(time) }
  if (is.null(t_max)) { t_max = max(time) }
  
  # prepare the arguments for the main optimization step
  prep_arg_list = list(Y         = Y,
                       Kt        = Kt,
                       time      = time,
                       t_min     = t_min,
                       t_max     = t_max,
                       periodic  = periodic,
                       seed      = seed,
                       subsample = subsample,
                       verbose   = verbose)
  prepared_args = do.call(bfpca_argPreparation, args = prep_arg_list)
  
  # arguments for the main optimization step
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
                   knots            = prepared_args$knots,
                   Theta_phi        = prepared_args$Theta_phi,
                   xi               = prepared_args$xi,
                   alpha_coefs      = prepared_args$alpha_coefs)
  
  # choose the number of FPCs based on the explained share of variance
  if (!is.null(npc_varExplained)) {
    prep_arg_list_vexp    = prep_arg_list
    prep_arg_list_vexp$Kt = 20
    prepared_args_vexp    = do.call(bfpca_argPreparation, args = prep_arg_list_vexp)
    
    arg_list_vexp     = arg_list
    arg_list_vexp$npc = 20
    arg_list_vexp$Kt  = 20
    arg_list_vexp$knots       = prepared_args_vexp$knots
    arg_list_vexp$Theta_phi   = prepared_args_vexp$Theta_phi
    arg_list_vexp$alpha_coefs = prepared_args_vexp$alpha_coefs
    
    list_npcChoice    = do.call(bfpca_optimization, args = arg_list_vexp)
    arg_list$npc      = list_npcChoice$npc
    arg_list$npc_varExplained = NULL
  }
  
  # main optimization step
  fpca_resList = do.call(bfpca_optimization, args = arg_list)
  
  if (!is.null(npc_varExplained)) { # use (approximated) shares of variance
    fpca_resList$evalues     = list_npcChoice$evalues
    fpca_resList$evalues_sum = list_npcChoice$evalues_sum
  }
  
  ret = list(
    fpca_type     = "variationalEM",
    t_vec         = fpca_resList$t_vec,
    knots         = prepared_args$knots, 
    alpha         = fpca_resList$Theta_phi_mean %*% fpca_resList$alpha_coefs,
    mu            = fpca_resList$Theta_phi_mean %*% fpca_resList$alpha_coefs, # return this to be consistent with refund.shiny
    efunctions    = fpca_resList$efunctions, 
    evalues       = fpca_resList$evalues,
    evalues_sum   = fpca_resList$evalues_sum,
    npc           = fpca_resList$npc,
    scores        = fpca_resList$scores,
    subject_coefs = fpca_resList$subject_coef,
    Yhat          = fpca_resList$fittedVals, 
    Y             = Y, 
    family        = "binomial",
    error         = fpca_resList$error
  )

  class(ret) = "fpca" 
  return(ret)
} 


#' Internal main preparation function for bfpca
#' 
#' @inheritParams bfpca
#' @param Y,time,t_min,t_max Internal objects created in \code{bfpca}.
#' 
#' @importFrom splines bs
#' @importFrom pbs pbs
#' 
#' @return List with elements \code{knots}, \code{Theta_phi}, \code{xi},
#' \code{alpha_coefs}.
bfpca_argPreparation <- function(Y, Kt, time, t_min, t_max, periodic, seed,
                                 subsample, verbose) {
  
  if (periodic) {
    # if periodic, then we want more global knots, because the resulting object from pbs 
    # only has (knots+intercept) columns.
    knots     = quantile(time, probs = seq(0, 1, length = Kt + 1))[-c(1, Kt + 1)]
    Theta_phi = pbs(c(t_min, t_max, time), knots = knots, intercept = TRUE)[-(1:2),]
  } else {
    # if not periodic, then we want fewer global knots, because the resulting object from bs
    # has (knots+degree+intercept) columns, and degree is set to 3 by default.
    knots     = quantile(time, probs = seq(0, 1, length = Kt - 2))[-c(1, Kt - 2)]
    Theta_phi = bs(c(t_min, t_max, time), knots = knots, intercept = TRUE)[-(1:2),]
  }
  
  ## initialize all your vectors
  set.seed(seed)
  xi          = matrix(rnorm(dim(Y)[1]), ncol = 1) * 0.5
  
  nrows_basis = nrow(Theta_phi)
  if (subsample && nrows_basis > 1000000) {
    if (verbose > 2) {
      message("bfpca: Running Sub-sampling")
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
    message("bfpca: running GLM")
  }    
  if (requireNamespace("fastglm", quietly = TRUE)) {
    glm_obj = fastglm::fastglm(y = Y$value[subsampling_index], x = Theta_phi[subsampling_index,], family = "binomial", method=2)
  } else {
    glm_obj = glm(Y$value[subsampling_index] ~ 0 + Theta_phi[subsampling_index,], family = "binomial",
                  control = list(trace = verbose > 3))
  }
  rm(subsampling_index)
  if (verbose > 2) {
    message("bfpca: GLM finished")
  }    
  alpha_coefs = coef(glm_obj)
  alpha_coefs = matrix(alpha_coefs, Kt, 1)
  
  
  return(list(knots       = knots,
              Theta_phi   = Theta_phi,
              xi          = xi,
              alpha_coefs = alpha_coefs))
}



#' Internal main optimization for bfpca
#' 
#' Main optimization function for \code{bfpca}. If \code{npc_varExplained}
#' is specified, the function simply returns a list with elements \code{npc}
#' (chosen number of FPCs), \code{evalues} (estimated variances of the first 'npc'
#' FPCs) and \code{evalues_sum} (sum of the estimated variances of the first 20
#' FPCs, as approximation of the overall variance).
#' 
#' @inheritParams bfpca
#' @param Y,rows,I,knots,Theta_phi,xi,alpha_coefs Internal objects created in
#' \code{bfpca}.
#' 
#' @importFrom splines bs
#' @importFrom pbs pbs
#' @importFrom stats rnorm
#' 
#' @return list with elements \code{t_vec}, \code{Theta_phi_mean}, \code{alpha_coefs},
#' \code{efunctions}, \code{evalues}, \code{evalues_sum}, \code{scores},
#' \code{subject_coef}, \code{fittedVals}, \code{error}. See documentation of
#' \code{\link{fpca_gauss}} for details.
bfpca_optimization <- function(npc, npc_varExplained = NULL, Kt, maxiter, print.iter,
                               seed, periodic, error_thresh, verbose,
                               Y, rows, I, knots, Theta_phi, xi, alpha_coefs) {
  
  curr_iter = 1
  error     = rep(NA, maxiter)
  error[1]  = 100.0
  
  set.seed(seed)
  psi_coefs = matrix(rnorm(Kt * npc), Kt, npc) * 0.5
  
  temp_alpha_coefs = alpha_coefs
  temp_psi_coefs   = psi_coefs
  
  phi_a  = list(NA, I)
  phi_b  = matrix(0, nrow = Kt * (npc+1), ncol = I)
  scores = matrix(NA, I, npc)
  
  while (curr_iter < maxiter && error[curr_iter] > error_thresh) {
    
    if (print.iter) {
      message("current iteration: ", curr_iter)
      message("current error: ", error[curr_iter])
    }
    
    for (i in 1:I) {
      
      subject_rows = rows$first_row[i]:rows$last_row[i]
      
      Yi           = Y$value[subject_rows]
      Theta_i      = Theta_phi[subject_rows, ] 
      xi_i         = xi[subject_rows,]
      Theta_i_quad = squareTheta(xi_i, Theta_i)
      
      # posterior scores
      mlist = expectedScores(Yi, temp_alpha_coefs, temp_psi_coefs, Theta_i, Theta_i_quad)
      
      Ci = mlist$Ci
      mi = mlist$mi
      mm = Ci + tcrossprod(mi)
      
      # variational parameter xi
      xi[ subject_rows, 1] = expectedXi(Theta_i, temp_alpha_coefs, mi, temp_psi_coefs, Ci)
      xi_i = xi[subject_rows, ]
      
      Theta_i_quad = squareTheta(xi_i, Theta_i)
      
      # estimate phi by maximum likelihood
      si = rbind(mi, 1)
      ss = cbind(rbind(mm, t(mi)), si)
      
      phi_a[[i]] = 2.0 * kronecker(Theta_i_quad, ss)  
      phi_b[,i]  = t((Yi - 0.5) %*% kronecker(Theta_i, t(si)) )
      
      scores[i,] = mi
      
    } # end loop over I
    
    phi_a_sum = Reduce("+", phi_a)
    
    phi_vec = -solve(phi_a_sum) %*% rowSums(phi_b)
    phi_mat = matrix(phi_vec, nrow = Kt, ncol = npc + 1, byrow = TRUE)
    
    alpha_coefs = phi_mat[, npc+1]
    psi_coefs   = phi_mat[, 1:npc, drop = FALSE]
    
    ## calculate error
    curr_iter        = curr_iter + 1
    error[curr_iter] = sum((psi_coefs-temp_psi_coefs)^2) + sum((alpha_coefs-temp_alpha_coefs)^2)
    
    temp_psi_coefs   = psi_coefs
    temp_alpha_coefs = alpha_coefs
    
  } ## end while loop
  
  if (curr_iter < maxiter) {
    if (verbose > 2) {
      message("BFPCA converged.")
    }
  } else {
    warning("BFPCA convergence not reached. Try increasing maxiter.")
  }
  
  ## evaluate mean and eigenfunctions on a grid of length 100
  t_vec = seq(min(Y$index), max(Y$index), length.out = 100)
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
    return(list(npc         = npc,
                evalues     = evalues[1:npc],
                evalues_sum = evalues_sum))
  }
  
  # orthogonalize scores
  d_diag     = if (length(psi_svd$d) == 1) { matrix(psi_svd$d) } else { diag(psi_svd$d) }
  scores     = scores_unorthogonal %*% psi_svd$v %*% d_diag
  
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
              alpha_coefs    = alpha_coefs,
              npc            = npc,
              efunctions     = efunctions,
              evalues        = evalues,
              evalues_sum    = evalues_sum,
              scores         = scores,
              subject_coef   = subject_coef,
              fittedVals     = fittedVals,
              error          = error[!is.na(error)]))
  
}
