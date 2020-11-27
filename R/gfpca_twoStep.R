#' Generalized functional principal component analysis
#' 
#' Function for applying FPCA to different exponential family distributions.
#' Used in the FPCA step for registering functional data,
#' called by \code{\link{register_fpca}} when \code{GFPCA_type = "two-step"}. \cr \cr
#' The method implements the `two-step approach` of Gertheiss et al. (2017)
#' and is based on the approach of Hall et al. (2008) to estimate functional
#' principal components. \cr \cr
#' This function is an adaptation of the implementation of Jan
#' Gertheiss for Gertheiss et al. (2017), with focus on higher (RAM) efficiency
#' for large data settings.
#' 
#' @param family One of \code{c("gaussian","binomial","gamma")}. Defaults to
#' \code{"gaussian"}.
#' @param index_relevantDigits Positive integer \code{>= 2}, stating the number
#' of relevant digits to which the index grid should be rounded. Coarsening the
#' index grid is necessary since otherwise the covariance surface matrix
#' explodes in size in the presence of too many unique index values (which is
#' always the case after some registration step). Defaults to 4. Set to
#' \code{NULL} to prevent rounding.
#' @param estimation_accuracy One of \code{c("high","low")}. When set to \code{"low"},
#' the mixed model estimation step in \code{lme4} is performed with lower
#' accuracy, reducing computation time. Defaults to \code{"high"}.
#' @param start_params Optional start values for gamm4.
#' @param periodic Only contained for full consistency with \code{fpca_gauss}
#' and \code{bfpca}. If TRUE, returns the knots vector for periodic b-spline
#' basis functions. Defaults to FALSE. This parameter does not change the
#' results of the two-step GFPCA.
#' @param ... Additional arguments passed to \code{\link{cov_hall}}.
#' @inheritParams register_fpca
#' @inheritParams fpca_gauss
#' 
#' @return An object of class \code{fpca} containing:
#' \item{t_vec}{Time vector over which the mean \code{mu} was evaluated.
#' The resolution is can be specified by setting \code{index_relevantDigits}.}
#' \item{knots}{Cutpoints for B-spline basis used to rebuild \code{alpha}.}
#' \item{efunctions}{\eqn{D \times npc} matrix of estimated FPC basis functions.}
#' \item{evalues}{Estimated variance of the FPC scores.}
#' \item{evalues_sum}{Sum of all eigenvalues, before restricting the \code{evalues}
#' vector to the first \code{npc} elements.}
#' \item{npc}{number of FPCs.}
#' \item{scores}{\eqn{I \times npc} matrix of estimated FPC scores.}
#' \item{alpha}{Estimated population-level mean.}
#' \item{mu}{Estimated population-level mean. Same value as \code{alpha} but included for compatibility
#' with \code{refund.shiny} package.}
#' \item{subject_coefs}{Always \code{NA} but included for full consistency
#' with \code{fpca_gauss} and \code{bfpca}.} 
#' \item{Yhat}{FPC approximation of subject-specific means.}
#' \item{Y}{The observed data.}
#' \item{family}{\code{binomial}, for compatibility with \code{refund.shiny} package.}
#' \item{gamm4_theta}{Estimated parameters of the mixed model.}
#' @export
#' 
#' @references Gertheiss, J., Goldsmith, J., & Staicu, A. M. (2017). A note on
#' modeling sparse exponential-family functional response curves.
#' \emph{Computational statistics & data analysis}, 105, 46--52.
#' 
#' Hall, P., Müller, H. G., & Yao, F. (2008). Modelling sparse
#' generalized longitudinal observations with latent Gaussian processes.
#' \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)},
#' 70(4), 703--723.
#' 
#' @author Alexander Bauer \email{alexander.bauer@@stat.uni-muenchen.de},
#' based on work of Jan Gertheiss
#' 
#' @importFrom gamm4 gamm4
#' @importFrom stats formula
#' @import lme4 mgcv
#' 
#' @examples
#' data(growth_incomplete)
#' 
#' fpca_obj = gfpca_twoStep(Y = growth_incomplete, npc = 2, family = "gaussian")
#' plot_fpca(fpca_obj)
#' 
gfpca_twoStep = function (Y, family = "gaussian", npc = 1, Kt = 8,
                          t_min = NULL, t_max = NULL,
                          row_obj = NULL, index_relevantDigits = 4L,
                          estimation_accuracy = "high", start_params = NULL,
                          periodic = FALSE,
                          ...) {
  
  if (family == "gamma" & any(Y$value <= 0))
    stop("family = 'gamma' can only be applied to strictly positive data.")
  
  # clean data
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
  if (!is.null(index_relevantDigits) && index_relevantDigits < 2) {
    stop("'relevant_digits' must be a positive integer >= 2.")
  }
  
  # coarsen the time index to better handle bigger data
  time_orig = Y$index
  if (!is.null(index_relevantDigits)) {
    Y$index   = coarsen_index(Y$index, index_relevantDigits)
  }
  
  time = Y$index
  if (is.null(t_min)) { t_min = min(time) }
  if (is.null(t_max)) { t_max = max(time) }
  output_index = sort(unique(time))
  
  # define knots vector. Not used here, but returned for full consistency with fpca_gauss
  if (periodic) {
    # if periodic, then we want more global knots, because the resulting object from pbs 
    # only has (knots+intercept) columns.
    knots = quantile(time_orig, probs = seq(0, 1, length = Kt + 1))[-c(1, Kt + 1)]
  } else {
    # if not periodic, then we want fewer global knots, because the resulting object from bs
    # has (knots+degree+intercept) columns, and degree is set to 3 by default.
    knots = quantile(time_orig, probs = seq(0, 1, length = Kt - 2))[-c(1, Kt - 2)]
  }
  
  Y.vec       = Y['value'][[1]]
  input_index = Y['index'][[1]]
  id.vec      = Y['id'][[1]]
  ids         = as.character(unique(id.vec))
  
  D     = length(output_index)
  I     = length(unique(id.vec))
  Y.obs = matrix(NA, nrow = I, ncol = D)
  
  for (i in 1:I) {
    Yi     = Y.vec[id.vec == ids[i]]
    ti     = input_index[id.vec == ids[i]]
    indexi = sapply(ti, function(t) which(output_index == t))
    Y.obs[i,indexi] = Yi
  }
  
  dat = data.frame(
    value = as.vector(t(Y.obs)),
    id    = rep(1:I, rep(D,I)),
    index = rep(output_index, I)
  )
  
  ### obtain eigenfunctions and eigenvalues
  # estimate the covariance matrix of the latent Gaussian process
  hmy_cov = cov_hall(Y, index_evalGrid = output_index, Kt = Kt,
                     family = family, diag_epsilon = 0.01, ...)
  
  # spectral decomposition
  eigen_HMY  = eigen(hmy_cov)
  fit.lambda = eigen_HMY$values
  fit.phi    = eigen_HMY$vectors
  
  # remove negative eigenvalues
  if (any(fit.lambda[1:npc] < 0)) {
    warning("The first 'npc' eigenvalues contained negative eigenvalues. These were removed.")
    wp         = which(fit.lambda > 0)
    fit.lambda = fit.lambda[wp]
    fit.phi    = fit.phi[,wp]
  }
  
  efunctions  = fit.phi[,1:npc, drop = FALSE]
  evalues     = fit.lambda[1:npc]
  evalues_sum = sum(fit.lambda[fit.lambda > 0])
  
  # prepare data for mixed model estimation
  for (i in 1:npc) {
    dat = cbind(dat, rep(efunctions[,i], I))
  }
  names(dat)[4:(4 + npc - 1)] = c(paste0("psi", 1:npc))
  dat_fit = dat[!is.na(dat$value),]
  
  # define the model family
  if (family == "gamma") {
    family_mgcv = mgcv::Tweedie(p = 2, link = "log")
  } else {
    family_mgcv = family
  }
  
  # mixed model
  random.structure = paste(paste0("psi", 1:npc), collapse = "+")
  random.formula   = stats::formula(paste("~(0+", random.structure, "|| id)"))
  if (estimation_accuracy == "high") {
    model = gamm4::gamm4(value ~ s(index, bs = "ps", k = Kt),
                         family = family_mgcv,
                         data   = dat_fit,
                         random = random.formula,
                         start  = start_params)
    
  } else if (estimation_accuracy == "low") {
    # Note: The following control arguments are commented out since they
    #       require a coming update of the gamm4 package.
    #       When the gamm4 update goes public:
    #       1) Uncomment the arguments
    #       2) Document the exact arguments for the optimization routine in more detail
    model = gamm4::gamm4(value ~ s(index, bs = "ps", k = Kt),
                         family  = family_mgcv,
                         data    = dat_fit,
                         random  = random.formula,
                         start   = start_params)
                         # control = lme4::glmerControl(optimizer     = "nloptwrap",
                         #                              calc.derivs   = FALSE,
                         #                              nAGQ0initStep = FALSE,
                         #                              boundary.tol  = 0,
                         #                              tolPwrss      = 1e-3,
                         #                              optCtrl = list(algorithm = "NLOPT_LN_BOBYQA", # default in lme4
                         #                                             xtol_abs  = 1e-04,             # default: 1e-08
                         #                                             ftol_abs  = 1e-04,             # default: 1e-08
                         #                                             maxeval   = 1e+02)),           # default: 1e+05
                         # nAGQ = 0)
  }
  
  gamm4_theta = attributes(model$mer)$theta
  
  scores       = as.matrix(coef(model$mer)$id[,npc:1])
  Z.gamm.fpca  = scores %*% t(efunctions[,1:npc])
  alpha        = as.vector(predict.gam(model$gam, newdata = data.frame(index = output_index)))
  alpha_matrix = matrix(rep(alpha, I), nrow = I, byrow = TRUE)
  Yhat         = alpha_matrix + Z.gamm.fpca
  
  ## format output
  fittedVals = data.frame(
    id    = rep(1:I, each = D),
    index = rep(output_index, I),
    value = as.vector(t(Yhat))
  )
  
  ret = list(t_vec         = output_index,
             knots         = knots,
             alpha         = matrix(alpha, ncol = 1), # return matrix for consistency with fpca_gauss()
             mu            = matrix(alpha, ncol = 1),
             efunctions    = efunctions,
             evalues       = evalues,
             evalues_sum   = evalues_sum,
             npc           = npc,
             scores        = scores,
             subject_coefs = NA,
             Yhat          = fittedVals,
             Y             = Y,
             family        = family,
             gamm4_theta   = gamm4_theta) # current estimates of gamm4
  class(ret) = "fpca"
  
  return(ret)
}