#' Covariance estimation after Hall et al. (2008)
#' 
#' Internal function for the estimation of the covariance matrix of the latent
#' process using the approach of Hall et al. (2008). Used in the
#' two-step GFPCA approach implemented in \code{\link{gfpca_twoStep}}. \cr \cr
#' This function is an adaptation of the implementation of Jan
#' Gertheiss and Ana-Maria Staicu for Gertheiss et al. (2017), with focus on
#' higher (RAM) efficiency for large data settings.
#' 
#' @param index_evalGrid Grid for the evaluation of the covariance structure.
#' @param Kt Number of P-spline basis functions for the estimation of the
#' marginal mean. Defaults to 8.
#' @param Kc Number of marginal P-spline basis functions for smoothing the
#' covariance surface. Defaults to 8.
#' @param diag_epsilon Small constant to which diagonal elements of the
#' covariance matrix are set if they are smaller. Defaults to 0.01.
#' @inheritParams gfpca_twoStep
#' 
#' @return Covariance matrix with dimension \code{time_evalGrid x time_evalGrid}.
#' 
#' @author Alexander Bauer \email{alexander.bauer@@stat.uni-muenchen.de},
#' based on work of Jan Gertheiss and Ana-Maria Staicu
#' 
#' @references Hall, P., Müller, H. G., & Yao, F. (2008). Modelling sparse
#' generalized longitudinal observations with latent Gaussian processes.
#' \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)},
#' 70(4), 703--723.
#' 
#' Gertheiss, J., Goldsmith, J., & Staicu, A. M. (2017). A note on
#' modeling sparse exponential-family functional response curves.
#' \emph{Computational statistics & data analysis}, 105, 46--52.
#' 
#' @import mgcv
#' 
#' @examples
#' data(growth_incomplete)
#' 
#' index_grid = c(1.25, seq(from = 2, to = 18, by = 1))
#' cov_matrix = registr:::cov_hall(growth_incomplete, index_evalGrid = index_grid)
#' 
cov_hall = function(Y, index_evalGrid, Kt = 8, Kc = 8, family = "gaussian",
                    diag_epsilon = 0.01){
  
  family_mgcv = family
  
  ids  = as.character(unique(Y$id))
  grid = unique(Y$index)
  D    = length(grid)
  I    = length(unique(Y$id))
  # empty matrix of measurements
  Y_miss = matrix(NA, nrow = I, ncol = D)
  
  # estimate the marginal mean mu(t) of the LGP X(t) (before applying the response function)
  out = mgcv::gam(value ~ s(index, bs = "ps", k = Kt),
                  family = family_mgcv, data = Y)
  mu  = as.vector(mgcv::predict.gam(out, newdata = data.frame(index = index_evalGrid)))
  
  # fill the Y_miss measurement matrix
  for (i in 1:I) {
    Yi     = Y$value[Y$id == ids[i]]
    ti     = Y$index[Y$id == ids[i]]
    indexi = unlist(sapply(ti, function(t) which(grid == t)))
    Y_miss[i,indexi] = Yi
  }
  
  # estimate the covariance surface beta(s,t) of g(X(t))
  combi_mat = expand.grid(index1 = 1:D, index2 = 1:D)
  Yi_2mom   = matrix(sapply(1:nrow(combi_mat), function(i) { 
    mean(Y_miss[,combi_mat$index1[i]] * Y_miss[,combi_mat$index2[i]], na.rm = TRUE)
  }), ncol = D)
  diag(Yi_2mom) = NA # omit the diagonal as of its discontinuuos structure
  
  row.vec = rep(grid, each = D) # set up row variable for bivariate smoothing
  col.vec = rep(grid, D)        # set up column variable for bivariate smoothing
  
  # smooth the covariance surface
  model_smoothed = mgcv::gam(as.vector(Yi_2mom) ~ te(row.vec, col.vec, k = Kc, bs = "ps"))
  Yi_2mom_sm     = matrix(
    mgcv::predict.gam(model_smoothed,
                      newdata = data.frame(row.vec = rep(index_evalGrid, each  = length(index_evalGrid)),
                                           col.vec = rep(index_evalGrid, times = length(index_evalGrid)))),
    ncol = length(index_evalGrid))
  Yi_2mom_sm = (Yi_2mom_sm + t(Yi_2mom_sm)) / 2 # ensure symmetry
  
  # estimate the marginal mean alpha(t) of g(X(t))
  Y_miss_mean = colMeans(Y_miss, na.rm = TRUE)
  Y.mean_sm = as.vector(mgcv::predict.gam(
    mgcv::gam(Y_miss_mean ~ s(grid, bs = "ps", k = Kt), method = "REML"),
    newdata = data.frame(grid = index_evalGrid)
  ))
  
  ## estimate tau(s,t) as final estimate for the covariance surface
  # calculate the numerator
  Yi.cov_sm = Yi_2mom_sm - (matrix(Y.mean_sm, ncol = 1) %*% matrix(Y.mean_sm, nrow = 1))
  
  # divide the numerator by the denominator, dependent on the first derivative
  # of the link function
  if (family == "gaussian") {
    Zi.cov_sm = diag(1 / mu) %*% Yi.cov_sm %*% diag(1 / mu)
  } else if (family == "binomial") {
    Zi.cov_sm = diag(1 / deriv.inv.logit(mu)) %*% Yi.cov_sm %*% diag(1 / deriv.inv.logit(mu))
  }
  
  # ensure a proper diagonal of the covariance surface
  ddd = diag(Zi.cov_sm)
  diag(Zi.cov_sm) = ifelse(ddd < diag_epsilon, diag_epsilon, ddd)
  
  return(Zi.cov_sm)
}
