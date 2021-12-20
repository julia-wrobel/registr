#' Covariance estimation after Hall et al. (2008)
#' 
#' Internal function for the estimation of the covariance matrix of the latent
#' process using the approach of Hall et al. (2008). Used in the
#' two-step GFPCA approach implemented in \code{\link{gfpca_twoStep}}. \cr \cr
#' This function is an adaptation of the implementation of Jan
#' Gertheiss and Ana-Maria Staicu for Gertheiss et al. (2017), with focus on
#' higher (RAM) efficiency for large data settings.
#' 
#' The implementation deviates from the algorithm described in Hall (2008) in
#' one crucial step -- we compute the crossproducts of \emph{centered}
#' observations and smooth the surface of these crossproducts directly instead
#' of computing and smoothing the surface of crossproducts of uncentered
#' observations and subsequently subtracting the (crossproducts of the) mean
#' function. The former seems to yield smoother eigenfunctions and 
#' fewer non-positive-definite covariance estimates.
#' 
#' If the data \code{Y} or the crossproduct matrix contain more than
#' \code{100,000} rows or elements, the estimation of the marginal mean or
#' the smoothing step of the covariance matrix are performed by
#' using the discretization-based estimation algorithm in \code{\link[mgcv]{bam}}
#' rather than the \code{\link[mgcv]{gam}} estimation algorithm.
#' 
#' @param index_evalGrid Grid for the evaluation of the covariance structure.
#' @param Kt Number of P-spline basis functions for the estimation of the
#' marginal mean. Defaults to 25.
#' @param Kc Number of marginal P-spline basis functions for smoothing the
#' covariance surface. Defaults to 10.
#' @param diag_epsilon Small constant to which diagonal elements of the
#' covariance matrix are set if they are smaller. Defaults to 0.01.
#' @param make_pd Indicator if positive (semi-)definiteness of the returned
#' latent covariance should be ensured via \code{Matrix::near_PD()}. Defaults to
#' TRUE.
#' @inheritParams gfpca_twoStep
#' 
#' @return Covariance matrix with dimension \code{time_evalGrid x time_evalGrid}.
#' 
#' @author Alexander Bauer \email{alexander.bauer@@stat.uni-muenchen.de} and 
#' Fabian Scheipl, based on work of Jan Gertheiss and Ana-Maria Staicu
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
#' @importFrom Matrix nearPD
#' @importFrom tidyr nest unnest
#' @importFrom dplyr select mutate filter n_distinct
#' @importFrom purrr map cross_df
#' @importFrom stats na.omit predict
#' @examples
#' data(growth_incomplete)
#' 
#' index_grid = c(1.25, seq(from = 2, to = 18, by = 1))
#' cov_matrix = registr:::cov_hall(growth_incomplete, index_evalGrid = index_grid)
#' 
cov_hall = function(Y, index_evalGrid, Kt = 25, Kc = 10, family = "gaussian",
  diag_epsilon = 0.01, make_pd = TRUE){
  
  if (family == "gamma") {
    # let the data start at 1.01 to make the marginal cov estimation more stable
    if (min(Y$value) < 1.01)
      Y$value = Y$value - min(Y$value) + 1.01
  }

  # define the model family
  if (family == "gamma") {
    family_mgcv = mgcv::Tweedie(p = 2, link = "log")
  } else {
    family_mgcv = family
  }
  
  this_gam = if (nrow(Y) > 1e5) {
    function(...) mgcv::bam(..., discrete = length(index_evalGrid))
  } else mgcv::gam
  # estimate the marginal mean mu(t) of the LGP X(t) (before applying the response function)
  model_mean = this_gam(value ~ s(index, bs = "ps", k = Kt), family = family_mgcv, data = Y)
  
  # create crossproducts of centered data
  Y$centered = Y$value - predict(model_mean, type = "response")
  # not at all sure about this cutoff here...
  Y_crossprods = if (dplyr::n_distinct(Y$index) > .5 * nrow(Y)) {
    crossprods_irregular(Y)
  } else {
    crossprods_regular(Y)
  }

  # smooth the covariance surface
  that_gam = if (nrow(Y_crossprods) > 1e5) {
    function(...) mgcv::bam(..., discrete = length(index_evalGrid))
  } else mgcv::gam
  model_smoothed = that_gam(cross ~ te(i1, i2, k = Kc, bs = "ps"), 
    data = Y_crossprods)
  Yi_2mom_sm     = matrix(
    predict(model_smoothed,
      newdata = data.frame(i1 = rep(index_evalGrid, each  = length(index_evalGrid)),
        i2 = rep(index_evalGrid, times = length(index_evalGrid)))),
    ncol = length(index_evalGrid))
  Yi_2mom_sm = (Yi_2mom_sm + t(Yi_2mom_sm)) / 2 # ensure symmetry
  
  # divide the numerator by the denominator, dependent on the first derivative
  # of the link function
  if (family == "gaussian") { # identity link
    Zi.cov_sm = Yi_2mom_sm # first derivative of identity function = 1
  } else {
    mu = as.vector(
      predict.gam(model_mean, type = "link", 
        newdata = data.frame(index = index_evalGrid)))
  }  
  if (family == "binomial") { # logit link
    Zi.cov_sm = diag(1 / deriv.inv.logit(mu)) %*% Yi_2mom_sm %*% 
      diag(1 / deriv.inv.logit(mu))
  } 
  if (family %in% c("gamma","poisson")) { # log link
    stopifnot(model_mean$family$link == "log")
    Zi.cov_sm = diag(mu) %*% Yi_2mom_sm %*% diag(mu)
  }
  
  # ensure positive diagonal of the covariance surface
  ddd = diag(Zi.cov_sm)
  diag(Zi.cov_sm) = ifelse(ddd < diag_epsilon, diag_epsilon, ddd)
  
  if (make_pd) {
    Zi.cov_sm = as.matrix(Matrix::nearPD(Zi.cov_sm, do2eigen = TRUE)$mat)
  }  
  Zi.cov_sm
}


#' Crossproduct computation for highly irregular grids
#' 
#' Compute the crossproduct in a fast way for highly irregular grids
#' (index values are mostly unique).
#' Only used internally in \code{cov_hall()}.
#' 
#' @param Y Dataframe with the centered observations.
#' Should have values id, centered, index.
crossprods_irregular = function(Y) {
  
  # some NULL variable definitions to appease CRAN package checks regarding the use of ggplot2
  crossprods = NULL
  
  
  Y_crossprods = select(Y, .data$id, .data$index, .data$centered) %>% 
    nest(data = c(.data$index, .data$centered)) %>% 
    mutate(crossprods = map(.data$data, ~ {
      cross_df(list(i1 = .x$index, i2 = .x$index)) %>% 
        mutate(cross = 
            cross_df(list(v1 = .x$centered, v2 = .x$centered)) %>% 
            {.[, 1] * .[, 2]} %>% unlist
        )
    })) %>% select(.data$id, .data$crossprods) %>% unnest(cols = crossprods) %>% 
    # drop diagonal
    filter(.data$i1 != .data$i2)
}


#' Crossproduct computation for mostly regular grids
#' 
#' Compute the crossproduct in a fast way for mostly regular grids
#' (index values are mostly *not* unique).
#' Only used internally in \code{cov_hall()}.
#' 
#' @inheritParams crossprods_irregular
crossprods_regular = function(Y) {
  ids  = as.character(unique(Y$id))
  grid = sort(unique(Y$index))
  D    = length(grid)
  I    = length(ids)
  
  # empty matrix of measurements
  Y_miss = matrix(NA, nrow = I, ncol = D)
  # fill the Y_miss measurement matrix
  for (i in 1:I) {
    ti     = Y$index[Y$id == ids[i]]
    indexi = unlist(sapply(ti, function(t) which(grid == t)))
    Y_miss[i, indexi] = Y$centered[Y$id == ids[i]]
  }
  
  # estimate the covariance surface beta(s,t) of g(X(t))
  combi_mat = expand.grid(index1 = 1:D, index2 = 1:D)
  Yi_2mom   = matrix(sapply(1:nrow(combi_mat), function(i) { 
    mean(Y_miss[,combi_mat$index1[i]] * Y_miss[,combi_mat$index2[i]], na.rm = TRUE)
  }), ncol = D)
  diag(Yi_2mom) = NA # omit the diagonal
  na.omit(data.frame(
    i1 = rep(grid, each = D), 
    i2 = rep(grid, D), times = D,
    cross = as.vector(Yi_2mom)))
}