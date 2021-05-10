registr_helper = function(obj = NULL, Y = NULL, Kt = 8, Kh = 4, family = "gaussian", gradient = TRUE,
                          incompleteness = NULL, lambda_inc = NULL,
                          Y_template = NULL,
                          beta = NULL, t_min = NULL, t_max = NULL, row_obj = NULL,
                          periodic = FALSE, warping = "nonparametric",
                          gamma_scales = NULL, cores = 1L,  subsample = TRUE,
                          verbose = 1,
                          ...){ 
  if (!is.null(incompleteness)) {
    if (warping != "nonparametric") {
      stop("The functionality for incomplete curves is only available for 'warping = 'nonparametric''")
    }
    if (!(incompleteness %in% c("leading","trailing","full"))) {
      stop("'incompleteness' must be either 'leading', 'trailing' or 'full'.")
    }
    if ((is.null(lambda_inc) || (lambda_inc < 0))) {
      stop("For incomplete curves the penalization parameter 'lambda_inc' has to be set to some nonnegative value.")
    }
  }
  
  if (is.null(Y)) { 
    Y = obj$Y
  }
  
  if (family == "gamma" & any(Y$value <= 0)) {
    stop("family = 'gamma' can only be applied to strictly positive data.")
  } else if (family == "poisson" & any(Y$value < 0)) {
    stop("family = 'poisson' can only be applied to nonnegative data.")
  }
  
  if(is.null(obj)) { 
    if(warping == "nonparametric"){
      Y$tstar = Y$index
      
    } else if(warping == "piecewise_linear2"){
      # scale time to 0, 1 for parametric warping
      if(is.null(Y$index_scaled)){
        stop("For warping = piecewise_linear2, need an index_scaled variable that ranges from 0 to 1.")
      }
      Y$tstar = Y$index_scaled
    }
  }
  
  if (is.null(row_obj)) {
    data = data_clean(Y)
    Y    = data$Y
    rows = data$Y_rows
    I    = data$I
  } else{
    rows = row_obj
    I    = nrow(rows)
  }
  
  if (Kh < 4) {
    stop("Kh must be greater than or equal to 4.")
  }
  if (Kt < 4) {
    stop("Kt must be greater than or equal to 4.")
  }
  
  if(!(warping %in% c("nonparametric", "piecewise_linear2"))){
    stop("warping argument can only take values of nonparametric or piecewise_linear2")
  }
  
  if (gradient & !(family %in% c("gaussian","binomial"))) {
    warning("gradient = TRUE is only available for families 'gaussian' and 'binomial'. Setting gradient = FALSE.")
    gradient = FALSE
  } else if (gradient & periodic){
    warning("gradient = TRUE is only available for periodic = FALSE. Setting gradient = FALSE.")
    gradient = FALSE
  } else if (gradient & warping != "nonparametric"){
    warning("gradient = TRUE is only available for warping = nonparametric. Setting gradient = FALSE.")
    gradient = FALSE
  }
  
  if (is.null(t_min)) { t_min = min(Y$tstar) }
  if (is.null(t_max)) { t_max = max(Y$tstar) }
  stopifnot(!is.na(t_min),
            !is.na(t_max))
  
  if (!is.null(obj)) { # template function = GFPCA representation
    global_knots = obj$knots
    mean_coefs   = obj$subject_coefs
    
  } else { # template function = mean of all curves
    
    if (!is.null(Y_template)) {
      if (verbose > 2) {
        message("Registr: Extracting Y_template")
      }
      if (!all(c("id", "index", "value") %in% names(Y_template))) {
        stop("Y_template must have variables 'id', 'index', and 'value'.")
      } else if (!identical(range(Y_template$index), range(Y$index))) {
        stop("The range of 'index' must be equal for Y_template and Y.")
      }
      Y_template$tstar = Y_template$index
      mean_dat         = Y_template
    } else {
      mean_dat = Y
    }
    
    if (verbose > 2) {
      message("Registr: Getting nnots and basis functions")
    }
    if (periodic) {
      # if periodic, then we want more global knots, because the resulting object from pbs 
      # only has (knots+intercept) columns.
      global_knots = quantile(mean_dat$tstar, probs = seq(0, 1, length = Kt+1))[-c(1, Kt+1)]
      mean_basis   = pbs(c(t_min, t_max, mean_dat$tstar), knots = global_knots, intercept = TRUE)[-(1:2),]
      
    } else {
      # if not periodic, then we want fewer global knots, because the resulting object from bs
      # has (knots+degree+intercept) columns, and degree is set to 3 by default.
      global_knots = quantile(mean_dat$tstar, probs = seq(0, 1, length = Kt - 2))[-c(1, Kt - 2)]
      mean_basis   =  bs(c(t_min, t_max, mean_dat$tstar), knots = global_knots, intercept = TRUE)[-(1:2),]
    } 
    
    if (family == "gamma") {
      mean_family = stats::Gamma(link = "log")
    } else {
      mean_family = family
    }
    nrows_basis = nrow(mean_basis)
    # if greater than 10M, subsample
    if (nrows_basis > 10000000 && subsample) {
      if (verbose > 2) {
        message("Registr: Running Sub-sampling")
      }       
      uids = unique(mean_dat$id)
      avg_rows_per_id = nrows_basis / length(uids)
      size = round(10000000 / avg_rows_per_id)
      ids = sample(uids, size = size, replace = FALSE)
      rm(uids)
      subsampling_index = which(mean_dat$id %in% ids)
      rm(ids)
      mean_dat = mean_dat[subsampling_index, ]
      mean_basis = mean_basis[subsampling_index, ]
      rm(subsampling_index)
    }
    if (verbose > 2) {
      message("Registr: Running GLM")
    }   
    if (requireNamespace("fastglm", quietly = TRUE)) {
      mean_coefs = fastglm::fastglm(
        x = mean_basis, y = mean_dat$value,
        family = mean_family, method=2)
      mean_coefs = coef(mean_coefs)
    } else {
      mean_coefs = coef(glm(mean_dat$value ~ 0 + mean_basis, family = mean_family,
                            control = list(trace = verbose > 3)))
    }
    rm(mean_basis)
    rm(mean_dat)
  }
  
  ### Calculate warping functions  
  args = list(obj            = obj,            Y            = Y,
              Kt             = Kt,             Kh           = Kh,
              family         = family,         gradient     = gradient,
              incompleteness = incompleteness, lambda_inc   = lambda_inc,
              beta           = beta,           t_min        = t_min,
              t_max          = t_max,          rows         = rows,
              periodic       = periodic,       warping      = warping,
              global_knots   = global_knots,   mean_coefs   = mean_coefs,
              gamma_scales   = gamma_scales,
              ...)

  # main function call
  ids = Y$id
  Y = split(Y, ids)
  if (!is.null(beta)) {
    beta_list = apply(beta, 2, list)
    beta_list = lapply(beta_list, unlist)
    Y = mapply(function(y, b) {
      attr(y, "beta") = b
      y
    }, Y, beta_list, SIMPLIFY = FALSE)
    rm(beta_list)
  }
  if (!is.null(mean_coefs)) {
    if (is.matrix(mean_coefs)) {
      mc_list = apply(mean_coefs, 2, list)
      mc_list = lapply(mc_list, unlist)
    } else {
      mc_list = lapply(1:length(Y), function(blah) {
        mean_coefs
      })
    }
    Y = mapply(function(y, b) {
      attr(y, "mean_coefs") = b
      y
    }, Y, mc_list, SIMPLIFY = FALSE)
    rm(mc_list)
  }  
  args$Y = Y
  
  return(args)
}

#' Register (in)complete curves from exponential family
#' 
#' Function used in the registration step of an FPCA-based approach for 
#' registering exponential-family, potentially incomplete functional data,
#' called by \code{\link{register_fpca}}. 
#' This method uses constrained optimization to estimate spline 
#' coefficients for warping functions, where the objective function for optimization comes from 
#' maximizing the EF likelihood subject to monotonicity constraints on the warping functions. 
#' You have to either specify \code{obj}, which is a fpca 
#' object from an earlier step, or \code{Y}, a dataframe in long format with variables 
#' \code{id}, \code{index}, and \code{value} to indicate subject IDs, times, and observations, 
#' respectively. \cr \cr
#' Warping functions by default are forced to start and end on the diagonal to be
#' domain-preserving. This behavior can be changed by setting
#' \code{incompleteness} to some other value than NULL and a reasonable \code{lambda_inc} value.
#' For further details see the accompanying vignette. \cr \cr
#' By specifying \code{cores > 1} the registration call can be parallelized.
#' 
#' The template function for the registration is defined by argument \code{obj}
#' or \code{Y_template}, depending on if \code{obj} is NULL or not, respectively.
#' 
#' @param obj Current estimate of FPC object. 
#' Can be NULL only if Y argument is selected.
#' @param Y Dataframe. Should have values id, value, index.
#' @param Kt Number of B-spline basis functions used to estimate mean functions. Default is 8.
#' @param Kh Number of B-spline basis functions used to estimate warping functions \emph{h}. Default is 4.
#' @param family One of \code{c("gaussian","binomial","gamma","poisson")}. Defaults to
#' \code{"gaussian"}.
#' @param gradient If \code{TRUE}, uses analytic gradient to calculate derivative. 
#' If \code{FALSE}, calculates gradient numerically. Not available for families
#' \code{"gamma","poisson"}.
#' @param incompleteness Optional specification of incompleteness structure.
#' One of \code{c("leading","trailing","full")}, specifying that incompleteness
#' is present only in the initial measurements, only in the trailing measurements, or
#' in both, respectively. For details see the accompanying vignette.
#' Defaults to NULL, i.e. no incompleteness structure.
#' Can only be set when \code{warping = "nonparametric"}.
#' @param lambda_inc Penalization parameter to control the amount of
#' overall dilation of the domain.
#' The higher this lambda, the more the registered domains are forced to have the
#' same length as the observed domains.
#' Only used if \code{incompleteness} is not NULL.
#' @param Y_template Optional dataframe with the same structure as \code{Y}.
#' Only used if \code{obj} is NULL. If \code{Y_template} is NULL,
#' curves are registered to the overall mean of all curves in \code{Y} as template function.
#' If \code{Y_template} is specified, the template function is taken as the mean
#' of all curves in \code{Y_template}. Default is NULL.
#' @param beta Current estimates for beta for each subject. Default is NULL. 
#' @param t_min Minimum value to be evaluated on the time domain.
#' if `NULL`, taken to be minimum observed value.
#' @param t_max Maximum value to be evaluated on the time domain. 
#' if `NULL`, taken to be maximum observed value.
#' @param row_obj If NULL, the function cleans the data and calculates row indices. 
#' Keep this NULL if you are using standalone \code{registr} function.
#' @param periodic If \code{TRUE}, uses periodic b-spline basis functions. Default is \code{FALSE}.
#' @param warping If \code{nonparametric} (default), inverse warping functions are estimated nonparametrically. 
#' If \code{piecewise_linear2} they follow a piecewise linear function with 2 knots.
#' @param gamma_scales Only used for \code{family = "gamma"}.
#' Vector with one entry for each subject, containing the current estimate for the scale parameter of its
#' gamma distribution. Default is NULL, which sets the starting value for the scale parameter to 1.5.
#' @param cores Number of cores to be used. If \code{cores > 1}, the registration
#' call is parallelized by using \code{parallel::mclapply} (for Unix-based
#' systems) or \code{parallel::parLapply} (for Windows). Defaults to 1,
#' no parallelized call.
#' @param verbose Can be set to integers between 0 and 4 to control the level of
#' detail of the printed diagnostic messages. Higher numbers lead to more detailed
#' messages. Defaults to 1.
#' @param ... additional arguments passed to or from other functions
#' @param subsample if the number of rows of the data is greater than 
#' 10 million rows, the `id` values are subsampled to get the mean coefficients.
#' 
#' @return An list containing:
#' \item{Y}{The observed data. The variables \code{index} and \code{index_scaled}
#' contain the new estimated time domain.}
#' \item{loss}{Value of the loss function after registraton.}
#' \item{hinv_innerKnots}{List of inner knots for setting up the spline bases
#' for the inverse warping functions. Only contains \code{NULL} values for
#' \code{Kh <= 4}.}
#' \item{hinv_beta}{Matrix of B-spline basis coefficients used to construct
#' subject-specific inverse warping functions. See examples on how to
#' reconstruct a warping function based on \code{hinv_innerKnots} and
#' \code{hinv_beta}.}
#' 
#' @author Julia Wrobel \email{julia.wrobel@@cuanschutz.edu},
#' Erin McDonnell \email{eim2117@@cumc.columbia.edu},
#' Alexander Bauer \email{alexander.bauer@@stat.uni-muenchen.de}
#' @export
#' 
#' @importFrom stats glm coef Gamma quantile pbeta
#' @importFrom splines bs
#' @importFrom pbs pbs
#' @importFrom parallel mclapply makePSOCKcluster clusterExport clusterEvalQ parLapply stopCluster
#' 
#' @examples
#' ### complete binomial curves
#' Y = simulate_unregistered_curves()
#' register_step = registr(obj = NULL, Y = Y, Kt = 6, Kh = 4, family = "binomial", 
#'                         gradient = TRUE)
#' \donttest{
#' ### incomplete Gaussian curves
#' data(growth_incomplete)
#' 
#' # Force the warping functions to start and end on the diagonal to preserve the domain
#' register_step2a = registr(obj = NULL, Y = growth_incomplete, Kt = 6, Kh = 4,
#'                           family = "gaussian", gradient = TRUE,
#'                           incompleteness = NULL)
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   library(ggplot2)
#' 
#'   ggplot(register_step2a$Y, aes(x = tstar, y = index, group = id)) +
#'     geom_line(alpha = 0.2) +
#'     ggtitle("Estimated warping functions")
#'   ggplot(register_step2a$Y, aes(x = index, y = value, group = id)) +
#'     geom_line(alpha = 0.2) +
#'     ggtitle("Registered curves")
#' }
#'   
#' # Example for how to recreate an estimated inverse warping function given
#' # the output of registr(). Focus on id "boy01".
#' id         = "boy01"
#' index_obsRange_i = range(growth_incomplete$index[growth_incomplete$id == id])
#' index      = seq(min(index_obsRange_i), max(index_obsRange_i), length.out = 100)
#' # (note that 'index' must contain both the observed min and max in index_obsRange_i)
#' Theta_h_i  = splines::bs(index, knots = register_step2a$hinv_innerKnots[[id]], intercept = TRUE)
#' index_reg  = as.vector(Theta_h_i %*% register_step2a$hinv_beta[,id])
#' warp_dat_i = data.frame(index_observed   = index,
#'                         index_registered = index_reg)
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   ggplot(warp_dat_i, aes(x = index_observed, y = index_registered)) + geom_line() +
#'     ggtitle("Extracted warping function for id 'boy01'")
#' }
#' 
#' # Allow the warping functions to not start / end on the diagonal.
#' # The higher lambda_inc, the more the starting points and endpoints are
#' # forced towards the diagonal.
#' register_step2b = registr(obj = NULL, Y = growth_incomplete, Kt = 6, Kh = 4,
#'                           family = "gaussian", gradient = TRUE,
#'                           incompleteness = "full", lambda_inc = 1)
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   ggplot(register_step2b$Y, aes(x = tstar, y = index, group = id)) +
#'     geom_line(alpha = 0.2) +
#'     ggtitle("Estimated warping functions")
#'   ggplot(register_step2b$Y, aes(x = index, y = value, group = id)) +
#'     geom_line(alpha = 0.2) +
#'     ggtitle("Registered curves")
#' }
#' 
#' # Define the template function only over a subset of the curves
#' # (even though not very reasonable in this example)
#' template_ids    = c("girl12","girl13","girl14")
#' Y_template      = growth_incomplete[growth_incomplete$id %in% template_ids,]
#' register_step2c = registr(obj = NULL, Y = growth_incomplete, Kt = 6, Kh = 4,
#'                           family = "gaussian", gradient = TRUE,
#'                           Y_template = Y_template,
#'                           incompleteness = "full", lambda_inc = 1)
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   ggplot(register_step2c$Y, aes(x = index, y = value, group = id)) +
#'     geom_line(alpha = 0.2) +
#'     ggtitle("Registered curves")
#' }
#' 
#' }
#' 
registr = function(obj = NULL, Y = NULL, Kt = 8, Kh = 4, family = "gaussian", gradient = TRUE,
                   incompleteness = NULL, lambda_inc = NULL,
                   Y_template = NULL,
                   beta = NULL, t_min = NULL, t_max = NULL, row_obj = NULL,
                   periodic = FALSE, warping = "nonparametric",
                   gamma_scales = NULL, cores = 1L,  subsample = TRUE,
                   verbose = 1,
                   ...){
  
  
  args = registr_helper(
    obj = obj, 
    Y = Y, 
    Kt = Kt, 
    Kh = Kh, 
    family = family, 
    gradient = gradient,
    incompleteness = incompleteness, 
    lambda_inc = lambda_inc,
    Y_template = Y_template,
    beta = beta, 
    t_min = t_min,
    t_max = t_max, 
    row_obj = row_obj,
    periodic = periodic, 
    warping = warping,
    gamma_scales = gamma_scales, 
    cores = cores,  
    subsample = subsample,
    verbose = verbose,
    ...)
  
  mean_coefs = args$mean_coefs
  Y = args$Y
  rows = args$rows
  args$rows = NULL
  beta = args$beta
  tstar = Y$tstar
  family = args$family
  t_min = args$t_min
  t_max = args$t_max
  
  args$Y = NULL
  
  run_one_curve = function(r) {
    args$Y = r
    args$beta = attr(r, "beta")
    args$mean_coefs = attr(r, "mean_coefs")
    args$verbose = verbose
    do.call(registr_oneCurve, args = args)
  }
  if (cores == 1) { # serial call
    if (requireNamespace("pbapply", quietly = TRUE) && verbose > 1) {
      results_list = pbapply::pblapply(Y, run_one_curve)
    } else {
      results_list = lapply(Y, run_one_curve)
    }
    
  } else if (.Platform$OS.type == "unix") { # parallelized call on Unix-based systems
    results_list = parallel::mclapply(Y, run_one_curve, mc.cores = cores)
    
  } else { # parallelized call on Windows
    local_cluster = parallel::makePSOCKcluster(rep("localhost", cores)) # set up cluster
    # export functions and packages to the cluster
    parallel::clusterExport(cl = local_cluster, c("args"), envir=environment())
    parallel::clusterEvalQ(cl = local_cluster, c(library(registr), library(stats)))
    
    results_list = parallel::parLapply(
      cl  = local_cluster,
      X = Y, fun = run_one_curve)    
    
    parallel::stopCluster(cl = local_cluster) # close cluster
  }
  if (!is.null(beta)) {
    Y = lapply(Y, function(x) {
      attr(x, "beta") = NULL
      attr(x, "mean_coefs") = NULL
      x
    })
  }
  Y = dplyr::bind_rows(Y)
  res = gather_results_list(results_list, Y, rows, t_max, family)
  
  return(res) 
} 

gather_results_list = function(results_list, Y, rows, t_max, family) {
  
  tstar = Y$tstar
  # gather the results
  hinv_innerKnots        = lapply(results_list, function(x) x$hinv_innerKnots)
  names(hinv_innerKnots) = rows$id
  hinv_beta              = sapply(results_list, function(x) x$hinv_beta)
  colnames(hinv_beta)    = rows$id
  t_hat                  = unlist(sapply(results_list, function(x) as.vector(x$t_hat), simplify = FALSE), use.names = FALSE)
  loss_subjects          = unlist(sapply(results_list, function(x) as.vector(x$loss),  simplify = FALSE), use.names = FALSE)
  
  Y$index        = t_hat
  Y$index_scaled = t_hat / t_max
  Y$tstar        = tstar
  
  res = list(Y               = Y,
             loss            = sum(loss_subjects),
             hinv_innerKnots = hinv_innerKnots,
             hinv_beta       = hinv_beta)
  if (family == "gamma") {
    gamma_scales    = unlist(sapply(results_list, function(x) 
      as.vector(x$gamma_scale),  simplify = FALSE))
    res$gamma_scales = gamma_scales
  }
  return(res)
}


#' Internal function to register one curve
#' 
#' This internal function is only to be used from within \code{registr}.
#' It performs the main optimization step with \code{constrOptim} for the
#' registration of one curve.
#' 
#' @inheritParams registr
#' @param global_knots knots for the basis/splines, passed to [pbs::pbs()] 
#' or [stats::bs()]
#' @param mean_coefs Mean coefficients for the mean of all curves or 
#' GFPCA based.  May extract from `obj` object
#' @param just_return_list Do not use.  For developers only
#' 
#' @return An list containing:
#' \item{hinv_innerKnots}{Inner knots for setting up the spline basis
#' for the inverse warping function.}
#' \item{hinv_beta}{Estimated B-spline basis coefficients used to construct
#' subject-specific inverse warping functions.}
#' \item{t_hat}{Vector of registered time domain.}
#' \item{loss}{Loss of the optimal solution.}
#' 
#' @author Julia Wrobel \email{julia.wrobel@@cuanschutz.edu},
#' Erin McDonnell \email{eim2117@@cumc.columbia.edu},
#' Alexander Bauer \email{alexander.bauer@@stat.uni-muenchen.de}
#' 
#' @importFrom pbs pbs
#' @importFrom splines bs
#' @importFrom stats constrOptim Gamma
#' 
registr_oneCurve = function(
  obj = NULL,
  Y = NULL, 
  Kt = 8, 
  Kh = 4, 
  family = "gaussian",
  gradient = TRUE,
  incompleteness = NULL, 
  lambda_inc = NULL,
  beta = NULL, 
  t_min = NULL, 
  t_max = NULL,
  periodic = FALSE, 
  warping = "nonparametric",
  gamma_scales = NULL,
  global_knots = NULL,
  mean_coefs = NULL,
  ...,
  verbose = 1,
  just_return_list = FALSE) {
  
  t_range_i     = range(Y$tstar)
  t_min_i       = t_range_i[1]
  t_max_i       = t_range_i[2]
  Y_i           = Y
  #! needs to change
  tstar_cropped = Y$tstar
  rm(Y)
  D_i           = nrow(Y_i)
  
  loss_h_family = family
  
  # spline basis on the curve-specific time interval.
  # Pre-calculate the knots (for setting up the curve-specific Theta_h spline
  # bases in registr_oneCurve) similarly to splines::bs to make the final splines::bs call faster.
  degree          = 3 # cubic splines
  n_innerKnots    = Kh - (1 + degree)
  if (n_innerKnots <= 0) {
    hinv_innerKnots_i = NULL
  } else if (n_innerKnots > 0) {
    p_quantiles = seq.int(from = 0, to = 1, length.out = n_innerKnots + 2)[-c(1,n_innerKnots + 2)]
    hinv_innerKnots_i = quantile(tstar_cropped, probs = p_quantiles)
  }
  if (verbose > 2) {
    message("Getting spline basis")
  }
  tstar_bs_i      = c(Y_i$tstar, range(tstar_cropped))
  Theta_h_i       = splines::bs(tstar_bs_i, knots = hinv_innerKnots_i, intercept = TRUE)
  Theta_h_i       = Theta_h_i[1:nrow(Y_i),]
  
  # start parameters
  if (is.null(beta)) { # newly initialize start parameters
    t_vec_i     = Y_i$tstar
    beta_full_i = initial_params(warping = warping,
                                 K       = Theta_h_i,
                                 t_vec   = t_vec_i)
    
  } else { # use the last estimates (from the joint approach) as start parameters
    beta_full_i = beta
  }
  # account the beta vector for potential incompleteness constraints
  if (warping == "nonparametric") {
    if (is.null(incompleteness)) { # initial and final parameters are fixed
      beta_i = beta_full_i[-c(1, length(beta_full_i))]
    } else if (incompleteness == "leading") { # final parameter is fixed
      beta_i = beta_full_i[-length(beta_full_i)]
    } else if (incompleteness == "trailing") { # initial parameter is fixed
      beta_i = beta_full_i[-1]
    } else if (incompleteness == "full") { # no parameter is fixed
      beta_i = beta_full_i
    }
  } else if (warping == "piecewise_linear2") {
    beta_i = beta_full_i
  }
  
  
  if (verbose > 2) {
    message("Getting mean coefficients")
  }
  if (is.null(obj)) { # template function = mean of all curves
    mean_coefs_i = mean_coefs
    
  } else { # template function = current GFPCA representation
    
    if (obj$fpca_type == "variationalEM") { # GFPCA based on fpca_gauss or bfpca
      # In this case, the FPCA is explicitly based on the spline basis 'mean_basis'
      mean_coefs_i = mean_coefs
      
      if (family %in% c("gamma","poisson")) {
        # special family only for loss_h since the variationalEM FPCA for Gamma/Poisson data
        # is performed with Gaussian family, which raises the need to handle
        # some FPCA results in loss_h a little differently.
        loss_h_family = paste0(family,"-varEM")
      }
      
    } else { # GFPCA based on gfpca_twoStep
      # In this case, the FPCA is not based on the spline basis 'mean_basis'.
      # Accordingly, smooth over the GFPCA representation of the i'th function
      # using 'mean_basis'.
      mean_dat_i = obj$Yhat[obj$Yhat$id == Y_i$id[1],]
      
      if (periodic) {
        mean_basis = pbs::pbs(c(t_min, t_max, mean_dat_i$index),
                              knots = global_knots, intercept = TRUE)[-(1:2),]
      } else {
        mean_basis = splines::bs(c(t_min, t_max, mean_dat_i$index),
                                 knots = global_knots, intercept = TRUE)[-(1:2),]
      } 
      
      if (family == "gamma") {
        mean_family      = stats::Gamma(link = "log")
        mean_dat_i$value = exp(mean_dat_i$value)
        # set very small positive values to some lowest threshold to prevent numerical problems
        mean_dat_i$value[mean_dat_i$value > 0 & mean_dat_i$value < 1e-4] = 1e-4
        
      } else {
        mean_family = "gaussian"
      }
      mean_coefs_i = coef(glm(value ~ 0 + mean_basis, data = mean_dat_i, family = mean_family))
    }
  }
  
  if (gradient) { 
    gradf = loss_h_gradient 
  } else { 
    gradf = NULL 
  }
  
  # optimization constraints
  if (is.null(incompleteness)) { # warping functions start and end on the diagonal
    if (warping == "nonparametric") {
      const_i = constraints(Kh, t_min_i, t_max_i, warping = warping)
      ui_i    = const_i$ui
      ui_i    = ui_i[-nrow(ui_i), -ncol(ui_i)]
      ci_i    = const_i$ci
      ci_i    = ci_i[-(length(ci_i) - 1)]
    } else if (warping == "piecewise_linear2") {
      const_i = constraints(Kh - 1, t_min, t_max_i, warping = warping)
      ui_i    = const_i$ui
      ci_i    = const_i$ci
    }
    
  } else { # warping functions not necessarily start and/or end on the diagonal
    dim_const   = Kh + ifelse(incompleteness == "full", 1, 0)
    t_min_const = ifelse(incompleteness == "trailing", t_min_i, t_min)
    t_max_const = ifelse(incompleteness == "leading",  t_max_i, t_max)
    
    const_i     = constraints(dim_const, t_min_const, t_max_const)
    ui_i        = const_i$ui
    ci_i        = const_i$ci
  }
  
  if (family == "gamma") { # add the scale parameter to be optimized as last element of beta_i
    ui_i   = cbind(ui_i, 0)
    ui_i   = rbind(ui_i, c(rep(0, ncol(ui_i) - 1), 1))
    ci_i   = append(ci_i, 0)
    scale  = ifelse(!is.null(gamma_scales), gamma_scales, 1.5)
    beta_i = append(beta_i, scale)
  }
  
  
  # workaround: sometimes constrOptim states that the starting values are not
  # inside the feasible region. This is only a numerical error, so let's simply
  # substract a minor value from ci
  # (source: 
  # https://stackoverflow.com/questions/50472525/constroptim-in-r-init-val-is-not-in-the-interior-of-the-feasible-region-error
  # )
  ci_i = ci_i - 1e-6
  
  # when an analytic gradient is used, constrOptim sometimes leads to beta
  # values slightly outside the possible domain, or to slightly nonmonotone beta
  # values that don't fulfill the constraints.
  # Correct these slight inconsistencies to ensure proper beta values.
  if (warping != "piecewise_linear2") {
    
    if (family == "gamma") { # remove scale parameter as last element
      scale = tail(beta_i, 1)
      beta_i = beta_i[1:(length(beta_i) - 1)]
    }
    if (verbose > 2) {
      message("Ensuring proper beta")
    }
    beta_i = ensure_proper_beta(beta  = beta_i,
                                t_min = ifelse(is.null(incompleteness) || incompleteness == "trailing",
                                               t_min_i, t_min),
                                t_max = ifelse(is.null(incompleteness) || incompleteness == "leading",
                                               t_max_i, t_max))
    
    if (family == "gamma") # add scale parameter again as last element
      beta_i = c(beta_i, scale)
  }
  
  if (verbose > 2) {
    message("Running optimization")
  }
  out_args = list(theta          = beta_i,
                  f              = loss_h,
                  grad           = gradf,
                  ui             = ui_i,
                  ci             = ci_i,
                  Y              = Y_i$value, 
                  Theta_h        = Theta_h_i,
                  mean_coefs     = mean_coefs_i, 
                  knots          = global_knots, 
                  family         = loss_h_family,
                  incompleteness = incompleteness,
                  lambda_inc     = lambda_inc,
                  t_min          = t_min,
                  t_max          = t_max,
                  t_min_curve    = t_min_i,
                  t_max_curve    = t_max_i,
                  periodic       = periodic,
                  Kt             = Kt,
                  warping        = warping,
                  ...)
  if (just_return_list) {
    return(out_args)
  }
  # main registration step	
  beta_optim = do.call(constrOptim, args = out_args)
  
  beta_inner = beta_optim$par
  
  if (family == "gamma") {
    scale      = tail(beta_inner, 1)
    beta_inner = beta_inner[1:(length(beta_inner)-1)]
  }
  
  if (warping == "nonparametric") {
    # create full beta vector
    if (is.null(incompleteness)) { # initial and final parameters are fixed
      beta_full_i = c(t_min_i, beta_inner, t_max_i)
    } else if (incompleteness == "leading") { # final parameter is fixed
      beta_full_i = c(beta_inner, t_max_i)
    } else if (incompleteness == "trailing") { # initial parameter is fixed
      beta_full_i = c(t_min_i, beta_inner)
    } else if (incompleteness == "full") { # no parameter is fixed
      beta_full_i = beta_inner
    }
    #t_hat = as.vector(cbind(1, Theta_h_i) %*% beta_full_i)
    t_hat = as.vector(Theta_h_i %*% beta_full_i)
    
  } else if (warping == "piecewise_linear2") {
    beta_full_i = beta_inner
    t_hat       = piecewise_linear2_hinv(grid = seq(0, t_max_i, length.out = D_i),
                                         knot_locations = beta_inner)
  }
  
  res = list(hinv_innerKnots = hinv_innerKnots_i,
             hinv_beta       = beta_full_i,
             t_hat           = t_hat,
             loss            = beta_optim$value)
  if (family == "gamma")
    res$gamma_scale = scale
  
  return(res)
  
}
