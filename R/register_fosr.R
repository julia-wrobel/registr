#' Covariate-Adjusted Registration Using VL-FoSR
#'
#' Alternates VL-FoSR fits with majorization-minimization (MM) or
#' constrained optimization (CO) to estimate subject-specific inverse
#' warping functions for binomial functional data.
#'
#' @param Y Long-format data containing \code{id}, \code{index},
#' \code{value}, and subject-level covariates.
#' @inheritParams registr_MM
#' @param max_iterations Maximum number of alternating registration iterations.
#' @param npc Number of functional principal components.
#' @param fosr_formula Fixed-effect formula used in VL-FoSR.
#' @param Kt_final Number of time B-spline basis functions for the final fit.
#' @param method_registration Registration method, \code{"MM"} or \code{"CO"}.
#' @param fosr_covariate_model Covariate model, one of \code{"linear"},
#' \code{"additive"}, or \code{"mixed"}.
#' @param fosr_M_spline_x Number of B-spline basis functions for additive
#' covariate effects.
#' @param fosr_additive_cols Names of additive covariates.
#' @param fosr_linear_cols Names of linear covariates.
#' @param verbose Verbosity level.
#' @param cores Number of C++ threads for VL-FoSR or workers for CO;
#' MM registration remains serial.
#' @param ... Additional arguments passed to the registration optimizer.
#' @return An object of class \code{registration} containing:
#' \item{fpca_obj}{Final VL-FoSR fit.}
#' \item{Y}{Input data with original \code{tstar} and estimated \code{t_hat}.}
#' \item{time_warps}{Registered time vectors from each iteration.}
#' \item{loss, delta}{Registration losses and alternating convergence values.}
#' \item{nu}{Interior inverse warping spline coefficients by subject.}
#' \item{family, method_registration}{Model family and registration method.}
#' The object also retains the covariate-model configuration.
#' @author Linlin Wu \email{linlin.wu@@emory.edu},
#' Jeff Goldsmith \email{ajg2202@@cumc.columbia.edu},
#' Julia Wrobel \email{julia.wrobel@@emory.edu}
#' @export
register_fosr <- function(
    Y, Kt = 8, Kh = 4, family = "binomial",
    max_iterations = 10, npc = 1, x_cols = NULL,
    fosr_formula = value ~ 1, Kt_final = Kt,
    method_registration = c("MM", "CO"),
    fosr_covariate_model = c("linear", "additive", "mixed"),
    fosr_M_spline_x = 5, fosr_additive_cols = NULL,
    fosr_linear_cols = NULL, verbose = 0, cores = 1L, ...
) {
  method_registration <- match.arg(method_registration)
  cores <- as.integer(cores)
  if (family != "binomial") stop("Only the binomial family is implemented.")
  fosr_covariate_model <- match.arg(fosr_covariate_model)
  x_cols <- unique(as.character(x_cols))
  formula_cols <- all.vars(stats::delete.response(stats::terms(fosr_formula)))
  if (!setequal(x_cols, formula_cols)) {
    stop("x_cols must match the covariates in fosr_formula.")
  }

  spec <- normalize_covariate_spec(
    covariate_names = x_cols,
    covariate_model = fosr_covariate_model,
    additive_cols = fosr_additive_cols,
    linear_cols = fosr_linear_cols
  )
  fosr_additive_cols <- spec$additive_cols
  fosr_linear_cols <- spec$linear_cols

  dat <- registr::data_clean(Y)
  Y <- dat$Y
  rows <- dat$Y_rows
  Y$tstar <- Y$index

  time_warps <- vector("list", max_iterations + 2)
  time_warps[[1]] <- Y$index
  loss_vec <- rep(NA_real_, max_iterations + 1)

  # Initial registration: both methods use the SAME covariate-adjusted
  # fixed-effect templates, but the MM path stays identical to the old code.
  if (verbose > 0) message("Initial registration (", method_registration, ")...")
  if (method_registration == "MM") {
    registr_step <- registr_MM(
      Y = Y, Kt = Kt, Kh = Kh, family = family,
      row_obj = rows, x_cols = x_cols,
      covariate_model = fosr_covariate_model,
      M_spline_x = fosr_M_spline_x,
      additive_cols = fosr_additive_cols,
      linear_cols = fosr_linear_cols, ...
    )
  } else {
    # Julia's registr() ordinarily fits a pooled initial template when obj
    # is NULL. Supply our covariate-adjusted template explicitly instead.
    init_template <- build_initial_template_by_covariates(
      Y_raw = Y, Kt = Kt, x_cols = x_cols,
      covariate_model = fosr_covariate_model,
      M_spline_x = fosr_M_spline_x,
      additive_cols = fosr_additive_cols,
      linear_cols = fosr_linear_cols
    )
    init_obj <- list(
      Y = Y,
      knots = init_template$knots,
      subject_coefs = init_template$subject_coefs,
      fpca_type = "variationalEM"
    )
    registr_step <- registr::registr(
      obj = init_obj, Y = Y, Kt = Kt, Kh = Kh,
      family = family, row_obj = rows, cores = cores,
      verbose = verbose, ...
    )
  }
  time_warps[[2]] <- registr_step$Y$index
  loss_vec[1] <- registr_step$loss

  # Alternate VL-FoSR fits and MM updates.
  iter <- 1
  delta <- rep(NA_real_, max_iterations)
  delta[1] <- 1e2

  while (iter < max_iterations && delta[iter] > 1e-4) {
    if (verbose > 0) message("Iteration: ", iter)

    fpca_step <- fosr(
      data = registr_step$Y, formula = fosr_formula,
      M_spline = Kt, npc = npc,
      maxiter = 50, error_thresh = 1e-4,
      print.iter = (verbose > 1), seed = 1988 + iter,
      row_obj = rows, covariate_model = fosr_covariate_model,
      M_spline_x = fosr_M_spline_x,
      additive_cols = fosr_additive_cols,
      linear_cols = fosr_linear_cols, cores = cores
    )

    # Reuse the preceding warping coefficients as the optimizer's start.
    if (method_registration == "MM") {
      registr_step <- registr_MM(
        obj = fpca_step, Kt = Kt, Kh = Kh,
        family = family, row_obj = rows,
        nu = registr_step$nu, x_cols = x_cols, ...
      )
    } else {
      # CO expects the full Kh x I coefficient matrix (including endpoints).
      registr_step <- registr::registr(
        obj = fpca_step, Kt = Kt, Kh = Kh,
        family = family, row_obj = rows,
        beta = registr_step$hinv_beta, cores = cores,
        verbose = verbose, ...
      )
    }

    time_warps[[iter + 2]] <- registr_step$Y$index
    loss_vec[iter + 1] <- registr_step$loss
    delta[iter + 1] <- mean(
      (time_warps[[iter + 2]] - time_warps[[iter + 1]])^2,
      na.rm = TRUE
    )
    iter <- iter + 1
  }

  # Final VL-FoSR refit at the estimated registered times.
  if (verbose > 0) message("Final FoSR...")
  fpca_final <- fosr(
    data = registr_step$Y, formula = fosr_formula,
    M_spline = Kt_final, npc = npc,
    maxiter = 50, error_thresh = 1e-4,
    print.iter = (verbose > 1), seed = 1988 + iter,
    row_obj = rows, covariate_model = fosr_covariate_model,
    M_spline_x = fosr_M_spline_x,
    additive_cols = fosr_additive_cols,
    linear_cols = fosr_linear_cols, cores = cores
  )

  Y$t_hat <- registr_step$Y$index
  # Keep the original MM output; CO returns full warping coefficients.
  if (method_registration == "MM") {
    nu_df <- as.data.frame(t(registr_step$nu))
  } else {
    beta_full <- registr_step$hinv_beta
    beta_full <- beta_full[, as.character(rows$id), drop = FALSE]
    nu_df <- as.data.frame(t(beta_full[2:(Kh - 1), , drop = FALSE]))
  }
  nu_df$id <- rows$id

  ret <- list(
    fpca_obj = fpca_final, Y = Y,
    time_warps = Filter(Negate(is.null), time_warps),
    loss = loss_vec[!is.na(loss_vec)], family = family,
    nu = nu_df, delta = delta[!is.na(delta)],
    fosr_covariate_model = fosr_covariate_model,
    fosr_M_spline_x = fosr_M_spline_x,
    fosr_additive_cols = fosr_additive_cols,
    fosr_linear_cols = fosr_linear_cols,
    method_registration = method_registration
  )
  class(ret) <- "registration"
  ret
}

