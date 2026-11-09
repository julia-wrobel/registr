#' Construct Initial Covariate-Adjusted Templates
#'
#' Uses a binomial ridge fit to construct subject-specific fixed-effect
#' time-spline templates for MM or CO registration.
#'
#' @inheritParams registr_MM
#' @inheritParams fosr
#' @param Y_raw Long-format input data containing \code{id}, \code{tstar},
#' \code{value}, and the selected covariates.
#' @return A list with \code{knots} (interior time-spline knots) and
#' \code{subject_coefs} (subject-specific template coefficients).
build_initial_template_by_covariates <- function(
    Y_raw, Kt, x_cols = NULL,
    covariate_model = c("linear", "additive", "mixed"),
    M_spline_x = 5, additive_cols = NULL,
    linear_cols = NULL, lambda_init = 0
) {
  covariate_model <- match.arg(covariate_model)
  x_cols <- unique(as.character(x_cols))
  id_vec <- as.character(Y_raw$id)
  id_order <- sort(unique(id_vec))
  first_index <- match(id_order, id_vec)
  subject_idx <- match(id_vec, id_order)

  tstar_all <- as.numeric(Y_raw$tstar)
  t_min <- min(tstar_all)
  t_max <- max(tstar_all)
  global_knots <- stats::quantile(
    tstar_all, probs = seq(0, 1, length.out = Kt - 2)
  )[-c(1, Kt - 2)]

  B_full <- splines::bs(
    c(t_min, t_max, tstar_all),
    knots = global_knots, intercept = TRUE,
    Boundary.knots = c(t_min, t_max)
  )[-(1:2), , drop = FALSE]
  B_sp <- Matrix::Matrix(B_full, sparse = TRUE)

  # Use the same centered additive features as VL-FoSR.
  features <- prepare_subject_features(
    subject_df = Y_raw[first_index, , drop = FALSE],
    covariate_names = x_cols,
    covariate_model = covariate_model,
    M_spline_x = M_spline_x,
    additive_cols = additive_cols,
    linear_cols = linear_cols
  )
  fixed_features <- features$fixed_features_by_id
  p_fixed <- ncol(fixed_features)

  if (p_fixed > 0) {
    F_obs <- fixed_features[subject_idx, , drop = FALSE]
    blocks <- lapply(seq_len(p_fixed), function(j) {
      Matrix::Diagonal(x = as.numeric(F_obs[, j])) %*% B_sp
    })
    design_sp <- do.call(cbind, c(list(B_sp), blocks))
  } else {
    design_sp <- B_sp
  }

  fit <- glmnet::glmnet(
    x = Matrix::drop0(design_sp), y = Y_raw$value,
    family = "binomial", alpha = 0, lambda = lambda_init,
    standardize = FALSE, intercept = FALSE,
    thresh = 1e-6, maxit = 500
  )
  coefs <- drop(as.matrix(stats::coef(fit, s = lambda_init))[-1, , drop = FALSE])
  alpha0 <- coefs[seq_len(Kt)]

  if (p_fixed > 0) {
    alpha_fixed <- matrix(coefs[-seq_len(Kt)], nrow = Kt)
    colnames(alpha_fixed) <- colnames(fixed_features)
    subject_coefs <- matrix(alpha0, Kt, length(id_order)) +
      alpha_fixed %*% t(fixed_features)
  } else {
    alpha_fixed <- matrix(numeric(0), Kt, 0)
    subject_coefs <- matrix(alpha0, Kt, length(id_order))
  }
  colnames(subject_coefs) <- id_order

  list(
    knots = global_knots,
    subject_coefs = subject_coefs
  )
}

#' Register Binomial Curves by Majorization-Minimization
#'
#' Estimates monotone subject-specific inverse warping functions using a
#' quadratic majorization step and backtracking line search. Registration
#' is serial and accepts warm-start coefficients.
#'
#' @param obj Optional VL-FoSR object with \code{Y}, \code{knots}, and
#' \code{subject_coefs}. If \code{NULL}, initial templates are fitted.
#' @param Y Long-format data frame with \code{id}, \code{index}, \code{value},
#' and selected covariates; required if \code{obj} is \code{NULL}.
#' @param Kt Number of cubic time B-spline basis functions for initial templates.
#' @param Kh Number of B-spline basis functions for inverse warping (at least 4).
#' @param family Distribution family; only \code{"binomial"} is supported.
#' @param nu Optional \eqn{(Kh-2) \times I} matrix of starting interior
#' warping coefficients.
#' @param row_obj Optional subject row locations returned by
#' \code{data_clean()}.
#' @param x_cols Names of subject-level covariates for initial templates.
#' @inheritParams fosr
#' @param mm_max_iter Maximum number of MM iterations per subject.
#' @param mm_tol Relative coefficient-change tolerance.
#' @param mm_ridge Positive ridge term in the quadratic subproblem.
#' @param mm_max_backtrack Maximum number of backtracking steps.
#' @param mm_armijo_eps Required decrease in Bernoulli loss to accept a step.
#' @param mm_shrink Multiplicative backtracking step-size factor.
#' @return A list containing \code{Y} with \code{index} replaced by
#' registered time, \code{loss} (sum of subject Bernoulli losses),
#' and \code{nu} (matrix of interior warping coefficients).
#' @author Linlin Wu \email{linlin.wu@@emory.edu},
#' Jeff Goldsmith \email{ajg2202@@cumc.columbia.edu},
#' Julia Wrobel \email{julia.wrobel@@emory.edu}
#' @export
registr_MM <- function(
    obj = NULL, Y = NULL, Kt = 8, Kh = 4,
    family = "binomial", nu = NULL, row_obj = NULL,
    x_cols = NULL, covariate_model = c("linear", "additive", "mixed"),
    M_spline_x = 5, additive_cols = NULL, linear_cols = NULL,
    mm_max_iter = 50, mm_tol = 1e-6, mm_ridge = 1e-4,
    mm_max_backtrack = 8, mm_armijo_eps = 1e-6, mm_shrink = 0.5
) {
  if (family != "binomial") stop("Only the binomial family is implemented.")
  if (Kh < 4) stop("Kh must be at least 4.")
  covariate_model <- match.arg(covariate_model)

  if (is.null(Y)) Y <- obj$Y
  if (is.null(obj)) Y$tstar <- Y$index

  if (is.null(row_obj)) {
    dat <- data_clean(Y)
    Y <- dat$Y
    rows <- dat$Y_rows
  } else {
    rows <- row_obj
  }
  I <- nrow(rows)

  t_min0 <- min(Y$tstar)
  t_max0 <- max(Y$tstar)

  if (is.null(obj)) {
    template <- build_initial_template_by_covariates(
      Y_raw = Y, Kt = Kt, x_cols = x_cols,
      covariate_model = covariate_model,
      M_spline_x = M_spline_x,
      additive_cols = additive_cols, linear_cols = linear_cols
    )
    global_knots <- template$knots
    subject_template_coefs <- template$subject_coefs
  } else {
    global_knots <- obj$knots
    subject_template_coefs <- obj$subject_coefs
    if (is.null(colnames(subject_template_coefs))) {
      colnames(subject_template_coefs) <- as.character(rows$id)
    }
  }
  subject_template_coefs <- subject_template_coefs[
    , as.character(rows$id), drop = FALSE
  ]

  # Cubic template spline knots for fast predictor evaluation.
  ord <- 4
  template_knots_full <- c(
    rep(t_min0, ord), global_knots, rep(t_max0, ord)
  )

  t_hat <- rep(NA_real_, nrow(Y))
  loss_subjects <- rep(NA_real_, I)
  nu_new <- matrix(NA_real_, Kh - 2, I)
  loss_trace <- vector("list", I)

  for (i in seq_len(I)) {
    idx <- rows$first_row[i]:rows$last_row[i]
    Yi <- Y$value[idx]
    tstar_i <- Y$tstar[idx]
    t_min_i <- min(tstar_i)
    t_max_i <- max(tstar_i)

    # Warping knots and design matrix may differ by subject.
    n_inner_knots <- Kh - 4
    inner_knots <- NULL
    if (n_inner_knots > 0) {
      probs <- seq(0, 1, length.out = n_inner_knots + 2)
      inner_knots <- stats::quantile(tstar_i, probs = probs[-c(1, length(probs))])
    }

    Theta_full <- splines::bs(
      c(tstar_i, t_min_i, t_max_i),
      knots = inner_knots, intercept = TRUE,
      Boundary.knots = c(t_min_i, t_max_i)
    )[seq_along(tstar_i), , drop = FALSE]
    Theta_inner <- Theta_full[, 2:(Kh - 1), drop = FALSE]

    # Fix endpoint coefficients and constrain the interior coefficients.
    constrs <- constraints(Kh, t_min_i, t_max_i, warping = "nonparametric")
    ui <- constrs$ui[-nrow(constrs$ui), -ncol(constrs$ui), drop = FALSE]
    ci <- constrs$ci[-(length(constrs$ci) - 1)]

    nu0 <- seq(t_min_i, t_max_i, length.out = Kh)[-c(1, Kh)]
    nu_inner <- if (is.null(nu)) nu0 else as.numeric(nu[, i])
    base_t <- Theta_full[, 1] * t_min_i + Theta_full[, Kh] * t_max_i
    rebuild_t <- function(coefs) drop(base_t + Theta_inner %*% coefs)

    mean_coefs_i <- as.numeric(subject_template_coefs[, as.character(rows$id[i])])
    loss_fn <- function(coefs) {
      loss_h_binomial_nonparametric_cpp(
        Y = Yi, Theta_h = Theta_full,
        mean_coefs = mean_coefs_i,
        template_knots_full = template_knots_full,
        beta_inner = coefs,
        t_min_curve = t_min_i, t_max_curve = t_max_i,
        inc_mode = 0, lambda_inc = 0
      )
    }

    loss_cur <- loss_fn(nu_inner)
    trace_i <- rep(NA_real_, mm_max_iter)

    for (mm in seq_len(mm_max_iter)) {
      t_cur <- pmin(pmax(rebuild_t(nu_inner), t_min_i), t_max_i)
      q <- eval_mm_quantities_cpp(
        t_vec = t_cur, knots_full = template_knots_full,
        alpha = mean_coefs_i, ord = ord, eps = 1e-8
      )

      w <- as.numeric(q$w)
      H <- 0.25 * crossprod(
        Theta_inner, Theta_inner * as.numeric(q$W2)
      ) + diag(mm_ridge, ncol(Theta_inner))
      g <- as.numeric(crossprod(
        Theta_inner, w * (as.numeric(q$mu) - Yi)
      ))

      b_shift <- as.numeric(ci - ui %*% nu_inner)
      sol <- quadprog::solve.QP(
        Dmat = H, dvec = -g,
        Amat = t(ui), bvec = b_shift
      )
      z <- as.numeric(sol$solution)

      # Accept a step only if it decreases the original Bernoulli loss.
      accepted <- FALSE
      step <- 1
      if (is.finite(loss_cur)) {
        loss_try <- loss_fn(nu_inner + step * z)
        if (is.finite(loss_try) && loss_try <= loss_cur - mm_armijo_eps) {
          accepted <- TRUE
        } else {
          for (bt in seq_len(mm_max_backtrack)) {
            step <- step * mm_shrink
            if (step < 1e-6 || bt >= mm_max_backtrack) break
            loss_try <- loss_fn(nu_inner + step * z)
            if (is.finite(loss_try) && loss_try <= loss_cur - mm_armijo_eps) {
              accepted <- TRUE
              break
            }
          }
        }
      }

      nu_next <- if (accepted) nu_inner + step * z else nu_inner
      if (accepted) loss_cur <- loss_try
      trace_i[mm] <- loss_cur
      relchg <- sqrt(sum((nu_next - nu_inner)^2)) /
        (sqrt(sum(nu_inner^2)) + 1e-12)
      nu_inner <- nu_next
      if (relchg < mm_tol) break
    }

    nu_new[, i] <- nu_inner
    t_hat[idx] <- rebuild_t(nu_inner)
    loss_subjects[i] <- loss_cur
    loss_trace[[i]] <- trace_i[seq_len(mm)]
  }

  colnames(nu_new) <- as.character(rows$id)
  Y$index <- t_hat
  list(
    Y = Y,
    loss = sum(loss_subjects),
    nu = nu_new
  )
}

