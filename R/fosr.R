#' Construct a Centered Covariate Spline Basis
#'
#' Builds and centers a cubic B-spline basis for an additive covariate.
#'
#' @param x Numeric vector of observed covariate values.
#' @param M_spline_x Number of B-spline basis functions for the covariate.
#' @param degree_x Degree of the covariate B-splines.
#' @return A list containing the centered and raw bases (\code{basis},
#' \code{raw_basis}), their column means (\code{center}), spline knots,
#' boundary knots, degree, basis dimension, and covariate range and mean.
make_x_spline_basis <- function(x, M_spline_x = 5, degree_x = 3) {
  x <- as.numeric(x)
  if (any(!is.finite(x)) || length(unique(x)) < 4) {
    stop("An additive covariate needs at least four distinct finite values.")
  }
  Bx_obj <- splines::bs(
    x, df = M_spline_x, degree = degree_x, intercept = FALSE
  )
  Bx_raw <- unclass(Bx_obj)
  center <- colMeans(Bx_raw)
  list(
    basis = sweep(Bx_raw, 2, center, "-"),
    raw_basis = Bx_raw,
    center = center,
    knots = attr(Bx_obj, "knots"),
    Boundary.knots = attr(Bx_obj, "Boundary.knots"),
    degree = degree_x,
    M_spline_x = ncol(Bx_raw),
    x_range = range(x),
    x_mean = mean(x)
  )
}

#' Evaluate a Centered Covariate Spline Basis
#'
#' Evaluates a previously fitted additive-covariate basis at new values,
#' using the original centering constants.
#'
#' @inheritParams make_x_spline_basis
#' @param surface_info Spline metadata returned by
#' \code{prepare_subject_features()}.
#' @return A matrix of centered B-spline basis values.
eval_x_spline_basis <- function(x, surface_info) {
  Bx <- splines::bs(
    as.numeric(x),
    knots = surface_info$knots,
    degree = surface_info$degree,
    intercept = FALSE,
    Boundary.knots = surface_info$Boundary.knots
  )
  sweep(unclass(Bx), 2, surface_info$center, "-")
}

#' Normalize the Covariate Specification
#'
#' Checks that covariates are consistently assigned to linear and additive
#' effects under the selected model.
#'
#' @param covariate_names Character vector of covariate names.
#' @param covariate_model One of \code{"linear"}, \code{"additive"}, or
#' \code{"mixed"}.
#' @param additive_cols Character vector of additive covariate names.
#' @param linear_cols Character vector of linear covariate names.
#' @return A list with \code{covariate_model}, \code{additive_cols}, and
#' \code{linear_cols}.
normalize_covariate_spec <- function(
    covariate_names,
    covariate_model = c("linear", "additive", "mixed"),
    additive_cols = NULL,
    linear_cols = NULL
) {
  covariate_model <- match.arg(covariate_model)
  covariate_names <- unique(as.character(covariate_names))
  additive_cols <- unique(as.character(additive_cols))
  linear_cols <- unique(as.character(linear_cols))

  if (covariate_model == "linear") {
    if (length(additive_cols)) {
      stop("Linear mode cannot have additive covariates.")
    }
    if (!length(linear_cols)) linear_cols <- covariate_names

  } else if (covariate_model == "additive") {
    if (length(linear_cols)) {
      stop("Additive mode cannot have linear covariates.")
    }
    if (!length(additive_cols)) additive_cols <- covariate_names

  } else if (!length(additive_cols) || !length(linear_cols)) {
    stop("Mixed mode needs both additive_cols and linear_cols.")
  }

  if (covariate_model == "additive" && !length(additive_cols)) {
    stop("Additive mode needs at least one covariate.")
  }

  if (length(intersect(additive_cols, linear_cols)) ||
      !setequal(c(additive_cols, linear_cols), covariate_names)) {
    stop("additive_cols and linear_cols must partition the formula covariates.")
  }

  list(
    covariate_model = covariate_model,
    additive_cols = additive_cols,
    linear_cols = linear_cols
  )
}

#' Resolve Additive Spline Basis Dimensions
#'
#' Matches user-supplied basis dimensions to additive covariates.
#'
#' @param M_spline_x A scalar, a vector with one entry per additive
#' covariate, or a named vector of basis dimensions.
#' @inheritParams normalize_covariate_spec
#' @return A named integer vector of basis dimensions for additive covariates.
resolve_additive_basis_sizes <- function(M_spline_x, additive_cols) {
  if (!length(additive_cols)) {
    return(setNames(integer(0), character(0)))
  }
  if (!is.null(names(M_spline_x))) {
    sizes <- M_spline_x[additive_cols]
  } else if (length(M_spline_x) == 1) {
    sizes <- rep(M_spline_x, length(additive_cols))
  } else if (length(M_spline_x) == length(additive_cols)) {
    sizes <- M_spline_x
  } else {
    stop("M_spline_x must be a scalar or one value per additive covariate.")
  }
  if (anyNA(sizes) || any(!is.finite(sizes)) ||
      any(sizes < 3) || any(sizes != as.integer(sizes))) {
    stop("Each additive basis size must be an integer >= 3.")
  }
  setNames(as.integer(sizes), additive_cols)
}

#' Prepare Subject-Level Covariate Features
#'
#' Constructs the fixed-effect design matrix and stores the spline metadata
#' needed to evaluate additive effects.
#'
#' @param subject_df Data frame with one row per subject and the selected
#' covariates.
#' @inheritParams normalize_covariate_spec
#' @inheritParams resolve_additive_basis_sizes
#' @return A list containing the subject-level matrices \code{Q_by_id},
#' \code{X_by_id}, and \code{fixed_features_by_id}; feature names;
#' covariate specifications; and additive spline metadata
#' (\code{surface_infos} and, for one additive covariate,
#' \code{surface_info}).
prepare_subject_features <- function(
    subject_df,
    covariate_names,
    covariate_model = c("linear", "additive", "mixed"),
    M_spline_x = 5,
    additive_cols = NULL,
    linear_cols = NULL
) {
  spec <- normalize_covariate_spec(
    covariate_names, covariate_model, additive_cols, linear_cols
  )
  additive_cols <- spec$additive_cols
  linear_cols <- spec$linear_cols
  cols <- c(linear_cols, additive_cols)
  I <- nrow(subject_df)

  # All covariates are finite numeric variables.
  if (length(cols) && (!all(cols %in% names(subject_df)) ||
                       !all(vapply(subject_df[cols], function(x) {
                         is.numeric(x) && all(is.finite(x))
                       }, logical(1))))) {
    stop("All specified covariates must exist and contain finite numeric values.")
  }

  sizes <- resolve_additive_basis_sizes(M_spline_x, additive_cols)
  additive_blocks <- vector("list", length(additive_cols))
  surface_infos <- vector("list", length(additive_cols))
  names(additive_blocks) <- additive_cols
  names(surface_infos) <- additive_cols
  # Center each additive basis separately.
  for (col in additive_cols) {
    bx <- make_x_spline_basis(subject_df[[col]], sizes[[col]])
    feature_names <- paste0(col, "_basis", seq_len(ncol(bx$basis)))
    additive_blocks[[col]] <- bx$basis
    colnames(additive_blocks[[col]]) <- feature_names
    surface_infos[[col]] <- c(
      list(
        x_name = col,
        linear_cols = linear_cols,
        centered = TRUE,
        additive_feature_names = feature_names,
        linear_feature_names = linear_cols
      ),
      bx[c(
        "center", "knots", "Boundary.knots",
        "degree", "M_spline_x", "x_range", "x_mean"
      )]
    )
  }
  linear_mat <- if (length(linear_cols)) {
    as.matrix(subject_df[, linear_cols, drop = FALSE])
  } else {
    matrix(numeric(0), nrow = I, ncol = 0)
  }
  fixed_features <- if (length(cols)) {
    do.call(cbind, c(list(linear_mat), additive_blocks))
  } else {
    matrix(numeric(0), nrow = I, ncol = 0)
  }
  storage.mode(fixed_features) <- "double"
  Q_by_id <- cbind(`(Intercept)` = rep(1, I), fixed_features)
  X_by_id <- if (length(cols)) {
    as.matrix(subject_df[, cols, drop = FALSE])
  } else {
    matrix(numeric(0), nrow = I, ncol = 0)
  }
  list(
    Q_by_id = Q_by_id,
    X_by_id = X_by_id,
    fixed_features_by_id = fixed_features,
    fixed_feature_names = colnames(fixed_features),
    covariate_names = colnames(X_by_id),
    covariate_model = spec$covariate_model,
    additive_cols = additive_cols,
    linear_cols = linear_cols,
    M_spline_x_by_additive = sizes,
    surface_infos = surface_infos,
    surface_info = if (length(additive_cols) == 1) {
      surface_infos[[1]]
    } else {
      NULL
    }
  )
}

#' Prepare Inputs for Variational Logistic FoSR
#'
#' Constructs time and covariate bases and ridge-based initial coefficients
#' for \code{fosr()}.
#'
#' @inheritParams fosr
#' @param rows Subject row locations returned by \code{data_clean()}.
#' @return A list containing \code{B_full}, subject covariate designs,
#' spline knots and boundary knots, initial variational parameters
#' (\code{zeta}), initial coefficients (\code{alpha_init}), and
#' information about the initialization subsample.
bfpca_cov_argPreparation <- function(
    data,
    M_spline,
    formula,
    rows,
    max_init_rows = 1e6,
    covariate_model = c("linear", "additive", "mixed"),
    M_spline_x = 5,
    additive_cols = NULL,
    linear_cols = NULL,
    lambda_init = 1e-4
) {
  covariate_model <- match.arg(covariate_model)
  fml <- if (inherits(formula, "formula")) {
    formula
  } else {
    stats::as.formula(formula)
  }
  terms_obj <- stats::terms(fml)
  response_var <- all.vars(fml)[1]
  N <- nrow(data)
  I <- nrow(rows)
  B_full_obj <- splines::bs(
    data$index, df = M_spline, degree = 3, intercept = TRUE
  )
  B_full <- unclass(B_full_obj)
  storage.mode(B_full) <- "double"
  subject_df <- data[rows$first_row, , drop = FALSE]
  rhs_vars <- all.vars(stats::delete.response(terms_obj))
  term_labels <- attr(terms_obj, "term.labels")
  if (length(term_labels) != length(rhs_vars) ||
      !setequal(term_labels, rhs_vars)) {
    stop("Use a simple formula, such as value ~ age + income + BMI.")
  }
  features <- prepare_subject_features(
    subject_df, rhs_vars, covariate_model,
    M_spline_x, additive_cols, linear_cols
  )
  Q_by_id <- features$Q_by_id
  p_fixed <- ncol(Q_by_id) - 1

  # Use a subset of subjects for large-data initialization.
  subject_n <- rows$last_row - rows$first_row + 1
  if (N <= max_init_rows) {
    init_subject_idx <- seq_len(I)
  } else {
    n_keep <- max(1, min(I, floor(max_init_rows / mean(subject_n))))
    init_subject_idx <- unique(
      as.integer(round(seq(1, I, length.out = n_keep)))
    )
  }
  init_lengths <- subject_n[init_subject_idx]
  init_rows <- unlist(
    Map(
      seq.int,
      rows$first_row[init_subject_idx],
      rows$last_row[init_subject_idx]
    ),
    use.names = FALSE
  )
  B_init_sp <- Matrix::Matrix(
    B_full[init_rows, , drop = FALSE], sparse = TRUE
  )
  if (p_fixed > 0) {
    subject_rep <- rep(init_subject_idx, times = init_lengths)
    F_init <- features$fixed_features_by_id[subject_rep, , drop = FALSE]
    fixed_blocks <- lapply(seq_len(p_fixed), function(j) {
      Matrix::Diagonal(x = as.numeric(F_init[, j])) %*% B_init_sp
    })
    design_init <- do.call(cbind, c(list(B_init_sp), fixed_blocks))
  } else {
    design_init <- B_init_sp
  }
  fit <- glmnet::glmnet(
    x = Matrix::drop0(design_init),
    y = data[[response_var]][init_rows],
    family = "binomial",
    alpha = 0,
    lambda = lambda_init,
    standardize = FALSE,
    intercept = FALSE,
    thresh = 1e-6,
    maxit = 200
  )
  coefs_ridge <- drop(
    as.matrix(stats::coef(fit, s = lambda_init))[-1, , drop = FALSE]
  )
  Q <- ncol(Q_by_id)
  alpha_mat_init <- matrix(
    coefs_ridge, nrow = M_spline, ncol = Q
  )
  colnames(alpha_mat_init) <- colnames(Q_by_id)
  fixed_init <- if (p_fixed > 0) {
    split(
      coefs_ridge[-seq_len(M_spline)],
      rep(seq_len(p_fixed), each = M_spline)
    )
  } else {
    list()
  }

  list(
    B_full = B_full,
    Q_by_id = Q_by_id,
    X_by_id = features$X_by_id,
    fixed_features_by_id = features$fixed_features_by_id,
    fixed_feature_names = features$fixed_feature_names,
    covariate_names = features$covariate_names,
    covariate_model = covariate_model,
    additive_cols = features$additive_cols,
    linear_cols = features$linear_cols,
    surface_info = features$surface_info,
    surface_infos = features$surface_infos,
    M_spline_x_by_additive = features$M_spline_x_by_additive,
    knots = attr(B_full_obj, "knots"),
    Boundary.knots = attr(B_full_obj, "Boundary.knots"),
    zeta = rnorm(N) * 0.5,
    alpha_init = list(
      baseline = coefs_ridge[seq_len(M_spline)],
      covariates = fixed_init,
      alpha_mat = alpha_mat_init
    ),
    init_n_subjects = length(init_subject_idx),
    init_n_rows = length(init_rows)
  )
}

#' Optimize the Variational Logistic FoSR Model
#'
#' Runs the variational EM updates for fixed and random functional effects.
#'
#' @inheritParams fosr
#' @param df Long-format input data.
#' @param B_full Observation-level time B-spline design matrix.
#' @param Q_by_id Subject-level fixed-effect design matrix.
#' @param zeta_init Initial observation-level variational parameters.
#' @param alpha_init Initial spline coefficients from
#'   \code{bfpca_cov_argPreparation()}.
#' @param rows Subject row locations returned by \code{data_clean()}.
#' @param Y_vec Optional response vector; defaults to \code{df$value}.
#' @return A list with \code{alpha_mat} (fixed-effect coefficients),
#'   \code{psi_coefs} (functional random-effect coefficients),
#'   \code{scores} (subject scores), and \code{error} (EM convergence values).
bfpca_cov_optimize <- function(
    df,
    B_full,
    Q_by_id,
    M_spline,
    zeta_init,
    alpha_init,
    rows,
    npc = 1,
    maxiter = 50,
    error_thresh = 1e-3,
    print.iter = FALSE,
    Y_vec = NULL,
    cores = 1L
) {
  if (!is.numeric(cores) || length(cores) != 1L ||
      !is.finite(cores) || cores < 1L || cores != as.integer(cores)) {
    stop("cores must be a positive integer.")
  }
  cores <- as.integer(cores)
  I <- nrow(rows)
  first0 <- as.integer(rows$first_row - 1)
  last0 <- as.integer(rows$last_row - 1)
  if (is.null(Y_vec)) Y_vec <- df$value
  Y_vec <- as.numeric(Y_vec)
  zeta <- pmax(as.numeric(zeta_init), 1e-3)
  Q <- ncol(Q_by_id)
  M <- M_spline
  alpha_mat <- alpha_init$alpha_mat
  set.seed(2025)
  psi_coefs <- matrix(
    rnorm(M * npc, sd = 0.1), nrow = M, ncol = npc
  )
  scores <- matrix(0, nrow = I, ncol = npc)
  error_vec <- rep(NA_real_, maxiter)
  prev_alpha <- as.numeric(alpha_mat)
  prev_psi <- psi_coefs
  iter_done <- maxiter

  # Constant within one FoSR fit.
  BtYhalf <- subject_BtYhalf_cpp(
    B_full, Y_vec, first0, last0, cores = cores
  )

  for (iter in seq_len(maxiter)) {
    if (print.iter) cat("FoSR iter", iter, "...\n")

    em_out <- em_iteration_structured_cpp(
      B_full = B_full,
      Y = Y_vec,
      zeta_in = zeta,
      Q_by_id = Q_by_id,
      alpha_mat = alpha_mat,
      psi_coefs = psi_coefs,
      BtYhalf = BtYhalf,
      first0 = first0,
      last0 = last0,
      cores = cores
    )

    scores <- as.matrix(em_out$scores)
    zeta <- pmax(as.numeric(em_out$zeta), 1e-3)

    coef_mat <- matrix(
      as.numeric(em_out$phi_vec),
      nrow = M,
      ncol = Q + npc
    )

    alpha_mat <- coef_mat[, seq_len(Q), drop = FALSE]
    colnames(alpha_mat) <- colnames(Q_by_id)
    psi_coefs <- coef_mat[, Q + seq_len(npc), drop = FALSE]

    alpha_coefs <- as.numeric(alpha_mat)

    if (iter > 1) {
      da <- sum((alpha_coefs - prev_alpha)^2, na.rm = TRUE)
      db <- sum((psi_coefs - prev_psi)^2, na.rm = TRUE)
      error_vec[iter] <- da + db

      if (print.iter) cat("  error =", error_vec[iter], "\n")

      if (error_vec[iter] < error_thresh) {
        iter_done <- iter
        break
      }
    } else {
      error_vec[iter] <- Inf
    }

    prev_alpha <- alpha_coefs
    prev_psi <- psi_coefs
  }

  list(
    alpha_mat = alpha_mat,
    psi_coefs = psi_coefs,
    scores = scores,
    error = error_vec[seq_len(iter_done)]
  )
}

#' Variational Logistic Function-on-Scalar Regression
#'
#' Estimates time-varying covariate effects and subject-specific functional
#' random effects for binomial curves by variational EM. Covariate effects
#' can be linear, additive, or mixed.
#'
#' @param data Long-format data frame containing \code{id}, \code{index},
#' \code{value}, and covariates appearing in \code{formula}.
#' @param formula Formula for main-effect scalar covariates, such as
#' \code{value ~ age + BMI}. Interactions are not supported.
#' @param M_spline Number of cubic B-spline basis functions for time.
#' @param npc Number of functional principal components.
#' @param maxiter Maximum number of variational EM iterations.
#' @param error_thresh Convergence tolerance for variational EM.
#' @param print.iter Logical; print EM iteration information.
#' @param seed Random seed for initialization.
#' @param row_obj Optional subject row locations returned by
#' \code{data_clean()}.
#' @param max_init_rows Maximum number of observations used for ridge
#' initialization.
#' @param covariate_model One of \code{"linear"}, \code{"additive"}, or
#' \code{"mixed"}.
#' @param M_spline_x Number of B-spline basis functions for additive
#' covariate effects; a scalar or one value per additive covariate.
#' @param additive_cols Names of covariates modeled with additive effects.
#' @param linear_cols Names of covariates modeled with linear effects.
#' @param lambda_init Ridge penalty used to initialize the fixed effects.
#' @param cores Number of C++ threads when OpenMP is supported; otherwise
#' the kernels run serially.
#' @return An object of class \code{fpca} containing:
#'   \item{fpca_type}{Variational EM fit identifier.}
#'   \item{t_vec, knots, Boundary.knots}{Time grid and spline information.}
#'   \item{alpha, a0_coefs, aJ_coefs}{Population mean and fixed-effect
#'   spline coefficients.}
#'   \item{efunctions, evalues, npc, scores}{Estimated principal
#'   components and subject scores.}
#'   \item{subject_coefs}{Subject-specific time-spline coefficients used
#'   for registration and prediction.}
#'   \item{Yhat, Y, family, error}{Fitted linear predictors, input data,
#'   family, and EM convergence history.}
#'   \item{surface_infos, surface_coefs_list, linear_coefs, linear_means}{
#'   Information used to evaluate additive and linear covariate effects.}
#'   Other elements retain the covariate specification and initialization
#'   information.
#' @author Linlin Wu \email{linlin.wu@@emory.edu},
#'   Jeff Goldsmith \email{ajg2202@@cumc.columbia.edu},
#'   Julia Wrobel \email{julia.wrobel@@emory.edu}
#' @export
fosr <- function(
    data,
    formula,
    M_spline,
    npc = 1,
    maxiter = 50,
    error_thresh = 1e-4,
    print.iter = FALSE,
    seed = 1988,
    row_obj = NULL,
    max_init_rows = 1e6,
    covariate_model = c("linear", "additive", "mixed"),
    M_spline_x = 5,
    additive_cols = NULL,
    linear_cols = NULL,
    lambda_init = 1e-4,
    cores = 1L
) {
  set.seed(seed)
  covariate_model <- match.arg(covariate_model)
  fml <- if (inherits(formula, "formula")) {
    formula
  } else {
    stats::as.formula(formula)
  }

  resp_name <- all.vars(fml)[1]
  df <- as.data.frame(data)

  # Use value as the response column name.
  if (!("value" %in% names(df)) && !(resp_name %in% names(df))) {
    stop("Response variable not found in data.")
  }
  if (!("value" %in% names(df))) {
    names(df)[match(resp_name, names(df))] <- "value"
  }
  fml[[2]] <- as.name("value")

  if (is.null(row_obj)) {
    dat <- registr::data_clean(df)
    df <- dat$Y
    rows <- dat$Y_rows
  } else {
    rows <- as.data.frame(row_obj)
  }

  first0 <- as.integer(rows$first_row - 1)
  last0 <- as.integer(rows$last_row - 1)

  prep <- bfpca_cov_argPreparation(
    data = df,
    M_spline = M_spline,
    formula = fml,
    rows = rows,
    max_init_rows = max_init_rows,
    covariate_model = covariate_model,
    M_spline_x = M_spline_x,
    additive_cols = additive_cols,
    linear_cols = linear_cols,
    lambda_init = lambda_init
  )

  opt <- bfpca_cov_optimize(
    df = df,
    B_full = prep$B_full,
    Q_by_id = prep$Q_by_id,
    M_spline = M_spline,
    zeta_init = prep$zeta,
    alpha_init = prep$alpha_init,
    rows = rows,
    npc = npc,
    maxiter = maxiter,
    error_thresh = error_thresh,
    print.iter = print.iter,
    Y_vec = df$value,
    cores = cores
  )

  alpha_mat <- opt$alpha_mat
  psi_coefs <- opt$psi_coefs
  scores_unorth <- opt$scores

  M <- M_spline
  Q <- ncol(prep$Q_by_id)
  t_vec <- seq(min(df$index), max(df$index), length.out = 100)

  Theta_phi_mean <- splines::bs(
    t_vec,
    knots = prep$knots,
    degree = 3,
    intercept = TRUE,
    Boundary.knots = prep$Boundary.knots
  )

  a0 <- alpha_mat[, 1]
  aJ <- alpha_mat[, -1, drop = FALSE]
  colnames(aJ) <- prep$fixed_feature_names

  qbar <- matrix(colMeans(prep$Q_by_id), ncol = 1)
  rownames(qbar) <- colnames(prep$Q_by_id)
  alpha_coefs_pop <- as.numeric(alpha_mat %*% qbar)

  # Fixed and random contributions to subject-specific templates.
  subject_coefs <- alpha_mat %*% t(prep$Q_by_id) +
    psi_coefs %*% t(scores_unorth)
  colnames(subject_coefs) <- as.character(rows$id)

  fits <- fitted_values_by_subject_cpp(
    B_full = prep$B_full,
    subject_coefs = subject_coefs,
    first0 = first0,
    last0 = last0,
    cores = cores
  )

  fittedVals <- data.frame(
    id = df$id,
    index = df$index,
    value = as.numeric(fits)
  )

  # Orthogonalize random-effect functions.
  psi_svd <- svd(Theta_phi_mean %*% psi_coefs)
  efunctions <- psi_svd$u
  evalues <- psi_svd$d^2

  d_diag <- if (length(psi_svd$d) == 1) {
    matrix(psi_svd$d, 1, 1)
  } else {
    diag(psi_svd$d)
  }

  scores_orth <- scores_unorth %*% psi_svd$v %*% d_diag

  xbar <- if (ncol(prep$X_by_id) > 0) {
    xb <- colMeans(prep$X_by_id)
    matrix(xb, ncol = 1, dimnames = list(names(xb), NULL))
  } else {
    matrix(0, nrow = 0, ncol = 1)
  }

  surface_coefs_list <- lapply(
    prep$surface_infos,
    function(info) {
      aJ[, info$additive_feature_names, drop = FALSE]
    }
  )

  surface_coefs <- if (length(surface_coefs_list) == 1) {
    surface_coefs_list[[1]]
  } else if (length(surface_coefs_list)) {
    surface_coefs_list
  } else {
    NULL
  }

  linear_coefs <- if (length(prep$linear_cols)) {
    aJ[, prep$linear_cols, drop = FALSE]
  } else {
    NULL
  }

  linear_means <- if (length(prep$linear_cols)) {
    colMeans(prep$X_by_id[, prep$linear_cols, drop = FALSE])
  } else {
    numeric(0)
  }

  M_spline_x_out <- prep$M_spline_x_by_additive
  if (!length(M_spline_x_out)) {
    M_spline_x_out <- NULL
  } else if (length(M_spline_x_out) == 1) {
    M_spline_x_out <- unname(M_spline_x_out[[1]])
  }

  ret <- list(
    fpca_type = "variationalEM",
    t_vec = t_vec,
    knots = prep$knots,
    Boundary.knots = prep$Boundary.knots,
    alpha = as.numeric(Theta_phi_mean %*% alpha_coefs_pop),
    efunctions = efunctions,
    evalues = evalues,
    npc = npc,
    scores = scores_orth,
    subject_coefs = subject_coefs,
    Yhat = fittedVals,
    Y = df,
    family = "binomial",
    error = opt$error,
    covariate_model = prep$covariate_model,
    additive_cols = prep$additive_cols,
    linear_cols = prep$linear_cols,
    M_spline_x = M_spline_x_out,
    M_spline_x_by_additive = prep$M_spline_x_by_additive,
    a0_coefs = a0,
    aJ_coefs = aJ,
    surface_coefs = surface_coefs,
    surface_info = prep$surface_info,
    surface_coefs_list = surface_coefs_list,
    surface_infos = prep$surface_infos,
    linear_coefs = linear_coefs,
    linear_means = linear_means,

    covariate_names = prep$covariate_names,
    fixed_feature_names = prep$fixed_feature_names,

    alpha_mat = alpha_mat,
    Xbar = xbar,
    Qbar = qbar,
    init_n_subjects = prep$init_n_subjects,
    init_n_rows = prep$init_n_rows
  )

  class(ret) <- "fpca"
  ret
}

#' Predict a Fixed-Effect Surface from a VL-FoSR Fit
#'
#' Evaluates the fitted linear predictor over time and one selected additive
#' covariate while holding the other covariates at supplied values or their
#' sample means.
#'
#' @param object A fitted object returned by \code{fosr()}.
#' @param t Numeric vector of time points; defaults to \code{object$t_vec}.
#' @param x Values of the additive covariate to vary; defaults to a grid
#' over its observed range.
#' @param linear_values Optional named values for the linear covariates.
#' @param additive_values Optional named values for additive covariates
#' held fixed.
#' @param varying_additive Name of the additive covariate to vary; by
#' default, the first additive covariate.
#' @return A data frame with \code{t}, \code{x}, \code{surface_effect},
#' \code{other_additive_adjustment}, \code{linear_adjustment},
#' \code{eta_fixed} (on the logit scale), and \code{x_name}.
predict_fosr_surface <- function(
    object,
    t = object$t_vec,
    x = NULL,
    linear_values = NULL,
    additive_values = NULL,
    varying_additive = NULL
) {
  infos <- object$surface_infos

  if (!length(infos)) {
    stop("The model has no additive covariates.")
  }
  if (is.null(varying_additive)) {
    varying_additive <- names(infos)[1]
  }

  chosen <- infos[[varying_additive]]
  if (is.null(chosen)) {
    stop("Unknown varying_additive.")
  }

  if (is.null(x)) {
    x <- seq(
      chosen$x_range[1],
      chosen$x_range[2],
      length.out = 100
    )
  }
  x <- as.numeric(x)

  # Hold other additive covariates at their sample means.
  held <- vapply(
    infos, function(info) info$x_mean, numeric(1)
  )

  if (length(additive_values)) {
    if (is.null(names(additive_values)) ||
        !all(names(additive_values) %in% names(infos))) {
      stop("additive_values must be named by additive covariates.")
    }

    held[names(additive_values)] <- as.numeric(additive_values)
  }
  Bt <- splines::bs(
    t,
    knots = object$knots,
    degree = 3,
    intercept = TRUE,
    Boundary.knots = object$Boundary.knots
  )
  Bx <- eval_x_spline_basis(x, chosen)
  chosen_effect <- Bt %*%
    object$surface_coefs_list[[varying_additive]] %*%
    t(Bx)
  other_effect_t <- rep(0, length(t))
  for (name in setdiff(names(infos), varying_additive)) {
    bx_held <- eval_x_spline_basis(
      held[[name]], infos[[name]]
    )

    other_effect_t <- other_effect_t + as.numeric(
      Bt %*% object$surface_coefs_list[[name]] %*% t(bx_held)
    )
  }
  # Hold linear covariates at their sample means.
  linear_effect_t <- rep(0, length(t))
  if (length(object$linear_cols)) {
    if (is.null(linear_values)) {
      linear_values <- object$linear_means
    }
    if (is.null(names(linear_values))) {
      names(linear_values) <- object$linear_cols
    }
    if (!all(object$linear_cols %in% names(linear_values))) {
      stop("Missing linear_values for some linear covariates.")
    }
    linear_effect_t <- as.numeric(
      Bt %*% object$linear_coefs %*%
        matrix(
          as.numeric(linear_values[object$linear_cols]),
          ncol = 1
        )
    )
  }
  baseline_t <- as.numeric(Bt %*% object$a0_coefs) +
    other_effect_t + linear_effect_t
  total <- sweep(chosen_effect, 1, baseline_t, "+")
  out <- expand.grid(t = t, x = x)
  out$surface_effect <- as.numeric(chosen_effect)
  out$other_additive_adjustment <- rep(
    other_effect_t, times = length(x)
  )
  out$linear_adjustment <- rep(
    linear_effect_t, times = length(x)
  )
  out$eta_fixed <- as.numeric(total)
  out$x_name <- varying_additive
  out
}
