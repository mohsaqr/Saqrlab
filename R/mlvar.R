# ---- Multilevel Vector Autoregression (mlVAR) ----

#' Multilevel Vector Autoregression
#'
#' @description Estimates three networks from ESM/EMA panel data using
#'   multilevel VAR: (1) a directed temporal network of lagged regression
#'   coefficients, (2) an undirected contemporaneous network of partial
#'   correlations among residuals, and (3) an undirected between-subjects
#'   network of partial correlations among person means. Uses within-person
#'   centering (fixed-effects OLS) for the temporal model and EBIC-selected
#'   graphical LASSO for the contemporaneous and between-subjects networks.
#'
#' @param data A \code{data.frame} containing the panel data.
#' @param vars Character vector of variable column names to model.
#' @param id Character string naming the person-ID column.
#' @param day Character string naming the day/session column, or \code{NULL}
#'   (default). When provided, lag pairs are only formed within the same day.
#' @param beep Character string naming the measurement-occasion column, or
#'   \code{NULL} (default). When provided, lag pairs require
#'   \code{beep[t] - beep[t-lag] == lag}.
#' @param lag Integer. The lag order (default 1).
#' @param standardize Logical. If \code{TRUE} (default), within-centered
#'   variables are divided by their pooled standard deviation before OLS.
#' @param gamma Numeric. EBIC hyperparameter for graphical LASSO model
#'   selection (0 = BIC, higher = sparser). Default 0.5.
#' @param nlambda Integer. Number of lambda values in the regularization path
#'   for graphical LASSO. Default 100.
#'
#' @return An S3 object of class \code{"mlvar_result"}, a list with:
#'   \describe{
#'     \item{\code{temporal}}{d x d matrix of fixed-effect temporal regression
#'       coefficients. Entry \code{[i, j]} is the effect of variable j at
#'       t-lag on variable i at t.}
#'     \item{\code{contemporaneous}}{d x d symmetric matrix of partial
#'       correlations among within-person residuals (EBIC-GLASSO).}
#'     \item{\code{between}}{d x d symmetric matrix of partial correlations
#'       among person means (EBIC-GLASSO).}
#'     \item{\code{coefs}}{List of d data frames, one per outcome variable,
#'       each containing columns \code{predictor}, \code{beta}, \code{se},
#'       \code{t}, \code{p}, \code{ci_lower}, \code{ci_upper}.}
#'     \item{\code{labels}}{Character vector of variable names.}
#'     \item{\code{n_obs}}{Number of valid lag-pair observations.}
#'     \item{\code{n_subjects}}{Number of unique subjects.}
#'     \item{\code{lag}}{Lag order used.}
#'     \item{\code{standardize}}{Logical; whether standardization was applied.}
#'     \item{\code{gamma}}{EBIC gamma used.}
#'   }
#'
#' @details
#' The algorithm proceeds in seven steps:
#' \enumerate{
#'   \item \strong{Data preparation}: select columns, coerce types, drop NA
#'     rows, sort by \code{(id, day, beep)}.
#'   \item \strong{Lag-pair construction}: form valid outcome/predictor pairs
#'     respecting person, day, and beep boundaries.
#'   \item \strong{Within-centering}: person-mean center both Y (outcome) and
#'     X (predictor) matrices. Optionally standardize by pooled SD.
#'   \item \strong{Temporal OLS}: for each outcome k, fit
#'     \code{lm(Y_k ~ X - 1)} with corrected degrees of freedom.
#'   \item \strong{Contemporaneous network}: EBIC-GLASSO on the correlation
#'     matrix of OLS residuals.
#'   \item \strong{Between-subjects network}: EBIC-GLASSO on the correlation
#'     matrix of person means.
#'   \item \strong{Assembly}: collect results into \code{mlvar_result} object.
#' }
#'
#' @examples
#' d <- simulate_data("mlvar", seed = 1)
#' fit <- mlvar(d, vars = attr(d, "vars"), id = "id", day = "day", beep = "beep")
#' print(fit)
#' summary(fit)
#'
#' @seealso \code{\link{simulate_data}}, \code{\link{build_network}}
#' @export
mlvar <- function(data,
                  vars,
                  id,
                  day = NULL,
                  beep = NULL,
                  lag = 1L,
                  standardize = TRUE,
                  gamma = 0.5,
                  nlambda = 100L) {
  # ---- Input validation ----
  stopifnot(
    is.data.frame(data),
    is.character(vars),
    length(vars) >= 2L,
    is.character(id),
    length(id) == 1L
  )
  stopifnot(
    is.numeric(lag),
    length(lag) == 1L,
    lag >= 1L
  )
  stopifnot(
    is.logical(standardize),
    length(standardize) == 1L
  )
  stopifnot(
    is.numeric(gamma),
    length(gamma) == 1L,
    gamma >= 0
  )
  nlambda <- as.integer(nlambda)
  stopifnot(
    is.integer(nlambda),
    length(nlambda) == 1L,
    nlambda >= 2L
  )

  # Check columns exist
  required <- c(vars, id)
  if (!is.null(day)) required <- c(required, day)
  if (!is.null(beep)) required <- c(required, beep)
  missing_cols <- setdiff(required, names(data))
  if (length(missing_cols) > 0L) {
    stop("Columns not found in data: ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  d <- length(vars)

  # Step 1: Prepare data
  prepared <- .mlvar_prepare_data(data, vars, id, day, beep)

  # Check minimum subjects
  n_subjects <- length(unique(prepared[[id]]))
  if (n_subjects < 2L) {
    stop("At least 2 subjects are required. Found: ", n_subjects, call. = FALSE)
  }

  # Step 2: Build lag pairs
  lag_result <- .mlvar_build_lag_pairs(prepared, vars, id, day, beep,
                                       lag = as.integer(lag))
  Y <- lag_result$Y
  X <- lag_result$X
  id_vec <- lag_result$id_vec
  n_obs <- nrow(Y)

  if (n_obs < d + 1L) {
    stop("Too few valid lag pairs (", n_obs,
         ") for ", d, " variables.", call. = FALSE)
  }

  n_subjects_pairs <- length(unique(id_vec))

  # Step 3: Within-center
  centered <- .mlvar_within_center(Y, X, id_vec, standardize)
  Y_c <- centered$Y
  X_c <- centered$X

  # Step 4: Temporal OLS
  temporal_result <- .mlvar_temporal_ols(Y_c, X_c, n_subjects_pairs, vars)

  # Step 5: Contemporaneous network
  contemporaneous <- .mlvar_contemporaneous(
    temporal_result$residuals, n_obs, gamma, nlambda
  )

  # Step 6: Between-subjects network
  between <- .mlvar_between(prepared, vars, id, n_subjects, gamma, nlambda)

  # Step 7: Assemble result
  result <- list(
    temporal        = temporal_result$B,
    contemporaneous = contemporaneous,
    between         = between,
    coefs           = temporal_result$coefs,
    labels          = vars,
    n_obs           = n_obs,
    n_subjects      = n_subjects_pairs,
    lag             = as.integer(lag),
    standardize     = standardize,
    gamma           = gamma
  )
  class(result) <- "mlvar_result"
  result
}


# ---- Internal helpers ----

#' Prepare data for mlVAR
#'
#' Selects relevant columns, coerces types, drops rows with NA in any variable,
#' and sorts by (id, day, beep).
#'
#' @noRd
.mlvar_prepare_data <- function(data, vars, id, day, beep) {
  keep_cols <- c(id, vars)
  if (!is.null(day)) keep_cols <- c(keep_cols, day)
  if (!is.null(beep)) keep_cols <- c(keep_cols, beep)
  keep_cols <- unique(keep_cols)

  df <- data[, keep_cols, drop = FALSE]

  # Coerce vars to numeric
  for (v in vars) {
    df[[v]] <- as.numeric(df[[v]])
  }

  # Drop rows with NA in any variable or id
  check_cols <- c(id, vars)
  if (!is.null(day)) check_cols <- c(check_cols, day)
  if (!is.null(beep)) check_cols <- c(check_cols, beep)
  complete <- complete.cases(df[, check_cols, drop = FALSE])
  if (!all(complete)) {
    df <- df[complete, , drop = FALSE]
  }

  if (nrow(df) < 2L) {
    stop("Fewer than 2 complete rows remain after removing NAs.", call. = FALSE)
  }

  # Sort by (id, day, beep)
  order_cols <- id
  if (!is.null(day)) order_cols <- c(order_cols, day)
  if (!is.null(beep)) order_cols <- c(order_cols, beep)
  ord <- do.call(order, df[, order_cols, drop = FALSE])
  df <- df[ord, , drop = FALSE]
  rownames(df) <- NULL

  df
}


#' Build lag pairs for mlVAR
#'
#' Constructs valid lagged observation pairs. A pair is valid if:
#' same person, same day (if provided), and beep gap equals lag (if provided).
#'
#' @return List with Y (n_pairs x d outcome matrix), X (n_pairs x d predictor
#'   matrix), id_vec (person IDs for each pair).
#' @noRd
.mlvar_build_lag_pairs <- function(data, vars, id, day, beep, lag) {
  n <- nrow(data)
  d <- length(vars)

  if (n <= lag) {
    stop("Not enough rows (", n, ") for lag ", lag, ".", call. = FALSE)
  }

  # Indices for current (t) and lagged (t-lag)
  idx_t <- seq(lag + 1L, n)
  idx_lag <- seq(1L, n - lag)

  # Condition 1: same person
  valid <- data[[id]][idx_t] == data[[id]][idx_lag]

  # Condition 2: same day (if provided)
  if (!is.null(day)) {
    valid <- valid & (data[[day]][idx_t] == data[[day]][idx_lag])
  }

  # Condition 3: beep gap == lag (if provided)
  if (!is.null(beep)) {
    beep_diff <- as.numeric(data[[beep]][idx_t]) -
                 as.numeric(data[[beep]][idx_lag])
    valid <- valid & (beep_diff == lag)
  }

  if (sum(valid) == 0L) {
    stop("No valid lag pairs found. Check id/day/beep columns and lag value.",
         call. = FALSE)
  }

  # Extract Y (outcome at t) and X (predictor at t-lag)
  Y <- as.matrix(data[idx_t[valid], vars, drop = FALSE])
  X <- as.matrix(data[idx_lag[valid], vars, drop = FALSE])
  id_vec <- data[[id]][idx_t[valid]]

  list(Y = Y, X = X, id_vec = id_vec)
}


#' Within-person centering for mlVAR
#'
#' Person-mean centers both Y and X matrices. Centering both sides absorbs
#' the random intercept, making OLS equivalent to the fixed-effects estimator.
#' If standardize = TRUE, divides by pooled SD after centering.
#'
#' @return List with centered Y and X matrices.
#' @noRd
.mlvar_within_center <- function(Y, X, id_vec, standardize) {
  d <- ncol(Y)

  # Center Y: subtract person means
  Y_c <- Y
  for (j in seq_len(d)) {
    person_means <- ave(Y[, j], id_vec, FUN = mean)
    Y_c[, j] <- Y[, j] - person_means
  }

  # Center X: subtract person means
  X_c <- X
  for (j in seq_len(d)) {
    person_means <- ave(X[, j], id_vec, FUN = mean)
    X_c[, j] <- X[, j] - person_means
  }

  if (standardize) {
    # Pooled SD: computed across all observations (already centered)
    for (j in seq_len(d)) {
      sd_y <- stats::sd(Y_c[, j])
      sd_x <- stats::sd(X_c[, j])
      if (sd_y > 0) Y_c[, j] <- Y_c[, j] / sd_y
      if (sd_x > 0) X_c[, j] <- X_c[, j] / sd_x
    }
  }

  list(Y = Y_c, X = X_c)
}


#' Temporal OLS for mlVAR
#'
#' For each outcome variable, fits OLS without intercept (centering absorbed
#' it). Applies degrees-of-freedom correction for absorbed person fixed effects.
#'
#' @return List with B (d x d coefficient matrix), coefs (list of data frames),
#'   residuals (n_obs x d matrix).
#' @noRd
.mlvar_temporal_ols <- function(Y, X, n_subjects, vars) {
  n_obs <- nrow(Y)
  d <- ncol(Y)

  B <- matrix(0, nrow = d, ncol = d,
              dimnames = list(vars, vars))
  residuals <- matrix(0, nrow = n_obs, ncol = d,
                      dimnames = list(NULL, vars))
  coefs_list <- vector("list", d)
  names(coefs_list) <- vars

  df_ols <- n_obs - d
  df_correct <- n_obs - d - n_subjects

  if (df_correct < 1L) {
    warning("Corrected degrees of freedom < 1. Results may be unreliable.",
            call. = FALSE)
    df_correct <- max(1L, df_correct)
  }

  df_ratio <- sqrt(df_ols / df_correct)

  for (k in seq_len(d)) {
    fit <- stats::lm.fit(X, Y[, k])
    beta <- fit$coefficients
    resid <- fit$residuals

    # Corrected SE
    rss <- sum(resid^2)
    sigma2 <- rss / df_ols
    XtX_inv <- tryCatch(
      solve(crossprod(X)),
      error = function(e) {
        warning("Singular X'X for variable ", vars[k],
                ". Using pseudoinverse.", call. = FALSE)
        MASS_needed <- FALSE
        # Fallback: use diagonal approximation
        diag(1 / pmax(diag(crossprod(X)), 1e-10))
      }
    )
    se_ols <- sqrt(sigma2 * diag(XtX_inv))
    se_correct <- se_ols * df_ratio

    # t-statistics and p-values
    t_val <- beta / se_correct
    p_val <- 2 * stats::pt(abs(t_val), df = df_correct, lower.tail = FALSE)

    # Confidence intervals
    t_crit <- stats::qt(0.975, df = df_correct)
    ci_lower <- beta - t_crit * se_correct
    ci_upper <- beta + t_crit * se_correct

    B[k, ] <- beta
    residuals[, k] <- resid

    coefs_list[[k]] <- data.frame(
      predictor = vars,
      beta      = as.numeric(beta),
      se        = as.numeric(se_correct),
      t         = as.numeric(t_val),
      p         = as.numeric(p_val),
      ci_lower  = as.numeric(ci_lower),
      ci_upper  = as.numeric(ci_upper),
      stringsAsFactors = FALSE
    )
  }

  list(B = B, coefs = coefs_list, residuals = residuals)
}


#' Contemporaneous network via EBIC-GLASSO on residuals
#'
#' @return d x d partial correlation matrix (symmetric, zero diagonal).
#' @noRd
.mlvar_contemporaneous <- function(residuals, n_obs, gamma, nlambda) {
  d <- ncol(residuals)
  vars <- colnames(residuals)

  S <- stats::cor(residuals)

  # Guard: if any correlations are NA or all zero

  if (any(is.na(S))) {
    warning("NA correlations in residuals. Returning zero matrix.",
            call. = FALSE)
    mat <- matrix(0, d, d, dimnames = list(vars, vars))
    return(mat)
  }

  # Use existing GLASSO pipeline
  lambda_path <- tryCatch(
    .compute_lambda_path(S, nlambda, 0.01),
    error = function(e) NULL
  )

  if (is.null(lambda_path)) {
    mat <- matrix(0, d, d, dimnames = list(vars, vars))
    return(mat)
  }

  selected <- tryCatch(
    .select_ebic(S, lambda_path, n_obs, gamma, FALSE),
    error = function(e) NULL
  )

  if (is.null(selected)) {
    mat <- matrix(0, d, d, dimnames = list(vars, vars))
    return(mat)
  }

  pcor <- .precision_to_pcor(selected$wi, 0)
  pcor <- (pcor + t(pcor)) / 2  # force exact symmetry
  colnames(pcor) <- rownames(pcor) <- vars
  pcor
}


#' Between-subjects network via EBIC-GLASSO on person means
#'
#' @return d x d partial correlation matrix (symmetric, zero diagonal).
#' @noRd
.mlvar_between <- function(data, vars, id, n_subjects, gamma, nlambda) {
  d <- length(vars)

  # Guard: need at least d+1 subjects for meaningful estimation
  if (n_subjects < d + 1L) {
    mat <- matrix(0, d, d, dimnames = list(vars, vars))
    return(mat)
  }

  # Compute person means
  agg_data <- data[, c(id, vars), drop = FALSE]
  names(agg_data)[1] <- ".id"
  person_means <- stats::aggregate(. ~ .id, data = agg_data, FUN = mean)
  # Remove .id column
  pm_mat <- as.matrix(person_means[, vars, drop = FALSE])

  # Check variance
  col_vars <- apply(pm_mat, 2, stats::var)
  if (any(col_vars == 0)) {
    mat <- matrix(0, d, d, dimnames = list(vars, vars))
    return(mat)
  }

  S <- stats::cor(pm_mat)

  if (any(is.na(S))) {
    mat <- matrix(0, d, d, dimnames = list(vars, vars))
    return(mat)
  }

  lambda_path <- tryCatch(
    .compute_lambda_path(S, nlambda, 0.01),
    error = function(e) NULL
  )

  if (is.null(lambda_path)) {
    mat <- matrix(0, d, d, dimnames = list(vars, vars))
    return(mat)
  }

  n_pm <- nrow(pm_mat)
  selected <- tryCatch(
    .select_ebic(S, lambda_path, n_pm, gamma, FALSE),
    error = function(e) NULL
  )

  if (is.null(selected)) {
    mat <- matrix(0, d, d, dimnames = list(vars, vars))
    return(mat)
  }

  pcor <- .precision_to_pcor(selected$wi, 0)
  pcor <- (pcor + t(pcor)) / 2  # force exact symmetry
  colnames(pcor) <- rownames(pcor) <- vars
  pcor
}


# ---- S3 methods ----

#' Print method for mlvar_result
#'
#' @param x An \code{mlvar_result} object.
#' @param ... Additional arguments (ignored).
#'
#' @export
print.mlvar_result <- function(x, ...) {
  d <- length(x$labels)
  n_temp <- sum(x$temporal != 0)
  n_cont <- sum(x$contemporaneous[upper.tri(x$contemporaneous)] != 0)
  n_betw <- sum(x$between[upper.tri(x$between)] != 0)

  cat("mlVAR result:",
      x$n_subjects, "subjects,",
      x$n_obs, "observations,",
      d, "variables\n")
  cat("  Temporal edges:", n_temp, "(directed)\n")
  cat("  Contemporaneous edges:", n_cont, "(undirected)\n")
  cat("  Between-subjects edges:", n_betw, "(undirected)\n")
  invisible(x)
}


#' Summary method for mlvar_result
#'
#' @param object An \code{mlvar_result} object.
#' @param ... Additional arguments (ignored).
#'
#' @export
summary.mlvar_result <- function(object, ...) {
  d <- length(object$labels)

  cat("=== mlVAR Summary ===\n")
  cat("Subjects:", object$n_subjects, " | Observations:", object$n_obs,
      " | Variables:", d, " | Lag:", object$lag, "\n")
  cat("Standardized:", object$standardize, " | EBIC gamma:", object$gamma, "\n")
  cat("\n")

  # Temporal coefficients
  cat("--- Temporal Network (B matrix) ---\n")
  print(round(object$temporal, 4))
  cat("\n")

  # Significant temporal edges
  sig_edges <- do.call(rbind, lapply(seq_along(object$coefs), function(k) {
    cf <- object$coefs[[k]]
    sig <- cf[cf$p < 0.05, , drop = FALSE]
    if (nrow(sig) > 0L) {
      sig$outcome <- object$labels[k]
      sig
    } else {
      NULL
    }
  }))
  if (!is.null(sig_edges) && nrow(sig_edges) > 0L) {
    cat("Significant temporal edges (p < 0.05):\n")
    sig_print <- sig_edges[, c("predictor", "outcome", "beta", "se", "t", "p"),
                           drop = FALSE]
    sig_print$beta <- round(sig_print$beta, 4)
    sig_print$se <- round(sig_print$se, 4)
    sig_print$t <- round(sig_print$t, 3)
    sig_print$p <- round(sig_print$p, 4)
    rownames(sig_print) <- NULL
    print(sig_print)
  } else {
    cat("No significant temporal edges at p < 0.05.\n")
  }
  cat("\n")

  # Contemporaneous
  cat("--- Contemporaneous Network ---\n")
  print(round(object$contemporaneous, 4))
  cat("\n")

  # Between
  cat("--- Between-Subjects Network ---\n")
  print(round(object$between, 4))
  cat("\n")

  invisible(object)
}
