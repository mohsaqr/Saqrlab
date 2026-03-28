# ===========================================================================
# simulate_longitudinal — panel / ESM data for VAR, mlVAR, growth models
# ===========================================================================

#' Simulate Longitudinal Panel Data
#'
#' @description Generate synthetic multilevel time-series data suitable for
#'   testing longitudinal models: multilevel VAR (mlVAR), random-intercept
#'   cross-lagged panel models (RI-CLPM), latent growth curves, and similar.
#'
#'   The data-generating process is a VAR(1) model with three layers:
#'   \describe{
#'     \item{Temporal}{Autoregressive and cross-lagged effects via matrix
#'       \code{B}. At each time point, \code{y(t) = mu_i + B (y(t-1) - mu_i)
#'       + innovation}.}
#'     \item{Contemporaneous}{Within-time correlations among innovations, given
#'       by the Cholesky of \code{contemporaneous}.}
#'     \item{Between-person}{Person-specific means \code{mu_i} drawn from
#'       \code{N(grand_means, between)}.}
#'   }
#'
#' @param n Integer. Number of subjects (persons). Default 50.
#' @param tp Integer. Number of time points per subject. Default 50.
#' @param vars Integer or character vector. If integer, the number of
#'   variables (named \code{V1}…\code{Vp}). If character, used as variable
#'   names directly. Default 4.
#' @param temporal Numeric matrix (\code{p x p}). The VAR(1) coefficient matrix
#'   \code{B}. Element \code{[i,j]} is the effect of variable \code{j} at
#'   \code{t-1} on variable \code{i} at \code{t}. Eigenvalues should have
#'   modulus < 1 for stationarity. Default: diagonal 0.3 (autoregressive only,
#'   no cross-lags).
#' @param contemporaneous Numeric matrix (\code{p x p}). Correlation matrix of
#'   the innovation terms (within-time, within-person). Must be symmetric and
#'   positive-definite. Default: identity (uncorrelated innovations).
#' @param between Numeric matrix (\code{p x p}). Covariance matrix of the
#'   person-specific means. Controls between-person variability. Default:
#'   identity.
#' @param grand_means Numeric vector of length \code{p}. Population-level means
#'   for each variable. Default: all zeros.
#' @param innovation_sd Numeric scalar or vector of length \code{p}. Standard
#'   deviation(s) of the innovation terms (before applying contemporaneous
#'   correlation). Default: 1.
#' @param beeps_per_day Integer or NULL. If non-NULL, adds \code{day} and
#'   \code{beep} columns (ESM/EMA structure). \code{tp} must be divisible by
#'   \code{beeps_per_day}. Default: NULL (single \code{time} column).
#' @param ar_range Numeric vector of length 2. When \code{temporal} is NULL,
#'   diagonal (autoregressive) coefficients are drawn uniformly from this range.
#'   Default: \code{c(0.2, 0.5)}.
#' @param cross_range Numeric vector of length 2. When \code{temporal} is NULL,
#'   off-diagonal (cross-lag) coefficients are drawn uniformly from this range
#'   (with random sign). Default: \code{c(0.05, 0.2)}.
#' @param n_cross Integer. When \code{temporal} is NULL, the number of non-zero
#'   cross-lagged effects to include. Default: \code{min(3, p*(p-1))}.
#' @param complexity Character. Edge-case injection, same semantics as
#'   \code{\link{simulate_data}}. One of \code{"clean"} (default),
#'   \code{"auto"}, or a character vector of specific cases (e.g.
#'   \code{c("na", "outliers")}). Supported cases: \code{"na"},
#'   \code{"outliers"}, \code{"heavy_tailed"}, \code{"heteroscedastic"},
#'   \code{"tiny_n"}.
#' @param seed Integer or NULL. Random seed.
#'
#' @return A \code{\link{saqr_sim}} object with:
#'   \describe{
#'     \item{\code{$data}}{data.frame with columns \code{id}, \code{time}
#'       (or \code{day}/\code{beep}), and \code{V1}…\code{Vp}.}
#'     \item{\code{$params}}{list with \code{temporal} (B matrix),
#'       \code{contemporaneous} (innovation correlation),
#'       \code{between} (person-mean covariance), \code{grand_means},
#'       \code{innovation_sd}, \code{n}, \code{tp}, \code{var_names}.}
#'   }
#'
#' @details
#' The temporal matrix \code{B} must be stationary (all eigenvalues inside the
#' unit circle). If an auto-generated or user-supplied \code{B} violates this,
#' it is rescaled by \code{0.95 / max(Mod(eigen(B)$values))}.
#'
#' For ESM-style data with \code{beeps_per_day}, the first observation of each
#' day is treated as a new starting point (no carry-over from previous day's
#' last beep), matching the assumption in most ESM-VAR software.
#'
#' @examples
#' # Basic: 30 subjects, 60 time points, 3 variables
#' r <- simulate_longitudinal(n = 30, tp = 60, vars = 3, seed = 42)
#' head(r$data)
#' r$params$temporal  # the true VAR(1) matrix
#'
#' # Explicit temporal structure
#' B <- matrix(c(0.4,  0.0, 0.1,
#'               0.2,  0.3, 0.0,
#'               0.0, -0.1, 0.5), nrow = 3, byrow = TRUE)
#' r2 <- simulate_longitudinal(n = 50, tp = 100, vars = 3, temporal = B, seed = 1)
#'
#' # ESM structure with days and beeps
#' r3 <- simulate_longitudinal(n = 40, tp = 70, vars = 4,
#'                              beeps_per_day = 7, seed = 7)
#' head(r3$data)  # has day + beep columns
#'
#' # With edge-case injection
#' r4 <- simulate_longitudinal(n = 30, tp = 50, vars = 3,
#'                              complexity = c("na", "outliers"), seed = 5)
#'
#' @export
simulate_longitudinal <- function(
    n               = 50L,
    tp              = 50L,
    vars            = 4L,
    temporal        = NULL,
    contemporaneous = NULL,
    between         = NULL,
    grand_means     = NULL,
    innovation_sd   = 1,
    beeps_per_day   = NULL,
    ar_range        = c(0.2, 0.5),
    cross_range     = c(0.05, 0.2),
    n_cross         = NULL,
    complexity      = "clean",
    seed            = NULL
) {
  # --- capture call for reproducibility ---
  mc <- match.call()

  # --- seed ---
  used_seed <- seed
  if (!is.null(seed)) {
    set.seed(as.integer(seed))
  } else {
    used_seed <- sample.int(.Machine$integer.max, 1L)
    set.seed(used_seed)
  }

  # --- resolve vars ---
  if (is.numeric(vars) && length(vars) == 1L) {
    p <- as.integer(vars)
    var_names <- paste0("V", seq_len(p))
  } else if (is.character(vars)) {
    var_names <- vars
    p <- length(var_names)
  } else {
    stop("`vars` must be a single integer or a character vector of variable names")
  }
  stopifnot(p >= 2L)

  # --- validate n, tp ---
  n  <- as.integer(n)
  tp <- as.integer(tp)
  stopifnot(n >= 2L, tp >= 3L)

  # --- beeps_per_day ---
  if (!is.null(beeps_per_day)) {
    beeps_per_day <- as.integer(beeps_per_day)
    if (tp %% beeps_per_day != 0L) {
      stop("`tp` (", tp, ") must be divisible by `beeps_per_day` (", beeps_per_day, ")")
    }
    n_days <- tp %/% beeps_per_day
  }

  # --- complexity: resolve early for tiny_n override ---
  cases <- .resolve_longitudinal_complexity(complexity)
  if ("tiny_n" %in% cases) {
    n  <- min(n, 5L)
    tp <- min(tp, 10L)
    if (!is.null(beeps_per_day)) {
      tp <- beeps_per_day * max(1L, tp %/% beeps_per_day)
      n_days <- tp %/% beeps_per_day
    }
  }

  # --- temporal matrix B ---
  if (is.null(temporal)) {
    temporal <- .generate_var_matrix(p, ar_range, cross_range, n_cross)
  } else {
    stopifnot(is.matrix(temporal), is.numeric(temporal),
              nrow(temporal) == p, ncol(temporal) == p)
  }
  colnames(temporal) <- rownames(temporal) <- var_names
  temporal <- .ensure_stationary(temporal)

  # --- contemporaneous correlation ---
  if (is.null(contemporaneous)) {
    contemporaneous <- diag(p)
  } else {
    stopifnot(is.matrix(contemporaneous), is.numeric(contemporaneous),
              nrow(contemporaneous) == p, ncol(contemporaneous) == p,
              isSymmetric(contemporaneous))
  }
  colnames(contemporaneous) <- rownames(contemporaneous) <- var_names

  # --- between covariance ---
  if (is.null(between)) {
    between <- diag(p)
  } else {
    stopifnot(is.matrix(between), is.numeric(between),
              nrow(between) == p, ncol(between) == p,
              isSymmetric(between))
  }
  colnames(between) <- rownames(between) <- var_names

  # --- grand means ---
  if (is.null(grand_means)) {
    grand_means <- rep(0, p)
  } else {
    stopifnot(is.numeric(grand_means), length(grand_means) == p)
  }
  names(grand_means) <- var_names

  # --- innovation sd ---
  if (length(innovation_sd) == 1L) {
    innovation_sd <- rep(innovation_sd, p)
  }
  stopifnot(is.numeric(innovation_sd), length(innovation_sd) == p,
            all(innovation_sd > 0))
  names(innovation_sd) <- var_names

  # --- Cholesky of innovation covariance ---
  innov_cov <- diag(innovation_sd) %*% contemporaneous %*% diag(innovation_sd)
  innov_cov <- .ensure_pd(innov_cov)
  L_innov <- chol(innov_cov)  # upper Cholesky

  # --- person means ---
  between_pd <- .ensure_pd(between)
  L_between <- chol(between_pd)
  person_means <- matrix(stats::rnorm(n * p), nrow = n, ncol = p) %*% L_between
  person_means <- sweep(person_means, 2L, grand_means, "+")

  # --- generate VAR(1) data per person ---
  heavy <- "heavy_tailed" %in% cases
  hetero <- "heteroscedastic" %in% cases

  all_rows <- vector("list", n)
  for (subj in seq_len(n)) {
    mu <- person_means[subj, ]
    subj_data <- matrix(NA_real_, nrow = tp, ncol = p)

    # First observation
    subj_data[1L, ] <- mu + as.vector(.draw_innovations(1L, p, L_innov, heavy))

    for (t in seq(2L, tp)) {
      # Day break: reset carry-over if ESM and new day
      is_day_break <- !is.null(beeps_per_day) && ((t - 1L) %% beeps_per_day == 0L)

      innov <- .draw_innovations(1L, p, L_innov, heavy)

      # Heteroscedastic: scale innovations by subject index
      if (hetero) {
        subj_scale <- 0.5 + 1.5 * (subj - 1L) / max(n - 1L, 1L)
        innov <- innov * subj_scale
      }

      if (is_day_break) {
        subj_data[t, ] <- mu + as.vector(innov)
      } else {
        subj_data[t, ] <- mu +
          as.vector(temporal %*% (subj_data[t - 1L, ] - mu)) +
          as.vector(innov)
      }
    }

    all_rows[[subj]] <- subj_data
  }

  # --- assemble data.frame ---
  big_mat <- do.call(rbind, all_rows)
  colnames(big_mat) <- var_names

  df <- data.frame(
    id   = rep(seq_len(n), each = tp),
    time = rep(seq_len(tp), times = n)
  )

  if (!is.null(beeps_per_day)) {
    df$day  <- rep(rep(seq_len(n_days), each = beeps_per_day), times = n)
    df$beep <- rep(rep(seq_len(beeps_per_day), times = n_days), times = n)
    df$time <- NULL
  }

  df[var_names] <- big_mat

  # --- post-hoc complexity injection ---
  df <- .inject_longitudinal_complexity(df, cases, var_names)

  # --- return saqr_sim ---
  saqr_sim(
    data   = df,
    params = list(
      temporal        = temporal,
      contemporaneous = contemporaneous,
      between         = between,
      grand_means     = grand_means,
      innovation_sd   = innovation_sd,
      n               = n,
      tp              = tp,
      var_names       = var_names
    ),
    type = "longitudinal",
    seed = used_seed,
    call = mc
  )
}


# ===========================================================================
# Internal helpers
# ===========================================================================

#' Generate a random VAR(1) matrix with controlled sparsity
#' @keywords internal
#' @noRd
.generate_var_matrix <- function(p, ar_range, cross_range, n_cross) {
  B <- diag(stats::runif(p, ar_range[1L], ar_range[2L]))

  # Off-diagonal cross-lags
  off_diag <- which(!diag(TRUE, p), arr.ind = TRUE)
  max_cross <- nrow(off_diag)
  if (is.null(n_cross)) n_cross <- min(3L, max_cross)
  n_cross <- min(n_cross, max_cross)

  if (n_cross > 0L) {
    chosen <- off_diag[sample.int(max_cross, n_cross), , drop = FALSE]
    vals <- stats::runif(n_cross, cross_range[1L], cross_range[2L]) *
      sample(c(-1L, 1L), n_cross, replace = TRUE)
    for (k in seq_len(n_cross)) {
      B[chosen[k, 1L], chosen[k, 2L]] <- vals[k]
    }
  }
  B
}


#' Rescale B so all eigenvalues are inside the unit circle
#' @keywords internal
#' @noRd
.ensure_stationary <- function(B) {
  ev <- eigen(B, only.values = TRUE)$values
  max_mod <- max(Mod(ev))
  if (max_mod >= 1) {
    B <- B * (0.95 / max_mod)
  }
  B
}


#' Nearest positive-definite matrix (eigenvalue clamp)
#' @keywords internal
#' @noRd
.ensure_pd <- function(M) {
  eig <- eigen(M, symmetric = TRUE)
  vals <- pmax(eig$values, 1e-6)
  out <- eig$vectors %*% diag(vals) %*% t(eig$vectors)
  (out + t(out)) / 2
}


#' Draw innovations — normal or heavy-tailed
#' @keywords internal
#' @noRd
.draw_innovations <- function(n_obs, p, L_chol, heavy = FALSE) {
  if (heavy) {
    df_t <- sample(3L:6L, 1L)
    Z <- matrix(stats::rt(n_obs * p, df = df_t), nrow = n_obs, ncol = p)
  } else {
    Z <- matrix(stats::rnorm(n_obs * p), nrow = n_obs, ncol = p)
  }
  Z %*% L_chol
}


#' Resolve complexity parameter for longitudinal data
#' @keywords internal
#' @noRd
.resolve_longitudinal_complexity <- function(complexity) {
  valid <- c("na", "outliers", "heavy_tailed", "heteroscedastic", "tiny_n")
  if (identical(complexity, "clean")) return(character(0L))
  if (identical(complexity, "auto")) {
    n_cases <- sample(0L:2L, 1L, prob = c(0.4, 0.4, 0.2))
    if (n_cases == 0L) return(character(0L))
    return(sample(valid, n_cases))
  }
  bad <- setdiff(complexity, valid)
  if (length(bad) > 0L) {
    stop("Unknown complexity cases: ", paste(bad, collapse = ", "),
         ". Valid: ", paste(valid, collapse = ", "))
  }
  complexity
}


#' Inject edge cases into assembled data.frame
#' @keywords internal
#' @noRd
.inject_longitudinal_complexity <- function(df, cases, var_names) {
  if (length(cases) == 0L) return(df)

  n_total <- nrow(df)

  if ("na" %in% cases) {
    n_na <- max(1L, as.integer(n_total * 0.03))
    na_idx <- sample.int(n_total, n_na)
    na_col <- sample(var_names, n_na, replace = TRUE)
    for (k in seq_along(na_idx)) {
      df[na_idx[k], na_col[k]] <- NA_real_
    }
  }

  if ("outliers" %in% cases) {
    n_out <- max(1L, as.integer(n_total * 0.02))
    out_idx <- sample.int(n_total, n_out)
    out_col <- sample(var_names, n_out, replace = TRUE)
    for (k in seq_along(out_idx)) {
      df[out_idx[k], out_col[k]] <- df[out_idx[k], out_col[k]] +
        sample(c(-1, 1), 1L) * stats::runif(1, 5, 10) * stats::sd(df[[out_col[k]]], na.rm = TRUE)
    }
  }

  df
}
