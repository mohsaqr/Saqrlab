#' Simulate Latent Profile Analysis Data
#'
#' @description Generate synthetic continuous-indicator data with known latent
#'   profile structure. All distributional parameters are explicitly specified
#'   so that parameter recovery can be verified in tests.
#'
#' @param means Numeric matrix of dimension \code{n_vars x n_profiles}. Each
#'   column is the vector of variable means for that profile.
#' @param sds Numeric matrix of dimension \code{n_vars x n_profiles}, or a
#'   scalar / vector recycled to that shape. Standard deviations within each
#'   profile for each variable. All values must be positive.
#' @param props Numeric vector of length \code{n_profiles}. Mixing
#'   proportions. Need not sum to 1 — normalised internally.
#' @param n Positive integer. Total sample size.
#' @param seed Integer or NULL. Random seed for reproducibility.
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{\code{data}}{data.frame with columns \code{y1}…\code{yp} and
#'       \code{true_profile} (integer 1…K).}
#'     \item{\code{params}}{list with \code{means}, \code{sds} (full matrix),
#'       and \code{props} (as supplied, un-normalised).}
#'   }
#'
#' @examples
#' means <- matrix(c(0, 0, 10, 10), nrow = 2, ncol = 2)
#' r <- simulate_lpa(means = means, sds = 0.5, props = c(0.5, 0.5), n = 200, seed = 1)
#' r$data        # data.frame
#' r$params      # ground-truth parameters
#'
#' @export
simulate_lpa <- function(means, sds, props, n, seed = NULL) {
  stopifnot(is.matrix(means), is.numeric(means))
  p <- nrow(means)
  K <- ncol(means)
  if (length(props) != K) {
    stop("`props` must have one entry per profile (ncol of `means`)")
  }
  stopifnot(
    is.numeric(props), all(is.finite(props)), all(props >= 0), sum(props) > 0
  )
  stopifnot(is.numeric(n), length(n) == 1L, n >= 1L)

  sds_mat <- matrix(sds, nrow = p, ncol = K)
  stopifnot(is.numeric(sds_mat), all(sds_mat > 0))

  if (!is.null(seed)) set.seed(as.integer(seed))

  norm_props   <- props / sum(props)
  class_assign <- sample(K, n, replace = TRUE, prob = norm_props)

  X <- vapply(seq_len(p), function(i) {
    stats::rnorm(n,
                 mean = means[i, class_assign],
                 sd   = sds_mat[i, class_assign])
  }, numeric(n))

  df              <- as.data.frame(X)
  colnames(df)    <- paste0("y", seq_len(p))
  df$true_profile <- class_assign

  list(
    data   = df,
    params = list(means = means, sds = sds_mat, props = props)
  )
}


#' Simulate Latent Class Analysis Data
#'
#' @description Generate synthetic binary-indicator data with known latent
#'   class structure. All item response probabilities are explicitly specified.
#'
#' @param item_probs Numeric matrix of dimension \code{n_items x n_classes}.
#'   Each column gives the probability that each item equals 1 for that class.
#'   All values must be in \eqn{[0, 1]}.
#' @param class_probs Numeric vector of length \code{n_classes}. Class mixing
#'   proportions. Need not sum to 1 — normalised internally.
#' @param n Positive integer. Total sample size.
#' @param seed Integer or NULL. Random seed for reproducibility.
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{\code{data}}{data.frame with binary columns \code{item1}…\code{itemm}
#'       and integer column \code{true_class} (labels 1…K).}
#'     \item{\code{params}}{list with \code{item_probs} and \code{class_probs}
#'       (as supplied, un-normalised).}
#'   }
#'
#' @examples
#' item_probs <- matrix(c(0.9, 0.1, 0.9, 0.1,
#'                        0.1, 0.9, 0.1, 0.9), nrow = 4, ncol = 2)
#' r <- simulate_lca(item_probs = item_probs, class_probs = c(0.5, 0.5),
#'                   n = 200, seed = 1)
#' r$data        # data.frame
#' r$params      # ground-truth parameters
#'
#' @export
simulate_lca <- function(item_probs, class_probs, n, seed = NULL) {
  stopifnot(is.matrix(item_probs), is.numeric(item_probs))
  if (any(item_probs < 0 | item_probs > 1, na.rm = TRUE)) {
    stop("`item_probs` must contain values in [0, 1]")
  }
  m <- nrow(item_probs)
  K <- ncol(item_probs)
  if (length(class_probs) != K) {
    stop("`class_probs` must have one entry per class (ncol of `item_probs`)")
  }
  stopifnot(
    is.numeric(class_probs), all(is.finite(class_probs)),
    all(class_probs >= 0), sum(class_probs) > 0
  )
  stopifnot(is.numeric(n), length(n) == 1L, n >= 1L)

  if (!is.null(seed)) set.seed(as.integer(seed))

  norm_probs   <- class_probs / sum(class_probs)
  class_assign <- sample(K, n, replace = TRUE, prob = norm_probs)

  X <- vapply(seq_len(m), function(j) {
    stats::rbinom(n, size = 1L, prob = item_probs[j, class_assign])
  }, integer(n))

  df            <- as.data.frame(X)
  colnames(df)  <- paste0("item", seq_len(m))
  df$true_class <- class_assign

  list(
    data   = df,
    params = list(item_probs = item_probs, class_probs = class_probs)
  )
}


#' Simulate Linear Regression Data with Known Coefficients
#'
#' @description Generate a dataset suitable for linear regression with
#'   fully specified ground-truth coefficients and predictor standard
#'   deviations. Designed so that \code{lm(y ~ ., data = r$data)} recovers
#'   the true coefficients at large \code{n}.
#'
#' @param coefs Named numeric vector. May include \code{"(Intercept)"}. All
#'   other names become predictor column names. Predictors are generated as
#'   independent \eqn{N(0, \sigma_j)} variables.
#' @param predictor_sds Named numeric vector. Names must exactly match the
#'   non-intercept names in \code{coefs}. All values must be positive.
#' @param error_sd Positive numeric. Standard deviation of the residuals.
#' @param n Positive integer. Sample size.
#' @param seed Integer or NULL. Random seed.
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{\code{data}}{data.frame with column \code{y} and one column per
#'       predictor.}
#'     \item{\code{params}}{list with \code{coefs}, \code{predictor_sds}, and
#'       \code{error_sd}.}
#'   }
#'
#' @examples
#' coefs <- c("(Intercept)" = 2, x1 = 3, x2 = -1)
#' r <- simulate_regression(coefs = coefs,
#'                           predictor_sds = c(x1 = 1, x2 = 1),
#'                           error_sd = 0.5, n = 500, seed = 42)
#' coef(lm(y ~ ., data = r$data))  # should be close to c(2, 3, -1)
#'
#' @export
simulate_regression <- function(coefs, predictor_sds, error_sd, n, seed = NULL) {
  stopifnot(is.numeric(coefs), length(coefs) >= 1L, !is.null(names(coefs)))
  if (error_sd <= 0 || !is.finite(error_sd)) {
    stop("`error_sd` must be a positive finite number")
  }
  stopifnot(is.numeric(n), length(n) == 1L, n >= 1L)

  pred_names <- setdiff(names(coefs), "(Intercept)")

  if (!is.numeric(predictor_sds) ||
      !identical(sort(names(predictor_sds)), sort(pred_names))) {
    stop(
      "`predictor_sds` must be a named numeric vector with names matching ",
      "the non-intercept entries of `coefs`. Expected: ",
      paste(pred_names, collapse = ", ")
    )
  }
  if (any(predictor_sds <= 0 | !is.finite(predictor_sds))) {
    stop("`predictor_sds` must all be positive and finite")
  }

  if (!is.null(seed)) set.seed(as.integer(seed))

  intercept <- if ("(Intercept)" %in% names(coefs)) coefs[["(Intercept)"]] else 0

  X <- vapply(pred_names, function(nm) {
    stats::rnorm(n, mean = 0, sd = predictor_sds[[nm]])
  }, numeric(n))
  if (length(pred_names) == 1L) X <- matrix(X, ncol = 1L)

  slopes <- coefs[pred_names]
  y      <- intercept + drop(X %*% slopes) + stats::rnorm(n, 0, error_sd)

  df           <- as.data.frame(X)
  colnames(df) <- pred_names
  df$y         <- y
  df           <- df[, c("y", pred_names), drop = FALSE]

  list(
    data   = df,
    params = list(coefs = coefs, predictor_sds = predictor_sds, error_sd = error_sd)
  )
}


#' Simulate Factor Analysis Data with Known Parameters
#'
#' @description Generate multivariate normal data from an explicit factor model
#'   for factor analysis parameter recovery testing. The implied covariance
#'   matrix \eqn{\Sigma = \Lambda \Phi \Lambda' + \Psi} is fully specified and
#'   stored in the return value for direct comparison with estimated covariances.
#'
#' @param loadings Numeric matrix of dimension \code{p x m} (observed variables
#'   x factors). Must be a matrix (not a vector). Zero entries indicate no
#'   loading; use \code{matrix(val, nrow = p, ncol = 1)} for single-factor
#'   models.
#' @param phi Numeric matrix of dimension \code{m x m}. Factor correlation
#'   matrix. Must be symmetric and positive-definite. Defaults to
#'   \code{diag(m)} (orthogonal factors).
#' @param psi Numeric vector of length \code{p}. Unique variances. All values
#'   must be positive. Defaults to
#'   \code{1 - diag(loadings \%*\% phi \%*\% t(loadings))} (communality
#'   complement). An error is raised if any auto-computed value is \eqn{\le 0}
#'   (over-factored model).
#' @param n Positive integer. Sample size.
#' @param seed Integer or NULL. Random seed.
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{\code{data}}{data.frame with columns \code{y1}…\code{yp}.}
#'     \item{\code{params}}{list with \code{loadings}, \code{phi}, \code{psi},
#'       and \code{sigma_implied} (the \eqn{p x p} model-implied covariance
#'       matrix).}
#'   }
#'
#' @examples
#' # Two orthogonal factors, 6 indicators (3 per factor)
#' loadings <- matrix(c(0.8, 0.7, 0.6, 0,   0,   0,
#'                      0,   0,   0,   0.8, 0.7, 0.6),
#'                    nrow = 6, ncol = 2)  # byrow = FALSE (default)
#' r <- simulate_fa(loadings = loadings, n = 500, seed = 1)
#' r$data                  # data.frame: y1..y6
#' r$params$sigma_implied  # implied covariance matrix
#'
#' # Oblique model
#' phi <- matrix(c(1, 0.4, 0.4, 1), nrow = 2)
#' r2  <- simulate_fa(loadings = loadings, phi = phi, n = 500, seed = 1)
#'
#' @export
simulate_fa <- function(loadings, phi = NULL, psi = NULL, n, seed = NULL) {
  if (!is.matrix(loadings) || !is.numeric(loadings)) {
    stop("`loadings` must be a numeric matrix — for single-factor, use matrix(val, nrow = p, ncol = 1)")
  }
  stopifnot(!anyNA(loadings), nrow(loadings) >= 1L, ncol(loadings) >= 1L)
  p <- nrow(loadings)
  m <- ncol(loadings)

  # --- phi ---
  if (is.null(phi)) {
    phi <- diag(m)
  } else {
    if (!is.matrix(phi) || !is.numeric(phi) || nrow(phi) != m || ncol(phi) != m) {
      stop("`phi` must be a numeric ", m, " x ", m, " matrix (one row/col per factor)")
    }
    if (!isSymmetric(phi)) {
      stop("`phi` must be symmetric")
    }
    pd_check <- tryCatch(chol(phi), error = function(e) e)
    if (inherits(pd_check, "error")) {
      stop("`phi` must be positive-definite: ", pd_check$message)
    }
  }

  # --- psi ---
  communalities <- diag(loadings %*% phi %*% t(loadings))
  if (is.null(psi)) {
    psi <- 1 - communalities
    bad <- which(psi <= 0)
    if (length(bad) > 0L) {
      stop(
        "`psi` auto-computed values <= 0 for variable(s): ",
        paste(bad, collapse = ", "),
        ". Reduce loadings (communality >= 1 indicates an over-factored model)."
      )
    }
  } else {
    if (!is.numeric(psi) || length(psi) != p) {
      stop("`psi` must be a numeric vector of length p (", p, ")")
    }
    bad <- which(psi <= 0)
    if (length(bad) > 0L) {
      stop("`psi` values must all be positive; non-positive at index: ",
           paste(bad, collapse = ", "))
    }
  }

  # --- implied covariance (raw eigenvalue clamp, no diagonal normalisation) ---
  Sigma <- loadings %*% phi %*% t(loadings) + diag(psi)
  eig   <- eigen(Sigma, symmetric = TRUE)
  Sigma <- eig$vectors %*% diag(pmax(eig$values, 1e-6)) %*% t(eig$vectors)
  Sigma <- (Sigma + t(Sigma)) / 2   # enforce exact symmetry

  # --- generate data via Cholesky ---
  if (!is.null(seed)) set.seed(as.integer(seed))
  L <- chol(Sigma)   # upper Cholesky: Sigma = t(L) %*% L
  Z <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
  X <- Z %*% L

  df           <- as.data.frame(X)
  colnames(df) <- paste0("y", seq_len(p))

  list(
    data   = df,
    params = list(
      loadings      = loadings,
      phi           = phi,
      psi           = psi,
      sigma_implied = Sigma
    )
  )
}
