# ===========================================================================
# Multilevel (MLM) and latent growth-curve simulation
# All return saqr_sim with $data and $params
# ===========================================================================


#' Simulate Two-Level Multilevel (Hierarchical) Data
#'
#' @description Generate data for a two-level model where level-1 units are
#'   nested within clusters. The outcome follows
#'   \deqn{y_{ij} = (\beta_0 + u_{0j}) + \sum_k \beta_k x_{kij}
#'         + (u_{1j} x_{1ij}) + e_{ij},}
#'   where the first predictor carries the fixed \code{slope}, cluster random
#'   intercepts \eqn{u_{0j} \sim N(0, \tau_{00})} with
#'   \eqn{\tau_{00} = \frac{icc}{1 - icc}\,\sigma_e^2}, and (optionally) a
#'   cluster random slope \eqn{u_{1j} \sim N(0, \mathrm{slope\_sd}^2)} on the
#'   first predictor. Designed so that the intraclass correlation and fixed
#'   slope are recovered at large \code{n_clusters}.
#'
#' @param n_clusters Integer. Number of level-2 clusters. Default: 30.
#' @param cluster_size Integer scalar or integer vector of length
#'   \code{n_clusters}. Number of level-1 units per cluster. A vector gives an
#'   unbalanced design. Default: 20.
#' @param intercept Numeric. Fixed intercept \eqn{\beta_0}. Default: 0.
#' @param slope Numeric. Fixed slope \eqn{\beta_1} on the first predictor
#'   (\code{x1}). Default: 0.5.
#' @param residual_sd Positive numeric. Level-1 residual standard deviation
#'   \eqn{\sigma_e}. Default: 1.
#' @param icc Numeric in \code{[0, 1)}. Target intraclass correlation; used to
#'   derive the random-intercept variance
#'   \eqn{\tau_{00} = \frac{icc}{1 - icc}\,\sigma_e^2}. Default: 0.1.
#' @param n_predictors Integer >= 1. Number of level-1 predictors
#'   (\code{x1}..\code{xK}). The first carries \code{slope}; the remaining
#'   predictors get fixed coefficients drawn from \code{N(0, 1)}. Default: 1.
#' @param random_slope Logical. If \code{TRUE}, add a cluster-level random
#'   slope on the first predictor. Default: \code{FALSE}.
#' @param slope_sd Non-negative numeric. Standard deviation of the random
#'   slope \eqn{u_{1j}} (used only when \code{random_slope = TRUE}).
#'   Default: 0.
#' @param seed Integer or NULL. Random seed.
#'
#' @return A \code{\link{saqr_sim}} object with:
#'   \describe{
#'     \item{\code{$data}}{data.frame with columns \code{cluster_id} (factor),
#'       \code{y} (numeric), and \code{x1}..\code{xK} (numeric predictors).}
#'     \item{\code{$params}}{list with \code{intercept}, \code{slope},
#'       \code{betas} (full fixed-effect vector), \code{icc}, \code{tau00},
#'       \code{residual_sd}, \code{slope_sd}, \code{random_slope},
#'       \code{cluster_size}, \code{cluster_intercepts} (true \eqn{u_{0j}}),
#'       and \code{cluster_slopes} (true \eqn{u_{1j}}, \code{NULL} unless
#'       \code{random_slope}).}
#'   }
#'
#' @examples
#' r <- simulate_mlm(n_clusters = 40, cluster_size = 25, slope = 0.8,
#'                   icc = 0.2, seed = 1)
#' print(r)
#' r$params$tau00
#'
#' # Unbalanced clusters with a random slope
#' r2 <- simulate_mlm(n_clusters = 30, cluster_size = sample(10:30, 30, TRUE),
#'                    random_slope = TRUE, slope_sd = 0.3, seed = 42)
#' summary(r2)
#'
#' @export
simulate_mlm <- function(n_clusters = 30, cluster_size = 20, intercept = 0,
                         slope = 0.5, residual_sd = 1, icc = 0.1,
                         n_predictors = 1, random_slope = FALSE,
                         slope_sd = 0, seed = NULL) {
  stopifnot(
    is.numeric(n_clusters), length(n_clusters) == 1L, n_clusters >= 2L,
    is.numeric(cluster_size), all(cluster_size >= 1L),
    length(cluster_size) == 1L || length(cluster_size) == n_clusters,
    is.numeric(intercept), length(intercept) == 1L,
    is.numeric(slope), length(slope) == 1L,
    is.numeric(residual_sd), length(residual_sd) == 1L, residual_sd > 0,
    is.numeric(icc), length(icc) == 1L, icc >= 0, icc < 1,
    is.numeric(n_predictors), length(n_predictors) == 1L, n_predictors >= 1L,
    is.logical(random_slope), length(random_slope) == 1L,
    is.numeric(slope_sd), length(slope_sd) == 1L, slope_sd >= 0
  )
  n_clusters <- as.integer(n_clusters)
  n_predictors <- as.integer(n_predictors)

  if (!is.null(seed)) set.seed(as.integer(seed))

  # Per-cluster level-1 sizes (scalar recycled, or unbalanced vector)
  sizes <- if (length(cluster_size) == 1L) {
    rep(as.integer(cluster_size), n_clusters)
  } else {
    as.integer(cluster_size)
  }
  n_total <- sum(sizes)

  # Fixed-effect vector: first predictor carries `slope`, rest drawn N(0, 1)
  betas <- c(slope, if (n_predictors > 1L) {
    stats::rnorm(n_predictors - 1L)
  } else {
    numeric(0)
  })
  names(betas) <- paste0("x", seq_len(n_predictors))

  # Random-intercept variance derived from ICC
  tau00 <- icc / (1 - icc) * residual_sd^2

  # Cluster-level random effects
  u0 <- stats::rnorm(n_clusters, mean = 0, sd = sqrt(tau00))
  u1 <- if (random_slope) {
    stats::rnorm(n_clusters, mean = 0, sd = slope_sd)
  } else {
    rep(0, n_clusters)
  }

  cluster_index <- rep(seq_len(n_clusters), times = sizes)

  # Level-1 predictor matrix (n_total x n_predictors)
  X <- matrix(stats::rnorm(n_total * n_predictors), nrow = n_total,
              ncol = n_predictors)
  colnames(X) <- names(betas)

  fixed_part <- intercept + drop(X %*% betas)
  random_intercept <- u0[cluster_index]
  random_slope_part <- u1[cluster_index] * X[, 1L]
  residual <- stats::rnorm(n_total, mean = 0, sd = residual_sd)

  y <- fixed_part + random_intercept + random_slope_part + residual

  df <- data.frame(
    cluster_id = factor(cluster_index),
    y = y
  )
  df <- cbind(df, as.data.frame(X))

  saqr_sim(
    data   = df,
    params = list(
      intercept          = intercept,
      slope              = slope,
      betas              = betas,
      icc                = icc,
      tau00              = tau00,
      residual_sd        = residual_sd,
      slope_sd           = slope_sd,
      random_slope       = random_slope,
      cluster_size       = stats::setNames(sizes, paste0("cluster_",
                                                         seq_len(n_clusters))),
      cluster_intercepts = stats::setNames(u0, paste0("cluster_",
                                                      seq_len(n_clusters))),
      cluster_slopes     = if (random_slope) {
        stats::setNames(u1, paste0("cluster_", seq_len(n_clusters)))
      } else {
        NULL
      }
    ),
    type = "mlm",
    seed = seed,
    call = match.call()
  )
}


#' Simulate Latent Growth-Curve Data
#'
#' @description Generate longitudinal data from a linear latent growth-curve
#'   model. Each subject draws a random intercept and slope from a bivariate
#'   normal with specified means, standard deviations, and correlation, then
#'   \deqn{y_{it} = b_{0i} + b_{1i} \, t + e_{it},}
#'   with \code{time} coded \code{0..(n_time - 1)}. The result is returned in
#'   long format. Designed so that the mean per-subject slope recovers
#'   \code{slope_mean} and the spread of per-subject intercepts recovers
#'   \code{intercept_sd} at large \code{n}.
#'
#' @param n Integer. Number of subjects. Default: 200.
#' @param n_time Integer >= 2. Number of measurement occasions per subject.
#'   Default: 5.
#' @param intercept_mean Numeric. Mean of subject intercepts \eqn{b_{0i}}.
#'   Default: 0.
#' @param slope_mean Numeric. Mean of subject slopes \eqn{b_{1i}}. Default: 1.
#' @param intercept_sd Positive numeric. Standard deviation of subject
#'   intercepts. Default: 1.
#' @param slope_sd Positive numeric. Standard deviation of subject slopes.
#'   Default: 0.5.
#' @param intercept_slope_cor Numeric in \code{[-1, 1]}. Correlation between
#'   subject intercepts and slopes. Default: 0.
#' @param residual_sd Positive numeric. Level-1 (occasion) residual standard
#'   deviation \eqn{\sigma_e}. Default: 1.
#' @param seed Integer or NULL. Random seed.
#'
#' @return A \code{\link{saqr_sim}} object with:
#'   \describe{
#'     \item{\code{$data}}{long-format data.frame with columns \code{subject}
#'       (factor), \code{time} (numeric, \code{0..n_time-1}), and \code{y}
#'       (numeric).}
#'     \item{\code{$params}}{list with \code{means} (intercept/slope),
#'       \code{sds} (intercept/slope), \code{correlation}, \code{residual_sd},
#'       \code{n_time}, \code{subject_intercepts} (true \eqn{b_{0i}}), and
#'       \code{subject_slopes} (true \eqn{b_{1i}}).}
#'   }
#'
#' @examples
#' r <- simulate_growth(n = 300, n_time = 6, slope_mean = 2,
#'                      intercept_sd = 1.5, seed = 1)
#' print(r)
#' head(r)
#'
#' # Correlated intercepts and slopes
#' r2 <- simulate_growth(n = 250, intercept_slope_cor = 0.5, seed = 42)
#' r2$params$correlation
#'
#' @export
simulate_growth <- function(n = 200, n_time = 5, intercept_mean = 0,
                            slope_mean = 1, intercept_sd = 1, slope_sd = 0.5,
                            intercept_slope_cor = 0, residual_sd = 1,
                            seed = NULL) {
  stopifnot(
    is.numeric(n), length(n) == 1L, n >= 2L,
    is.numeric(n_time), length(n_time) == 1L, n_time >= 2L,
    is.numeric(intercept_mean), length(intercept_mean) == 1L,
    is.numeric(slope_mean), length(slope_mean) == 1L,
    is.numeric(intercept_sd), length(intercept_sd) == 1L, intercept_sd > 0,
    is.numeric(slope_sd), length(slope_sd) == 1L, slope_sd > 0,
    is.numeric(intercept_slope_cor), length(intercept_slope_cor) == 1L,
    intercept_slope_cor >= -1, intercept_slope_cor <= 1,
    is.numeric(residual_sd), length(residual_sd) == 1L, residual_sd > 0
  )
  n <- as.integer(n)
  n_time <- as.integer(n_time)

  if (!is.null(seed)) set.seed(as.integer(seed))

  # Bivariate normal (intercept, slope) via Cholesky transform of i.i.d.
  # standard normals -- no loop.
  means <- c(intercept = intercept_mean, slope = slope_mean)
  sds   <- c(intercept = intercept_sd, slope = slope_sd)
  cov_is <- intercept_slope_cor * intercept_sd * slope_sd
  sigma <- matrix(c(intercept_sd^2, cov_is,
                    cov_is, slope_sd^2), nrow = 2L)
  L <- chol(sigma)
  Z <- matrix(stats::rnorm(n * 2L), nrow = n, ncol = 2L)
  re <- sweep(Z %*% L, 2L, means, "+")
  b0 <- re[, 1L]
  b1 <- re[, 2L]

  # Long format: time varies fastest within subject.
  time_grid <- 0:(n_time - 1L)
  subject_index <- rep(seq_len(n), each = n_time)
  time_long <- rep(time_grid, times = n)

  y <- b0[subject_index] + b1[subject_index] * time_long +
    stats::rnorm(n * n_time, mean = 0, sd = residual_sd)

  df <- data.frame(
    subject = factor(subject_index),
    time    = time_long,
    y       = y
  )

  saqr_sim(
    data   = df,
    params = list(
      means              = means,
      sds                = sds,
      correlation        = intercept_slope_cor,
      residual_sd        = residual_sd,
      n_time             = n_time,
      subject_intercepts = stats::setNames(b0, paste0("s", seq_len(n))),
      subject_slopes     = stats::setNames(b1, paste0("s", seq_len(n)))
    ),
    type = "growth",
    seed = seed,
    call = match.call()
  )
}
