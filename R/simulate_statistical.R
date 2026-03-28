# ===========================================================================
# Explicit-parameter simulation for common statistical tests
# All return saqr_sim with $data and $params
# ===========================================================================


#' Simulate Two-Group Comparison Data (t-test)
#'
#' @description Generate data for a two-group comparison with fully specified
#'   ground-truth means and standard deviations. Designed so that
#'   \code{t.test(score ~ group, data = r$data)} recovers the true difference
#'   at large \code{n}.
#'
#' @param n_a Integer. Sample size for group A.
#' @param n_b Integer. Sample size for group B.
#' @param mean_a Numeric. Population mean for group A.
#' @param mean_b Numeric. Population mean for group B.
#' @param sd_a Positive numeric. Standard deviation for group A. Default: 1.
#' @param sd_b Positive numeric. Standard deviation for group B. Default:
#'   same as \code{sd_a}.
#' @param labels Character vector of length 2. Group labels. Default:
#'   \code{c("A", "B")}.
#' @param seed Integer or NULL. Random seed.
#'
#' @return A \code{\link{saqr_sim}} object with:
#'   \describe{
#'     \item{\code{$data}}{data.frame with columns \code{group} (factor) and
#'       \code{score} (numeric).}
#'     \item{\code{$params}}{list with \code{mean_a}, \code{mean_b}, \code{sd_a},
#'       \code{sd_b}, \code{n_a}, \code{n_b}, \code{cohens_d} (true Cohen's d
#'       based on pooled SD).}
#'   }
#'
#' @examples
#' r <- simulate_ttest(n_a = 50, n_b = 50, mean_a = 100, mean_b = 105, seed = 1)
#' t.test(score ~ group, data = r$data)
#' r$params$cohens_d
#'
#' @export
simulate_ttest <- function(n_a, n_b, mean_a, mean_b,
                           sd_a = 1, sd_b = sd_a,
                           labels = c("A", "B"),
                           seed = NULL) {
  stopifnot(
    is.numeric(n_a), length(n_a) == 1L, n_a >= 2L,
    is.numeric(n_b), length(n_b) == 1L, n_b >= 2L,
    is.numeric(mean_a), length(mean_a) == 1L,
    is.numeric(mean_b), length(mean_b) == 1L,
    is.numeric(sd_a), length(sd_a) == 1L, sd_a > 0,
    is.numeric(sd_b), length(sd_b) == 1L, sd_b > 0,
    is.character(labels), length(labels) == 2L
  )
  n_a <- as.integer(n_a)
  n_b <- as.integer(n_b)

  if (!is.null(seed)) set.seed(as.integer(seed))

  score_a <- stats::rnorm(n_a, mean = mean_a, sd = sd_a)
  score_b <- stats::rnorm(n_b, mean = mean_b, sd = sd_b)

  pooled_sd <- sqrt(((n_a - 1L) * sd_a^2 + (n_b - 1L) * sd_b^2) /
                      (n_a + n_b - 2L))

  df <- data.frame(
    group = factor(c(rep(labels[1L], n_a), rep(labels[2L], n_b)),
                   levels = labels),
    score = c(score_a, score_b)
  )

  saqr_sim(
    data   = df,
    params = list(
      mean_a   = mean_a,
      mean_b   = mean_b,
      sd_a     = sd_a,
      sd_b     = sd_b,
      n_a      = n_a,
      n_b      = n_b,
      cohens_d = (mean_b - mean_a) / pooled_sd
    ),
    type = "ttest",
    seed = seed,
    call = match.call()
  )
}


#' Simulate Multi-Group Comparison Data (ANOVA)
#'
#' @description Generate data for a one-way ANOVA with fully specified
#'   group means and standard deviations. Designed so that
#'   \code{aov(score ~ group, data = r$data)} recovers the true group effects
#'   at large \code{n}.
#'
#' @param n Integer or integer vector. If a single value, the per-group sample
#'   size (equal groups). If a vector, must have length equal to
#'   \code{length(means)} giving per-group sizes.
#' @param means Numeric vector. Population mean for each group.
#' @param sds Numeric scalar or vector. Standard deviation(s) for each group.
#'   Scalar is recycled. Default: 1.
#' @param labels Character vector or NULL. Group labels. Default:
#'   \code{paste0("G", seq_along(means))}.
#' @param seed Integer or NULL. Random seed.
#'
#' @return A \code{\link{saqr_sim}} object with:
#'   \describe{
#'     \item{\code{$data}}{data.frame with columns \code{group} (factor) and
#'       \code{score} (numeric).}
#'     \item{\code{$params}}{list with \code{means}, \code{sds} (as vector),
#'       \code{n} (per-group sizes), \code{labels}, \code{eta_squared} (true
#'       population \eqn{\eta^2}).}
#'   }
#'
#' @examples
#' r <- simulate_anova(n = 30, means = c(10, 12, 15), seed = 1)
#' summary(aov(score ~ group, data = r$data))
#' r$params$eta_squared
#'
#' # Unequal groups and heteroscedastic
#' r2 <- simulate_anova(n = c(50, 30, 20), means = c(5, 5, 10),
#'                       sds = c(1, 2, 3), seed = 42)
#'
#' @export
simulate_anova <- function(n, means, sds = 1, labels = NULL, seed = NULL) {
  k <- length(means)
  stopifnot(is.numeric(means), k >= 2L)

  # Resolve n to per-group vector
  if (length(n) == 1L) {
    stopifnot(is.numeric(n), n >= 2L)
    n_per <- rep(as.integer(n), k)
  } else {
    stopifnot(is.numeric(n), length(n) == k, all(n >= 2L))
    n_per <- as.integer(n)
  }

  # Resolve sds
  if (length(sds) == 1L) sds <- rep(sds, k)
  stopifnot(is.numeric(sds), length(sds) == k, all(sds > 0))

  # Labels
  if (is.null(labels)) labels <- paste0("G", seq_len(k))
  stopifnot(is.character(labels), length(labels) == k)

  if (!is.null(seed)) set.seed(as.integer(seed))

  scores <- unlist(
    mapply(function(m, sz, s) stats::rnorm(sz, mean = m, sd = s),
           means, n_per, sds, SIMPLIFY = FALSE)
  )
  groups <- rep(labels, times = n_per)

  # True population eta-squared
  grand_mean <- sum(n_per * means) / sum(n_per)
  ss_between <- sum(n_per * (means - grand_mean)^2)
  ss_within  <- sum((n_per - 1L) * sds^2)
  eta_sq     <- ss_between / (ss_between + ss_within)

  df <- data.frame(
    group = factor(groups, levels = labels),
    score = scores
  )

  saqr_sim(
    data   = df,
    params = list(
      means       = stats::setNames(means, labels),
      sds         = stats::setNames(sds, labels),
      n           = stats::setNames(n_per, labels),
      labels      = labels,
      eta_squared = eta_sq
    ),
    type = "anova",
    seed = seed,
    call = match.call()
  )
}


#' Simulate Correlated Multivariate Data
#'
#' @description Generate multivariate normal data from an explicit correlation
#'   (or covariance) matrix. Designed so that \code{cor(r$data)} recovers
#'   \code{sigma} at large \code{n}.
#'
#' @param n Integer. Sample size.
#' @param sigma Numeric matrix. Either a correlation matrix (all diagonal = 1)
#'   or a covariance matrix. Must be symmetric and positive-definite (or
#'   near-PD; corrected internally).
#' @param means Numeric vector or NULL. Population means for each variable.
#'   Default: all zeros.
#' @param var_names Character vector or NULL. Variable names. Default:
#'   \code{paste0("x", seq_len(ncol(sigma)))}.
#' @param seed Integer or NULL. Random seed.
#'
#' @return A \code{\link{saqr_sim}} object with:
#'   \describe{
#'     \item{\code{$data}}{data.frame with columns \code{x1}…\code{xp}.}
#'     \item{\code{$params}}{list with \code{sigma} (the input matrix),
#'       \code{means}, \code{is_correlation} (logical).}
#'   }
#'
#' @examples
#' # Correlation matrix
#' R <- matrix(c(1, 0.6, 0.3,
#'               0.6, 1, 0.5,
#'               0.3, 0.5, 1), nrow = 3)
#' r <- simulate_correlation(n = 200, sigma = R, seed = 1)
#' cor(r$data)           # should approximate R
#' r$params$sigma        # the true matrix
#'
#' # With means and custom names
#' r2 <- simulate_correlation(n = 500, sigma = R,
#'                             means = c(10, 20, 30),
#'                             var_names = c("IQ", "GPA", "income"), seed = 42)
#'
#' @export
simulate_correlation <- function(n, sigma, means = NULL, var_names = NULL,
                                 seed = NULL) {
  stopifnot(is.numeric(n), length(n) == 1L, n >= 3L)
  stopifnot(is.matrix(sigma), is.numeric(sigma), isSymmetric(sigma))
  p <- ncol(sigma)
  stopifnot(p >= 2L)

  is_cor <- all(abs(diag(sigma) - 1) < 1e-10)

  if (is.null(means)) means <- rep(0, p)
  stopifnot(is.numeric(means), length(means) == p)

  if (is.null(var_names)) var_names <- paste0("x", seq_len(p))
  stopifnot(is.character(var_names), length(var_names) == p)

  n <- as.integer(n)
  if (!is.null(seed)) set.seed(as.integer(seed))

  # Ensure PD
  eig <- eigen(sigma, symmetric = TRUE)
  vals <- pmax(eig$values, 1e-6)
  sigma_pd <- eig$vectors %*% diag(vals) %*% t(eig$vectors)
  sigma_pd <- (sigma_pd + t(sigma_pd)) / 2

  L <- chol(sigma_pd)
  Z <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
  X <- sweep(Z %*% L, 2L, means, "+")

  df <- as.data.frame(X)
  colnames(df) <- var_names

  saqr_sim(
    data   = df,
    params = list(
      sigma          = sigma,
      means          = stats::setNames(means, var_names),
      is_correlation = is_cor
    ),
    type = "correlation",
    seed = seed,
    call = match.call()
  )
}


#' Simulate Cluster Data with Known Centers
#'
#' @description Generate multivariate data from a mixture of Gaussians with
#'   fully specified cluster centers and standard deviations. Designed for
#'   testing clustering algorithms where the ground truth is known.
#'
#' @param n Integer or integer vector. If a single value, the total sample
#'   size (allocated by \code{props}). If a vector of length \code{k}, the
#'   per-cluster sizes.
#' @param centers Numeric matrix (\code{k x d}). Each row is a cluster
#'   centroid. \code{k} = number of clusters, \code{d} = number of dimensions.
#' @param sds Numeric scalar, vector of length \code{k}, or matrix
#'   (\code{k x d}). Standard deviations per cluster (and optionally per
#'   dimension). Scalar is recycled. Default: 1.
#' @param props Numeric vector of length \code{k} or NULL. Mixing proportions
#'   when \code{n} is a single value. Normalised internally. Default: equal
#'   mixing.
#' @param seed Integer or NULL. Random seed.
#'
#' @return A \code{\link{saqr_sim}} object with:
#'   \describe{
#'     \item{\code{$data}}{data.frame with columns \code{x1}…\code{xd} and
#'       integer column \code{true_cluster} (1…k).}
#'     \item{\code{$params}}{list with \code{centers} (matrix),
#'       \code{sds} (k x d matrix), \code{props} (normalised), \code{n}
#'       (per-cluster sizes).}
#'   }
#'
#' @examples
#' centers <- matrix(c(0, 0,
#'                     5, 5,
#'                     10, 0), nrow = 3, byrow = TRUE)
#' r <- simulate_clusters(n = 300, centers = centers, seed = 1)
#' plot(x2 ~ x1, data = r$data, col = r$data$true_cluster)
#' r$params$centers
#'
#' # Per-cluster sizes and SDs
#' r2 <- simulate_clusters(n = c(100, 50, 150), centers = centers,
#'                          sds = c(0.5, 1.0, 2.0), seed = 42)
#'
#' @export
simulate_clusters <- function(n, centers, sds = 1, props = NULL,
                              seed = NULL) {
  stopifnot(is.matrix(centers), is.numeric(centers))
  k <- nrow(centers)
  d <- ncol(centers)
  stopifnot(k >= 2L, d >= 1L)

  # Resolve per-cluster sizes
  if (length(n) == 1L) {
    n_total <- as.integer(n)
    stopifnot(n_total >= k)
    if (is.null(props)) props <- rep(1, k)
    stopifnot(is.numeric(props), length(props) == k, all(props >= 0),
              sum(props) > 0)
    norm_props <- props / sum(props)
    raw_sizes <- round(norm_props * n_total)
    raw_sizes <- pmax(raw_sizes, 1L)
    raw_sizes[k] <- n_total - sum(raw_sizes[-k])
    n_per <- as.integer(raw_sizes)
  } else {
    stopifnot(is.numeric(n), length(n) == k, all(n >= 1L))
    n_per <- as.integer(n)
    n_total <- sum(n_per)
    if (is.null(props)) {
      norm_props <- n_per / n_total
    } else {
      norm_props <- props / sum(props)
    }
  }

  # Resolve sds to k x d matrix
  if (is.matrix(sds)) {
    stopifnot(nrow(sds) == k, ncol(sds) == d, all(sds > 0))
    sd_mat <- sds
  } else if (length(sds) == k) {
    stopifnot(all(sds > 0))
    sd_mat <- matrix(rep(sds, each = d), nrow = k, ncol = d, byrow = TRUE)
  } else if (length(sds) == 1L) {
    stopifnot(sds > 0)
    sd_mat <- matrix(sds, nrow = k, ncol = d)
  } else {
    stop("`sds` must be a scalar, vector of length k, or k x d matrix")
  }

  if (!is.null(seed)) set.seed(as.integer(seed))

  data_list <- lapply(seq_len(k), function(cl) {
    sweep(
      matrix(stats::rnorm(n_per[cl] * d), nrow = n_per[cl], ncol = d) *
        rep(sd_mat[cl, ], each = n_per[cl]),
      2L, centers[cl, ], "+"
    )
  })

  X <- do.call(rbind, data_list)
  true_cluster <- rep(seq_len(k), times = n_per)

  df <- as.data.frame(X)
  colnames(df) <- paste0("x", seq_len(d))
  df$true_cluster <- true_cluster

  saqr_sim(
    data   = df,
    params = list(
      centers = centers,
      sds     = sd_mat,
      props   = norm_props,
      n       = stats::setNames(n_per, paste0("cluster_", seq_len(k)))
    ),
    type = "clusters",
    seed = seed,
    call = match.call()
  )
}


#' Simulate Prediction/Regression Data with Known Coefficients
#'
#' @description Generate a regression dataset with both continuous and
#'   categorical predictors, a non-linear term, and known ground-truth
#'   coefficients. A richer version of \code{\link{simulate_regression}} for
#'   testing prediction workflows.
#'
#' @param n Integer. Sample size.
#' @param coefs Named numeric vector. Must include names matching continuous
#'   predictors (e.g. \code{x1}, \code{x2}). May include
#'   \code{"(Intercept)"}. Non-linear terms are NOT automatically generated;
#'   this specifies the linear part.
#' @param cat_levels Named list. Each element is a character vector of factor
#'   levels for a categorical predictor. Default: \code{NULL} (no categorical
#'   predictors).
#' @param cat_effects Named list of numeric vectors. Effect of each level
#'   (same length and names as \code{cat_levels}). Default: random effects.
#' @param error_sd Positive numeric. Residual standard deviation. Default: 1.
#' @param predictor_means Named numeric vector or NULL. Means for continuous
#'   predictors. Default: all zeros.
#' @param predictor_sds Named numeric vector or NULL. SDs for continuous
#'   predictors. Default: all ones.
#' @param seed Integer or NULL. Random seed.
#'
#' @return A \code{\link{saqr_sim}} object with:
#'   \describe{
#'     \item{\code{$data}}{data.frame with columns \code{y}, continuous
#'       predictors, and categorical predictors (factors).}
#'     \item{\code{$params}}{list with \code{coefs}, \code{cat_effects},
#'       \code{error_sd}, \code{predictor_means}, \code{predictor_sds},
#'       \code{r_squared} (population \eqn{R^2}).}
#'   }
#'
#' @examples
#' r <- simulate_prediction(
#'   n = 200,
#'   coefs = c("(Intercept)" = 5, x1 = 2, x2 = -1),
#'   cat_levels = list(treatment = c("control", "drug_a", "drug_b")),
#'   cat_effects = list(treatment = c(0, 3, 5)),
#'   error_sd = 2, seed = 42
#' )
#' summary(lm(y ~ ., data = r$data))
#' r$params$r_squared
#'
#' @export
simulate_prediction <- function(n, coefs, cat_levels = NULL,
                                cat_effects = NULL, error_sd = 1,
                                predictor_means = NULL,
                                predictor_sds = NULL,
                                seed = NULL) {
  stopifnot(is.numeric(n), length(n) == 1L, n >= 5L)
  n <- as.integer(n)
  stopifnot(is.numeric(coefs), !is.null(names(coefs)))
  stopifnot(is.numeric(error_sd), length(error_sd) == 1L, error_sd > 0)

  pred_names <- setdiff(names(coefs), "(Intercept)")
  p <- length(pred_names)

  # Predictor means/sds
  if (is.null(predictor_means)) {
    predictor_means <- stats::setNames(rep(0, p), pred_names)
  }
  if (is.null(predictor_sds)) {
    predictor_sds <- stats::setNames(rep(1, p), pred_names)
  }
  stopifnot(length(predictor_means) == p, length(predictor_sds) == p,
            all(predictor_sds > 0))

  if (!is.null(seed)) set.seed(as.integer(seed))

  intercept <- if ("(Intercept)" %in% names(coefs)) coefs[["(Intercept)"]] else 0

  # Continuous predictors
  X <- vapply(pred_names, function(nm) {
    stats::rnorm(n, mean = predictor_means[[nm]], sd = predictor_sds[[nm]])
  }, numeric(n))
  if (p == 1L) X <- matrix(X, ncol = 1L, dimnames = list(NULL, pred_names))

  linear <- intercept + drop(X %*% coefs[pred_names])

  # Categorical predictors
  df <- as.data.frame(X)
  colnames(df) <- pred_names

  if (!is.null(cat_levels)) {
    stopifnot(is.list(cat_levels))
    if (is.null(cat_effects)) {
      cat_effects <- lapply(cat_levels, function(lvls) {
        stats::rnorm(length(lvls), sd = 2)
      })
    }
    stopifnot(is.list(cat_effects),
              identical(names(cat_effects), names(cat_levels)))

    for (cat_nm in names(cat_levels)) {
      lvls <- cat_levels[[cat_nm]]
      effs <- cat_effects[[cat_nm]]
      stopifnot(length(effs) == length(lvls))
      assignment <- sample(lvls, n, replace = TRUE)
      df[[cat_nm]] <- factor(assignment, levels = lvls)
      linear <- linear + effs[match(assignment, lvls)]
    }
  }

  noise <- stats::rnorm(n, sd = error_sd)
  signal_var <- stats::var(linear)
  r_squared <- signal_var / (signal_var + error_sd^2)

  df$y <- linear + noise
  df <- df[, c("y", setdiff(names(df), "y")), drop = FALSE]

  saqr_sim(
    data   = df,
    params = list(
      coefs           = coefs,
      cat_effects     = cat_effects,
      error_sd        = error_sd,
      predictor_means = predictor_means,
      predictor_sds   = predictor_sds,
      r_squared       = r_squared
    ),
    type = "prediction",
    seed = seed,
    call = match.call()
  )
}
