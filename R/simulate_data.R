#' Simulate Ready-to-Use Statistical Datasets
#'
#' @description Generate synthetic datasets suitable for common statistical
#'   analyses. Each dataset type contains real signal so the intended analysis
#'   works out of the box. Different seeds produce structurally different
#'   datasets (varying n, effect sizes, number of groups/variables).
#'
#' @param type Character string specifying the dataset type. One of:
#'   \describe{
#'     \item{`"ttest"`}{Two-group comparison. Columns: `group` (factor),
#'       `score` (numeric). Use with `t.test(score ~ group, data = d)`.}
#'     \item{`"anova"`}{Multi-group comparison (3--5 groups). Columns: `group`
#'       (factor), `score` (numeric). Use with `aov(score ~ group, data = d)`.}
#'     \item{`"correlation"`}{Correlated variables (4--7 columns). All numeric
#'       columns `x1`--`xp`. Use with `cor(d)` or `pairs(d)`.}
#'     \item{`"clusters"`}{Cluster structure (3--5 clusters, 2--5 dimensions).
#'       Numeric columns `x1`--`xd` plus integer `true_cluster`. Use with
#'       `kmeans(d[, -ncol(d)], centers = max(d$true_cluster))`.}
#'     \item{`"factor_analysis"`}{Latent factor structure (2--4 factors,
#'       6--16 observed variables). All numeric columns `x1`--`xp`. Use with
#'       `factanal(d, factors = attr(d, "n_factors"))` or `prcomp(d)`.
#'       Attributes `n_factors` and `loadings` store the true generating
#'       parameters.}
#'     \item{`"prediction"`}{Regression dataset. Columns: `y` (numeric),
#'       `x1`--`x3` (numeric linear), `x4` (numeric nonlinear), `cat1` and
#'       `cat2` (factors). Use with `lm(y ~ ., data = d)`.}
#'     \item{`"mlvar"`}{Multilevel VAR panel data. Columns: `id` (integer
#'       person ID), `day` (integer), `beep` (integer), `V1`--`Vd` (numeric).
#'       Use with `mlvar(d, vars = attr(d, "vars"), id = "id", day = "day",
#'       beep = "beep")`. Attributes `true_temporal` (d x d matrix),
#'       `true_contemporaneous` (d x d matrix), and `vars` (variable names)
#'       store the generating parameters.}
#'   }
#' @param seed Integer or NULL. Random seed for reproducibility. When non-NULL,
#'   the seed also determines dataset characteristics (sample size, effect
#'   sizes, number of groups/variables). Default: NULL.
#' @param ... Optional overrides for structural parameters. Any parameter not
#'   specified is drawn randomly (controlled by seed). Supported parameters
#'   by type:
#'   \describe{
#'     \item{All types}{`n` (integer, 10--500): number of rows.}
#'     \item{`"ttest"`}{`effect_size` (numeric, 0.3--1.2): Cohen's d.}
#'     \item{`"anova"`}{`n_groups` (integer, 3--5): number of groups.}
#'     \item{`"correlation"`}{`n_vars` (integer, 4--7): number of variables.}
#'     \item{`"clusters"`}{`n_clusters` (integer, 3--5): number of clusters;
#'       `n_dims` (integer, 2--5): number of dimensions.}
#'     \item{`"factor_analysis"`}{`n_factors` (integer, 2--4): number of
#'       latent factors; `items_per_factor` (integer, 3--4): observed
#'       variables per factor.}
#'     \item{`"prediction"`}{`n_cat1` (integer, 2--4): levels for cat1;
#'       `n_cat2` (integer, 2--3): levels for cat2.}
#'     \item{`"mlvar"`}{`n_subjects` (integer, 20--80): number of persons;
#'       `d` (integer, 3--6): number of variables;
#'       `n_days` (integer, 3--7): days per person;
#'       `beeps_per_day` (integer, 5--15): beeps per day.}
#'   }
#'
#' @return A `data.frame` with attributes:
#'   \describe{
#'     \item{`type`}{The dataset type string.}
#'     \item{`info`}{A usage hint string describing how to analyze the data.}
#'   }
#'   For `"factor_analysis"` type, additional attributes `n_factors` (integer)
#'   and `loadings` (matrix) store the true generating parameters.
#'
#' @details
#' Datasets have 10--500 rows depending on the seed. Signal strength varies
#' so that the intended analysis produces significant results for most seeds,
#' but specific outcomes depend on the random draw.
#'
#' No additional packages are required --- all generation uses base R and
#' stats functions only.
#'
#' @examples
#' # Two-group t-test
#' d <- simulate_data("ttest", seed = 42)
#' t.test(score ~ group, data = d)
#'
#' # One-way ANOVA
#' d <- simulate_data("anova", seed = 7)
#' summary(aov(score ~ group, data = d))
#'
#' # Correlation matrix
#' d <- simulate_data("correlation", seed = 1)
#' round(cor(d), 2)
#'
#' # Cluster analysis
#' d <- simulate_data("clusters", seed = 99)
#' km <- kmeans(d[, -ncol(d)], centers = max(d$true_cluster))
#' table(km$cluster, d$true_cluster)
#'
#' # Factor analysis
#' d <- simulate_data("factor_analysis", seed = 12)
#' factanal(d, factors = attr(d, "n_factors"))
#'
#' # Regression
#' d <- simulate_data("prediction", seed = 5)
#' summary(lm(y ~ ., data = d))
#'
#' # mlVAR panel data
#' d <- simulate_data("mlvar", seed = 1)
#' fit <- mlvar(d, vars = attr(d, "vars"), id = "id", day = "day", beep = "beep")
#'
#' # Override structural parameters
#' d <- simulate_data("clusters", seed = 1, n = 200, n_clusters = 4)
#' d <- simulate_data("factor_analysis", seed = 1, n_factors = 3)
#' d <- simulate_data("anova", n_groups = 4, n = 100)
#'
#' @seealso [simulate_sequences()], [simulate_long_data()],
#'   [simulate_igraph()]
#' @export
simulate_data <- function(type, seed = NULL, ...) {
  valid_types <- c(
    "ttest", "anova", "correlation", "clusters",
    "factor_analysis", "prediction", "mlvar"
  )
  stopifnot(
    is.character(type),
    length(type) == 1L,
    type %in% valid_types
  )
  if (!is.null(seed)) {
    stopifnot(is.numeric(seed), length(seed) == 1L, !is.na(seed))
    set.seed(as.integer(seed))
  }

  dots <- list(...)
  generator <- switch(type,
    ttest           = .simulate_ttest,
    anova           = .simulate_anova,
    correlation     = .simulate_correlation,
    clusters        = .simulate_clusters,
    factor_analysis = .simulate_factor_analysis,
    prediction      = .simulate_prediction,
    mlvar           = .simulate_mlvar
  )
  result <- generator(dots)

  attr(result, "type") <- type
  info_map <- c(
    ttest           = "t.test(score ~ group, data = d)",
    anova           = "summary(aov(score ~ group, data = d))",
    correlation     = "cor(d); pairs(d)",
    clusters        = "kmeans(d[, -ncol(d)], centers = max(d$true_cluster))",
    factor_analysis = "factanal(d, factors = attr(d, 'n_factors'))",
    prediction      = "summary(lm(y ~ ., data = d))",
    mlvar           = "mlvar(d, vars = attr(d, 'vars'), id = 'id', day = 'day', beep = 'beep')"
  )
  attr(result, "info") <- info_map[[type]]
  result
}


# ---- Internal generators ----

.simulate_ttest <- function(opts = list()) {
  n <- opts$n %||% sample(10:500, 1L)
  balance <- stats::runif(1L, 0.35, 0.65)
  n_a <- max(3L, round(n * balance))
  n_b <- n - n_a
  cohens_d <- opts$effect_size %||% stats::runif(1L, 0.3, 1.2)
  sd_a <- stats::runif(1L, 0.8, 1.5)
  sd_b <- stats::runif(1L, 0.8, 1.5)
  pooled_sd <- sqrt((sd_a^2 + sd_b^2) / 2)
  mean_diff <- cohens_d * pooled_sd
  mean_a <- stats::runif(1L, 40, 60)
  mean_b <- mean_a + mean_diff

  score_a <- stats::rnorm(n_a, mean = mean_a, sd = sd_a)
  score_b <- stats::rnorm(n_b, mean = mean_b, sd = sd_b)

  df <- data.frame(
    group = factor(c(rep("A", n_a), rep("B", n_b))),
    score = c(score_a, score_b)
  )
  df[sample.int(nrow(df)), ]
}


.simulate_anova <- function(opts = list()) {
  k <- opts$n_groups %||% sample(3:5, 1L)
  n <- opts$n %||% sample(10:500, 1L)
  # generate mildly imbalanced group sizes, min 5 each
  raw_sizes <- stats::runif(k, 0.5, 1.5)
  sizes <- pmax(5L, round(raw_sizes / sum(raw_sizes) * n))
  # adjust last group to hit target n

  sizes[k] <- max(5L, n - sum(sizes[-k]))
  actual_n <- sum(sizes)

  spread <- stats::runif(1L, 3, 10)
  base_mean <- stats::runif(1L, 30, 70)
  group_means <- base_mean + seq(0, spread, length.out = k) +
    stats::rnorm(k, sd = spread * 0.1)
  within_sd <- stats::runif(1L, 2, 5)

  groups <- rep(paste0("G", seq_len(k)), times = sizes)
  scores <- unlist(
    mapply(
      function(m, sz) stats::rnorm(sz, mean = m, sd = within_sd),
      group_means, sizes,
      SIMPLIFY = FALSE
    )
  )

  df <- data.frame(
    group = factor(groups, levels = paste0("G", seq_len(k))),
    score = scores
  )
  df[sample.int(nrow(df)), ]
}


.simulate_correlation <- function(opts = list()) {
  p <- opts$n_vars %||% sample(4:7, 1L)
  n <- opts$n %||% sample(10:500, 1L)

  # Build correlation matrix
  R <- diag(p)
  # Plant 2-3 strong correlations
  n_strong <- sample(2:3, 1L)
  pairs_pool <- which(upper.tri(R), arr.ind = TRUE)
  chosen <- pairs_pool[sample.int(nrow(pairs_pool), min(n_strong, nrow(pairs_pool))), , drop = FALSE]
  vapply(seq_len(nrow(chosen)), function(idx) {
    i <- chosen[idx, 1L]
    j <- chosen[idx, 2L]
    r <- stats::runif(1L, 0.5, 0.85) * sample(c(-1, 1), 1L)
    R[i, j] <<- r
    R[j, i] <<- r
    0
  }, numeric(1L))

  # Plant 1-2 moderate correlations among remaining
  remaining <- which(upper.tri(R) & R == 0, arr.ind = TRUE)
  if (nrow(remaining) > 0L) {
    n_moderate <- sample(1:min(2, nrow(remaining)), 1L)
    chosen_mod <- remaining[sample.int(nrow(remaining), n_moderate), , drop = FALSE]
    vapply(seq_len(nrow(chosen_mod)), function(idx) {
      i <- chosen_mod[idx, 1L]
      j <- chosen_mod[idx, 2L]
      r <- stats::runif(1L, 0.25, 0.49) * sample(c(-1, 1), 1L)
      R[i, j] <<- r
      R[j, i] <<- r
      0
    }, numeric(1L))
  }

  R <- .nearest_pd(R)
  L <- chol(R)
  Z <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
  X <- Z %*% L

  # Scale each variable to its own mean and SD
  means <- stats::runif(p, 10, 100)
  sds <- stats::runif(p, 1, 20)
  X <- sweep(X, 2L, sds, `*`)
  X <- sweep(X, 2L, means, `+`)

  df <- as.data.frame(X)
  colnames(df) <- paste0("x", seq_len(p))
  df
}


.simulate_clusters <- function(opts = list()) {
  k <- opts$n_clusters %||% sample(3:5, 1L)
  d <- opts$n_dims %||% sample(2:5, 1L)
  n <- opts$n %||% sample(10:500, 1L)

  # Imbalanced cluster sizes
  raw_sizes <- stats::runif(k, 0.3, 1.5)
  sizes <- pmax(3L, round(raw_sizes / sum(raw_sizes) * n))
  sizes[k] <- max(3L, n - sum(sizes[-k]))

  # Generate well-separated centers
  center_spread <- stats::runif(1L, 4, 10)
  centers <- matrix(stats::runif(k * d, -center_spread, center_spread),
    nrow = k, ncol = d
  )
  # Push centers apart to ensure minimum separation
  for (iter in seq_len(20L)) {
    dists <- as.matrix(stats::dist(centers))
    diag(dists) <- Inf
    min_dist <- min(dists)
    if (min_dist > 2) break
    pair <- which(dists == min_dist, arr.ind = TRUE)[1L, ]
    direction <- centers[pair[1L], ] - centers[pair[2L], ]
    norm_dir <- direction / max(sqrt(sum(direction^2)), 1e-8)
    push <- (2.5 - min_dist) / 2
    centers[pair[1L], ] <- centers[pair[1L], ] + norm_dir * push
    centers[pair[2L], ] <- centers[pair[2L], ] - norm_dir * push
  }

  within_sds <- stats::runif(k, 0.5, 1.5)
  data_list <- mapply(
    function(center, sz, sd_val) {
      matrix(stats::rnorm(sz * d, mean = rep(center, each = sz), sd = sd_val),
        nrow = sz, ncol = d
      )
    },
    split(centers, row(centers)),  # list of center vectors per cluster
    sizes, within_sds,
    SIMPLIFY = FALSE
  )

  # mapply with split(centers, row(centers)) gives named list of center components

  # Rebuild properly: generate per-cluster
  data_list <- lapply(seq_len(k), function(cl) {
    center <- centers[cl, ]
    sz <- sizes[cl]
    sd_val <- within_sds[cl]
    sweep(
      matrix(stats::rnorm(sz * d, sd = sd_val), nrow = sz, ncol = d),
      2L, center, `+`
    )
  })

  X <- do.call(rbind, data_list)
  true_cluster <- rep(seq_len(k), times = sizes)

  df <- as.data.frame(X)
  colnames(df) <- paste0("x", seq_len(d))
  df$true_cluster <- true_cluster

  shuffle <- sample.int(nrow(df))
  df[shuffle, ]
}


.simulate_factor_analysis <- function(opts = list()) {
  n_factors <- opts$n_factors %||% sample(2:4, 1L)
  items_per <- opts$items_per_factor %||% sample(3:4, 1L)
  p <- n_factors * items_per
  n <- opts$n %||% sample(max(p + 10L, 50L):500, 1L)

  # Build loading matrix
  Lambda <- matrix(stats::runif(p * n_factors, -0.15, 0.15),
    nrow = p, ncol = n_factors
  )
  # Set primary loadings
  vapply(seq_len(n_factors), function(f) {
    rows <- ((f - 1L) * items_per + 1L):(f * items_per)
    Lambda[rows, f] <<- stats::runif(items_per, 0.5, 0.85)
    0
  }, numeric(1L))

  # Inter-factor correlation matrix
  Phi <- diag(n_factors)
  if (n_factors > 1L) {
    phi_pairs <- which(upper.tri(Phi), arr.ind = TRUE)
    vapply(seq_len(nrow(phi_pairs)), function(idx) {
      r <- stats::runif(1L, 0, 0.35)
      Phi[phi_pairs[idx, 1L], phi_pairs[idx, 2L]] <<- r
      Phi[phi_pairs[idx, 2L], phi_pairs[idx, 1L]] <<- r
      0
    }, numeric(1L))
  }

  # Model-implied covariance: Lambda %*% Phi %*% t(Lambda) + Psi
  Sigma <- Lambda %*% Phi %*% t(Lambda)
  communalities <- diag(Sigma)
  Psi <- diag(pmax(0.1, 1 - communalities))
  Sigma <- Sigma + Psi
  Sigma <- .nearest_pd(Sigma)

  L <- chol(Sigma)
  Z <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
  X <- Z %*% L

  # Shift to Likert-like range
  means <- stats::runif(p, 2.5, 5.5)
  X <- sweep(X, 2L, means, `+`)

  df <- as.data.frame(X)
  colnames(df) <- paste0("x", seq_len(p))
  attr(df, "n_factors") <- n_factors
  attr(df, "loadings") <- Lambda
  df
}


.simulate_prediction <- function(opts = list()) {
  n <- opts$n %||% sample(10:500, 1L)

  # Numeric predictors
  x1 <- stats::rnorm(n, mean = stats::runif(1L, 0, 10), sd = stats::runif(1L, 1, 5))
  x2 <- stats::rnorm(n, mean = stats::runif(1L, 0, 10), sd = stats::runif(1L, 1, 5))
  x3 <- stats::rnorm(n, mean = stats::runif(1L, 0, 10), sd = stats::runif(1L, 1, 5))
  x4 <- stats::rnorm(n, mean = stats::runif(1L, 0, 5), sd = stats::runif(1L, 1, 3))

  # Categorical predictors
  n_cat1 <- opts$n_cat1 %||% sample(2:4, 1L)
  n_cat2 <- opts$n_cat2 %||% sample(2:3, 1L)
  cat1 <- factor(sample(paste0("L", seq_len(n_cat1)), n, replace = TRUE))
  cat2 <- factor(sample(paste0("M", seq_len(n_cat2)), n, replace = TRUE))

  # Coefficients
  b1 <- stats::runif(1L, -3, 3)
  b2 <- stats::runif(1L, -3, 3)
  b3 <- stats::runif(1L, -3, 3) * sample(c(1, 2), 1L)  # one possibly amplified

  b4_quad <- stats::runif(1L, 0.3, 1.5) * sample(c(-1, 1), 1L)

  # Categorical effects
  cat1_effects <- stats::rnorm(n_cat1, sd = 2)
  cat2_effects <- stats::rnorm(n_cat2, sd = 1.5)

  linear <- b1 * x1 + b2 * x2 + b3 * x3 + b4_quad * x4^2 +
    cat1_effects[as.integer(cat1)] + cat2_effects[as.integer(cat2)]

  # Noise calibrated for R^2 ~ 0.4-0.8
  signal_sd <- stats::sd(linear)
  noise_ratio <- stats::runif(1L, 0.5, 1.2)  # noise relative to signal
  noise_sd <- signal_sd * noise_ratio
  y <- linear + stats::rnorm(n, sd = noise_sd)

  data.frame(
    y = y,
    x1 = x1, x2 = x2, x3 = x3, x4 = x4,
    cat1 = cat1, cat2 = cat2
  )
}


.simulate_mlvar <- function(opts = list()) {
  # Structural parameters (seed-driven or overridden)
  n_subjects <- opts$n_subjects %||% sample(20:80, 1L)
  d <- opts$d %||% sample(3:6, 1L)
  n_days <- opts$n_days %||% sample(3:7, 1L)
  beeps_per_day <- opts$beeps_per_day %||% sample(5:15, 1L)
  obs_per_subject <- n_days * beeps_per_day

  var_names <- paste0("V", seq_len(d))

  # ---- True temporal B matrix ----
  # Diagonal: autoregressive (0.2-0.5), off-diagonal: sparse +/- 0.3
  B <- diag(stats::runif(d, 0.2, 0.5))
  # Add 1-3 sparse off-diagonal effects
  off_diag <- which(!diag(TRUE, d), arr.ind = TRUE)
  n_off <- min(sample(1:3, 1L), nrow(off_diag))
  chosen_off <- off_diag[sample.int(nrow(off_diag), n_off), , drop = FALSE]
  vapply(seq_len(nrow(chosen_off)), function(idx) {
    B[chosen_off[idx, 1L], chosen_off[idx, 2L]] <<-
      stats::runif(1L, 0.1, 0.3) * sample(c(-1, 1), 1L)
    0
  }, numeric(1L))
  colnames(B) <- rownames(B) <- var_names

  # ---- True contemporaneous noise covariance ----
  # Sparse partial correlation structure for innovations
  noise_cor <- diag(d)
  n_noise_edges <- min(sample(1:2, 1L), d * (d - 1L) / 2L)
  noise_pairs <- which(upper.tri(noise_cor), arr.ind = TRUE)
  chosen_noise <- noise_pairs[
    sample.int(nrow(noise_pairs), n_noise_edges), , drop = FALSE
  ]
  vapply(seq_len(nrow(chosen_noise)), function(idx) {
    r <- stats::runif(1L, 0.2, 0.5) * sample(c(-1, 1), 1L)
    noise_cor[chosen_noise[idx, 1L], chosen_noise[idx, 2L]] <<- r
    noise_cor[chosen_noise[idx, 2L], chosen_noise[idx, 1L]] <<- r
    0
  }, numeric(1L))
  noise_cov <- .nearest_pd(noise_cor)
  L_noise <- chol(noise_cov)

  # ---- Person-level means ----
  person_means <- matrix(stats::rnorm(n_subjects * d, mean = 0, sd = 1.5),
                         nrow = n_subjects, ncol = d)

  # ---- Simulate VAR(1) per person ----
  all_rows <- vector("list", n_subjects)

  for (subj in seq_len(n_subjects)) {
    mu <- person_means[subj, ]
    # Initialize from person mean + noise
    y_prev <- mu + stats::rnorm(d, sd = 0.5)

    subj_data <- matrix(NA_real_, nrow = obs_per_subject, ncol = d)
    subj_data[1L, ] <- y_prev

    for (t in seq(2L, obs_per_subject)) {
      # VAR(1): y_t = mu + B %*% (y_{t-1} - mu) + noise
      innovation <- as.vector(stats::rnorm(d) %*% L_noise) * 0.5
      y_new <- mu + as.vector(B %*% (y_prev - mu)) + innovation
      subj_data[t, ] <- y_new
      y_prev <- y_new
    }

    # Build data frame rows
    day_vec <- rep(seq_len(n_days), each = beeps_per_day)
    beep_vec <- rep(seq_len(beeps_per_day), times = n_days)

    subj_df <- data.frame(
      id = rep(subj, obs_per_subject),
      day = day_vec,
      beep = beep_vec
    )
    subj_df[var_names] <- subj_data
    all_rows[[subj]] <- subj_df
  }

  df <- do.call(rbind, all_rows)
  rownames(df) <- NULL

  attr(df, "true_temporal") <- B
  attr(df, "true_contemporaneous") <- noise_cor
  attr(df, "vars") <- var_names
  df
}


# ---- Helper ----

#' Nearest positive-definite matrix projection
#'
#' Clamps eigenvalues to a minimum threshold and reconstructs the matrix
#' with unit diagonal. Used internally by correlation and factor analysis
#' generators.
#'
#' @param R A symmetric matrix.
#' @return A positive-definite symmetric matrix with unit diagonal (if input
#'   had unit diagonal) or normalized to unit diagonal.
#' @keywords internal
#' @noRd
.nearest_pd <- function(R) {
  eig <- eigen(R, symmetric = TRUE)
  vals <- pmax(eig$values, 1e-6)
  R_pd <- eig$vectors %*% diag(vals) %*% t(eig$vectors)
  # Normalize to unit diagonal
  d <- sqrt(diag(R_pd))
  R_pd <- R_pd / outer(d, d)
  # Force exact symmetry
  (R_pd + t(R_pd)) / 2
}
