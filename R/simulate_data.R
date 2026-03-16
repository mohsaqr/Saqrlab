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
#'       `score` (numeric).}
#'     \item{`"anova"`}{Multi-group comparison (3--5 groups). Columns: `group`
#'       (factor), `score` (numeric).}
#'     \item{`"correlation"`}{Correlated variables (4--7 columns). All numeric
#'       `x1`--`xp`.}
#'     \item{`"clusters"`}{Cluster structure. Numeric `x1`--`xd` plus integer
#'       `true_cluster`.}
#'     \item{`"factor_analysis"`}{Latent factor structure. All numeric
#'       `x1`--`xp`. Attributes `n_factors` and `loadings`.}
#'     \item{`"prediction"`}{Regression dataset. Columns: `y`, `x1`--`x4`,
#'       `cat1`, `cat2`.}
#'     \item{`"mlvar"`}{Multilevel VAR panel data. Columns: `id`, `day`,
#'       `beep`, `V1`--`Vd`. Attributes `true_temporal`, `true_contemporaneous`,
#'       `vars`.}
#'     \item{`"batch"`}{Wildcard: generates all 7 types. Returns a named list
#'       with `n_batch` (default 1000) datasets per type.}
#'   }
#' @param seed Integer or NULL. Random seed. Also determines structural
#'   parameters (n, effect sizes, etc.). Default: NULL.
#' @param n_batch Integer or NULL. When provided, returns a list of `n_batch`
#'   datasets instead of a single dataset. When `type = "batch"`, defaults to
#'   1000L. Each item has `seed`, `batch_id`, and `complexity` attributes.
#' @param complexity Character. Controls edge-case injection for stress testing.
#'   \describe{
#'     \item{`"clean"`}{(default for single datasets) No edge cases.}
#'     \item{`"auto"`}{(default in batch mode) Randomly injects 0--3 edge
#'       cases per dataset (seed-driven, fully reproducible).}
#'     \item{character vector}{Inject specific cases. One or more of:
#'       \code{"na"}, \code{"outliers"}, \code{"ties"}, \code{"duplicates"},
#'       \code{"constant_col"}, \code{"all_na_col"}, \code{"tiny_n"},
#'       \code{"heavy_tailed"}, \code{"heteroscedastic"},
#'       \code{"extreme_imbalance"}, \code{"multicollinear"}.}
#'   }
#' @param ... Optional overrides for structural parameters (n, n_groups, etc.).
#'
#' @return A `data.frame` (single dataset) or a list (when `n_batch` is
#'   non-NULL or `type = "batch"`). Each dataset has attributes `type`,
#'   `info`, and `complexity`. Batch items additionally have `seed` and
#'   `batch_id`.
#'
#' @examples
#' # Single clean dataset
#' d <- simulate_data("ttest", seed = 42)
#'
#' # With edge cases
#' d <- simulate_data("correlation", seed = 1, complexity = "auto")
#' d <- simulate_data("ttest", seed = 1, complexity = c("na", "outliers"))
#'
#' # Batch of 100 ttest datasets
#' batch <- simulate_data("ttest", seed = 1, n_batch = 100)
#'
#' # All types, 1000 datasets each
#' all_batches <- simulate_data("batch", seed = 1)
#'
#' @seealso [simulate_sequences()], [simulate_igraph()]
#' @export
simulate_data <- function(type, seed = NULL, complexity = "clean",
                          ..., n_batch = NULL) {
  valid_types <- c(.SIM_VALID_TYPES, "batch")
  stopifnot(
    is.character(type), length(type) == 1L, type %in% valid_types
  )
  if (!is.null(seed)) {
    stopifnot(is.numeric(seed), length(seed) == 1L, !is.na(seed))
  }
  if (!is.null(n_batch)) {
    stopifnot(is.numeric(n_batch), length(n_batch) == 1L,
              !is.na(n_batch), n_batch >= 1L)
    n_batch <- as.integer(n_batch)
  }
  if (!is.character(complexity) || length(complexity) == 0L) {
    stop("`complexity` must be 'clean', 'auto', or a character vector of edge case names")
  }
  if (!identical(complexity, "clean") && !identical(complexity, "auto")) {
    bad <- setdiff(complexity, .SIM_ALL_CASES)
    if (length(bad) > 0L) {
      stop(sprintf("Unknown complexity case(s): %s. Valid: %s",
                   paste(bad, collapse = ", "),
                   paste(.SIM_ALL_CASES, collapse = ", ")))
    }
  }

  base_seed <- if (!is.null(seed)) as.integer(seed) else 42L
  dots <- list(...)

  # ---- batch mode: type = "batch" ----
  if (type == "batch") {
    nb <- n_batch %||% 1000L
    eff_complexity <- if (identical(complexity, "clean")) "auto" else complexity
    result <- lapply(
      stats::setNames(.SIM_VALID_TYPES, .SIM_VALID_TYPES),
      function(tp) {
        type_idx <- which(.SIM_VALID_TYPES == tp)
        .sim_batch_type(tp, type_idx, base_seed, nb, eff_complexity, dots)
      }
    )
    class(result) <- c("sim_batch_full", "list")
    attr(result, "n_batch")    <- nb
    attr(result, "seed")       <- base_seed
    attr(result, "complexity") <- eff_complexity
    return(result)
  }

  # ---- batch mode: n_batch provided for a single type ----
  if (!is.null(n_batch)) {
    type_idx <- which(.SIM_VALID_TYPES == type)
    eff_complexity <- if (identical(complexity, "clean")) "auto" else complexity
    result <- .sim_batch_type(type, type_idx, base_seed, n_batch,
                              eff_complexity, dots)
    class(result) <- c("sim_batch_type", "list")
    attr(result, "type")       <- type
    attr(result, "n_batch")    <- n_batch
    attr(result, "seed")       <- base_seed
    attr(result, "complexity") <- eff_complexity
    return(result)
  }

  # ---- single dataset ----
  if (!is.null(seed)) set.seed(base_seed)
  cases <- if (identical(complexity, "auto")) {
    .sample_auto_complexity()
  } else if (identical(complexity, "clean")) {
    character(0L)
  } else {
    complexity
  }

  result <- .sim_one_dataset(type, cases, dots)
  attr(result, "type")       <- type
  attr(result, "info")       <- .SIM_INFO_MAP[[type]]
  attr(result, "complexity") <- cases
  result
}


# ===========================================================================
# Batch helpers
# ===========================================================================

.sim_batch_type <- function(type, type_idx, base_seed, n_batch,
                             complexity, dots) {
  lapply(seq_len(n_batch), function(i) {
    item_seed <- .sim_item_seed(base_seed, type_idx, i)
    set.seed(item_seed)
    cases <- if (identical(complexity, "auto")) {
      .sample_auto_complexity()
    } else if (identical(complexity, "clean")) {
      character(0L)
    } else {
      complexity
    }
    result <- .sim_one_dataset(type, cases, dots)
    attr(result, "type")       <- type
    attr(result, "info")       <- .SIM_INFO_MAP[[type]]
    attr(result, "complexity") <- cases
    attr(result, "seed")       <- item_seed
    attr(result, "batch_id")   <- i
    result
  })
}

# Deterministic per-item seed hash: same (base_seed, type_idx, batch_idx) always
# maps to the same seed regardless of n_batch.
.sim_item_seed <- function(base_seed, type_idx, batch_idx) {
  as.integer(
    (as.numeric(base_seed) * 7919 + (type_idx - 1L) * 997 + batch_idx) %%
      .Machine$integer.max
  )
}

.sim_one_dataset <- function(type, cases, dots) {
  opts <- c(dots, list(complexity = cases))
  generator <- switch(type,
    ttest           = .simulate_ttest,
    anova           = .simulate_anova,
    correlation     = .simulate_correlation,
    clusters        = .simulate_clusters,
    factor_analysis = .simulate_factor_analysis,
    prediction      = .simulate_prediction,
    mlvar           = .simulate_mlvar
  )
  result <- generator(opts)
  .inject_complexity_post(result, cases)
}


# ===========================================================================
# Auto-complexity sampler (Gumbel-max trick — weighted, no-replacement, vectorized)
# ===========================================================================

.sample_auto_complexity <- function() {
  n_cases <- sample(0L:3L, 1L, prob = c(0.25, 0.35, 0.25, 0.15))
  if (n_cases == 0L) return(character(0L))
  pool    <- c("na", "outliers", "ties", "tiny_n", "duplicates",
               "heavy_tailed", "heteroscedastic", "extreme_imbalance",
               "constant_col", "multicollinear", "all_na_col")
  weights <- c(0.18, 0.18, 0.14, 0.12, 0.10, 0.10, 0.07, 0.05, 0.04, 0.04, 0.02)
  scores  <- -log(stats::runif(length(pool))) / weights
  pool[order(scores)[seq_len(min(n_cases, length(pool)))]]
}


# ===========================================================================
# Post-injection edge cases (universal, applied after generation)
# ===========================================================================

.inject_complexity_post <- function(df, cases) {
  if (length(cases) == 0L) return(df)

  num_cols <- names(df)[vapply(df, is.numeric, logical(1L))]
  n <- nrow(df)

  # NA injection: 5-15% in 50-100% of numeric columns
  if ("na" %in% cases && length(num_cols) > 0L) {
    na_rate    <- stats::runif(1L, 0.05, 0.15)
    n_target   <- max(1L, round(length(num_cols) * stats::runif(1L, 0.5, 1.0)))
    target_cols <- sample(num_cols, min(n_target, length(num_cols)))
    df[target_cols] <- lapply(df[target_cols], function(v) {
      n_na <- max(1L, round(length(v) * na_rate))
      v[sample.int(length(v), n_na)] <- NA_real_
      v
    })
  }

  # Outliers: 1-5 extreme values (3-8 SD from mean)
  if ("outliers" %in% cases && length(num_cols) > 0L) {
    n_out    <- sample(1L:5L, 1L)
    out_cols <- sample(num_cols, n_out, replace = TRUE)
    out_rows <- sample.int(n, n_out, replace = TRUE)
    unique_cols <- unique(out_cols)
    modified <- lapply(stats::setNames(unique_cols, unique_cols), function(col) {
      v    <- df[[col]]
      rows <- out_rows[out_cols == col]
      vals <- v[!is.na(v)]
      if (length(vals) < 2L) return(v)
      m <- mean(vals)
      s <- stats::sd(vals)
      if (s == 0) return(v)
      v[rows] <- m + sample(c(-1L, 1L), length(rows), replace = TRUE) *
                     stats::runif(length(rows), 3, 8) * s
      v
    })
    df[unique_cols] <- modified
  }

  # Ties: discretize 50% of numeric columns to ~25-35% unique values
  if ("ties" %in% cases && length(num_cols) > 0L) {
    n_cols      <- max(1L, round(length(num_cols) * 0.5))
    target_cols <- sample(num_cols, n_cols)
    df[target_cols] <- lapply(df[target_cols], function(v) {
      valid <- v[!is.na(v)]
      if (length(valid) < 4L) return(v)
      n_bins <- max(3L, round(length(valid) * stats::runif(1L, 0.25, 0.35)))
      rng    <- range(valid)
      breaks <- seq(rng[1L], rng[2L], length.out = n_bins + 1L)
      mids   <- (breaks[-length(breaks)] + breaks[-1L]) / 2L
      idx    <- findInterval(v, breaks, rightmost.closed = TRUE)
      idx    <- pmin(pmax(idx, 1L), n_bins)
      ifelse(is.na(v), NA_real_, mids[idx])
    })
  }

  # Duplicates: append 10-25% duplicate rows
  if ("duplicates" %in% cases && n > 5L) {
    n_dup    <- max(1L, round(n * stats::runif(1L, 0.10, 0.25)))
    dup_rows <- df[sample.int(n, n_dup, replace = TRUE), ]
    df       <- rbind(df, dup_rows)
    rownames(df) <- NULL
  }

  # Constant column: set one secondary numeric col to its mean
  if ("constant_col" %in% cases && length(num_cols) >= 2L) {
    col       <- sample(num_cols[-1L], 1L)
    df[[col]] <- mean(df[[col]], na.rm = TRUE)
  }

  # All-NA column: set one secondary numeric col to NA
  if ("all_na_col" %in% cases && length(num_cols) >= 2L) {
    cands <- num_cols[-1L]
    if ("constant_col" %in% cases && length(num_cols) >= 3L) {
      cands <- num_cols[-(1:2)]
    }
    if (length(cands) > 0L) df[[sample(cands, 1L)]] <- NA_real_
  }

  df
}


# ===========================================================================
# Generators (each reads opts$complexity for generation-time edge cases)
# ===========================================================================

.simulate_ttest <- function(opts = list()) {
  complexity <- opts$complexity %||% character(0L)

  n <- if ("tiny_n" %in% complexity) {
    max(6L, sample(6L:20L, 1L))
  } else {
    opts$n %||% sample(30L:500L, 1L)
  }
  n <- opts$n %||% n

  balance <- if ("extreme_imbalance" %in% complexity) {
    stats::runif(1L, 0.03, 0.10)
  } else {
    stats::runif(1L, 0.35, 0.65)
  }
  n_a <- max(2L, min(n - 2L, round(n * balance)))
  n_b <- max(2L, n - n_a)

  cohens_d  <- opts$effect_size %||% stats::runif(1L, 0.3, 1.2)
  sd_a      <- stats::runif(1L, 0.8, 1.5)
  sd_b_base <- stats::runif(1L, 0.8, 1.5)
  sd_b      <- if ("heteroscedastic" %in% complexity) {
    sd_b_base * stats::runif(1L, 3, 8)
  } else {
    sd_b_base
  }

  pooled_sd <- sqrt((sd_a^2 + sd_b^2) / 2)
  mean_a    <- stats::runif(1L, 40, 60)
  mean_b    <- mean_a + cohens_d * pooled_sd

  if ("heavy_tailed" %in% complexity) {
    df_t    <- sample(2L:5L, 1L)
    score_a <- mean_a + sd_a * stats::rt(n_a, df = df_t)
    score_b <- mean_b + sd_b * stats::rt(n_b, df = df_t)
  } else {
    score_a <- stats::rnorm(n_a, mean = mean_a, sd = sd_a)
    score_b <- stats::rnorm(n_b, mean = mean_b, sd = sd_b)
  }

  df <- data.frame(
    group = factor(c(rep("A", n_a), rep("B", n_b))),
    score = c(score_a, score_b)
  )
  df[sample.int(nrow(df)), ]
}


.simulate_anova <- function(opts = list()) {
  complexity <- opts$complexity %||% character(0L)

  k <- opts$n_groups %||% sample(3L:5L, 1L)
  n <- if ("tiny_n" %in% complexity) {
    k * max(3L, sample(3L:5L, 1L))
  } else {
    opts$n %||% sample(30L:500L, 1L)
  }
  n <- opts$n %||% n

  # Group sizes: extreme imbalance gives one dominant group
  sizes <- if ("extreme_imbalance" %in% complexity) {
    dominant <- round(n * stats::runif(1L, 0.70, 0.90))
    remainder <- max(k - 1L, 1L)
    rest <- pmax(2L, round(rep((n - dominant) / remainder, remainder)))
    rest[remainder] <- max(2L, n - dominant - sum(rest[-remainder]))
    c(dominant, rest)
  } else {
    raw <- stats::runif(k, 0.5, 1.5)
    s   <- pmax(5L, round(raw / sum(raw) * n))
    s[k] <- max(5L, n - sum(s[-k]))
    s
  }
  sizes <- sizes[seq_len(k)]

  spread     <- stats::runif(1L, 3, 10)
  base_mean  <- stats::runif(1L, 30, 70)
  group_means <- base_mean + seq(0, spread, length.out = k) +
                 stats::rnorm(k, sd = spread * 0.1)

  if ("heteroscedastic" %in% complexity) {
    within_sds <- stats::runif(k, 1, 8)
  } else {
    within_sds <- rep(stats::runif(1L, 2, 5), k)
  }

  make_group <- function(m, sz, s) {
    if ("heavy_tailed" %in% complexity) {
      df_t <- sample(2L:5L, 1L)
      m + s * stats::rt(sz, df = df_t)
    } else {
      stats::rnorm(sz, mean = m, sd = s)
    }
  }

  scores <- unlist(
    mapply(make_group, group_means, sizes, within_sds, SIMPLIFY = FALSE)
  )
  groups <- rep(paste0("G", seq_len(k)), times = sizes)

  df <- data.frame(
    group = factor(groups, levels = paste0("G", seq_len(k))),
    score = scores
  )
  df[sample.int(nrow(df)), ]
}


.simulate_correlation <- function(opts = list()) {
  complexity <- opts$complexity %||% character(0L)

  p <- opts$n_vars %||% sample(4L:7L, 1L)
  n <- if ("tiny_n" %in% complexity) {
    p + sample(2L:10L, 1L)
  } else {
    opts$n %||% sample(30L:500L, 1L)
  }
  n <- opts$n %||% n

  # Build correlation matrix
  R <- diag(p)
  n_strong <- sample(2:3, 1L)
  pairs_pool <- which(upper.tri(R), arr.ind = TRUE)
  chosen <- pairs_pool[sample.int(nrow(pairs_pool),
                                  min(n_strong, nrow(pairs_pool))), , drop = FALSE]
  vapply(seq_len(nrow(chosen)), function(idx) {
    i <- chosen[idx, 1L]; j <- chosen[idx, 2L]
    r <- stats::runif(1L, 0.5, 0.85) * sample(c(-1, 1), 1L)
    R[i, j] <<- r; R[j, i] <<- r
    0
  }, numeric(1L))

  remaining <- which(upper.tri(R) & R == 0, arr.ind = TRUE)
  if (nrow(remaining) > 0L) {
    n_mod    <- sample(1:min(2, nrow(remaining)), 1L)
    chosen_m <- remaining[sample.int(nrow(remaining), n_mod), , drop = FALSE]
    vapply(seq_len(nrow(chosen_m)), function(idx) {
      i <- chosen_m[idx, 1L]; j <- chosen_m[idx, 2L]
      r <- stats::runif(1L, 0.25, 0.49) * sample(c(-1, 1), 1L)
      R[i, j] <<- r; R[j, i] <<- r
      0
    }, numeric(1L))
  }

  R <- .nearest_pd(R)
  L <- chol(R)

  if ("heavy_tailed" %in% complexity) {
    df_t <- sample(2L:5L, 1L)
    Z <- matrix(stats::rt(n * p, df = df_t), nrow = n, ncol = p)
  } else {
    Z <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
  }
  X <- Z %*% L

  means <- stats::runif(p, 10, 100)
  sds   <- stats::runif(p, 1, 20)
  X     <- sweep(sweep(X, 2L, sds, `*`), 2L, means, `+`)

  # Multicollinear: add extra column nearly identical to x1
  if ("multicollinear" %in% complexity) {
    mc_col <- X[, 1L] + stats::rnorm(n, sd = stats::sd(X[, 1L]) * 0.05)
    X <- cbind(X, mc_col)
    p <- p + 1L
  }

  df <- as.data.frame(X)
  colnames(df) <- paste0("x", seq_len(p))
  df
}


.simulate_clusters <- function(opts = list()) {
  complexity <- opts$complexity %||% character(0L)

  k <- opts$n_clusters %||% sample(3L:5L, 1L)
  d <- opts$n_dims    %||% sample(2L:5L, 1L)
  n <- if ("tiny_n" %in% complexity) {
    k * 3L
  } else {
    opts$n %||% sample(30L:500L, 1L)
  }
  n <- opts$n %||% n

  raw_sizes <- stats::runif(k, 0.3, 1.5)
  sizes     <- pmax(3L, round(raw_sizes / sum(raw_sizes) * n))
  sizes[k]  <- max(3L, n - sum(sizes[-k]))

  center_spread <- stats::runif(1L, 4, 10)
  centers <- matrix(stats::runif(k * d, -center_spread, center_spread),
                    nrow = k, ncol = d)
  for (iter in seq_len(20L)) {
    dists    <- as.matrix(stats::dist(centers)); diag(dists) <- Inf
    min_dist <- min(dists)
    if (min_dist > 2) break
    pair      <- which(dists == min_dist, arr.ind = TRUE)[1L, ]
    direction <- centers[pair[1L], ] - centers[pair[2L], ]
    norm_dir  <- direction / max(sqrt(sum(direction^2)), 1e-8)
    push      <- (2.5 - min_dist) / 2
    centers[pair[1L], ] <- centers[pair[1L], ] + norm_dir * push
    centers[pair[2L], ] <- centers[pair[2L], ] - norm_dir * push
  }

  within_sds <- stats::runif(k, 0.5, 1.5)
  data_list  <- lapply(seq_len(k), function(cl) {
    center  <- centers[cl, ]
    sz      <- sizes[cl]
    sd_val  <- within_sds[cl]
    if ("heavy_tailed" %in% complexity) {
      df_t <- sample(2L:5L, 1L)
      sweep(matrix(sd_val * stats::rt(sz * d, df = df_t), nrow = sz, ncol = d),
            2L, center, `+`)
    } else {
      sweep(matrix(stats::rnorm(sz * d, sd = sd_val), nrow = sz, ncol = d),
            2L, center, `+`)
    }
  })

  X            <- do.call(rbind, data_list)
  true_cluster <- rep(seq_len(k), times = sizes)

  df            <- as.data.frame(X)
  colnames(df)  <- paste0("x", seq_len(d))
  df$true_cluster <- true_cluster
  df[sample.int(nrow(df)), ]
}


.simulate_factor_analysis <- function(opts = list()) {
  complexity <- opts$complexity %||% character(0L)

  n_factors <- opts$n_factors      %||% sample(2L:4L, 1L)
  items_per <- opts$items_per_factor %||% sample(3L:4L, 1L)
  p         <- n_factors * items_per
  n <- if ("tiny_n" %in% complexity) {
    p + sample(3L:8L, 1L)
  } else {
    opts$n %||% sample(max(p + 10L, 50L):500L, 1L)
  }
  n <- opts$n %||% n

  Lambda <- matrix(stats::runif(p * n_factors, -0.15, 0.15),
                   nrow = p, ncol = n_factors)
  vapply(seq_len(n_factors), function(f) {
    rows <- ((f - 1L) * items_per + 1L):(f * items_per)
    Lambda[rows, f] <<- stats::runif(items_per, 0.5, 0.85)
    0
  }, numeric(1L))

  # Multicollinear: strong cross-loading for factor 1 items on factor 2
  if ("multicollinear" %in% complexity && n_factors >= 2L) {
    rows_f1 <- 1L:items_per
    Lambda[rows_f1, 2L] <- stats::runif(items_per, 0.4, 0.7)
  }

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

  Sigma        <- Lambda %*% Phi %*% t(Lambda)
  communalities <- diag(Sigma)
  Psi          <- diag(pmax(0.1, 1 - communalities))
  Sigma        <- .nearest_pd(Sigma + Psi)
  L            <- chol(Sigma)

  if ("heavy_tailed" %in% complexity) {
    df_t <- sample(2L:5L, 1L)
    Z <- matrix(stats::rt(n * p, df = df_t), nrow = n, ncol = p)
  } else {
    Z <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
  }
  X    <- Z %*% L
  means <- stats::runif(p, 2.5, 5.5)
  X    <- sweep(X, 2L, means, `+`)

  df <- as.data.frame(X)
  colnames(df) <- paste0("x", seq_len(p))
  attr(df, "n_factors") <- n_factors
  attr(df, "loadings")  <- Lambda
  df
}


.simulate_prediction <- function(opts = list()) {
  complexity <- opts$complexity %||% character(0L)

  n <- if ("tiny_n" %in% complexity) {
    sample(10L:25L, 1L)
  } else {
    opts$n %||% sample(30L:500L, 1L)
  }
  n <- opts$n %||% n

  mu_x  <- stats::runif(4L, 0, 10)
  sd_x  <- stats::runif(4L, 1, 5)

  if ("heavy_tailed" %in% complexity) {
    df_t <- sample(2L:5L, 1L)
    x1 <- mu_x[1L] + sd_x[1L] * stats::rt(n, df = df_t)
    x2 <- mu_x[2L] + sd_x[2L] * stats::rt(n, df = df_t)
    x3 <- mu_x[3L] + sd_x[3L] * stats::rt(n, df = df_t)
    x4 <- mu_x[4L] + sd_x[4L] * stats::rt(n, df = df_t)
  } else {
    x1 <- stats::rnorm(n, mean = mu_x[1L], sd = sd_x[1L])
    x2 <- stats::rnorm(n, mean = mu_x[2L], sd = sd_x[2L])
    x3 <- stats::rnorm(n, mean = mu_x[3L], sd = sd_x[3L])
    x4 <- stats::rnorm(n, mean = mu_x[4L], sd = sd_x[4L])
  }

  # Multicollinear: x2 becomes near-duplicate of x1
  if ("multicollinear" %in% complexity) {
    x2 <- x1 + stats::rnorm(n, sd = stats::sd(x1) * 0.05)
  }

  n_cat1 <- opts$n_cat1 %||% sample(2L:4L, 1L)
  n_cat2 <- opts$n_cat2 %||% sample(2L:3L, 1L)
  cat1   <- factor(sample(paste0("L", seq_len(n_cat1)), n, replace = TRUE))
  cat2   <- factor(sample(paste0("M", seq_len(n_cat2)), n, replace = TRUE))

  b1 <- stats::runif(1L, -3, 3)
  b2 <- stats::runif(1L, -3, 3)
  b3 <- stats::runif(1L, -3, 3) * sample(c(1, 2), 1L)
  b4 <- stats::runif(1L, 0.3, 1.5) * sample(c(-1, 1), 1L)
  cat1_eff <- stats::rnorm(n_cat1, sd = 2)
  cat2_eff <- stats::rnorm(n_cat2, sd = 1.5)

  linear    <- b1 * x1 + b2 * x2 + b3 * x3 + b4 * x4^2 +
               cat1_eff[as.integer(cat1)] + cat2_eff[as.integer(cat2)]
  signal_sd <- stats::sd(linear)
  noise_sd  <- signal_sd * stats::runif(1L, 0.5, 1.2)

  noise <- if ("heavy_tailed" %in% complexity) {
    df_t <- sample(2L:5L, 1L)
    stats::rt(n, df = df_t) * noise_sd
  } else {
    stats::rnorm(n, sd = noise_sd)
  }

  data.frame(y = linear + noise, x1 = x1, x2 = x2, x3 = x3, x4 = x4,
             cat1 = cat1, cat2 = cat2)
}


.simulate_mlvar <- function(opts = list()) {
  complexity   <- opts$complexity %||% character(0L)

  n_subjects   <- opts$n_subjects %||% sample(20L:80L, 1L)
  d            <- opts$d %||% sample(3L:6L, 1L)

  if ("tiny_n" %in% complexity) {
    n_subjects    <- opts$n_subjects %||% sample(3L:8L, 1L)
    n_days        <- opts$n_days %||% sample(2L:3L, 1L)
    beeps_per_day <- opts$beeps_per_day %||% sample(3L:5L, 1L)
  } else {
    n_days        <- opts$n_days %||% sample(3L:7L, 1L)
    beeps_per_day <- opts$beeps_per_day %||% sample(5L:15L, 1L)
  }
  obs_per_subject <- n_days * beeps_per_day
  var_names       <- paste0("V", seq_len(d))

  B <- diag(stats::runif(d, 0.2, 0.5))
  off_diag  <- which(!diag(TRUE, d), arr.ind = TRUE)
  n_off     <- min(sample(1L:3L, 1L), nrow(off_diag))
  chosen_off <- off_diag[sample.int(nrow(off_diag), n_off), , drop = FALSE]
  vapply(seq_len(nrow(chosen_off)), function(idx) {
    B[chosen_off[idx, 1L], chosen_off[idx, 2L]] <<-
      stats::runif(1L, 0.1, 0.3) * sample(c(-1, 1), 1L)
    0
  }, numeric(1L))
  colnames(B) <- rownames(B) <- var_names

  noise_cor <- diag(d)
  n_ne      <- min(sample(1L:2L, 1L), d * (d - 1L) / 2L)
  np        <- which(upper.tri(noise_cor), arr.ind = TRUE)
  c_noise   <- np[sample.int(nrow(np), n_ne), , drop = FALSE]
  vapply(seq_len(nrow(c_noise)), function(idx) {
    r <- stats::runif(1L, 0.2, 0.5) * sample(c(-1, 1), 1L)
    noise_cor[c_noise[idx, 1L], c_noise[idx, 2L]] <<- r
    noise_cor[c_noise[idx, 2L], c_noise[idx, 1L]] <<- r
    0
  }, numeric(1L))
  noise_cov <- .nearest_pd(noise_cor)
  L_noise   <- chol(noise_cov)

  person_means <- matrix(stats::rnorm(n_subjects * d, mean = 0, sd = 1.5),
                         nrow = n_subjects, ncol = d)

  all_rows <- vector("list", n_subjects)
  for (subj in seq_len(n_subjects)) {
    mu     <- person_means[subj, ]
    y_prev <- mu + stats::rnorm(d, sd = 0.5)
    subj_data <- matrix(NA_real_, nrow = obs_per_subject, ncol = d)
    subj_data[1L, ] <- y_prev
    for (t in seq(2L, obs_per_subject)) {
      if ("heavy_tailed" %in% complexity) {
        df_t <- sample(2L:5L, 1L)
        raw_innov <- stats::rt(d, df = df_t)
      } else {
        raw_innov <- stats::rnorm(d)
      }
      innovation <- as.vector(raw_innov %*% L_noise) * 0.5
      y_new      <- mu + as.vector(B %*% (y_prev - mu)) + innovation
      subj_data[t, ] <- y_new
      y_prev <- y_new
    }
    day_vec  <- rep(seq_len(n_days), each = beeps_per_day)
    beep_vec <- rep(seq_len(beeps_per_day), times = n_days)
    subj_df  <- data.frame(id = rep(subj, obs_per_subject),
                           day = day_vec, beep = beep_vec)
    subj_df[var_names] <- subj_data
    all_rows[[subj]] <- subj_df
  }

  df <- do.call(rbind, all_rows)
  rownames(df) <- NULL
  attr(df, "true_temporal")       <- B
  attr(df, "true_contemporaneous") <- noise_cor
  attr(df, "vars")                 <- var_names
  df
}


# ===========================================================================
# S3 print methods for batch objects
# ===========================================================================

#' @export
print.sim_batch_type <- function(x, ...) {
  cat(sprintf(
    "Batch of %d '%s' datasets (base_seed=%d, complexity='%s')\n",
    length(x), attr(x, "type"), attr(x, "seed"),
    paste(if (length(attr(x, "complexity")) == 0L) "auto"
          else attr(x, "complexity"), collapse = "+")
  ))
  invisible(x)
}

#' @export
print.sim_batch_full <- function(x, ...) {
  cat(sprintf(
    "Full dataset batch: %d datasets \u00d7 %d types = %d total (base_seed=%d)\n",
    attr(x, "n_batch"), length(x),
    attr(x, "n_batch") * length(x), attr(x, "seed")
  ))
  invisible(x)
}


# ===========================================================================
# Constants
# ===========================================================================

.SIM_VALID_TYPES <- c(
  "ttest", "anova", "correlation", "clusters",
  "factor_analysis", "prediction", "mlvar"
)

.SIM_INFO_MAP <- c(
  ttest           = "t.test(score ~ group, data = d)",
  anova           = "summary(aov(score ~ group, data = d))",
  correlation     = "cor(d); pairs(d)",
  clusters        = "kmeans(d[, -ncol(d)], centers = max(d$true_cluster))",
  factor_analysis = "factanal(d, factors = attr(d, 'n_factors'))",
  prediction      = "summary(lm(y ~ ., data = d))",
  mlvar           = "mlvar(d, vars = attr(d, 'vars'), id = 'id', day = 'day', beep = 'beep')"
)

.SIM_ALL_CASES <- c(
  "na", "outliers", "ties", "duplicates", "constant_col", "all_na_col",
  "tiny_n", "heavy_tailed", "heteroscedastic", "extreme_imbalance",
  "multicollinear"
)


# ===========================================================================
# Helper: null coalescing (defined here; utils.R moved to Nestimate)
# ===========================================================================

`%||%` <- function(a, b) if (!is.null(a)) a else b


# ===========================================================================
# Helper: nearest positive-definite matrix
# ===========================================================================

#' @keywords internal
#' @noRd
.nearest_pd <- function(R) {
  eig   <- eigen(R, symmetric = TRUE)
  vals  <- pmax(eig$values, 1e-6)
  R_pd  <- eig$vectors %*% diag(vals) %*% t(eig$vectors)
  d_vec <- sqrt(diag(R_pd))
  R_pd  <- R_pd / outer(d_vec, d_vec)
  (R_pd + t(R_pd)) / 2
}
