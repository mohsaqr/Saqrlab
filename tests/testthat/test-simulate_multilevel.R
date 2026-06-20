# Tests for simulate_mlm() and simulate_growth()

# ---------------------------------------------------------------------------
# simulate_mlm — structure
# ---------------------------------------------------------------------------

test_that("simulate_mlm returns a well-formed saqr_sim", {
  r <- simulate_mlm(n_clusters = 20, cluster_size = 15, n_predictors = 2,
                    seed = 1)
  expect_s3_class(r, "saqr_sim")
  expect_equal(r$type, "mlm")
  expect_true(all(c("cluster_id", "y", "x1", "x2") %in% names(r$data)))
  expect_s3_class(r$data$cluster_id, "factor")
  expect_equal(nrow(r$data), 20 * 15)
  expect_equal(nlevels(r$data$cluster_id), 20)

  # params structure
  expect_true(all(c("intercept", "slope", "betas", "icc", "tau00",
                    "residual_sd", "slope_sd", "cluster_intercepts") %in%
                    names(r$params)))
  expect_equal(unname(r$params$betas["x1"]), r$params$slope)
  expect_length(r$params$cluster_intercepts, 20)
  expect_null(r$params$cluster_slopes)
})

test_that("simulate_mlm tau00 follows the ICC formula", {
  r <- simulate_mlm(n_clusters = 10, icc = 0.25, residual_sd = 2, seed = 7)
  expect_equal(r$params$tau00, 0.25 / 0.75 * 4)
})

test_that("simulate_mlm supports unbalanced clusters", {
  sizes <- c(10, 20, 30, 5, 8)
  r <- simulate_mlm(n_clusters = 5, cluster_size = sizes, seed = 3)
  expect_equal(nrow(r$data), sum(sizes))
  expect_equal(unname(as.integer(table(r$data$cluster_id))), sizes)
  expect_equal(unname(r$params$cluster_size), sizes)
})

test_that("simulate_mlm random slope is recorded", {
  r <- simulate_mlm(n_clusters = 15, random_slope = TRUE, slope_sd = 0.4,
                    seed = 5)
  expect_true(r$params$random_slope)
  expect_length(r$params$cluster_slopes, 15)
})

# ---------------------------------------------------------------------------
# simulate_mlm — parameter recovery (base R only)
# ---------------------------------------------------------------------------

test_that("simulate_mlm recovers ICC via variance decomposition", {
  r <- simulate_mlm(n_clusters = 200, cluster_size = 40, icc = 0.30,
                    slope = 0.5, residual_sd = 1, seed = 11)
  d <- r$data
  # The design ICC is conditional on the fixed predictors:
  # tau00 / (tau00 + sigma^2). Residualise out the fixed slope first so the
  # predictor variance does not contaminate the decomposition, then apply the
  # balanced one-way ANOVA estimator (var of cluster means is inflated by
  # sigma^2/m, so back it out via the within-cluster mean square).
  d$resid <- stats::residuals(stats::lm(y ~ x1, data = d))
  m <- 40
  cluster_means <- tapply(d$resid, d$cluster_id, mean)
  grand_mean <- mean(d$resid)
  k <- length(cluster_means)
  ms_between <- m * sum((cluster_means - grand_mean)^2) / (k - 1)
  within_ss <- sum(tapply(d$resid, d$cluster_id,
                          function(z) sum((z - mean(z))^2)))
  ms_within <- within_ss / (nrow(d) - k)
  tau00_hat <- (ms_between - ms_within) / m
  emp_icc <- tau00_hat / (tau00_hat + ms_within)
  expect_equal(emp_icc, 0.30, tolerance = 0.06)
})

test_that("simulate_mlm recovers the fixed slope via pooled lm", {
  r <- simulate_mlm(n_clusters = 150, cluster_size = 40, slope = 0.7,
                    icc = 0.15, seed = 13)
  fit <- stats::lm(y ~ x1, data = r$data)
  expect_equal(unname(stats::coef(fit)["x1"]), 0.7, tolerance = 0.05)
})

test_that("simulate_mlm recovers fixed slope with lme4 if available", {
  skip_if_not_installed("lme4")
  r <- simulate_mlm(n_clusters = 120, cluster_size = 40, slope = 0.6,
                    icc = 0.2, seed = 17)
  fit <- lme4::lmer(y ~ x1 + (1 | cluster_id), data = r$data)
  fe <- lme4::fixef(fit)
  expect_equal(unname(fe["x1"]), 0.6, tolerance = 0.05)
})

# ---------------------------------------------------------------------------
# simulate_growth — structure
# ---------------------------------------------------------------------------

test_that("simulate_growth returns a long-format saqr_sim", {
  r <- simulate_growth(n = 100, n_time = 5, seed = 1)
  expect_s3_class(r, "saqr_sim")
  expect_equal(r$type, "growth")
  expect_equal(names(r$data), c("subject", "time", "y"))
  expect_s3_class(r$data$subject, "factor")
  expect_equal(nrow(r$data), 100 * 5)
  expect_equal(sort(unique(r$data$time)), 0:4)
  expect_equal(nlevels(r$data$subject), 100)

  expect_true(all(c("means", "sds", "correlation", "residual_sd",
                    "subject_intercepts", "subject_slopes") %in%
                    names(r$params)))
  expect_length(r$params$subject_intercepts, 100)
  expect_length(r$params$subject_slopes, 100)
})

# ---------------------------------------------------------------------------
# simulate_growth — parameter recovery (base R only)
# ---------------------------------------------------------------------------

test_that("simulate_growth recovers mean slope from per-subject OLS", {
  r <- simulate_growth(n = 800, n_time = 6, slope_mean = 2,
                       intercept_mean = 5, residual_sd = 1, seed = 21)
  d <- r$data
  per_subject_slope <- vapply(
    split(d, d$subject),
    function(s) unname(stats::coef(stats::lm(y ~ time, data = s))["time"]),
    numeric(1)
  )
  expect_equal(mean(per_subject_slope), 2, tolerance = 0.05)
})

test_that("simulate_growth recovers intercept sd at large n", {
  r <- simulate_growth(n = 2000, n_time = 5, intercept_sd = 2,
                       slope_sd = 0.5, residual_sd = 1, seed = 23)
  # True per-subject intercepts are stored as ground truth.
  expect_equal(stats::sd(r$params$subject_intercepts), 2, tolerance = 0.15)

  # Empirical recovery from data: intercept = predicted y at time 0.
  d <- r$data
  per_subject_intercept <- vapply(
    split(d, d$subject),
    function(s) unname(stats::coef(stats::lm(y ~ time, data = s))["(Intercept)"]),
    numeric(1)
  )
  expect_equal(stats::sd(per_subject_intercept), 2, tolerance = 0.25)
})

test_that("simulate_growth recovers intercept-slope correlation", {
  r <- simulate_growth(n = 3000, intercept_slope_cor = 0.5, intercept_sd = 1.5,
                       slope_sd = 0.8, seed = 29)
  emp_cor <- stats::cor(r$params$subject_intercepts, r$params$subject_slopes)
  expect_equal(emp_cor, 0.5, tolerance = 0.08)
})

# ---------------------------------------------------------------------------
# reproducibility
# ---------------------------------------------------------------------------

test_that("both simulators are reproducible by seed", {
  expect_equal(simulate_mlm(seed = 99)$data, simulate_mlm(seed = 99)$data)
  expect_equal(simulate_growth(seed = 99)$data,
               simulate_growth(seed = 99)$data)
})
