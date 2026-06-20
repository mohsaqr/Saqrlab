# Tests for simulate_survival() — Cox proportional-hazards generative model

test_that("simulate_survival returns a well-formed saqr_sim", {
  r <- simulate_survival(n = 200, n_covariates = 2, seed = 1)

  expect_s3_class(r, "saqr_sim")
  expect_equal(r$type, "survival")

  d <- as.data.frame(r)
  expect_true(all(c("time", "status", "X1", "X2") %in% names(d)))
  expect_equal(nrow(d), 200L)

  # status is binary, time strictly positive.
  expect_true(all(d$status %in% c(0L, 1L)))
  expect_true(all(d$time > 0))
  expect_false(anyNA(d$time))
})

test_that("params carry the true generating parameters", {
  r <- simulate_survival(n = 100, betas = c(0.7, -0.4, 0.2),
                         baseline = "weibull", lambda = 0.15, shape = 1.5,
                         seed = 7)
  p <- r$params

  expect_equal(unname(p$betas), c(0.7, -0.4, 0.2))
  expect_equal(p$baseline, "weibull")
  expect_equal(p$lambda, 0.15)
  expect_equal(p$shape, 1.5)
  expect_true(is.numeric(p$realized_censoring_rate))
  expect_true(p$realized_censoring_rate >= 0 && p$realized_censoring_rate <= 1)
})

test_that("realized censoring rate is close to the target", {
  r <- simulate_survival(n = 4000, n_covariates = 2, censoring_rate = 0.3,
                         seed = 11)
  expect_equal(r$params$realized_censoring_rate, 0.3, tolerance = 0.05)

  r2 <- simulate_survival(n = 4000, n_covariates = 2, censoring_rate = 0.5,
                          seed = 12)
  expect_equal(r2$params$realized_censoring_rate, 0.5, tolerance = 0.05)
})

test_that("zero censoring yields all events", {
  r <- simulate_survival(n = 300, censoring_rate = 0, seed = 3)
  d <- as.data.frame(r)
  expect_true(all(d$status == 1L))
  expect_equal(r$params$realized_censoring_rate, 0)
})

test_that("larger linear predictor implies shorter event times", {
  # Among events, the true linear predictor X %*% betas should correlate
  # NEGATIVELY with the observed event time (higher hazard -> earlier event).
  r <- simulate_survival(n = 5000, betas = c(1.0, -0.8), censoring_rate = 0.2,
                         seed = 21)
  d <- as.data.frame(r)
  betas <- r$params$betas

  x_mat <- as.matrix(d[, c("X1", "X2")])
  eta <- as.numeric(x_mat %*% betas)

  events <- d$status == 1L
  rho <- cor(eta[events], d$time[events], method = "spearman")
  expect_lt(rho, 0)
})

test_that("all three baselines produce valid positive times", {
  for (bl in c("weibull", "exponential", "gompertz")) {
    r <- simulate_survival(n = 300, baseline = bl, seed = 99)
    d <- as.data.frame(r)
    expect_true(all(d$time > 0), info = bl)
    expect_true(all(is.finite(d$time)), info = bl)
    expect_true(all(d$status %in% c(0L, 1L)), info = bl)
  }
})

test_that("binary covariates are coded 0/1", {
  r <- simulate_survival(n = 200, n_covariates = 2,
                         covariate_type = "binary", seed = 5)
  d <- as.data.frame(r)
  expect_true(all(d$X1 %in% c(0, 1)))
  expect_true(all(d$X2 %in% c(0, 1)))
})

test_that("results are reproducible with a fixed seed", {
  r1 <- simulate_survival(n = 250, n_covariates = 3, seed = 2024)
  r2 <- simulate_survival(n = 250, n_covariates = 3, seed = 2024)
  expect_equal(as.data.frame(r1), as.data.frame(r2))
  expect_equal(r1$params, r2$params)
})

test_that("input validation rejects bad arguments", {
  expect_error(simulate_survival(n = 2))
  expect_error(simulate_survival(n = 100, lambda = -1))
  expect_error(simulate_survival(n = 100, shape = 0))
  expect_error(simulate_survival(n = 100, censoring_rate = 1))
  expect_error(simulate_survival(n = 100, censoring_rate = -0.1))
})

test_that("coxph recovers betas when the survival package is available", {
  skip_if_not_installed("survival")

  betas <- c(0.6, -0.4)
  r <- simulate_survival(n = 6000, betas = betas, baseline = "weibull",
                         lambda = 0.1, shape = 1.2, censoring_rate = 0.25,
                         seed = 314)
  d <- as.data.frame(r)

  fit <- survival::coxph(
    survival::Surv(time, status) ~ X1 + X2, data = d
  )
  est <- unname(coef(fit))
  expect_equal(est, betas, tolerance = 0.12)
})
