# Tests for inject_missingness() — MCAR / MAR / MNAR injection transformer

make_df <- function(n = 5000, seed = 1) {
  set.seed(seed)
  data.frame(
    x1 = stats::rnorm(n),
    x2 = stats::rnorm(n),
    y  = stats::rnorm(n)
  )
}

test_that("MCAR realized fraction approximates prop at large n", {
  df <- make_df(n = 8000, seed = 10)
  out <- inject_missingness(df, mechanism = "MCAR", prop = 0.2, seed = 42)
  info <- attr(out, "missing_info")

  expect_equal(info$mechanism, "MCAR")
  # tolerance is relative; ~5% of 0.2 = +/- 0.01 absolute
  expect_equal(info$realized_prop, 0.2, tolerance = 0.05)
})

test_that("attribute is present and correctly structured", {
  df <- make_df(n = 1000, seed = 11)
  out <- inject_missingness(df, mechanism = "MCAR", prop = 0.15,
                            cols = c("x1", "y"), seed = 3)
  info <- attr(out, "missing_info")

  expect_type(info, "list")
  expect_named(info, c("mechanism", "prop", "cols", "n_missing",
                       "realized_prop", "indicator"))
  expect_equal(info$cols, c("x1", "y"))
  expect_true(is.matrix(info$indicator))
  expect_true(is.logical(info$indicator))
  expect_equal(dim(info$indicator), c(1000L, 2L))
  expect_equal(colnames(info$indicator), c("x1", "y"))
  # n_missing must match both the indicator and the actual NA cells
  expect_equal(info$n_missing, sum(info$indicator))
  expect_equal(info$n_missing, sum(is.na(out[c("x1", "y")])))
})

test_that("columns not in cols are untouched", {
  df <- make_df(n = 2000, seed = 12)
  out <- inject_missingness(df, mechanism = "MCAR", prop = 0.3,
                            cols = c("x1"), seed = 5)

  expect_false(anyNA(out$x2))
  expect_false(anyNA(out$y))
  expect_true(anyNA(out$x1))
})

test_that("MAR missingness correlates with predictor and differs from random", {
  df <- make_df(n = 6000, seed = 13)
  out <- inject_missingness(df, mechanism = "MAR", prop = 0.25,
                            cols = "y", predictor = "x1", seed = 99)
  info <- attr(out, "missing_info")

  expect_equal(info$mechanism, "MAR")
  # predictor must be excluded from targets
  expect_false("x1" %in% info$cols)
  # realized fraction calibrated near prop
  expect_equal(info$realized_prop, 0.25, tolerance = 0.03)

  # Missingness flag should correlate with the predictor values
  flag <- info$indicator[, "y"]
  pred <- df$x1
  rho <- stats::cor(as.numeric(flag), pred)
  expect_gt(abs(rho), 0.15)

  # The predictor's mean among missing vs observed should differ clearly
  mean_missing  <- mean(pred[flag])
  mean_observed <- mean(pred[!flag])
  expect_gt(mean_missing - mean_observed, 0.2)
})

test_that("MAR requires a predictor", {
  df <- make_df(n = 100, seed = 14)
  expect_error(inject_missingness(df, mechanism = "MAR", prop = 0.2))
})

test_that("MNAR missingness correlates with the underlying value", {
  df <- make_df(n = 6000, seed = 15)
  out <- inject_missingness(df, mechanism = "MNAR", prop = 0.25,
                            cols = "y", seed = 77)
  info <- attr(out, "missing_info")

  expect_equal(info$mechanism, "MNAR")
  expect_equal(info$realized_prop, 0.25, tolerance = 0.03)

  flag <- info$indicator[, "y"]
  underlying <- df$y
  # Larger underlying values should be more likely to be missing
  mean_missing  <- mean(underlying[flag])
  mean_observed <- mean(underlying[!flag])
  expect_gt(mean_missing - mean_observed, 0.2)

  rho <- stats::cor(as.numeric(flag), underlying)
  expect_gt(rho, 0.15)
})

test_that("reproducibility: same seed identical, different seed differs", {
  df <- make_df(n = 3000, seed = 16)

  a <- inject_missingness(df, mechanism = "MCAR", prop = 0.2, seed = 100)
  b <- inject_missingness(df, mechanism = "MCAR", prop = 0.2, seed = 100)
  c <- inject_missingness(df, mechanism = "MCAR", prop = 0.2, seed = 200)

  expect_identical(attr(a, "missing_info")$indicator,
                   attr(b, "missing_info")$indicator)
  expect_false(identical(attr(a, "missing_info")$indicator,
                         attr(c, "missing_info")$indicator))

  # Same holds for MAR
  a2 <- inject_missingness(df, mechanism = "MAR", prop = 0.2,
                           predictor = "x1", cols = "y", seed = 100)
  b2 <- inject_missingness(df, mechanism = "MAR", prop = 0.2,
                           predictor = "x1", cols = "y", seed = 100)
  expect_identical(attr(a2, "missing_info")$indicator,
                   attr(b2, "missing_info")$indicator)
})

test_that("cols default NULL targets all columns", {
  df <- make_df(n = 500, seed = 17)
  out <- inject_missingness(df, mechanism = "MCAR", prop = 0.2, seed = 1)
  expect_equal(attr(out, "missing_info")$cols, names(df))
})

test_that("prop = 0 injects nothing, prop = 1 injects everything", {
  df <- make_df(n = 300, seed = 18)
  none <- inject_missingness(df, mechanism = "MNAR", prop = 0, cols = "y",
                             seed = 1)
  all_ <- inject_missingness(df, mechanism = "MNAR", prop = 1, cols = "y",
                             seed = 1)
  expect_equal(attr(none, "missing_info")$n_missing, 0L)
  expect_false(anyNA(none$y))
  expect_true(all(is.na(all_$y)))
})

test_that("accepts a saqr_sim object and operates on its $data", {
  sim <- simulate_ttest(n_a = 200, n_b = 200, mean_a = 0, mean_b = 1, seed = 2)
  out <- inject_missingness(sim, mechanism = "MCAR", prop = 0.2,
                            cols = "score", seed = 9)
  expect_s3_class(out, "data.frame")
  expect_true(anyNA(out$score))
  expect_equal(attr(out, "missing_info")$cols, "score")
})

test_that("invalid prop is rejected", {
  df <- make_df(n = 50, seed = 19)
  expect_error(inject_missingness(df, mechanism = "MCAR", prop = 1.5))
  expect_error(inject_missingness(df, mechanism = "MCAR", prop = -0.1))
})
