# ===========================================================================
# simulate_longitudinal — basic structure
# ===========================================================================

test_that("simulate_longitudinal: returns saqr_sim with correct structure", {
  r <- simulate_longitudinal(n = 10, tp = 20, vars = 3, seed = 1)

  expect_s3_class(r, "saqr_sim")
  expect_equal(r$type, "longitudinal")
  expect_equal(r$seed, 1L)
  expect_s3_class(r$data, "data.frame")
})

test_that("simulate_longitudinal: data dimensions match n * tp", {
  r <- simulate_longitudinal(n = 15, tp = 30, vars = 4, seed = 42)
  expect_equal(nrow(r$data), 15L * 30L)
  expect_true(all(c("id", "time", "V1", "V2", "V3", "V4") %in% names(r$data)))
})

test_that("simulate_longitudinal: id and time columns are correct", {
  r <- simulate_longitudinal(n = 5, tp = 10, vars = 2, seed = 7)
  expect_equal(sort(unique(r$data$id)), 1:5)
  # Each subject has tp time points
  counts <- table(r$data$id)
  expect_true(all(counts == 10L))
  # Time runs 1..tp within each subject
  subj1 <- r$data[r$data$id == 1, ]
  expect_equal(subj1$time, 1:10)
})

test_that("simulate_longitudinal: custom variable names", {
  r <- simulate_longitudinal(n = 5, tp = 10,
                              vars = c("mood", "energy", "stress"), seed = 1)
  expect_true(all(c("mood", "energy", "stress") %in% names(r$data)))
  expect_equal(r$params$var_names, c("mood", "energy", "stress"))
})

test_that("simulate_longitudinal: seed reproducibility", {
  r1 <- simulate_longitudinal(n = 10, tp = 20, vars = 3, seed = 99)
  r2 <- simulate_longitudinal(n = 10, tp = 20, vars = 3, seed = 99)
  expect_identical(r1$data, r2$data)
  expect_identical(r1$params$temporal, r2$params$temporal)
})

test_that("simulate_longitudinal: different seeds produce different data", {
  r1 <- simulate_longitudinal(n = 10, tp = 20, vars = 3, seed = 1)
  r2 <- simulate_longitudinal(n = 10, tp = 20, vars = 3, seed = 2)
  expect_false(identical(r1$data$V1, r2$data$V1))
})


# ===========================================================================
# Parameter recovery (temporal structure)
# ===========================================================================

test_that("simulate_longitudinal: explicit temporal matrix is stored and used", {
  B <- matrix(c(0.4, 0.0,
                0.2, 0.3), nrow = 2, byrow = TRUE)
  r <- simulate_longitudinal(n = 10, tp = 50, vars = 2, temporal = B, seed = 1)
  expect_equal(unname(r$params$temporal), B, tolerance = 1e-10)
})

test_that("simulate_longitudinal: auto-generated temporal is stationary", {
  r <- simulate_longitudinal(n = 5, tp = 20, vars = 5, seed = 42)
  ev <- eigen(r$params$temporal, only.values = TRUE)$values
  expect_true(all(Mod(ev) < 1))
})

test_that("simulate_longitudinal: temporal coefficients approximately recoverable", {
  B <- matrix(c(0.4,  0.0,
                0.15, 0.3), nrow = 2, byrow = TRUE)
  r <- simulate_longitudinal(
    n = 200, tp = 100, vars = 2, temporal = B,
    between = diag(2) * 0.5, innovation_sd = 0.5, seed = 42
  )
  d <- r$data
  # Person-mean center
  d$V1c <- d$V1 - ave(d$V1, d$id, FUN = mean)
  d$V2c <- d$V2 - ave(d$V2, d$id, FUN = mean)
  # Lagged values
  d$V1_lag <- ave(d$V1c, d$id, FUN = function(x) c(NA, x[-length(x)]))
  d$V2_lag <- ave(d$V2c, d$id, FUN = function(x) c(NA, x[-length(x)]))
  d <- d[!is.na(d$V1_lag), ]
  # OLS recovery of row 1 of B: V1 ~ V1_lag + V2_lag
  fit <- lm(V1c ~ V1_lag + V2_lag - 1, data = d)
  recovered <- coef(fit)
  # Should be in the ballpark (not exact due to finite sample)
  expect_true(abs(recovered["V1_lag"] - 0.4) < 0.15)
  expect_true(abs(recovered["V2_lag"] - 0.0) < 0.15)
})


# ===========================================================================
# params structure
# ===========================================================================

test_that("simulate_longitudinal: params contains all expected fields", {
  r <- simulate_longitudinal(n = 5, tp = 10, vars = 3, seed = 1)
  expected <- c("temporal", "contemporaneous", "between", "grand_means",
                "innovation_sd", "n", "tp", "var_names")
  expect_true(all(expected %in% names(r$params)))
})

test_that("simulate_longitudinal: explicit contemporaneous stored correctly", {
  C <- matrix(c(1.0, 0.4, 0.4, 1.0), nrow = 2)
  r <- simulate_longitudinal(n = 5, tp = 10, vars = 2,
                              contemporaneous = C, seed = 1)
  expect_equal(unname(r$params$contemporaneous), C)
})

test_that("simulate_longitudinal: explicit between stored correctly", {
  Bw <- matrix(c(2.0, 0.5, 0.5, 1.5), nrow = 2)
  r <- simulate_longitudinal(n = 5, tp = 10, vars = 2, between = Bw, seed = 1)
  expect_equal(unname(r$params$between), Bw)
})

test_that("simulate_longitudinal: grand_means shift data", {
  r0 <- simulate_longitudinal(n = 50, tp = 50, vars = 2,
                               grand_means = c(0, 0), seed = 1)
  r5 <- simulate_longitudinal(n = 50, tp = 50, vars = 2,
                               grand_means = c(5, 10), seed = 1)
  # Mean of V1 should be ~5 higher in r5
  diff1 <- mean(r5$data$V1) - mean(r0$data$V1)
  expect_true(abs(diff1 - 5) < 1.5)
})


# ===========================================================================
# ESM / beeps_per_day
# ===========================================================================

test_that("simulate_longitudinal: beeps_per_day creates day/beep columns", {
  r <- simulate_longitudinal(n = 5, tp = 21, vars = 2,
                              beeps_per_day = 7, seed = 1)
  expect_true(all(c("day", "beep") %in% names(r$data)))
  expect_false("time" %in% names(r$data))
  # 3 days × 7 beeps
  expect_equal(max(r$data$day), 3L)
  expect_equal(max(r$data$beep), 7L)
})

test_that("simulate_longitudinal: beeps_per_day errors on non-divisible tp", {
  expect_error(
    simulate_longitudinal(n = 5, tp = 20, vars = 2, beeps_per_day = 7, seed = 1),
    "divisible"
  )
})


# ===========================================================================
# Complexity injection
# ===========================================================================

test_that("simulate_longitudinal: complexity='clean' produces no NAs", {
  r <- simulate_longitudinal(n = 10, tp = 20, vars = 3,
                              complexity = "clean", seed = 1)
  expect_false(anyNA(r$data))
})

test_that("simulate_longitudinal: complexity='na' injects NAs", {
  r <- simulate_longitudinal(n = 20, tp = 50, vars = 3,
                              complexity = "na", seed = 1)
  expect_true(anyNA(r$data[, c("V1", "V2", "V3")]))
})

test_that("simulate_longitudinal: complexity='outliers' produces extreme values", {
  r_clean <- simulate_longitudinal(n = 20, tp = 50, vars = 3,
                                    complexity = "clean", seed = 1)
  r_out   <- simulate_longitudinal(n = 20, tp = 50, vars = 3,
                                    complexity = "outliers", seed = 1)
  # Outlier version should have wider range
  expect_true(max(abs(r_out$data$V1), na.rm = TRUE) >
                max(abs(r_clean$data$V1), na.rm = TRUE))
})

test_that("simulate_longitudinal: invalid complexity errors", {
  expect_error(
    simulate_longitudinal(n = 5, tp = 10, vars = 2,
                           complexity = "nonsense", seed = 1),
    "Unknown complexity"
  )
})

test_that("simulate_longitudinal: complexity vector works", {
  r <- simulate_longitudinal(n = 10, tp = 30, vars = 2,
                              complexity = c("na", "outliers"), seed = 1)
  expect_true(anyNA(r$data[, c("V1", "V2")]))
})


# ===========================================================================
# Edge cases and validation
# ===========================================================================

test_that("simulate_longitudinal: errors on p < 2", {
  expect_error(simulate_longitudinal(n = 5, tp = 10, vars = 1, seed = 1))
})

test_that("simulate_longitudinal: errors on mismatched temporal dimensions", {
  B3 <- diag(0.3, 3)
  expect_error(
    simulate_longitudinal(n = 5, tp = 10, vars = 2, temporal = B3, seed = 1)
  )
})

test_that("simulate_longitudinal: non-stationary B is auto-rescaled", {
  B_bad <- matrix(c(1.2, 0.0, 0.0, 0.9), nrow = 2)
  r <- simulate_longitudinal(n = 5, tp = 20, vars = 2,
                              temporal = B_bad, seed = 1)
  ev <- eigen(r$params$temporal, only.values = TRUE)$values
  expect_true(all(Mod(ev) < 1))
})

test_that("simulate_longitudinal: NULL seed still stores a seed", {
  r <- simulate_longitudinal(n = 5, tp = 10, vars = 2, seed = NULL)
  expect_true(!is.null(r$seed))
  expect_true(is.numeric(r$seed))
})

test_that("simulate_longitudinal: saqr_sim methods work", {
  r <- simulate_longitudinal(n = 5, tp = 10, vars = 2, seed = 1)
  expect_equal(nrow(r), 50L)
  expect_equal(ncol(r), 4L)  # id, time, V1, V2
  expect_equal(nrow(head(r, 3)), 3L)
  expect_output(print(r), "saqr_sim")
  expect_identical(as.data.frame(r), r$data)
})
