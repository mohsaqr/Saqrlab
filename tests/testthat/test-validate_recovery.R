# Tests for validate_recovery() — parameter-recovery scoring

test_that("exact recovery yields all within_tol and zero rel_error", {
  truth <- c(a = 2, b = -3, c = 0.5)
  res <- validate_recovery(truth, estimates = truth)

  expect_s3_class(res, "recovery_result")
  expect_s3_class(res, "data.frame")
  expect_identical(nrow(res), 3L)
  expect_true(all(res$within_tol))
  expect_equal(res$abs_error, c(0, 0, 0))
  expect_equal(res$rel_error, c(0, 0, 0))
})

test_that("known off-by-X estimate gives correct errors and flags", {
  truth <- c(a = 10, b = 100)
  # a is off by 1 (rel 0.10), b off by 20 (rel 0.20)
  est <- c(a = 11, b = 120)
  res <- validate_recovery(truth, estimates = est, tolerance = 0.1,
                           relative = TRUE)

  expect_equal(res$abs_error, c(1, 20))
  expect_equal(res$rel_error, c(0.1, 0.2))
  # a: rel_error 0.1 <= 0.1 TRUE; b: 0.2 <= 0.1 FALSE
  expect_equal(res$within_tol, c(TRUE, FALSE))

  # absolute tolerance flips both depending on threshold
  res_abs <- validate_recovery(truth, estimates = est, tolerance = 5,
                               relative = FALSE)
  expect_equal(res_abs$within_tol, c(TRUE, FALSE))
})

test_that("true == 0 yields NA rel_error", {
  res <- validate_recovery(c(z = 0), estimates = c(z = 0.05),
                           tolerance = 0.1, relative = TRUE)
  expect_true(is.na(res$rel_error))
  expect_equal(res$abs_error, 0.05)
  # relative comparison against NA -> within_tol NA
  expect_true(is.na(res$within_tol))
})

test_that("vector-valued params flatten into name.subname rows", {
  fake_sim <- structure(
    list(
      data = data.frame(x = 1),
      params = list(
        error_sd = 2,
        betas = c(x1 = 0.5, x2 = -0.3),
        label = "skip-me"   # non-numeric: skipped
      ),
      type = "fake", seed = NULL
    ),
    class = c("saqr_sim", "list")
  )

  est <- c("error_sd" = 2, "betas.x1" = 0.5, "betas.x2" = -0.3)
  res <- validate_recovery(fake_sim, estimates = est)

  expect_setequal(res$parameter, c("error_sd", "betas.x1", "betas.x2"))
  expect_false("label" %in% res$parameter)
  expect_true(all(res$within_tol))
})

test_that("params argument restricts which names are compared", {
  truth <- c(a = 1, b = 2, c = 3)
  est <- c(a = 1, b = 2, c = 3)
  res <- validate_recovery(truth, estimates = est, params = c("a", "c"))
  expect_setequal(res$parameter, c("a", "c"))
  expect_identical(nrow(res), 2L)
})

test_that("end-to-end regression coefficients are recovered within tolerance", {
  rsim <- simulate_regression(
    coefs = c("(Intercept)" = 1, x1 = 2, x2 = -1.5),
    predictor_sds = c(x1 = 1, x2 = 1),
    error_sd = 1, n = 4000, seed = 7
  )
  fit <- lm(y ~ x1 + x2, data = rsim$data)
  est <- stats::setNames(coef(fit), paste0("coefs.", names(coef(fit))))

  res <- validate_recovery(rsim, estimates = est, tolerance = 0.1,
                           relative = TRUE)

  expect_setequal(res$parameter,
                  c("coefs.(Intercept)", "coefs.x1", "coefs.x2"))
  expect_true(all(res$within_tol))
})

test_that("summary returns the documented one-row data.frame", {
  truth <- c(a = 10, b = 100)
  est <- c(a = 11, b = 120)
  res <- validate_recovery(truth, estimates = est, tolerance = 0.1)
  s <- summary(res)

  expect_s3_class(s, "data.frame")
  expect_identical(nrow(s), 1L)
  expect_identical(
    names(s),
    c("n_params", "n_within_tol", "pct_within_tol",
      "mean_abs_error", "mean_rel_error", "max_abs_error")
  )
  expect_equal(s$n_params, 2L)
  expect_equal(s$n_within_tol, 1L)
  expect_equal(s$pct_within_tol, 50)
  expect_equal(s$mean_abs_error, mean(c(1, 20)))
  expect_equal(s$max_abs_error, 20)
})

test_that("print runs without error", {
  truth <- c(a = 10, b = 100)
  res <- validate_recovery(truth, estimates = c(a = 11, b = 120))
  expect_output(print(res), "Parameter recovery")
})

test_that("validation errors on bad inputs", {
  expect_error(validate_recovery(c(a = 1), estimates = c(1, 2)))  # unnamed
  expect_error(validate_recovery(c(a = 1), estimates = c(a = 1),
                                 tolerance = -1))                  # tol <= 0
  expect_error(validate_recovery(c(a = 1), estimates = c(z = 1)))  # no overlap
})
