# ===========================================================================
# simulate_lpa
# ===========================================================================

test_that("simulate_lpa: returns list with $data and $params", {
  means <- matrix(c(2, 5, 3, 7), nrow = 2, ncol = 2)
  r <- simulate_lpa(means = means, sds = 0.5, props = c(0.5, 0.5), n = 100, seed = 1)
  expect_type(r, "list")
  expect_true(all(c("data", "params") %in% names(r)))
})

test_that("simulate_lpa: $data is data.frame with correct columns", {
  means <- matrix(c(2, 5, 3, 7), nrow = 2, ncol = 2)
  r <- simulate_lpa(means = means, sds = 0.5, props = c(0.5, 0.5), n = 100, seed = 1)
  expect_s3_class(r$data, "data.frame")
  expect_true(all(c("y1", "y2", "true_profile") %in% names(r$data)))
  expect_equal(nrow(r$data), 100L)
})

test_that("simulate_lpa: true_profile contains only valid class labels", {
  means <- matrix(c(1, 4, 8), nrow = 1, ncol = 3)
  r <- simulate_lpa(means = means, sds = 0.3, props = c(0.2, 0.5, 0.3), n = 200, seed = 1)
  expect_true(all(r$data$true_profile %in% 1:3))
  expect_true(all(c(1L, 2L, 3L) %in% r$data$true_profile))
})

test_that("simulate_lpa: $params holds exact input parameters", {
  means <- matrix(c(0, 5, 0, 5), nrow = 2, ncol = 2)
  sds   <- matrix(c(1, 1, 1, 1), nrow = 2, ncol = 2)
  props <- c(0.4, 0.6)
  r  <- simulate_lpa(means = means, sds = sds, props = props, n = 50, seed = 1)
  expect_identical(r$params$means, means)
  expect_identical(r$params$props, props)
})

test_that("simulate_lpa: profile means are approximately recovered at large n", {
  means <- matrix(c(0, 0, 10, 10), nrow = 2, ncol = 2)
  r     <- simulate_lpa(means = means, sds = 0.5, props = c(0.5, 0.5), n = 2000, seed = 42)
  d     <- r$data
  grp_means <- vapply(1:2, function(k) {
    colMeans(d[d$true_profile == k, c("y1", "y2")])
  }, numeric(2))
  expect_true(abs(grp_means[1, 1] - 0)  < 0.5)
  expect_true(abs(grp_means[1, 2] - 10) < 0.5)
})

test_that("simulate_lpa: seed produces identical results", {
  means <- matrix(c(1, 6), nrow = 1, ncol = 2)
  r1 <- simulate_lpa(means = means, sds = 1, props = c(0.5, 0.5), n = 50, seed = 7)
  r2 <- simulate_lpa(means = means, sds = 1, props = c(0.5, 0.5), n = 50, seed = 7)
  expect_identical(r1, r2)
})

test_that("simulate_lpa: different seeds produce different data", {
  means <- matrix(c(1, 6), nrow = 1, ncol = 2)
  r1 <- simulate_lpa(means = means, sds = 1, props = c(0.5, 0.5), n = 50, seed = 1)
  r2 <- simulate_lpa(means = means, sds = 1, props = c(0.5, 0.5), n = 50, seed = 2)
  expect_false(identical(r1$data$y1, r2$data$y1))
})

test_that("simulate_lpa: props normalised internally", {
  means <- matrix(c(0, 10), nrow = 1, ncol = 2)
  expect_no_error(
    simulate_lpa(means = means, sds = 1, props = c(1, 3), n = 100, seed = 1)
  )
})

test_that("simulate_lpa: errors on dimension mismatch between means and props", {
  means <- matrix(c(0, 5, 0, 5), nrow = 2, ncol = 2)
  expect_error(
    simulate_lpa(means = means, sds = 1, props = c(0.3, 0.3, 0.4), n = 100),
    regexp = "props"
  )
})

test_that("simulate_lpa: matrix sds accepted", {
  means <- matrix(c(0, 5, 0, 5), nrow = 2, ncol = 2)
  sds   <- matrix(c(1, 2, 1, 2), nrow = 2, ncol = 2)
  expect_no_error(
    simulate_lpa(means = means, sds = sds, props = c(0.5, 0.5), n = 100, seed = 1)
  )
})


# ===========================================================================
# simulate_lca
# ===========================================================================

test_that("simulate_lca: returns list with $data and $params", {
  item_probs  <- matrix(c(0.9, 0.1, 0.9, 0.1,
                          0.1, 0.9, 0.1, 0.9), nrow = 4, ncol = 2)
  r <- simulate_lca(item_probs = item_probs, class_probs = c(0.5, 0.5),
                    n = 100, seed = 1)
  expect_type(r, "list")
  expect_true(all(c("data", "params") %in% names(r)))
})

test_that("simulate_lca: $data is data.frame with correct columns", {
  item_probs  <- matrix(c(0.9, 0.1, 0.9, 0.1,
                          0.1, 0.9, 0.1, 0.9), nrow = 4, ncol = 2)
  r <- simulate_lca(item_probs = item_probs, class_probs = c(0.5, 0.5),
                    n = 100, seed = 1)
  expect_s3_class(r$data, "data.frame")
  expect_true(all(c("item1", "item2", "item3", "item4", "true_class") %in% names(r$data)))
  expect_equal(nrow(r$data), 100L)
})

test_that("simulate_lca: all item columns are binary (0/1)", {
  item_probs  <- matrix(c(0.8, 0.2, 0.8, 0.2, 0.2, 0.8, 0.2, 0.8),
                        nrow = 4, ncol = 2)
  r <- simulate_lca(item_probs = item_probs, class_probs = c(0.5, 0.5),
                    n = 200, seed = 1)
  item_cols <- setdiff(names(r$data), "true_class")
  vals <- unlist(r$data[, item_cols])
  expect_true(all(vals %in% c(0L, 1L)))
})

test_that("simulate_lca: true_class contains only valid labels", {
  item_probs  <- matrix(c(0.9, 0.1, 0.5,
                          0.1, 0.9, 0.5), nrow = 2, ncol = 3)
  r <- simulate_lca(item_probs = item_probs, class_probs = c(0.4, 0.4, 0.2),
                    n = 300, seed = 1)
  expect_true(all(r$data$true_class %in% 1:3))
})

test_that("simulate_lca: $params holds exact input parameters", {
  item_probs  <- matrix(c(0.9, 0.1, 0.1, 0.9), nrow = 2, ncol = 2)
  class_probs <- c(0.6, 0.4)
  r  <- simulate_lca(item_probs = item_probs, class_probs = class_probs,
                     n = 100, seed = 1)
  expect_identical(r$params$item_probs,  item_probs)
  expect_identical(r$params$class_probs, class_probs)
})

test_that("simulate_lca: item probabilities approximately recovered at large n", {
  item_probs  <- matrix(c(0.9, 0.1, 0.1, 0.9), nrow = 2, ncol = 2)
  r  <- simulate_lca(item_probs = item_probs, class_probs = c(0.5, 0.5),
                     n = 3000, seed = 42)
  cl1 <- r$data[r$data$true_class == 1L, ]
  expect_true(abs(mean(cl1$item1) - 0.9) < 0.05)
  expect_true(abs(mean(cl1$item2) - 0.1) < 0.05)
})

test_that("simulate_lca: seed reproducibility", {
  item_probs  <- matrix(c(0.8, 0.2, 0.2, 0.8), nrow = 2, ncol = 2)
  r1 <- simulate_lca(item_probs = item_probs, class_probs = c(0.5, 0.5),
                     n = 50, seed = 5)
  r2 <- simulate_lca(item_probs = item_probs, class_probs = c(0.5, 0.5),
                     n = 50, seed = 5)
  expect_identical(r1, r2)
})

test_that("simulate_lca: errors on dimension mismatch", {
  item_probs  <- matrix(c(0.9, 0.1, 0.1, 0.9), nrow = 2, ncol = 2)
  expect_error(
    simulate_lca(item_probs = item_probs, class_probs = c(0.4, 0.3, 0.3),
                 n = 100),
    regexp = "class_probs"
  )
})

test_that("simulate_lca: errors on out-of-range item_probs", {
  item_probs <- matrix(c(1.2, 0.1, 0.1, 0.9), nrow = 2, ncol = 2)
  expect_error(
    simulate_lca(item_probs = item_probs, class_probs = c(0.5, 0.5), n = 100),
    regexp = "item_probs"
  )
})


# ===========================================================================
# simulate_regression
# ===========================================================================

test_that("simulate_regression: returns list with $data and $params", {
  coefs <- c("(Intercept)" = 1, x1 = 2, x2 = -1)
  r <- simulate_regression(coefs = coefs, predictor_sds = c(x1 = 1, x2 = 1),
                            error_sd = 0.5, n = 100, seed = 1)
  expect_type(r, "list")
  expect_true(all(c("data", "params") %in% names(r)))
})

test_that("simulate_regression: $data is data.frame with correct columns", {
  coefs <- c("(Intercept)" = 1, x1 = 2, x2 = -1)
  r <- simulate_regression(coefs = coefs, predictor_sds = c(x1 = 1, x2 = 1),
                            error_sd = 0.5, n = 100, seed = 1)
  expect_s3_class(r$data, "data.frame")
  expect_true(all(c("y", "x1", "x2") %in% names(r$data)))
  expect_equal(nrow(r$data), 100L)
})

test_that("simulate_regression: y is numeric, predictors are numeric", {
  coefs <- c("(Intercept)" = 0, x1 = 3)
  r <- simulate_regression(coefs = coefs, predictor_sds = c(x1 = 1),
                            error_sd = 1, n = 50, seed = 1)
  expect_true(is.numeric(r$data$y))
  expect_true(is.numeric(r$data$x1))
})

test_that("simulate_regression: $params holds exact input", {
  coefs         <- c("(Intercept)" = 5, x1 = 2, x2 = -3)
  predictor_sds <- c(x1 = 1, x2 = 2)
  r  <- simulate_regression(coefs = coefs, predictor_sds = predictor_sds,
                              error_sd = 1, n = 100, seed = 1)
  expect_identical(r$params$coefs,         coefs)
  expect_identical(r$params$predictor_sds, predictor_sds)
  expect_equal(r$params$error_sd, 1)
})

test_that("simulate_regression: lm() recovers intercept and slope at large n", {
  coefs <- c("(Intercept)" = 3, x1 = 5)
  r     <- simulate_regression(coefs = coefs, predictor_sds = c(x1 = 1),
                                error_sd = 0.1, n = 5000, seed = 42)
  fit   <- lm(y ~ x1, data = r$data)
  est   <- coef(fit)
  expect_true(abs(est["(Intercept)"] - 3) < 0.1)
  expect_true(abs(est["x1"] - 5) < 0.1)
})

test_that("simulate_regression: lm() recovers multiple coefficients at large n", {
  coefs <- c("(Intercept)" = 1, x1 = 2, x2 = -3, x3 = 0.5)
  r     <- simulate_regression(coefs = coefs,
                                predictor_sds = c(x1 = 1, x2 = 1, x3 = 1),
                                error_sd = 0.5, n = 3000, seed = 42)
  fit   <- lm(y ~ ., data = r$data)
  est   <- coef(fit)
  expect_true(abs(est["(Intercept)"] - 1)  < 0.15)
  expect_true(abs(est["x1"] - 2)           < 0.15)
  expect_true(abs(est["x2"] - (-3))        < 0.15)
  expect_true(abs(est["x3"] - 0.5)         < 0.15)
})

test_that("simulate_regression: seed reproducibility", {
  coefs <- c("(Intercept)" = 0, x1 = 1)
  r1 <- simulate_regression(coefs = coefs, predictor_sds = c(x1 = 1),
                              error_sd = 1, n = 50, seed = 3)
  r2 <- simulate_regression(coefs = coefs, predictor_sds = c(x1 = 1),
                              error_sd = 1, n = 50, seed = 3)
  expect_identical(r1, r2)
})

test_that("simulate_regression: different seeds produce different data", {
  coefs <- c("(Intercept)" = 0, x1 = 1)
  r1 <- simulate_regression(coefs = coefs, predictor_sds = c(x1 = 1),
                              error_sd = 1, n = 50, seed = 1)
  r2 <- simulate_regression(coefs = coefs, predictor_sds = c(x1 = 1),
                              error_sd = 1, n = 50, seed = 2)
  expect_false(identical(r1$data$y, r2$data$y))
})

test_that("simulate_regression: errors when predictor_sds names don't match coefs", {
  coefs <- c("(Intercept)" = 0, x1 = 1, x2 = 2)
  expect_error(
    simulate_regression(coefs = coefs, predictor_sds = c(x1 = 1),
                        error_sd = 1, n = 100),
    regexp = "predictor_sds"
  )
})

test_that("simulate_regression: errors on non-positive error_sd", {
  coefs <- c("(Intercept)" = 0, x1 = 1)
  expect_error(
    simulate_regression(coefs = coefs, predictor_sds = c(x1 = 1),
                        error_sd = -1, n = 100),
    regexp = "error_sd"
  )
})

test_that("simulate_regression: works without intercept in coefs", {
  coefs <- c(x1 = 2, x2 = -1)
  expect_no_error(
    simulate_regression(coefs = coefs, predictor_sds = c(x1 = 1, x2 = 1),
                        error_sd = 0.5, n = 100, seed = 1)
  )
})


# ===========================================================================
# simulate_fa
# ===========================================================================

test_that("simulate_fa: returns list with $data and $params", {
  loadings <- matrix(c(0.8, 0.7, 0.6, 0, 0, 0,
                       0,   0,   0,   0.8, 0.7, 0.6), nrow = 6, ncol = 2)
  r <- simulate_fa(loadings = loadings, n = 100, seed = 1)
  expect_type(r, "list")
  expect_true(all(c("data", "params") %in% names(r)))
})

test_that("simulate_fa: $data has columns y1..yp", {
  loadings <- matrix(c(0.8, 0.7, 0.6, 0, 0, 0,
                       0,   0,   0,   0.8, 0.7, 0.6), nrow = 6, ncol = 2)
  r <- simulate_fa(loadings = loadings, n = 100, seed = 1)
  expect_s3_class(r$data, "data.frame")
  expect_true(all(paste0("y", 1:6) %in% names(r$data)))
  expect_equal(nrow(r$data), 100L)
  expect_equal(ncol(r$data), 6L)
})

test_that("simulate_fa: $params contains loadings, phi, psi, sigma_implied", {
  loadings <- matrix(c(0.8, 0.7, 0.6, 0, 0, 0,
                       0,   0,   0,   0.8, 0.7, 0.6), nrow = 6, ncol = 2)
  r <- simulate_fa(loadings = loadings, n = 100, seed = 1)
  expect_true(all(c("loadings", "phi", "psi", "sigma_implied") %in% names(r$params)))
})

test_that("simulate_fa: $params$loadings is identical to input", {
  loadings <- matrix(c(0.8, 0.7, 0.6, 0, 0, 0,
                       0,   0,   0,   0.8, 0.7, 0.6), nrow = 6, ncol = 2)
  r <- simulate_fa(loadings = loadings, n = 100, seed = 1)
  expect_identical(r$params$loadings, loadings)
})

test_that("simulate_fa: phi defaults to identity matrix", {
  loadings <- matrix(c(0.7, 0.6, 0, 0, 0.7, 0.6), nrow = 3, ncol = 2)
  r <- simulate_fa(loadings = loadings, n = 100, seed = 1)
  expect_equal(r$params$phi, diag(2))
})

test_that("simulate_fa: psi is auto-computed and all positive", {
  loadings <- matrix(c(0.8, 0.7, 0.6, 0, 0, 0,
                       0,   0,   0,   0.8, 0.7, 0.6), nrow = 6, ncol = 2)
  r <- simulate_fa(loadings = loadings, n = 100, seed = 1)
  expect_true(all(r$params$psi > 0))
})

test_that("simulate_fa: psi equals 1 - communalities for orthogonal model", {
  loadings <- matrix(c(0.8, 0.6, 0, 0, 0.7, 0.5), nrow = 3, ncol = 2)
  r        <- simulate_fa(loadings = loadings, n = 100, seed = 1)
  expected_psi <- 1 - diag(loadings %*% diag(2) %*% t(loadings))
  expect_equal(r$params$psi, expected_psi)
})

test_that("simulate_fa: sigma_implied equals L*Phi*t(L) + diag(psi)", {
  loadings <- matrix(c(0.8, 0.7, 0.6, 0, 0, 0,
                       0,   0,   0,   0.8, 0.7, 0.6), nrow = 6, ncol = 2)
  r     <- simulate_fa(loadings = loadings, n = 100, seed = 1)
  Sigma <- loadings %*% r$params$phi %*% t(loadings) + diag(r$params$psi)
  expect_true(isTRUE(all.equal(r$params$sigma_implied, Sigma, tolerance = 1e-10)))
})

test_that("simulate_fa: sample covariance approximates sigma_implied at large n", {
  loadings <- matrix(c(0.8, 0.7, 0.6, 0, 0, 0,
                       0,   0,   0,   0.8, 0.7, 0.6), nrow = 6, ncol = 2)
  r <- simulate_fa(loadings = loadings, n = 5000, seed = 42)
  # Use absolute tolerance (max element-wise difference)
  expect_true(max(abs(cov(r$data) - r$params$sigma_implied)) < 0.1)
})

test_that("simulate_fa: seed produces identical results", {
  loadings <- matrix(c(0.8, 0.7, 0.6, 0, 0, 0,
                       0,   0,   0,   0.8, 0.7, 0.6), nrow = 6, ncol = 2)
  r1 <- simulate_fa(loadings = loadings, n = 100, seed = 9)
  r2 <- simulate_fa(loadings = loadings, n = 100, seed = 9)
  expect_identical(r1$data, r2$data)
  expect_true(isTRUE(all.equal(r1$params$sigma_implied, r2$params$sigma_implied)))
})

test_that("simulate_fa: different seeds produce different data", {
  loadings <- matrix(c(0.8, 0.7, 0.6, 0, 0, 0,
                       0,   0,   0,   0.8, 0.7, 0.6), nrow = 6, ncol = 2)
  r1 <- simulate_fa(loadings = loadings, n = 100, seed = 1)
  r2 <- simulate_fa(loadings = loadings, n = 100, seed = 2)
  expect_false(identical(r1$data, r2$data))
})

test_that("simulate_fa: oblique model accepted (phi != I)", {
  loadings <- matrix(c(0.8, 0.7, 0.6, 0, 0, 0,
                       0,   0,   0,   0.8, 0.7, 0.6), nrow = 6, ncol = 2)
  phi <- matrix(c(1, 0.4, 0.4, 1), nrow = 2)
  r   <- simulate_fa(loadings = loadings, phi = phi, n = 100, seed = 1)
  expect_identical(r$params$phi, phi)
})

test_that("simulate_fa: user-supplied psi respected", {
  loadings <- matrix(c(0.7, 0.6, 0.5, 0, 0, 0,
                       0,   0,   0,   0.7, 0.6, 0.5), nrow = 6, ncol = 2)
  psi <- rep(0.3, 6)
  r   <- simulate_fa(loadings = loadings, psi = psi, n = 100, seed = 1)
  expect_equal(r$params$psi, psi)
})

test_that("simulate_fa: single-factor model works (matrix input)", {
  loadings <- matrix(c(0.8, 0.7, 0.6), nrow = 3, ncol = 1)
  expect_no_error(simulate_fa(loadings = loadings, n = 100, seed = 1))
})

test_that("simulate_fa: errors on non-matrix loadings", {
  expect_error(
    simulate_fa(loadings = c(0.8, 0.7, 0.6), n = 100),
    regexp = "matrix"
  )
})

test_that("simulate_fa: errors on phi dimension mismatch", {
  loadings <- matrix(c(0.8, 0.7, 0.6, 0, 0, 0,
                       0,   0,   0,   0.8, 0.7, 0.6), nrow = 6, ncol = 2)
  phi_bad  <- diag(3)  # wrong: 3×3 for 2-factor model
  expect_error(
    simulate_fa(loadings = loadings, phi = phi_bad, n = 100),
    regexp = "phi"
  )
})

test_that("simulate_fa: errors on non-positive-definite phi", {
  loadings <- matrix(c(0.8, 0.7, 0.6, 0, 0, 0,
                       0,   0,   0,   0.8, 0.7, 0.6), nrow = 6, ncol = 2)
  phi_npd  <- matrix(c(1, 1, 1, 1), nrow = 2)  # singular
  expect_error(
    simulate_fa(loadings = loadings, phi = phi_npd, n = 100),
    regexp = "phi"
  )
})

test_that("simulate_fa: errors when user psi has non-positive values", {
  loadings <- matrix(c(0.8, 0.7, 0.6, 0, 0, 0,
                       0,   0,   0,   0.8, 0.7, 0.6), nrow = 6, ncol = 2)
  expect_error(
    simulate_fa(loadings = loadings, psi = c(0.5, 0.5, 0.5, 0.5, 0.5, -0.1),
                n = 100),
    regexp = "psi"
  )
})

test_that("simulate_fa: errors when auto-psi has non-positive values (over-factored)", {
  # Loadings > 1 causes communality > 1, psi < 0
  loadings <- matrix(c(1.5, 0.8), nrow = 2, ncol = 1)
  expect_error(
    simulate_fa(loadings = loadings, n = 100),
    regexp = "psi"
  )
})
