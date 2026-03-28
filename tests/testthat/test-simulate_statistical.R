# ===========================================================================
# simulate_ttest
# ===========================================================================

test_that("simulate_ttest: returns saqr_sim with correct structure", {
  r <- simulate_ttest(n_a = 30, n_b = 30, mean_a = 10, mean_b = 12, seed = 1)
  expect_s3_class(r, "saqr_sim")
  expect_equal(r$type, "ttest")
  expect_s3_class(r$data, "data.frame")
  expect_true(all(c("group", "score") %in% names(r$data)))
  expect_equal(nrow(r$data), 60L)
})

test_that("simulate_ttest: params contain all generating parameters", {
  r <- simulate_ttest(n_a = 50, n_b = 40, mean_a = 5, mean_b = 8,
                       sd_a = 2, sd_b = 3, seed = 1)
  expect_equal(r$params$mean_a, 5)
  expect_equal(r$params$mean_b, 8)
  expect_equal(r$params$sd_a, 2)
  expect_equal(r$params$sd_b, 3)
  expect_equal(r$params$n_a, 50L)
  expect_equal(r$params$n_b, 40L)
  expect_true(is.numeric(r$params$cohens_d))
})

test_that("simulate_ttest: group means approximately recoverable", {
  r <- simulate_ttest(n_a = 5000, n_b = 5000, mean_a = 50, mean_b = 55, seed = 42)
  d <- r$data
  obs_a <- mean(d$score[d$group == "A"])
  obs_b <- mean(d$score[d$group == "B"])
  expect_true(abs(obs_a - 50) < 0.5)
  expect_true(abs(obs_b - 55) < 0.5)
})

test_that("simulate_ttest: seed reproducibility", {
  r1 <- simulate_ttest(n_a = 20, n_b = 20, mean_a = 0, mean_b = 1, seed = 7)
  r2 <- simulate_ttest(n_a = 20, n_b = 20, mean_a = 0, mean_b = 1, seed = 7)
  expect_identical(r1$data, r2$data)
})

test_that("simulate_ttest: custom labels", {
  r <- simulate_ttest(n_a = 10, n_b = 10, mean_a = 0, mean_b = 1,
                       labels = c("control", "treatment"), seed = 1)
  expect_equal(levels(r$data$group), c("control", "treatment"))
})

test_that("simulate_ttest: sd_b defaults to sd_a", {
  r <- simulate_ttest(n_a = 20, n_b = 20, mean_a = 0, mean_b = 1,
                       sd_a = 3, seed = 1)
  expect_equal(r$params$sd_b, 3)
})

test_that("simulate_ttest: Cohen's d computed correctly", {
  r <- simulate_ttest(n_a = 100, n_b = 100, mean_a = 0, mean_b = 1,
                       sd_a = 1, sd_b = 1, seed = 1)
  expect_equal(r$params$cohens_d, 1.0, tolerance = 1e-10)
})

test_that("simulate_ttest: validation errors", {
  expect_error(simulate_ttest(n_a = 1, n_b = 10, mean_a = 0, mean_b = 1))
  expect_error(simulate_ttest(n_a = 10, n_b = 10, mean_a = 0, mean_b = 1, sd_a = -1))
})


# ===========================================================================
# simulate_anova
# ===========================================================================

test_that("simulate_anova: returns saqr_sim with correct structure", {
  r <- simulate_anova(n = 30, means = c(10, 12, 15), seed = 1)
  expect_s3_class(r, "saqr_sim")
  expect_equal(r$type, "anova")
  expect_true(all(c("group", "score") %in% names(r$data)))
  expect_equal(nrow(r$data), 90L)  # 30 per group × 3
  expect_equal(nlevels(r$data$group), 3L)
})

test_that("simulate_anova: params contain all generating parameters", {
  r <- simulate_anova(n = 20, means = c(5, 10, 15), sds = c(1, 2, 3), seed = 1)
  expect_equal(as.numeric(r$params$means), c(5, 10, 15))
  expect_equal(as.numeric(r$params$sds), c(1, 2, 3))
  expect_equal(as.numeric(r$params$n), c(20, 20, 20))
  expect_true(is.numeric(r$params$eta_squared))
  expect_true(r$params$eta_squared > 0 && r$params$eta_squared < 1)
})

test_that("simulate_anova: unequal group sizes", {
  r <- simulate_anova(n = c(50, 30, 20), means = c(0, 0, 5), seed = 1)
  expect_equal(nrow(r$data), 100L)
  tab <- table(r$data$group)
  expect_equal(as.integer(tab), c(50L, 30L, 20L))
})

test_that("simulate_anova: group means approximately recoverable", {
  r <- simulate_anova(n = 5000, means = c(10, 20, 30), sds = 1, seed = 42)
  obs <- tapply(r$data$score, r$data$group, mean)
  expect_true(abs(obs["G1"] - 10) < 0.2)
  expect_true(abs(obs["G2"] - 20) < 0.2)
  expect_true(abs(obs["G3"] - 30) < 0.2)
})

test_that("simulate_anova: custom labels", {
  r <- simulate_anova(n = 10, means = c(1, 2),
                       labels = c("placebo", "drug"), seed = 1)
  expect_equal(levels(r$data$group), c("placebo", "drug"))
})

test_that("simulate_anova: seed reproducibility", {
  r1 <- simulate_anova(n = 20, means = c(1, 5), seed = 99)
  r2 <- simulate_anova(n = 20, means = c(1, 5), seed = 99)
  expect_identical(r1$data, r2$data)
})

test_that("simulate_anova: validation errors", {
  expect_error(simulate_anova(n = 1, means = c(1, 2)))
  expect_error(simulate_anova(n = c(10, 10), means = c(1, 2, 3)))  # length mismatch
})


# ===========================================================================
# simulate_correlation
# ===========================================================================

test_that("simulate_correlation: returns saqr_sim with correct structure", {
  R <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
  r <- simulate_correlation(n = 100, sigma = R, seed = 1)
  expect_s3_class(r, "saqr_sim")
  expect_equal(r$type, "correlation")
  expect_equal(ncol(r$data), 2L)
  expect_equal(nrow(r$data), 100L)
  expect_equal(names(r$data), c("x1", "x2"))
})

test_that("simulate_correlation: params contain input matrix", {
  R <- matrix(c(1, 0.6, 0.3, 0.6, 1, 0.5, 0.3, 0.5, 1), nrow = 3)
  r <- simulate_correlation(n = 50, sigma = R, seed = 1)
  expect_equal(r$params$sigma, R)
  expect_true(r$params$is_correlation)
})

test_that("simulate_correlation: correlation approximately recoverable", {
  R <- matrix(c(1, 0.7, 0.7, 1), nrow = 2)
  r <- simulate_correlation(n = 10000, sigma = R, seed = 42)
  obs <- cor(r$data)
  expect_true(abs(obs[1, 2] - 0.7) < 0.05)
})

test_that("simulate_correlation: custom means shift data", {
  R <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
  r <- simulate_correlation(n = 5000, sigma = R, means = c(100, 200), seed = 1)
  expect_true(abs(mean(r$data$x1) - 100) < 2)
  expect_true(abs(mean(r$data$x2) - 200) < 2)
})

test_that("simulate_correlation: custom var_names", {
  R <- matrix(c(1, 0.3, 0.3, 1), nrow = 2)
  r <- simulate_correlation(n = 50, sigma = R,
                             var_names = c("IQ", "GPA"), seed = 1)
  expect_equal(names(r$data), c("IQ", "GPA"))
})

test_that("simulate_correlation: seed reproducibility", {
  R <- matrix(c(1, 0.4, 0.4, 1), nrow = 2)
  r1 <- simulate_correlation(n = 50, sigma = R, seed = 7)
  r2 <- simulate_correlation(n = 50, sigma = R, seed = 7)
  expect_identical(r1$data, r2$data)
})

test_that("simulate_correlation: covariance matrix (non-correlation) works", {
  S <- matrix(c(4, 1, 1, 9), nrow = 2)
  r <- simulate_correlation(n = 100, sigma = S, seed = 1)
  expect_false(r$params$is_correlation)
})


# ===========================================================================
# simulate_clusters
# ===========================================================================

test_that("simulate_clusters: returns saqr_sim with correct structure", {
  centers <- matrix(c(0, 0, 5, 5), nrow = 2, byrow = TRUE)
  r <- simulate_clusters(n = 100, centers = centers, seed = 1)
  expect_s3_class(r, "saqr_sim")
  expect_equal(r$type, "clusters")
  expect_true(all(c("x1", "x2", "true_cluster") %in% names(r$data)))
  expect_equal(nrow(r$data), 100L)
})

test_that("simulate_clusters: params contain centers and sds", {
  centers <- matrix(c(0, 0, 10, 10, 20, 20), nrow = 3, byrow = TRUE)
  r <- simulate_clusters(n = 90, centers = centers, sds = 0.5, seed = 1)
  expect_equal(r$params$centers, centers)
  expect_true(all(r$params$sds == 0.5))
  expect_equal(sum(r$params$n), 90L)
})

test_that("simulate_clusters: cluster centers approximately recoverable", {
  centers <- matrix(c(0, 0, 10, 10), nrow = 2, byrow = TRUE)
  r <- simulate_clusters(n = 10000, centers = centers, sds = 0.5, seed = 42)
  d <- r$data
  obs1 <- colMeans(d[d$true_cluster == 1, c("x1", "x2")])
  obs2 <- colMeans(d[d$true_cluster == 2, c("x1", "x2")])
  expect_true(all(abs(obs1 - c(0, 0)) < 0.2))
  expect_true(all(abs(obs2 - c(10, 10)) < 0.2))
})

test_that("simulate_clusters: per-cluster sizes via vector n", {
  centers <- matrix(c(0, 5, 10), nrow = 3, ncol = 1)
  r <- simulate_clusters(n = c(50, 30, 20), centers = centers, seed = 1)
  tab <- table(r$data$true_cluster)
  expect_equal(as.integer(tab), c(50L, 30L, 20L))
})

test_that("simulate_clusters: per-cluster sds as matrix", {
  centers <- matrix(c(0, 0, 5, 5), nrow = 2, byrow = TRUE)
  sd_mat <- matrix(c(0.5, 1.0, 2.0, 3.0), nrow = 2, byrow = TRUE)
  r <- simulate_clusters(n = 100, centers = centers, sds = sd_mat, seed = 1)
  expect_equal(r$params$sds, sd_mat)
})

test_that("simulate_clusters: seed reproducibility", {
  centers <- matrix(c(0, 10), nrow = 2, ncol = 1)
  r1 <- simulate_clusters(n = 50, centers = centers, seed = 42)
  r2 <- simulate_clusters(n = 50, centers = centers, seed = 42)
  expect_identical(r1$data, r2$data)
})

test_that("simulate_clusters: props controls allocation", {
  centers <- matrix(c(0, 10), nrow = 2, ncol = 1)
  r <- simulate_clusters(n = 1000, centers = centers,
                          props = c(0.9, 0.1), seed = 1)
  tab <- table(r$data$true_cluster)
  expect_true(tab[1] > tab[2] * 5)
})


# ===========================================================================
# simulate_prediction
# ===========================================================================

test_that("simulate_prediction: returns saqr_sim with correct structure", {
  r <- simulate_prediction(n = 100,
                            coefs = c("(Intercept)" = 5, x1 = 2, x2 = -1),
                            seed = 1)
  expect_s3_class(r, "saqr_sim")
  expect_equal(r$type, "prediction")
  expect_true(all(c("y", "x1", "x2") %in% names(r$data)))
  expect_equal(nrow(r$data), 100L)
})

test_that("simulate_prediction: params contain all inputs", {
  coefs <- c("(Intercept)" = 3, x1 = 1)
  r <- simulate_prediction(n = 50, coefs = coefs, error_sd = 2, seed = 1)
  expect_identical(r$params$coefs, coefs)
  expect_equal(r$params$error_sd, 2)
  expect_true(is.numeric(r$params$r_squared))
})

test_that("simulate_prediction: coefficients approximately recoverable", {
  coefs <- c("(Intercept)" = 5, x1 = 3, x2 = -2)
  r <- simulate_prediction(n = 10000, coefs = coefs, error_sd = 0.5, seed = 42)
  fit <- lm(y ~ x1 + x2, data = r$data)
  recovered <- coef(fit)
  expect_true(abs(recovered["(Intercept)"] - 5) < 0.2)
  expect_true(abs(recovered["x1"] - 3) < 0.1)
  expect_true(abs(recovered["x2"] - (-2)) < 0.1)
})

test_that("simulate_prediction: categorical predictors", {
  r <- simulate_prediction(
    n = 200,
    coefs = c("(Intercept)" = 0, x1 = 1),
    cat_levels = list(tx = c("A", "B", "C")),
    cat_effects = list(tx = c(0, 5, 10)),
    seed = 1
  )
  expect_s3_class(r$data$tx, "factor")
  expect_equal(levels(r$data$tx), c("A", "B", "C"))
  expect_equal(r$params$cat_effects$tx, c(0, 5, 10))
})

test_that("simulate_prediction: seed reproducibility", {
  coefs <- c(x1 = 1)
  r1 <- simulate_prediction(n = 50, coefs = coefs, seed = 7)
  r2 <- simulate_prediction(n = 50, coefs = coefs, seed = 7)
  expect_identical(r1$data, r2$data)
})

test_that("simulate_prediction: custom predictor means/sds", {
  coefs <- c(x1 = 1)
  r <- simulate_prediction(n = 5000, coefs = coefs,
                            predictor_means = c(x1 = 100),
                            predictor_sds = c(x1 = 10), seed = 1)
  expect_true(abs(mean(r$data$x1) - 100) < 1)
  expect_true(abs(stats::sd(r$data$x1) - 10) < 1)
})
