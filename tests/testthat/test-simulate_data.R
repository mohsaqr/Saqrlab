# ---- Input validation ----

test_that("simulate_data rejects invalid type", {
  expect_error(simulate_data("bogus"))
  expect_error(simulate_data(123))
  expect_error(simulate_data(c("ttest", "anova")))
})

test_that("simulate_data rejects invalid seed", {
  expect_error(simulate_data("ttest", seed = "abc"))
  expect_error(simulate_data("ttest", seed = NA))
  expect_error(simulate_data("ttest", seed = c(1, 2)))
})

test_that("simulate_data works with NULL seed", {
  d <- simulate_data("ttest", seed = NULL)
  expect_s3_class(d, "data.frame")
  expect_gt(nrow(d), 0L)
})

# ---- Per-type structure: ttest ----

test_that("ttest has correct columns and types", {
  d <- simulate_data("ttest", seed = 1)
  expect_true(all(c("group", "score") %in% names(d)))
  expect_s3_class(d$group, "factor")
  expect_true(is.numeric(d$score))
})

test_that("ttest has 10-500 rows with no NAs", {
  d <- simulate_data("ttest", seed = 2)
  expect_gte(nrow(d), 10L)
  expect_lte(nrow(d), 500L)
  expect_false(anyNA(d))
})

test_that("ttest has exactly two groups", {
  d <- simulate_data("ttest", seed = 3)
  expect_equal(nlevels(d$group), 2L)
  expect_true(all(levels(d$group) %in% c("A", "B")))
})

test_that("ttest is reproducible with same seed", {
  d1 <- simulate_data("ttest", seed = 10)
  d2 <- simulate_data("ttest", seed = 10)
  expect_identical(d1, d2)
})

test_that("ttest differs with different seeds", {
  d1 <- simulate_data("ttest", seed = 10)
  d2 <- simulate_data("ttest", seed = 11)
  expect_false(identical(d1, d2))
})

# ---- Per-type structure: anova ----

test_that("anova has correct columns and types", {
  d <- simulate_data("anova", seed = 1)
  expect_true(all(c("group", "score") %in% names(d)))
  expect_s3_class(d$group, "factor")
  expect_true(is.numeric(d$score))
})

test_that("anova has 10-500 rows with no NAs", {
  d <- simulate_data("anova", seed = 2)
  expect_gte(nrow(d), 10L)
  expect_lte(nrow(d), 500L)
  expect_false(anyNA(d))
})

test_that("anova has 3-5 groups", {
  d <- simulate_data("anova", seed = 3)
  expect_gte(nlevels(d$group), 3L)
  expect_lte(nlevels(d$group), 5L)
})

test_that("anova is reproducible with same seed", {
  d1 <- simulate_data("anova", seed = 20)
  d2 <- simulate_data("anova", seed = 20)
  expect_identical(d1, d2)
})

test_that("anova differs with different seeds", {
  d1 <- simulate_data("anova", seed = 20)
  d2 <- simulate_data("anova", seed = 21)
  expect_false(identical(d1, d2))
})

# ---- Per-type structure: correlation ----

test_that("correlation has correct columns and types", {
  d <- simulate_data("correlation", seed = 1)
  expect_true(all(grepl("^x\\d+$", names(d))))
  expect_true(all(vapply(d, is.numeric, logical(1L))))
  expect_gte(ncol(d), 4L)
  expect_lte(ncol(d), 7L)
})

test_that("correlation has 10-500 rows with no NAs", {
  d <- simulate_data("correlation", seed = 2)
  expect_gte(nrow(d), 10L)
  expect_lte(nrow(d), 500L)
  expect_false(anyNA(d))
})

test_that("correlation is reproducible with same seed", {
  d1 <- simulate_data("correlation", seed = 30)
  d2 <- simulate_data("correlation", seed = 30)
  expect_identical(d1, d2)
})

test_that("correlation differs with different seeds", {
  d1 <- simulate_data("correlation", seed = 30)
  d2 <- simulate_data("correlation", seed = 31)
  expect_false(identical(d1, d2))
})

# ---- Per-type structure: clusters ----

test_that("clusters has correct columns and types", {
  d <- simulate_data("clusters", seed = 1)
  numeric_cols <- grep("^x\\d+$", names(d), value = TRUE)
  expect_gte(length(numeric_cols), 2L)
  expect_lte(length(numeric_cols), 5L)
  expect_true("true_cluster" %in% names(d))
  expect_true(is.integer(d$true_cluster) || is.numeric(d$true_cluster))
  expect_true(all(vapply(d[numeric_cols], is.numeric, logical(1L))))
})

test_that("clusters has 10-500 rows with no NAs", {
  d <- simulate_data("clusters", seed = 2)
  expect_gte(nrow(d), 10L)
  expect_lte(nrow(d), 500L)
  expect_false(anyNA(d))
})

test_that("clusters has 3-5 true clusters", {
  d <- simulate_data("clusters", seed = 3)
  k <- length(unique(d$true_cluster))
  expect_gte(k, 3L)
  expect_lte(k, 5L)
})

test_that("clusters is reproducible with same seed", {
  d1 <- simulate_data("clusters", seed = 40)
  d2 <- simulate_data("clusters", seed = 40)
  expect_identical(d1, d2)
})

test_that("clusters differs with different seeds", {
  d1 <- simulate_data("clusters", seed = 40)
  d2 <- simulate_data("clusters", seed = 41)
  expect_false(identical(d1, d2))
})

# ---- Per-type structure: factor_analysis ----

test_that("factor_analysis has correct columns and types", {
  d <- simulate_data("factor_analysis", seed = 1)
  expect_true(all(grepl("^x\\d+$", names(d))))
  expect_true(all(vapply(d, is.numeric, logical(1L))))
  expect_gte(ncol(d), 6L)
  expect_lte(ncol(d), 16L)
})

test_that("factor_analysis has 50-500 rows with no NAs", {
  d <- simulate_data("factor_analysis", seed = 2)
  expect_gte(nrow(d), 50L)
  expect_lte(nrow(d), 500L)
  expect_false(anyNA(d))
})

test_that("factor_analysis is reproducible with same seed", {
  d1 <- simulate_data("factor_analysis", seed = 50)
  d2 <- simulate_data("factor_analysis", seed = 50)
  expect_identical(d1, d2)
})

test_that("factor_analysis differs with different seeds", {
  d1 <- simulate_data("factor_analysis", seed = 50)
  d2 <- simulate_data("factor_analysis", seed = 51)
  expect_false(identical(d1, d2))
})

# ---- Per-type structure: prediction ----

test_that("prediction has correct columns and types", {
  d <- simulate_data("prediction", seed = 1)
  expect_true(all(c("y", "x1", "x2", "x3", "x4", "cat1", "cat2") %in% names(d)))
  expect_true(is.numeric(d$y))
  expect_true(is.numeric(d$x1))
  expect_true(is.numeric(d$x4))
  expect_s3_class(d$cat1, "factor")
  expect_s3_class(d$cat2, "factor")
})

test_that("prediction has 10-500 rows with no NAs", {
  d <- simulate_data("prediction", seed = 2)
  expect_gte(nrow(d), 10L)
  expect_lte(nrow(d), 500L)
  expect_false(anyNA(d))
})

test_that("prediction categorical levels are in range", {
  d <- simulate_data("prediction", seed = 3)
  expect_gte(nlevels(d$cat1), 2L)
  expect_lte(nlevels(d$cat1), 4L)
  expect_gte(nlevels(d$cat2), 2L)
  expect_lte(nlevels(d$cat2), 3L)
})

test_that("prediction is reproducible with same seed", {
  d1 <- simulate_data("prediction", seed = 60)
  d2 <- simulate_data("prediction", seed = 60)
  expect_identical(d1, d2)
})

test_that("prediction differs with different seeds", {
  d1 <- simulate_data("prediction", seed = 60)
  d2 <- simulate_data("prediction", seed = 61)
  expect_false(identical(d1, d2))
})

# ---- Analysis correctness ----

test_that("ttest data yields significant t-test for multiple seeds", {
  seeds <- c(1, 42, 100, 200, 300)
  p_vals <- vapply(seeds, function(s) {
    d <- simulate_data("ttest", seed = s)
    t.test(score ~ group, data = d)$p.value
  }, numeric(1L))
  # At least 3 out of 5 should be significant

  expect_gte(sum(p_vals < 0.05), 3L)
})

test_that("anova data yields significant F-test for multiple seeds", {
  seeds <- c(1, 7, 50, 150, 250)
  p_vals <- vapply(seeds, function(s) {
    d <- simulate_data("anova", seed = s)
    summary(aov(score ~ group, data = d))[[1]][["Pr(>F)"]][1L]
  }, numeric(1L))
  expect_gte(sum(p_vals < 0.05), 3L)
})

test_that("correlation data has at least one strong correlation", {
  seeds <- c(1, 5, 20, 80, 150)
  has_strong <- vapply(seeds, function(s) {
    d <- simulate_data("correlation", seed = s)
    R <- cor(d)
    diag(R) <- 0
    max(abs(R)) > 0.4
  }, logical(1L))
  expect_gte(sum(has_strong), 3L)
})

test_that("clusters data allows reasonable kmeans recovery", {
  seeds <- c(1, 99, 42, 200, 300)
  recoveries <- vapply(seeds, function(s) {
    d <- simulate_data("clusters", seed = s)
    k <- max(d$true_cluster)
    features <- d[, grep("^x\\d+$", names(d)), drop = FALSE]
    km <- kmeans(features, centers = k, nstart = 10L)
    # Compute cluster recovery via contingency table purity
    tab <- table(km$cluster, d$true_cluster)
    sum(apply(tab, 1L, max)) / nrow(d)
  }, numeric(1L))
  # At least 3 out of 5 should achieve >0.6 purity
  expect_gte(sum(recoveries > 0.6), 3L)
})

test_that("factor_analysis data allows factanal to converge", {
  seeds <- c(1, 12, 50, 100, 200)
  converged <- vapply(seeds, function(s) {
    d <- simulate_data("factor_analysis", seed = s)
    nf <- attr(d, "n_factors")
    tryCatch(
      {
        factanal(d, factors = nf)
        TRUE
      },
      error = function(e) FALSE
    )
  }, logical(1L))
  expect_gte(sum(converged), 3L)
})

test_that("prediction data yields lm with R-squared > 0.2", {
  seeds <- c(1, 5, 30, 100, 250)
  r2s <- vapply(seeds, function(s) {
    d <- simulate_data("prediction", seed = s)
    fit <- lm(y ~ ., data = d)
    summary(fit)$r.squared
  }, numeric(1L))
  expect_gte(sum(r2s > 0.2), 3L)
})

# ---- Attributes ----

test_that("attr(d, 'type') is correct for all types", {
  types <- c("ttest", "anova", "correlation", "clusters",
    "factor_analysis", "prediction")
  vapply(types, function(tp) {
    d <- simulate_data(tp, seed = 1)
    expect_equal(attr(d, "type"), tp)
    TRUE
  }, logical(1L))
})

test_that("attr(d, 'info') exists for all types", {
  types <- c("ttest", "anova", "correlation", "clusters",
    "factor_analysis", "prediction")
  vapply(types, function(tp) {
    d <- simulate_data(tp, seed = 1)
    info <- attr(d, "info")
    expect_true(is.character(info))
    expect_gt(nchar(info), 0L)
    TRUE
  }, logical(1L))
})

test_that("factor_analysis has n_factors and loadings attributes", {
  d <- simulate_data("factor_analysis", seed = 12)
  nf <- attr(d, "n_factors")
  loadings <- attr(d, "loadings")
  expect_true(is.numeric(nf))
  expect_gte(nf, 2L)
  expect_lte(nf, 4L)
  expect_true(is.matrix(loadings))
  expect_equal(ncol(loadings), nf)
  expect_equal(nrow(loadings), ncol(d))
})

# ---- Parameter overrides via ... ----

test_that("n override works for all types", {
  types <- c("ttest", "anova", "correlation", "clusters",
    "factor_analysis", "prediction")
  vapply(types, function(tp) {
    d <- simulate_data(tp, seed = 1, n = 100)
    expect_equal(nrow(d), 100L)
    TRUE
  }, logical(1L))
})

test_that("ttest effect_size override works", {
  d1 <- simulate_data("ttest", seed = 1, n = 300, effect_size = 0.3)
  d2 <- simulate_data("ttest", seed = 1, n = 300, effect_size = 1.2)
  # Larger effect size should produce larger mean difference
  diff1 <- abs(diff(tapply(d1$score, d1$group, mean)))
  diff2 <- abs(diff(tapply(d2$score, d2$group, mean)))
  expect_gt(diff2, diff1)
})

test_that("anova n_groups override works", {
  d <- simulate_data("anova", seed = 1, n_groups = 5)
  expect_equal(nlevels(d$group), 5L)
})

test_that("correlation n_vars override works", {
  d <- simulate_data("correlation", seed = 1, n_vars = 6)
  expect_equal(ncol(d), 6L)
})

test_that("clusters n_clusters override works", {
  d <- simulate_data("clusters", seed = 1, n_clusters = 4)
  expect_equal(length(unique(d$true_cluster)), 4L)
})

test_that("clusters n_dims override works", {
  d <- simulate_data("clusters", seed = 1, n_dims = 3)
  expect_equal(sum(grepl("^x\\d+$", names(d))), 3L)
})

test_that("factor_analysis n_factors override works", {
  d <- simulate_data("factor_analysis", seed = 1, n_factors = 3)
  expect_equal(attr(d, "n_factors"), 3L)
  expect_equal(ncol(attr(d, "loadings")), 3L)
})

test_that("factor_analysis items_per_factor override works", {
  d <- simulate_data("factor_analysis", seed = 1, n_factors = 2, items_per_factor = 4)
  expect_equal(ncol(d), 8L)  # 2 * 4
})

test_that("prediction n_cat1 and n_cat2 overrides work", {
  d <- simulate_data("prediction", seed = 1, n_cat1 = 3, n_cat2 = 2)
  expect_equal(nlevels(d$cat1), 3L)
  expect_equal(nlevels(d$cat2), 2L)
})
