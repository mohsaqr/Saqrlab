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

# ---- complexity parameter validation ----

test_that("complexity rejects non-character", {
  expect_error(simulate_data("ttest", seed = 1, complexity = 123))
  expect_error(simulate_data("ttest", seed = 1, complexity = TRUE))
})

test_that("complexity rejects empty character vector", {
  expect_error(simulate_data("ttest", seed = 1, complexity = character(0)))
})

test_that("complexity rejects unknown case names", {
  expect_error(simulate_data("ttest", seed = 1, complexity = "not_a_case"))
  expect_error(simulate_data("ttest", seed = 1, complexity = c("na", "bogus")))
})

test_that("complexity = 'clean' is backward compatible with default", {
  d1 <- simulate_data("ttest", seed = 42)
  d2 <- simulate_data("ttest", seed = 42, complexity = "clean")
  expect_identical(d1, d2)
})

test_that("type = 'batch' is accepted as valid type", {
  expect_no_error(simulate_data("batch", seed = 1, n_batch = 2L))
})

# ---- complexity = "auto" single datasets ----

test_that("auto complexity returns data.frame with complexity attribute", {
  d <- simulate_data("ttest", seed = 7, complexity = "auto")
  expect_s3_class(d, "data.frame")
  expect_true(is.character(attr(d, "complexity")))
})

test_that("auto complexity is reproducible with same seed", {
  d1 <- simulate_data("ttest", seed = 99, complexity = "auto")
  d2 <- simulate_data("ttest", seed = 99, complexity = "auto")
  expect_identical(d1, d2)
  expect_identical(attr(d1, "complexity"), attr(d2, "complexity"))
})

test_that("auto complexity varies across seeds", {
  complexities <- vapply(1:10, function(s) {
    d <- simulate_data("ttest", seed = s, complexity = "auto")
    paste(attr(d, "complexity"), collapse = "+")
  }, character(1L))
  # Not all 10 should have the exact same complexity profile
  expect_gt(length(unique(complexities)), 1L)
})

test_that("clean dataset has empty complexity attribute", {
  d <- simulate_data("ttest", seed = 1, complexity = "clean")
  expect_equal(attr(d, "complexity"), character(0L))
})

# ---- individual edge cases ----

test_that("na complexity injects NAs into ttest data", {
  d <- simulate_data("ttest", seed = 5, complexity = "na")
  expect_true(anyNA(d$score))
})

test_that("na complexity injects NAs into correlation data", {
  d <- simulate_data("correlation", seed = 5, n = 100, complexity = "na")
  expect_true(anyNA(d))
})

test_that("outliers complexity injects extreme values into ttest data", {
  # Check that some value is more than 5 MAD from median (robust outlier check)
  found <- vapply(1:5, function(s) {
    d <- simulate_data("ttest", seed = s, n = 200, complexity = "outliers")
    med <- stats::median(d$score, na.rm = TRUE)
    mad <- stats::mad(d$score, na.rm = TRUE)
    any(abs(d$score - med) > 5 * mad, na.rm = TRUE)
  }, logical(1L))
  expect_gte(sum(found), 3L)
})

test_that("ties complexity creates tied values in ttest score", {
  # Ties means many repeated values: fewer unique than total rows
  found_ties <- vapply(1:5, function(s) {
    d <- simulate_data("ttest", seed = s, n = 100, complexity = "ties")
    length(unique(d$score)) < 0.5 * nrow(d)
  }, logical(1L))
  expect_gte(sum(found_ties), 3L)
})

test_that("duplicates complexity increases row count", {
  original <- simulate_data("ttest", seed = 10, n = 80, complexity = "clean")
  with_dups <- simulate_data("ttest", seed = 10, n = 80, complexity = "duplicates")
  expect_gt(nrow(with_dups), nrow(original))
})

test_that("constant_col complexity creates a zero-variance column in correlation data", {
  d <- simulate_data("correlation", seed = 3, n = 50, n_vars = 5, complexity = "constant_col")
  col_vars <- vapply(d, function(v) stats::var(v, na.rm = TRUE), numeric(1L))
  expect_true(any(col_vars == 0 | is.na(col_vars)))
})

test_that("tiny_n complexity produces small datasets", {
  found_small <- vapply(1:5, function(s) {
    d <- simulate_data("ttest", seed = s, complexity = "tiny_n")
    nrow(d) < 30L
  }, logical(1L))
  expect_gte(sum(found_small), 3L)
})

test_that("heavy_tailed complexity generates data without error", {
  expect_no_error(simulate_data("ttest", seed = 1, n = 50, complexity = "heavy_tailed"))
  expect_no_error(simulate_data("anova", seed = 1, n = 60, complexity = "heavy_tailed"))
  expect_no_error(simulate_data("correlation", seed = 1, n = 50, complexity = "heavy_tailed"))
})

test_that("heteroscedastic complexity produces unequal group variances in ttest", {
  found_hetero <- vapply(1:10, function(s) {
    d <- simulate_data("ttest", seed = s, n = 100, complexity = "heteroscedastic")
    vars <- tapply(d$score, d$group, stats::var)
    ratio <- max(vars) / min(vars)
    ratio > 4
  }, logical(1L))
  expect_gte(sum(found_hetero), 6L)
})

test_that("extreme_imbalance complexity creates very unequal group sizes in ttest", {
  found_imbal <- vapply(1:5, function(s) {
    d <- simulate_data("ttest", seed = s, n = 100, complexity = "extreme_imbalance")
    tab <- table(d$group)
    min(tab) < 0.15 * nrow(d)
  }, logical(1L))
  expect_gte(sum(found_imbal), 3L)
})

test_that("multicollinear complexity produces near-perfect correlation in correlation data", {
  found_mc <- vapply(1:5, function(s) {
    d <- simulate_data("correlation", seed = s, n = 100, complexity = "multicollinear")
    R <- cor(d)
    diag(R) <- 0
    max(abs(R)) > 0.90
  }, logical(1L))
  expect_gte(sum(found_mc), 4L)
})

test_that("complexity attribute records applied cases", {
  d <- simulate_data("ttest", seed = 1, complexity = c("na", "outliers"))
  expect_equal(sort(attr(d, "complexity")), c("na", "outliers"))
})

# ---- n_batch single type ----

test_that("n_batch returns a list of the correct length", {
  result <- simulate_data("ttest", seed = 1, n_batch = 10L)
  expect_type(result, "list")
  expect_length(result, 10L)
})

test_that("n_batch result has sim_batch_type class", {
  result <- simulate_data("ttest", seed = 1, n_batch = 5L)
  expect_s3_class(result, "sim_batch_type")
})

test_that("n_batch result attributes are correct", {
  result <- simulate_data("ttest", seed = 42, n_batch = 8L)
  expect_equal(attr(result, "type"), "ttest")
  expect_equal(attr(result, "n_batch"), 8L)
  expect_equal(attr(result, "seed"), 42L)
})

test_that("n_batch all elements are data.frames", {
  result <- simulate_data("ttest", seed = 1, n_batch = 5L)
  all_df <- vapply(result, is.data.frame, logical(1L))
  expect_true(all(all_df))
})

test_that("n_batch elements have batch_id and seed attributes", {
  result <- simulate_data("ttest", seed = 1, n_batch = 5L)
  ids  <- vapply(result, function(d) attr(d, "batch_id"), integer(1L))
  seeds <- vapply(result, function(d) attr(d, "seed"), integer(1L))
  expect_equal(ids, 1:5)
  expect_true(all(!is.na(seeds)))
})

test_that("n_batch is reproducible with same seed", {
  r1 <- simulate_data("ttest", seed = 77, n_batch = 3L)
  r2 <- simulate_data("ttest", seed = 77, n_batch = 3L)
  expect_identical(r1, r2)
})

test_that("n_batch differs with different seeds", {
  r1 <- simulate_data("ttest", seed = 10, n_batch = 3L)
  r2 <- simulate_data("ttest", seed = 20, n_batch = 3L)
  expect_false(identical(r1[[1]], r2[[1]]))
})

test_that("n_batch datasets vary across batch items", {
  result <- simulate_data("ttest", seed = 1, n_batch = 5L)
  # Not all datasets should be identical
  identical_pairs <- vapply(2:5, function(i) identical(result[[1]], result[[i]]), logical(1L))
  expect_false(all(identical_pairs))
})

# ---- type = "batch" full batch ----

test_that("type batch returns named list of all types", {
  result <- simulate_data("batch", seed = 1, n_batch = 3L)
  expect_type(result, "list")
  expect_equal(length(result), 7L)
  expected_names <- c("ttest", "anova", "correlation", "clusters",
                      "factor_analysis", "prediction", "mlvar")
  expect_equal(sort(names(result)), sort(expected_names))
})

test_that("type batch has sim_batch_full class", {
  result <- simulate_data("batch", seed = 1, n_batch = 3L)
  expect_s3_class(result, "sim_batch_full")
})

test_that("type batch n_batch attribute is correct", {
  result <- simulate_data("batch", seed = 1, n_batch = 4L)
  expect_equal(attr(result, "n_batch"), 4L)
})

test_that("type batch each type element has correct length", {
  result <- simulate_data("batch", seed = 1, n_batch = 5L)
  lengths <- vapply(result, length, integer(1L))
  expect_true(all(lengths == 5L))
})

test_that("type batch ttest element has correct structure", {
  result <- simulate_data("batch", seed = 1, n_batch = 3L)
  ttest_batch <- result[["ttest"]]
  expect_true(all(vapply(ttest_batch, function(d) {
    all(c("group", "score") %in% names(d))
  }, logical(1L))))
})

test_that("type batch is reproducible with same seed", {
  r1 <- simulate_data("batch", seed = 5, n_batch = 3L)
  r2 <- simulate_data("batch", seed = 5, n_batch = 3L)
  expect_identical(r1, r2)
})

# ---- seed consistency between batch types ----

test_that("ttest item 1 from n_batch matches item 1 from full batch", {
  single_batch <- simulate_data("ttest", seed = 42, n_batch = 3L)
  full_batch   <- simulate_data("batch", seed = 42, n_batch = 3L)
  # Strip batch_id attribute difference — seeds and data should be identical
  expect_identical(
    single_batch[[1]][, c("group", "score")],
    full_batch$ttest[[1]][, c("group", "score")]
  )
})

# ---- print methods ----

test_that("print.sim_batch_type works without error", {
  result <- simulate_data("ttest", seed = 1, n_batch = 3L)
  expect_output(print(result), regexp = "ttest")
})

test_that("print.sim_batch_full works without error", {
  result <- simulate_data("batch", seed = 1, n_batch = 3L)
  expect_output(print(result), regexp = "batch")
})
