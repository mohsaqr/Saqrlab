# ---- bootstrap_network() Tests ----

# Helper: generate wide sequence data
.make_boot_wide <- function(n = 50, t = 10, states = c("A", "B", "C"),
                            seed = 42) {
  set.seed(seed)
  mat <- matrix(sample(states, n * t, replace = TRUE), nrow = n, ncol = t)
  colnames(mat) <- paste0("T", seq_len(t))
  as.data.frame(mat, stringsAsFactors = FALSE)
}

# Helper: generate frequency-like data for association methods
.make_boot_assoc <- function(n = 100, p = 5, seed = 42) {
  set.seed(seed)
  mat <- matrix(rpois(n * p, lambda = 10), nrow = n, ncol = p)
  colnames(mat) <- paste0("state_", seq_len(p))
  as.data.frame(mat)
}


# ---- Basic functionality: transition methods ----

test_that("bootstrap_network works with method='relative'", {
  wide <- .make_boot_wide()
  boot <- bootstrap_network(wide, method = "relative", iter = 30L, seed = 1)

  expect_s3_class(boot, "saqr_bootstrap")
  expect_s3_class(boot$original, "saqr_network")
  expect_s3_class(boot$model, "saqr_network")
  expect_equal(boot$method, "relative")
  expect_equal(boot$iter, 30L)
  expect_true(is.matrix(boot$mean))
  expect_true(is.matrix(boot$sd))
  expect_true(is.matrix(boot$p_values))
  expect_true(is.matrix(boot$ci_lower))
  expect_true(is.matrix(boot$ci_upper))
  expect_true(is.matrix(boot$significant))
  expect_equal(nrow(boot$mean), 3)
  expect_equal(ncol(boot$mean), 3)
})

test_that("bootstrap_network works with method='frequency'", {
  wide <- .make_boot_wide()
  boot <- bootstrap_network(wide, method = "frequency", iter = 30L, seed = 1)

  expect_s3_class(boot, "saqr_bootstrap")
  expect_equal(boot$method, "frequency")
  expect_true(all(boot$mean >= 0))
})

test_that("bootstrap_network works with method='co_occurrence'", {
  wide <- .make_boot_wide()
  boot <- bootstrap_network(wide, method = "co_occurrence", iter = 30L,
                            seed = 1)

  expect_s3_class(boot, "saqr_bootstrap")
  expect_equal(boot$method, "co_occurrence")
  expect_false(boot$original$directed)
})


# ---- Basic functionality: association methods ----

test_that("bootstrap_network works with method='cor'", {
  df <- .make_boot_assoc()
  boot <- bootstrap_network(df, method = "cor", iter = 20L, seed = 1)

  expect_s3_class(boot, "saqr_bootstrap")
  expect_equal(boot$method, "cor")
  expect_false(boot$original$directed)
  expect_true(is.matrix(boot$mean))
  expect_equal(nrow(boot$mean), 5)
})

test_that("bootstrap_network works with method='pcor'", {
  df <- .make_boot_assoc(n = 80, p = 4)
  boot <- bootstrap_network(df, method = "pcor", iter = 20L, seed = 1)

  expect_s3_class(boot, "saqr_bootstrap")
  expect_equal(boot$method, "pcor")
})

test_that("bootstrap_network works with method='glasso'", {
  df <- .make_boot_assoc(n = 80, p = 4)
  boot <- bootstrap_network(df, method = "glasso", iter = 20L, seed = 1,
                            params = list(nlambda = 20L))

  expect_s3_class(boot, "saqr_bootstrap")
  expect_equal(boot$method, "glasso")
})


# ---- Inference ----

test_that("stability inference produces valid p-values", {
  wide <- .make_boot_wide()
  boot <- bootstrap_network(wide, method = "relative", iter = 50L,
                            seed = 42, inference = "stability")

  expect_true(all(boot$p_values >= 0))
  expect_true(all(boot$p_values <= 1))
  expect_equal(boot$inference, "stability")
  # CR bounds should exist
  expect_true(is.matrix(boot$cr_lower))
  expect_true(is.matrix(boot$cr_upper))
})

test_that("threshold inference produces valid p-values", {
  wide <- .make_boot_wide()
  boot <- bootstrap_network(wide, method = "relative", iter = 50L,
                            seed = 42, inference = "threshold")

  expect_true(all(boot$p_values >= 0))
  expect_true(all(boot$p_values <= 1))
  expect_equal(boot$inference, "threshold")
  # edge_threshold should be auto-set

  expect_true(is.numeric(boot$edge_threshold))
})


# ---- CI correctness ----

test_that("ci_lower <= ci_upper everywhere", {
  wide <- .make_boot_wide()
  boot <- bootstrap_network(wide, method = "relative", iter = 50L, seed = 1)

  expect_true(all(boot$ci_lower <= boot$ci_upper + 1e-10))
})


# ---- Summary ----

test_that("summary returns correct columns", {
  wide <- .make_boot_wide()
  boot <- bootstrap_network(wide, method = "relative", iter = 30L, seed = 1)

  s <- summary(boot)
  expect_true(is.data.frame(s))
  expect_true(all(c("from", "to", "weight", "mean", "sd", "p_value",
                     "sig", "ci_lower", "ci_upper") %in% names(s)))

  # Stability inference should include CR columns
  expect_true(all(c("cr_lower", "cr_upper") %in% names(s)))

  # All edges should have non-zero original weight
  expect_true(all(s$weight != 0))
  # For directed: from != to
  expect_true(all(s$from != s$to))
})

test_that("summary for undirected keeps only upper triangle", {
  df <- .make_boot_assoc()
  boot <- bootstrap_network(df, method = "cor", iter = 20L, seed = 1)
  s <- summary(boot)

  if (nrow(s) > 0) {
    expect_true(all(s$from < s$to))
  }
})


# ---- Pruned model ----

test_that("pruned model has <= edges of original", {
  wide <- .make_boot_wide()
  boot <- bootstrap_network(wide, method = "relative", iter = 50L, seed = 1)

  expect_true(boot$model$n_edges <= boot$original$n_edges)
  # Pruned matrix should be zero where not significant
  non_sig <- boot$p_values >= boot$ci_level
  expect_true(all(boot$model$matrix[non_sig] == 0))
})


# ---- Composability ----

test_that("params stored and original matches standalone estimate", {
  wide <- .make_boot_wide()
  params <- list(format = "wide")
  boot <- bootstrap_network(wide, method = "relative", params = params,
                            iter = 20L, seed = 1)

  expect_identical(boot$params, params)
  expect_identical(boot$original$params, params)

  # Original should match standalone call
  standalone <- estimate_network(wide, method = "relative", params = params)
  expect_equal(boot$original$matrix, standalone$matrix)
})


# ---- Reproducibility ----

test_that("same seed produces identical results", {
  wide <- .make_boot_wide()
  boot1 <- bootstrap_network(wide, method = "relative", iter = 30L, seed = 99)
  boot2 <- bootstrap_network(wide, method = "relative", iter = 30L, seed = 99)

  expect_equal(boot1$p_values, boot2$p_values)
  expect_equal(boot1$mean, boot2$mean)
  expect_equal(boot1$ci_lower, boot2$ci_lower)
})


# ---- Print ----

test_that("print.saqr_bootstrap produces expected output", {
  wide <- .make_boot_wide()
  boot <- bootstrap_network(wide, method = "relative", iter = 30L, seed = 1)

  out <- capture.output(print(boot))
  expect_true(any(grepl("Bootstrap:", out)))
  expect_true(any(grepl("Iterations:", out)))
  expect_true(any(grepl("Nodes:", out)))
  expect_true(any(grepl("Significant edges:", out)))
})


# ---- Custom estimator ----

test_that("bootstrap works with custom estimator", {
  # Register a simple custom estimator
  custom_fn <- function(data, ...) {
    numeric_cols <- vapply(data, is.numeric, logical(1))
    mat <- as.matrix(data[, numeric_cols, drop = FALSE])
    S <- cor(mat)
    diag(S) <- 0
    list(matrix = S, nodes = colnames(S), directed = FALSE)
  }
  register_estimator("test_custom_boot", custom_fn,
                     "Test custom for bootstrap", directed = FALSE)
  on.exit(remove_estimator("test_custom_boot"), add = TRUE)

  df <- .make_boot_assoc(n = 50, p = 4)
  boot <- bootstrap_network(df, method = "test_custom_boot", iter = 20L,
                            seed = 1)

  expect_s3_class(boot, "saqr_bootstrap")
  expect_equal(boot$method, "test_custom_boot")
})


# ---- Method aliases ----

test_that("bootstrap resolves method aliases", {
  wide <- .make_boot_wide()
  boot <- bootstrap_network(wide, method = "transition", iter = 20L, seed = 1)
  expect_equal(boot$method, "relative")

  boot2 <- bootstrap_network(wide, method = "counts", iter = 20L, seed = 1)
  expect_equal(boot2$method, "frequency")
})


# ---- Validation ----

test_that("invalid inputs error", {
  wide <- .make_boot_wide()
  expect_error(bootstrap_network(wide, iter = 1L), "iter")
  expect_error(bootstrap_network(wide, ci_level = 0), "ci_level")
  expect_error(bootstrap_network(wide, ci_level = 1), "ci_level")
  expect_error(bootstrap_network(wide, method = 123), "is.character")
  expect_error(bootstrap_network(wide, inference = "bad"),
               "'arg' should be one of")
})


# ---- Scaling and threshold pass-through ----

test_that("scaling is applied to bootstrap replicates", {
  wide <- .make_boot_wide()
  boot <- bootstrap_network(wide, method = "relative", scaling = "max",
                            iter = 20L, seed = 1)
  # Original should have max scaling applied
  expect_true(max(abs(boot$original$matrix)) <= 1 + 1e-10)
  # Mean bootstrap values should also be bounded
  expect_true(max(abs(boot$mean)) <= 1 + 1e-10)
})
