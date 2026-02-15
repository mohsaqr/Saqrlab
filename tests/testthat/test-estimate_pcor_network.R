# Helper: generate reproducible frequency-like data
.make_freq_data <- function(n = 100, p = 5, seed = 42) {
  set.seed(seed)
  mat <- matrix(rpois(n * p, lambda = 10), nrow = n, ncol = p)
  colnames(mat) <- paste0("state_", seq_len(p))
  df <- as.data.frame(mat)
  df$rid <- seq_len(n)
  df
}


# ---- Input validation ----

test_that("pcor_network errors on non-numeric data only", {
  df <- data.frame(a = letters[1:10], b = letters[10:1])
  expect_error(
    pcor_network(df),
    "At least 2 numeric columns"
  )
})

test_that("pcor_network errors on non-symmetric matrix", {
  m <- matrix(1:9, 3, 3)
  expect_error(
    pcor_network(m, n = 50),
    "symmetric"
  )
})

test_that("pcor_network errors when n missing for matrix input", {
  m <- diag(5)
  expect_error(
    pcor_network(m),
    "Sample size 'n' is required"
  )
})


# ---- Auto-cleaning ----

test_that("zero-variance columns are dropped with message", {
  df <- .make_freq_data()
  df$constant <- 5
  expect_message(
    net <- pcor_network(df, nlambda = 20L),
    "Dropping zero-variance"
  )
  expect_equal(net$p, 5)
  expect_false("constant" %in% colnames(net$pcor_matrix))
})

test_that("non-syntactic column names are dropped with message", {
  df <- .make_freq_data(n = 80, p = 4)
  df$`%` <- rpois(80, 2)
  df$`*` <- rpois(80, 3)
  expect_message(
    net <- pcor_network(df, nlambda = 20L),
    "non-syntactic"
  )
  expect_equal(net$p, 4)
  expect_false("%" %in% colnames(net$pcor_matrix))
})

test_that("all-NA columns are dropped with message", {
  df <- .make_freq_data(n = 80, p = 5)
  df$empty <- NA_real_
  expect_message(
    net <- pcor_network(df, nlambda = 20L),
    "all-NA"
  )
  expect_equal(net$p, 5)
})

test_that("rows with NA are dropped with message", {
  df <- .make_freq_data(n = 80, p = 5)
  df$state_1[1:3] <- NA
  expect_message(
    net <- pcor_network(df, nlambda = 20L),
    "rows with NA"
  )
  expect_equal(net$n, 77)
})


# ---- Basic functionality ----

test_that("pcor_network works with data frame input", {
  df <- .make_freq_data(n = 80, p = 6)
  net <- pcor_network(df, nlambda = 20L)

  expect_s3_class(net, "pcor_network")
  expect_equal(net$n, 80)
  expect_equal(net$p, 6)
  expect_true(is.matrix(net$pcor_matrix))
  expect_equal(nrow(net$pcor_matrix), 6)
  expect_equal(ncol(net$pcor_matrix), 6)
  # Diagonal should be zero
  expect_true(all(diag(net$pcor_matrix) == 0))
  # Should be symmetric
  expect_equal(net$pcor_matrix, t(net$pcor_matrix))
  # Edges data frame
  expect_true(is.data.frame(net$edges))
  expect_true(all(c("from", "to", "weight") %in% names(net$edges)))
  expect_equal(net$n_edges, nrow(net$edges))
  # EBIC path length matches nlambda
  expect_equal(length(net$ebic_path), 20)
  expect_equal(length(net$lambda_path), 20)
})

test_that("pcor_network works with correlation matrix input", {
  df <- .make_freq_data(n = 100, p = 5)
  num_cols <- setdiff(names(df), "rid")
  S <- cor(df[, num_cols])

  net <- pcor_network(S, n = 100, nlambda = 20L)

  expect_s3_class(net, "pcor_network")
  expect_equal(net$n, 100)
  expect_equal(net$p, 5)
})

test_that("pcor_network works with covariance matrix input", {
  df <- .make_freq_data(n = 100, p = 5)
  num_cols <- setdiff(names(df), "rid")
  C <- cov(df[, num_cols])

  net <- pcor_network(C, n = 100, input_type = "cov", nlambda = 20L)

  expect_s3_class(net, "pcor_network")
  expect_equal(net$p, 5)
})


# ---- id_col exclusion ----

test_that("id_col columns are excluded from analysis", {
  df <- .make_freq_data(n = 80, p = 5)
  df$subject_id <- seq_len(80)

  net <- pcor_network(df, id_col = "subject_id", nlambda = 20L)

  # subject_id and rid should be excluded -> 5 variables
  expect_equal(net$p, 5)
  expect_false("subject_id" %in% colnames(net$pcor_matrix))
  expect_false("rid" %in% colnames(net$pcor_matrix))
})


# ---- Gamma effects ----

test_that("higher gamma produces sparser or equal networks", {
  df <- .make_freq_data(n = 150, p = 7, seed = 123)

  net_low <- pcor_network(df, gamma = 0, nlambda = 50L)
  net_high <- pcor_network(df, gamma = 1, nlambda = 50L)

  # Higher gamma should produce same or fewer edges
  expect_true(net_high$n_edges <= net_low$n_edges)
})


# ---- S3 print method ----

test_that("print.pcor_network produces expected output", {
  df <- .make_freq_data(n = 80, p = 5)
  net <- pcor_network(df, nlambda = 20L)

  out <- capture.output(print(net))
  expect_true(any(grepl("Partial Correlation Network", out)))
  expect_true(any(grepl("Variables: 5", out)))
  expect_true(any(grepl("Sample size: 80", out)))
  expect_true(any(grepl("Gamma:", out)))
  expect_true(any(grepl("Lambda:", out)))
})

test_that("print.pcor_network returns invisible(x)", {
  df <- .make_freq_data(n = 80, p = 5)
  net <- pcor_network(df, nlambda = 20L)
  ret <- capture.output(result <- print(net))
  expect_identical(result, net)
})


# ---- Correlation method argument ----

test_that("cor_method argument is respected", {
  df <- .make_freq_data(n = 80, p = 5)

  net_pearson <- pcor_network(df, cor_method = "pearson", nlambda = 20L)
  net_spearman <- pcor_network(df, cor_method = "spearman", nlambda = 20L)

  # Correlation matrices should differ
  expect_false(identical(net_pearson$cor_matrix, net_spearman$cor_matrix))
})


# ---- Edge data frame correctness ----

test_that("edges match non-zero upper triangle of pcor_matrix", {
  df <- .make_freq_data(n = 100, p = 6)
  net <- pcor_network(df, nlambda = 30L)

  pcor <- net$pcor_matrix
  upper_nz <- which(upper.tri(pcor) & pcor != 0, arr.ind = TRUE)
  expect_equal(nrow(net$edges), nrow(upper_nz))

  # Weights should match matrix values
  for (i in seq_len(nrow(net$edges))) {
    r <- match(net$edges$from[i], colnames(pcor))
    c <- match(net$edges$to[i], colnames(pcor))
    expect_equal(net$edges$weight[i], pcor[r, c])
  }
})
