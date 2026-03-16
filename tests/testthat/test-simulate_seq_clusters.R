# Helper: build a valid 3-state transition matrix
.make_trans <- function(rows) {
  m <- matrix(unlist(rows), nrow = length(rows), byrow = TRUE)
  rownames(m) <- colnames(m) <- LETTERS[seq_len(nrow(m))]
  m
}

K3 <- list(
  cluster1 = .make_trans(list(c(0.7, 0.2, 0.1), c(0.1, 0.8, 0.1), c(0.2, 0.1, 0.7))),
  cluster2 = .make_trans(list(c(0.1, 0.1, 0.8), c(0.2, 0.6, 0.2), c(0.8, 0.1, 0.1))),
  cluster3 = .make_trans(list(c(0.5, 0.3, 0.2), c(0.3, 0.4, 0.3), c(0.2, 0.3, 0.5)))
)

# ===========================================================================
# Return structure
# ===========================================================================

test_that("simulate_seq_clusters: returns list with $data and $params", {
  r <- simulate_seq_clusters(trans_list = K3, n = 60, seed = 1)
  expect_type(r, "list")
  expect_true(all(c("data", "params") %in% names(r)))
})

test_that("simulate_seq_clusters: $data is data.frame with T1..Tseq_length + true_cluster", {
  r <- simulate_seq_clusters(trans_list = K3, n = 60, seq_length = 10, seed = 1)
  expect_s3_class(r$data, "data.frame")
  expect_true(all(paste0("T", 1:10) %in% names(r$data)))
  expect_true("true_cluster" %in% names(r$data))
  expect_equal(nrow(r$data), 60L)
})

test_that("simulate_seq_clusters: nrow equals n", {
  r <- simulate_seq_clusters(trans_list = K3, n = 150, seed = 1)
  expect_equal(nrow(r$data), 150L)
})

test_that("simulate_seq_clusters: true_cluster contains only valid labels 1..K", {
  r <- simulate_seq_clusters(trans_list = K3, n = 120, seed = 1)
  expect_true(all(r$data$true_cluster %in% 1:3))
})

test_that("simulate_seq_clusters: all K clusters represented", {
  r <- simulate_seq_clusters(trans_list = K3, n = 150, seed = 1)
  expect_true(all(1:3 %in% r$data$true_cluster))
})

test_that("simulate_seq_clusters: sequences contain only valid state names", {
  r  <- simulate_seq_clusters(trans_list = K3, n = 60, seq_length = 8, seed = 1)
  state_cols <- paste0("T", 1:8)
  vals <- na.omit(unlist(r$data[, state_cols]))
  expect_true(all(vals %in% c("A", "B", "C")))
})

# ===========================================================================
# $params
# ===========================================================================

test_that("simulate_seq_clusters: $params contains trans_list, props, init_probs", {
  r <- simulate_seq_clusters(trans_list = K3, n = 60, seed = 1)
  expect_true(all(c("trans_list", "props", "init_probs") %in% names(r$params)))
})

test_that("simulate_seq_clusters: $params$trans_list identical to input", {
  r <- simulate_seq_clusters(trans_list = K3, n = 60, seed = 1)
  expect_identical(r$params$trans_list, K3)
})

test_that("simulate_seq_clusters: props default to equal mixing", {
  r <- simulate_seq_clusters(trans_list = K3, n = 300, seed = 1)
  expect_equal(r$params$props, rep(1/3, 3))
})

test_that("simulate_seq_clusters: custom props normalised and stored", {
  r <- simulate_seq_clusters(trans_list = K3, props = c(1, 2, 1), n = 200, seed = 1)
  expect_equal(r$params$props, c(0.25, 0.5, 0.25))
})

# ===========================================================================
# Mixing proportions reflected in cluster counts at large n
# ===========================================================================

test_that("simulate_seq_clusters: cluster proportions approximate props at large n", {
  r   <- simulate_seq_clusters(trans_list = K3, props = c(0.5, 0.3, 0.2),
                                n = 3000, seed = 42)
  obs <- table(r$data$true_cluster) / 3000
  expect_true(abs(obs["1"] - 0.5) < 0.05)
  expect_true(abs(obs["2"] - 0.3) < 0.05)
  expect_true(abs(obs["3"] - 0.2) < 0.05)
})

# ===========================================================================
# Seed reproducibility
# ===========================================================================

test_that("simulate_seq_clusters: seed produces identical results", {
  r1 <- simulate_seq_clusters(trans_list = K3, n = 60, seed = 7)
  r2 <- simulate_seq_clusters(trans_list = K3, n = 60, seed = 7)
  expect_identical(r1, r2)
})

test_that("simulate_seq_clusters: different seeds produce different data", {
  r1 <- simulate_seq_clusters(trans_list = K3, n = 60, seed = 1)
  r2 <- simulate_seq_clusters(trans_list = K3, n = 60, seed = 2)
  expect_false(identical(r1$data$T1, r2$data$T1))
})

# ===========================================================================
# init_probs
# ===========================================================================

test_that("simulate_seq_clusters: shared init_probs vector accepted", {
  init <- c(A = 0.5, B = 0.3, C = 0.2)
  expect_no_error(
    simulate_seq_clusters(trans_list = K3, n = 60, init_probs = init, seed = 1)
  )
})

test_that("simulate_seq_clusters: per-cluster init_probs list accepted", {
  init <- list(c(A=0.7, B=0.2, C=0.1), c(A=0.1, B=0.2, C=0.7),
               c(A=0.3, B=0.4, C=0.3))
  expect_no_error(
    simulate_seq_clusters(trans_list = K3, n = 60, init_probs = init, seed = 1)
  )
})

test_that("simulate_seq_clusters: init_probs stored in params", {
  init <- c(A = 0.5, B = 0.3, C = 0.2)
  r    <- simulate_seq_clusters(trans_list = K3, n = 60, init_probs = init, seed = 1)
  expect_equal(r$params$init_probs[[1]], init)
})

# ===========================================================================
# Auto-generation (trans_list = NULL)
# ===========================================================================

test_that("simulate_seq_clusters: auto-generates when trans_list is NULL", {
  expect_no_error(
    simulate_seq_clusters(n = 60, n_clusters = 2, n_states = 5, seed = 1)
  )
})

test_that("simulate_seq_clusters: auto mode uses n_states states (default 10)", {
  r <- simulate_seq_clusters(n = 120, n_clusters = 3, seed = 1)
  state_cols <- grep("^T", names(r$data), value = TRUE)
  vals <- na.omit(unlist(r$data[, state_cols]))
  expect_equal(length(unique(vals)), 10L)
})

test_that("simulate_seq_clusters: auto mode params$trans_list has K matrices", {
  r <- simulate_seq_clusters(n = 60, n_clusters = 3, n_states = 5, seed = 1)
  expect_equal(length(r$params$trans_list), 3L)
  expect_equal(nrow(r$params$trans_list[[1]]), 5L)
})

test_that("simulate_seq_clusters: auto mode rows of each matrix sum to 1", {
  r <- simulate_seq_clusters(n = 60, n_clusters = 2, n_states = 4, seed = 1)
  row_sums <- vapply(r$params$trans_list, function(m) max(abs(rowSums(m) - 1)),
                     numeric(1))
  expect_true(all(row_sums < 1e-10))
})

# ===========================================================================
# Validation errors
# ===========================================================================

test_that("simulate_seq_clusters: errors on non-square matrix", {
  bad <- list(matrix(1:6, nrow = 2))
  expect_error(simulate_seq_clusters(trans_list = bad, n = 10),
               regexp = "square")
})

test_that("simulate_seq_clusters: errors when matrix rows don't sum to 1", {
  bad_mat <- matrix(c(0.5, 0.5, 0.5, 0.1, 0.9, 0.1, 0.3, 0.3, 0.3),
                    nrow = 3, byrow = TRUE)
  rownames(bad_mat) <- colnames(bad_mat) <- c("A", "B", "C")
  expect_error(
    simulate_seq_clusters(trans_list = list(bad_mat), n = 10),
    regexp = "sum to 1"
  )
})

test_that("simulate_seq_clusters: errors when matrices have different dimensions", {
  m2 <- matrix(c(0.5, 0.5, 0.5, 0.5), nrow = 2,
               dimnames = list(c("A","B"), c("A","B")))
  m3 <- K3[[1]]
  expect_error(
    simulate_seq_clusters(trans_list = list(m2, m3), n = 20),
    regexp = "same"
  )
})

test_that("simulate_seq_clusters: errors when trans_list matrices lack rownames", {
  bad <- matrix(c(0.7, 0.3, 0.4, 0.6), nrow = 2)
  expect_error(
    simulate_seq_clusters(trans_list = list(bad), n = 10),
    regexp = "rownames"
  )
})
