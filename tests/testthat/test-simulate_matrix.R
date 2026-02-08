test_that("simulate_matrix creates valid transition matrix", {
  mat <- simulate_matrix(n_nodes = 5, seed = 42)

  expect_true(is.matrix(mat))
  expect_equal(dim(mat), c(5, 5))
  expect_equal(rownames(mat), colnames(mat))

  # Rows should sum to 1 for transition matrix
  expect_equal(unname(rowSums(mat)), rep(1, 5), tolerance = 1e-10)
})

test_that("simulate_matrix respects matrix_type", {
  # Transition matrix
  trans <- simulate_matrix(n_nodes = 4, matrix_type = "transition", seed = 42)
  expect_equal(unname(rowSums(trans)), rep(1, 4), tolerance = 1e-10)

  # Frequency matrix (integers)
  freq <- simulate_matrix(n_nodes = 4, matrix_type = "frequency", seed = 42)
  expect_true(all(freq >= 0))

  # Co-occurrence matrix (symmetric)
  cooc <- simulate_matrix(n_nodes = 4, matrix_type = "co-occurrence", seed = 42)
  expect_true(isSymmetric(cooc))

  # Adjacency matrix
  adj <- simulate_matrix(n_nodes = 4, matrix_type = "adjacency", weighted = FALSE, seed = 42)
  expect_true(all(adj %in% c(0, 1)))
})

test_that("simulate_matrix uses learning state names", {
  mat <- simulate_matrix(n_nodes = 6, seed = 42)
  names <- rownames(mat)

  expect_equal(length(names), 6)
  expect_true(all(nchar(names) > 0))
})

test_that("simulate_matrix respects custom names", {
  custom <- c("A", "B", "C", "D")
  mat <- simulate_matrix(n_nodes = 4, names = custom, seed = 42)

  expect_equal(rownames(mat), custom)
  expect_equal(colnames(mat), custom)
})

test_that("simulate_matrix is reproducible with seed", {
  mat1 <- simulate_matrix(n_nodes = 5, seed = 42)
  mat2 <- simulate_matrix(n_nodes = 5, seed = 42)

  expect_equal(mat1, mat2)
})

test_that("simulate_matrix respects edge_prob", {
  sparse <- simulate_matrix(n_nodes = 10, edge_prob = 0.1, seed = 42)
  dense <- simulate_matrix(n_nodes = 10, edge_prob = 0.9, seed = 42)

  # Dense should have more non-zero entries
  expect_gt(sum(dense > 0), sum(sparse > 0))
})
