test_that("simulate_htna creates valid output", {
  result <- simulate_htna(n_nodes = 4, n_types = 3, seed = 42)

  expect_type(result, "list")
  expect_true("matrix" %in% names(result))
  expect_true("node_types" %in% names(result))
  expect_true("n_nodes_per_type" %in% names(result))
})

test_that("simulate_htna matrix has correct dimensions", {
  result <- simulate_htna(n_nodes = 5, n_types = 2, seed = 42)

  expected_size <- 5 * 2  # n_nodes * n_types
  expect_equal(dim(result$matrix), c(expected_size, expected_size))
})

test_that("simulate_htna creates transition matrix", {
  result <- simulate_htna(n_nodes = 3, n_types = 2, seed = 42)

  # Rows should sum to approximately 1 (some may be 0 for isolated nodes)
  row_sums <- unname(rowSums(result$matrix))
  # Each row should sum to 0 or 1
  expect_true(all(row_sums >= 0 & row_sums <= 1.0001))
})

test_that("simulate_htna node_types has correct structure", {
  result <- simulate_htna(n_nodes = 4, n_types = 3, seed = 42)

  expect_equal(length(result$node_types), 3)
  expect_equal(unname(result$n_nodes_per_type), c(4, 4, 4))

  # Each type should have 4 nodes
  for (type_nodes in result$node_types) {
    expect_equal(length(type_nodes), 4)
  }
})

test_that("simulate_htna respects custom type_names", {
  result <- simulate_htna(
    n_nodes = 3,
    n_types = 2,
    type_names = c("Individual", "Group"),
    seed = 42
  )

  expect_equal(names(result$node_types), c("Individual", "Group"))
})

test_that("simulate_htna is reproducible with seed", {
  result1 <- simulate_htna(n_nodes = 4, n_types = 2, seed = 42)
  result2 <- simulate_htna(n_nodes = 4, n_types = 2, seed = 42)

  expect_equal(result1$matrix, result2$matrix)
  expect_equal(result1$node_types, result2$node_types)
})
