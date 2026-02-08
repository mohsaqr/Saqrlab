test_that("simulate_edge_list creates valid data frame", {
  edges <- simulate_edge_list(n_nodes = 10, n_edges = 20, seed = 42)

  expect_true(is.data.frame(edges))
  expect_true("source" %in% names(edges))
  expect_true("target" %in% names(edges))
  expect_true("weight" %in% names(edges))
  expect_true("class" %in% names(edges))
})

test_that("simulate_edge_list respects n_edges", {
  edges <- simulate_edge_list(n_nodes = 20, n_edges = 30, seed = 42)

  expect_equal(nrow(edges), 30)
})

test_that("simulate_edge_list respects weight_range", {
  edges <- simulate_edge_list(
    n_nodes = 10,
    n_edges = 20,
    weight_range = c(5, 10),
    seed = 42
  )

  expect_true(all(edges$weight >= 5))
  expect_true(all(edges$weight <= 10))
})

test_that("simulate_edge_list respects n_classes", {
  edges <- simulate_edge_list(
    n_nodes = 20,
    n_edges = 50,
    n_classes = 5,
    seed = 42
  )

  expect_true(all(edges$class >= 1 & edges$class <= 5))
})

test_that("simulate_edge_list with custom names works", {
  custom <- paste0("Person", 1:10)
  edges <- simulate_edge_list(
    n_nodes = 10,
    n_edges = 20,
    names = custom,
    seed = 42
  )

  all_nodes <- unique(c(edges$source, edges$target))
  expect_true(all(all_nodes %in% custom))
})

test_that("simulate_edge_list directed vs undirected", {
  edges_dir <- simulate_edge_list(n_nodes = 10, n_edges = 30, directed = TRUE, seed = 42)
  edges_undir <- simulate_edge_list(n_nodes = 10, n_edges = 30, directed = FALSE, seed = 42)

  expect_true(is.data.frame(edges_dir))
  expect_true(is.data.frame(edges_undir))
})

test_that("simulate_edge_list is reproducible with seed", {
  edges1 <- simulate_edge_list(n_nodes = 10, n_edges = 20, seed = 42)
  edges2 <- simulate_edge_list(n_nodes = 10, n_edges = 20, seed = 42)

  expect_equal(edges1, edges2)
})

test_that("simulate_edge_list handles edge density", {
  edges <- simulate_edge_list(n_nodes = 20, edge_density = 5, seed = 42)

  # Should have approximately n_nodes * edge_density edges
  expect_equal(nrow(edges), 100)
})
