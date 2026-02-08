test_that("simulate_igraph creates valid igraph object", {

  g <- simulate_igraph(n = 10, seed = 42)

  expect_s3_class(g, "igraph")
  expect_equal(igraph::vcount(g), 10)
  expect_true(all(nchar(igraph::V(g)$name) > 0))
})

test_that("simulate_igraph respects model parameter", {
  # Erdos-Renyi

g_er <- simulate_igraph(n = 20, model = "er", p = 0.3, seed = 42)
  expect_s3_class(g_er, "igraph")

  # Barabasi-Albert
  g_ba <- simulate_igraph(n = 20, model = "ba", seed = 42)
  expect_s3_class(g_ba, "igraph")

  # Watts-Strogatz
  g_ws <- simulate_igraph(n = 20, model = "ws", seed = 42)
  expect_s3_class(g_ws, "igraph")

  # Stochastic Block Model
  g_sbm <- simulate_igraph(n = 21, model = "sbm", blocks = 3, seed = 42)
  expect_s3_class(g_sbm, "igraph")
  expect_true("block" %in% igraph::vertex_attr_names(g_sbm))
})

test_that("simulate_igraph uses human names by default", {
  g <- simulate_igraph(n = 10, seed = 42)
  names <- igraph::V(g)$name

  # Names should not be V1, V2, etc.
  expect_false(any(grepl("^V\\d+$", names)))
})

test_that("simulate_igraph can use learning states", {
  g <- simulate_igraph(n = 8, name_source = "states", seed = 42)
  names <- igraph::V(g)$name

  # Should have valid state names
  expect_equal(length(names), 8)
  expect_true(all(nchar(names) > 0))
})

test_that("simulate_igraph respects regions parameter", {
  g <- simulate_igraph(n = 10, regions = "arab", seed = 42)
  names <- igraph::V(g)$name

  # All names should come from arab region
  arab_names <- GLOBAL_NAMES$arab
  expect_true(all(names %in% arab_names))
})

test_that("simulate_igraph adds weights when requested", {
  g <- simulate_igraph(n = 10, model = "er", p = 0.5, weighted = TRUE, seed = 42)

  if (igraph::ecount(g) > 0) {
    expect_true("weight" %in% igraph::edge_attr_names(g))
    weights <- igraph::E(g)$weight
    expect_true(all(weights >= 0.1 & weights <= 1.0))
  }
})

test_that("simulate_igraph with custom names works", {
  custom <- c("Alice", "Bob", "Carol", "Dave", "Eve")
  g <- simulate_igraph(n = 5, names = custom, seed = 42)

  expect_equal(igraph::V(g)$name, custom)
})

test_that("simulate_igraph with NULL n produces random size", {
  set.seed(42)
  g1 <- simulate_igraph(n = NULL, seed = 123)
  g2 <- simulate_igraph(n = NULL, seed = 456)

  # Both should be in 20-50 range
  expect_true(igraph::vcount(g1) >= 20 && igraph::vcount(g1) <= 50)
  expect_true(igraph::vcount(g2) >= 20 && igraph::vcount(g2) <= 50)
})

test_that("simulate_igraph is reproducible with seed", {
  g1 <- simulate_igraph(n = 15, seed = 42)
  g2 <- simulate_igraph(n = 15, seed = 42)

  expect_equal(igraph::V(g1)$name, igraph::V(g2)$name)
  expect_equal(igraph::ecount(g1), igraph::ecount(g2))
})
