test_that("simulate_network creates valid network object", {
  skip_if_not_installed("network")

  net <- simulate_network(n = 10, seed = 42)

  expect_s3_class(net, "network")
  expect_equal(network::network.size(net), 10)
})

test_that("simulate_network has vertex names", {
  skip_if_not_installed("network")

  net <- simulate_network(n = 10, seed = 42)
  names <- network::network.vertex.names(net)

  expect_equal(length(names), 10)
  expect_true(all(nchar(names) > 0))
})

test_that("simulate_network respects model parameter", {
  skip_if_not_installed("network")

  net_er <- simulate_network(n = 15, model = "er", seed = 42)
  net_ba <- simulate_network(n = 15, model = "ba", seed = 42)

  expect_s3_class(net_er, "network")
  expect_s3_class(net_ba, "network")
})

test_that("simulate_network can use learning states", {
  skip_if_not_installed("network")

  net <- simulate_network(n = 8, name_source = "states", seed = 42)
  names <- network::network.vertex.names(net)

  expect_equal(length(names), 8)
})

test_that("simulate_network is reproducible with seed", {
  skip_if_not_installed("network")

  net1 <- simulate_network(n = 10, seed = 42)
  net2 <- simulate_network(n = 10, seed = 42)

  expect_equal(
    network::network.vertex.names(net1),
    network::network.vertex.names(net2)
  )
})
