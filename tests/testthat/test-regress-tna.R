# Characterization / regression tests for TNA & sequence simulators.
# Assertions reflect ACTUAL observed behavior (inspected via str/class/names),
# not assumptions. All calls use small sizes + explicit seeds for speed and
# determinism. See "BUGS FOUND" notes inline for functions that error on basic use.

# ---------------------------------------------------------------------------
# generate_probabilities
# ---------------------------------------------------------------------------
test_that("generate_probabilities: smoke + structure", {
  p <- generate_probabilities(n_states = 4, seed = 1)
  expect_type(p, "list")
  expect_named(p, c("initial_probs", "transition_probs", "state_names"))

  expect_length(p$initial_probs, 4)
  expect_equal(unname(names(p$initial_probs)), LETTERS[1:4])
  expect_equal(sum(p$initial_probs), 1, tolerance = 1e-8)

  expect_true(is.matrix(p$transition_probs))
  expect_equal(dim(p$transition_probs), c(4, 4))
  expect_equal(unname(rowSums(p$transition_probs)), rep(1, 4), tolerance = 1e-8)
  expect_equal(dimnames(p$transition_probs), list(LETTERS[1:4], LETTERS[1:4]))

  expect_equal(p$state_names, LETTERS[1:4])
})

test_that("generate_probabilities: reproducible with same seed", {
  expect_identical(
    generate_probabilities(n_states = 4, seed = 1),
    generate_probabilities(n_states = 4, seed = 1)
  )
})

test_that("generate_probabilities: custom state names honored", {
  p <- generate_probabilities(n_states = 3, states = c("x", "y", "z"), seed = 2)
  expect_equal(p$state_names, c("x", "y", "z"))
  expect_equal(unname(names(p$initial_probs)), c("x", "y", "z"))
})

# ---------------------------------------------------------------------------
# simulate_tna_matrix / generate_tna_matrix (alias)
# ---------------------------------------------------------------------------
test_that("simulate_tna_matrix: smoke + structure", {
  m <- simulate_tna_matrix(nodes_per_group = 2, n_groups = 3, seed = 1, verbose = FALSE)
  expect_type(m, "list")
  expect_named(m, c("matrix", "node_types"))

  expect_true(is.matrix(m$matrix))
  # nodes_per_group * n_groups total nodes
  expect_equal(dim(m$matrix), c(6, 6))

  expect_type(m$node_types, "list")
  expect_length(m$node_types, 3)
  # each group lists its nodes_per_group node labels
  expect_true(all(vapply(m$node_types, length, integer(1)) == 2))
})

test_that("simulate_tna_matrix: reproducible with same seed", {
  expect_identical(
    simulate_tna_matrix(nodes_per_group = 2, n_groups = 3, seed = 1, verbose = FALSE),
    simulate_tna_matrix(nodes_per_group = 2, n_groups = 3, seed = 1, verbose = FALSE)
  )
})

test_that("generate_tna_matrix: alias matches simulate_tna_matrix", {
  expect_identical(
    generate_tna_matrix(nodes_per_group = 2, n_groups = 3, seed = 1, verbose = FALSE),
    simulate_tna_matrix(nodes_per_group = 2, n_groups = 3, seed = 1, verbose = FALSE)
  )
})

# ---------------------------------------------------------------------------
# simulate_tna_network
# ---------------------------------------------------------------------------
test_that("simulate_tna_network: smoke + structure", {
  skip_if_not_installed("seqHMM")
  skip_if_not_installed("tna")
  n <- simulate_tna_network(n_states = 4, n_sequences = 30, seq_length = 10, seed = 1)
  expect_s3_class(n, "tna")
  expect_named(n, c("weights", "inits", "labels", "data"))
  expect_true(is.matrix(n$weights))
  expect_equal(dim(n$weights), c(4, 4))
  expect_length(n$labels, 4)
})

test_that("simulate_tna_network: reproducible with same seed", {
  skip_if_not_installed("seqHMM")
  skip_if_not_installed("tna")
  a <- simulate_tna_network(n_states = 4, n_sequences = 30, seq_length = 10, seed = 1)
  b <- simulate_tna_network(n_states = 4, n_sequences = 30, seq_length = 10, seed = 1)
  expect_identical(a$weights, b$weights)
  expect_identical(a$labels, b$labels)
})

# ---------------------------------------------------------------------------
# simulate_tna_networks / generate_tna_networks (alias)
# ---------------------------------------------------------------------------
test_that("simulate_tna_networks: smoke + structure", {
  skip_if_not_installed("seqHMM")
  skip_if_not_installed("tna")
  nn <- simulate_tna_networks(
    n_networks = 2, n_states = 4, n_sequences = 30, seq_length = 10,
    seed = 1, verbose = FALSE
  )
  expect_type(nn, "list")
  expect_length(nn, 2)
  expect_named(nn, c("network_1", "network_2"))
  expect_s3_class(nn$network_1, "tna")
  expect_s3_class(nn$network_2, "tna")
  expect_equal(dim(nn$network_1$weights), c(4, 4))
})

test_that("simulate_tna_networks: reproducible with same seed", {
  skip_if_not_installed("seqHMM")
  skip_if_not_installed("tna")
  a <- simulate_tna_networks(n_networks = 2, n_states = 4, n_sequences = 30,
                             seq_length = 10, seed = 1, verbose = FALSE)
  b <- simulate_tna_networks(n_networks = 2, n_states = 4, n_sequences = 30,
                             seq_length = 10, seed = 1, verbose = FALSE)
  expect_identical(a$network_1$weights, b$network_1$weights)
  expect_identical(a$network_2$weights, b$network_2$weights)
})

test_that("generate_tna_networks: alias matches simulate_tna_networks", {
  skip_if_not_installed("seqHMM")
  skip_if_not_installed("tna")
  a <- generate_tna_networks(n_networks = 1, n_states = 4, n_sequences = 30,
                             seq_length = 10, seed = 1, verbose = FALSE)
  b <- simulate_tna_networks(n_networks = 1, n_states = 4, n_sequences = 30,
                             seq_length = 10, seed = 1, verbose = FALSE)
  expect_identical(a$network_1$weights, b$network_1$weights)
})

# ---------------------------------------------------------------------------
# simulate_tna_datasets / generate_tna_datasets + generate_sequence_data (aliases)
# ---------------------------------------------------------------------------
test_that("simulate_tna_datasets: smoke + structure", {
  skip_if_not_installed("seqHMM")
  skip_if_not_installed("tna")
  d <- simulate_tna_datasets(
    n_datasets = 2, n_states = 4, n_sequences = 20, seq_length = 8,
    seed = 1, verbose = FALSE
  )
  expect_type(d, "list")
  expect_length(d, 2)
  expect_named(d, c("dataset_1", "dataset_2"))

  one <- d[[1]]
  expect_named(one, c("model", "params", "transition_probs", "initial_probs", "sequences"))
  expect_s3_class(one$model, "tna")
  expect_true(is.matrix(one$transition_probs))
  expect_equal(dim(one$transition_probs), c(4, 4))
  expect_length(one$initial_probs, 4)
  expect_s3_class(one$sequences, "data.frame")
  expect_equal(dim(one$sequences), c(20, 8))
})

test_that("simulate_tna_datasets: reproducible with same seed", {
  skip_if_not_installed("seqHMM")
  skip_if_not_installed("tna")
  a <- simulate_tna_datasets(n_datasets = 2, n_states = 4, n_sequences = 20,
                             seq_length = 8, seed = 1, verbose = FALSE)
  b <- simulate_tna_datasets(n_datasets = 2, n_states = 4, n_sequences = 20,
                             seq_length = 8, seed = 1, verbose = FALSE)
  expect_identical(a$dataset_1$sequences, b$dataset_1$sequences)
  expect_identical(a$dataset_1$transition_probs, b$dataset_1$transition_probs)
})

test_that("generate_sequence_data: alias of simulate_tna_datasets", {
  skip_if_not_installed("seqHMM")
  skip_if_not_installed("tna")
  a <- generate_sequence_data(n_datasets = 1, n_states = 4, n_sequences = 20,
                              seq_length = 8, seed = 1, verbose = FALSE)
  b <- simulate_tna_datasets(n_datasets = 1, n_states = 4, n_sequences = 20,
                             seq_length = 8, seed = 1, verbose = FALSE)
  expect_named(a, "dataset_1")
  expect_identical(a$dataset_1$sequences, b$dataset_1$sequences)
})

# ---------------------------------------------------------------------------
# simulate_sequences_advanced
# ---------------------------------------------------------------------------
test_that("simulate_sequences_advanced: smoke + structure", {
  s <- simulate_sequences_advanced(
    n_sequences = 10, seq_length = 8, n_states = 4, seed = 1
  )
  expect_s3_class(s, "data.frame")
  expect_equal(dim(s), c(10, 8))
  expect_equal(colnames(s), paste0("V", 1:8))
  # state values are stored as character labels
  expect_type(s[[1]], "character")
})

test_that("simulate_sequences_advanced: reproducible with same seed", {
  a <- simulate_sequences_advanced(n_sequences = 10, seq_length = 8, n_states = 4, seed = 1)
  b <- simulate_sequences_advanced(n_sequences = 10, seq_length = 8, n_states = 4, seed = 1)
  expect_identical(a, b)
})

# ---------------------------------------------------------------------------
# sample_tna
# ---------------------------------------------------------------------------
test_that("sample_tna: smoke + structure (refits a tna model)", {
  skip_if_not_installed("seqHMM")
  skip_if_not_installed("tna")
  m <- simulate_tna_network(n_states = 4, n_sequences = 40, seq_length = 10, seed = 2)
  st <- sample_tna(m, sampling_percent = 0.5)
  expect_s3_class(st, "tna")
  expect_named(st, c("weights", "inits", "labels", "data"))
})

test_that("sample_tna: reproducible under fixed RNG state", {
  skip_if_not_installed("seqHMM")
  skip_if_not_installed("tna")
  m <- simulate_tna_network(n_states = 4, n_sequences = 40, seq_length = 10, seed = 2)
  set.seed(7); a <- sample_tna(m, 0.5)
  set.seed(7); b <- sample_tna(m, 0.5)
  expect_identical(a$weights, b$weights)
})

test_that("sample_tna: validates inputs", {
  skip_if_not_installed("seqHMM")
  skip_if_not_installed("tna")
  m <- simulate_tna_network(n_states = 4, n_sequences = 40, seq_length = 10, seed = 2)
  expect_error(sample_tna(m, sampling_percent = 0), "between 0 and 1")
  expect_error(sample_tna(m, sampling_percent = 1.5), "between 0 and 1")
  expect_error(sample_tna(list(1), 0.5), "TNA object or a data frame")
})

# ---------------------------------------------------------------------------
# simulate_group_tna_networks / generate_group_tna_networks (alias)
#
# Previously broken (called an undefined helper long_to_wide()), now FIXED.
# Characterization tests pin the now-working behavior: a "group_tna" object
# (a list of fitted "tna" models, one per group). Reproducibility verified:
# same seed -> identical, different seed -> different.
# ---------------------------------------------------------------------------
test_that("simulate_group_tna_networks: runs and returns a group_tna object", {
  skip_if_not_installed("tna")
  m <- simulate_group_tna_networks(n_groups = 2, n_actors = 5, n_states = 3,
                                   seq_length_range = c(5, 8), seed = 1,
                                   verbose = FALSE)
  expect_s3_class(m, "group_tna")
  expect_length(m, 2L)
  # Each group element is a fitted tna model with the expected components.
  expect_s3_class(m[[1]], "tna")
  expect_s3_class(m[[2]], "tna")
  expect_named(m[[1]], c("weights", "inits", "labels", "data"))
  expect_true(is.matrix(m[[1]]$weights))
})

test_that("simulate_group_tna_networks: is reproducible by seed", {
  skip_if_not_installed("tna")
  a <- simulate_group_tna_networks(n_groups = 2, n_actors = 5, n_states = 3,
                                   seq_length_range = c(5, 8), seed = 1,
                                   verbose = FALSE)
  b <- simulate_group_tna_networks(n_groups = 2, n_actors = 5, n_states = 3,
                                   seq_length_range = c(5, 8), seed = 1,
                                   verbose = FALSE)
  d <- simulate_group_tna_networks(n_groups = 2, n_actors = 5, n_states = 3,
                                   seq_length_range = c(5, 8), seed = 99,
                                   verbose = FALSE)
  expect_identical(a, b)        # same seed -> identical
  expect_false(identical(a, d)) # different seed -> different
})

test_that("generate_group_tna_networks: alias matches simulate_group_tna_networks", {
  skip_if_not_installed("tna")
  m2 <- generate_group_tna_networks(n_groups = 2, n_actors = 5, n_states = 3,
                                    seq_length_range = c(5, 8), seed = 2,
                                    verbose = FALSE)
  expect_s3_class(m2, "group_tna")
  expect_length(m2, 2L)
  expect_s3_class(m2[[1]], "tna")
  # Alias is the same function: same args + seed -> identical result.
  via_main <- simulate_group_tna_networks(n_groups = 2, n_actors = 5,
                                          n_states = 3, seq_length_range = c(5, 8),
                                          seed = 2, verbose = FALSE)
  expect_identical(m2, via_main)
})
