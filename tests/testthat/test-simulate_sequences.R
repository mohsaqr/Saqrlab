test_that("simulate_sequences creates valid output", {
  seq <- simulate_sequences(
    n_sequences = 10,
    seq_length = 15,
    n_states = 5,
    seed = 42
  )

  expect_true(is.data.frame(seq))
  expect_equal(nrow(seq), 10)
  expect_equal(ncol(seq), 15)
})

test_that("simulate_sequences uses correct number of states", {
  seq <- simulate_sequences(
    n_sequences = 20,
    seq_length = 20,
    n_states = 4,
    seed = 42
  )

  unique_states <- unique(unlist(seq))
  expect_lte(length(unique_states), 4)
})

test_that("simulate_sequences respects custom transition matrix", {
  trans_mat <- matrix(c(
    0.8, 0.1, 0.1,
    0.1, 0.8, 0.1,
    0.1, 0.1, 0.8
  ), nrow = 3, byrow = TRUE)
  rownames(trans_mat) <- colnames(trans_mat) <- c("A", "B", "C")
  init_probs <- c(A = 0.4, B = 0.3, C = 0.3)

  seq <- simulate_sequences(
    trans_matrix = trans_mat,
    init_probs = init_probs,
    n_sequences = 50,
    seq_length = 20,
    seed = 42
  )

  unique_states <- unique(unlist(seq))
  expect_true(all(unique_states %in% c("A", "B", "C")))
})

test_that("simulate_sequences with include_params returns list", {
  result <- simulate_sequences(
    n_sequences = 10,
    seq_length = 15,
    n_states = 4,
    include_params = TRUE,
    seed = 42
  )

  expect_type(result, "list")
  expect_true("sequences" %in% names(result))
  expect_true("transition_matrix" %in% names(result))
  expect_true("initial_probabilities" %in% names(result))
})

test_that("simulate_sequences can add NAs", {
  seq <- simulate_sequences(
    n_sequences = 20,
    seq_length = 15,
    n_states = 4,
    include_na = TRUE,
    na_range = c(1, 5),
    seed = 42
  )

  # At least some sequences should have NAs
  has_na <- apply(seq, 1, function(x) any(is.na(x)))
  expect_true(any(has_na))
})

test_that("simulate_sequences is reproducible with seed", {
  seq1 <- simulate_sequences(n_sequences = 10, seq_length = 15, n_states = 5, seed = 42)
  seq2 <- simulate_sequences(n_sequences = 10, seq_length = 15, n_states = 5, seed = 42)

  expect_equal(seq1, seq2)
})

test_that("simulate_sequences respects categories parameter", {
  seq <- simulate_sequences(
    n_sequences = 10,
    seq_length = 15,
    n_states = 5,
    categories = "metacognitive",
    seed = 42
  )

  unique_states <- unique(unlist(seq))
  metacog_states <- LEARNING_STATES$metacognitive

  expect_true(all(unique_states %in% metacog_states))
})
