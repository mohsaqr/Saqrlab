# Tests for simulate_hmm() — hidden Markov model sequence generator

test_that("simulate_hmm returns a well-formed saqr_sim", {
  r <- simulate_hmm(n_sequences = 20, seq_length = 15, seed = 1)

  expect_s3_class(r, "saqr_sim")
  expect_equal(r$type, "hmm")

  d <- as.data.frame(r)
  expect_true(all(c("sequence_id", "time", "symbol") %in% names(d)))
  expect_equal(nrow(d), 20L * 15L)
  expect_equal(sort(unique(d$sequence_id)), 1:20)
  expect_equal(sort(unique(d$time)), 1:15)
})

test_that("symbols stay within 1..n_symbols", {
  r <- simulate_hmm(n_sequences = 30, seq_length = 20, n_states = 2,
                    n_symbols = 4, seed = 2)
  d <- as.data.frame(r)
  expect_true(all(d$symbol >= 1L))
  expect_true(all(d$symbol <= 4L))
})

test_that("params carry trans, emission, init and hidden paths", {
  r <- simulate_hmm(n_sequences = 10, seq_length = 12, n_states = 3,
                    n_symbols = 3, seed = 3)
  p <- r$params

  expect_equal(dim(p$trans), c(3L, 3L))
  expect_equal(dim(p$emission), c(3L, 3L))
  expect_length(p$init, 3L)
  expect_equal(p$n_states, 3L)
  expect_equal(p$n_symbols, 3L)

  # hidden_paths matrix: one row per sequence, latent states in 1..n_states.
  expect_equal(dim(p$hidden_paths), c(10L, 12L))
  expect_true(all(p$hidden_paths >= 1L & p$hidden_paths <= 3L))

  # Auto-generated matrices are row-stochastic.
  expect_equal(unname(rowSums(p$trans)), rep(1, 3), tolerance = 1e-8)
  expect_equal(unname(rowSums(p$emission)), rep(1, 3), tolerance = 1e-8)
})

test_that("empirical transition frequencies recover trans", {
  # Large data: count state->state transitions from the TRUE hidden paths
  # and compare row-normalised counts against the generating trans matrix.
  tr <- matrix(c(0.7, 0.2, 0.1,
                 0.15, 0.75, 0.1,
                 0.1, 0.2, 0.7), nrow = 3, byrow = TRUE)
  em <- matrix(c(0.7, 0.2, 0.1,
                 0.1, 0.7, 0.2,
                 0.2, 0.1, 0.7), nrow = 3, byrow = TRUE)
  r <- simulate_hmm(n_sequences = 600, seq_length = 200, n_states = 3,
                    n_symbols = 3, trans = tr, emission = em, seed = 101)

  paths <- r$params$hidden_paths
  from <- as.integer(paths[, -ncol(paths)])
  to   <- as.integer(paths[, -1L])

  counts <- table(factor(from, levels = 1:3), factor(to, levels = 1:3))
  emp_trans <- counts / rowSums(counts)

  expect_equal(unname(as.matrix(emp_trans)), tr, tolerance = 0.03)
})

test_that("empirical emission frequencies recover emission", {
  tr <- matrix(c(0.8, 0.2,
                 0.2, 0.8), nrow = 2, byrow = TRUE)
  em <- matrix(c(0.6, 0.3, 0.1,
                 0.1, 0.3, 0.6), nrow = 2, byrow = TRUE)
  r <- simulate_hmm(n_sequences = 600, seq_length = 200, n_states = 2,
                    n_symbols = 3, trans = tr, emission = em, seed = 202)

  d <- as.data.frame(r)
  states <- as.integer(t(r$params$hidden_paths))  # match long-format order
  # Long data is ordered by sequence_id then time; t() of the path matrix
  # flattens row-by-row (sequence-major) to the same order.
  symbols <- d$symbol

  counts <- table(factor(states, levels = 1:2),
                  factor(symbols, levels = 1:3))
  emp_em <- counts / rowSums(counts)

  expect_equal(unname(as.matrix(emp_em)), em, tolerance = 0.03)
})

test_that("empirical initial distribution recovers init", {
  init <- c(0.7, 0.3)
  r <- simulate_hmm(n_sequences = 5000, seq_length = 5, n_states = 2,
                    n_symbols = 2, init = init, seed = 303)
  first_states <- r$params$hidden_paths[, 1L]
  emp_init <- as.numeric(table(factor(first_states, levels = 1:2))) /
    length(first_states)
  expect_equal(emp_init, init, tolerance = 0.03)
})

test_that("results are reproducible with a fixed seed", {
  r1 <- simulate_hmm(n_sequences = 40, seq_length = 25, n_states = 3,
                     n_symbols = 4, seed = 555)
  r2 <- simulate_hmm(n_sequences = 40, seq_length = 25, n_states = 3,
                     n_symbols = 4, seed = 555)
  expect_equal(as.data.frame(r1), as.data.frame(r2))
  expect_equal(r1$params$hidden_paths, r2$params$hidden_paths)
})

test_that("explicit non-stochastic matrices are rejected", {
  bad_trans <- matrix(c(0.5, 0.2,
                        0.3, 0.3), nrow = 2, byrow = TRUE)
  expect_error(
    simulate_hmm(n_sequences = 5, seq_length = 5, n_states = 2,
                 n_symbols = 2, trans = bad_trans)
  )
})

test_that("seq_length of 1 produces single-step sequences", {
  r <- simulate_hmm(n_sequences = 50, seq_length = 1, n_states = 2,
                    n_symbols = 2, seed = 9)
  d <- as.data.frame(r)
  expect_equal(nrow(d), 50L)
  expect_equal(dim(r$params$hidden_paths), c(50L, 1L))
})
