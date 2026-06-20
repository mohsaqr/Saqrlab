# Characterization / regression tests for network-comparison & model-fitting
# functions. Behaviour here is pinned to ACTUAL observed output (probed against
# the live package), not assumed semantics.
#
# NOTE: the internal helper `extract_transition_matrix()` is now implemented, so
# compare_networks(), compare_centralities(), compare_edge_recovery() and its
# alias calculate_edge_recovery() are exercised below with real characterization
# tests. These functions require models whose states match (identical dimnames /
# shared alphabet); otherwise they legitimately error with "no common states".

# ---- shared tiny fixtures -------------------------------------------------

make_seqs <- function(seed = 1) {
  simulate_sequences(n_sequences = 30, seq_length = 12, n_states = 4, seed = seed)
}

# ---- fit_network_model ----------------------------------------------------

test_that("fit_network_model: smoke + structure for every model type", {
  seqs <- make_seqs()
  for (mt in c("tna", "ftna", "ctna", "atna")) {
    m <- fit_network_model(seqs, mt)
    expect_s3_class(m, "tna")
    expect_true(all(c("weights", "inits", "labels") %in% names(m)))
    expect_true(is.matrix(m$weights))
    expect_equal(nrow(m$weights), ncol(m$weights))
  }
})

test_that("fit_network_model: determinism (sequences fixed -> identical weights)", {
  seqs <- make_seqs()
  a <- fit_network_model(seqs, "tna")
  b <- fit_network_model(seqs, "tna")
  expect_equal(a$weights, b$weights)
})

# ---- comparison fixtures (shared-state weight-matrix models) ---------------

# Plain list models with identical state dimnames exercise the list path.
make_weight_model <- function(seed, states = c("A", "B", "C")) {
  set.seed(seed)
  n <- length(states)
  list(weights = matrix(stats::runif(n * n), n, n,
                        dimnames = list(states, states)))
}

# ---- compare_networks -----------------------------------------------------

test_that("compare_networks: smoke + structure (shared-state list models)", {
  r <- compare_networks(make_weight_model(1), make_weight_model(2))
  expect_s3_class(r, "tna_comparison")
  expect_named(r, c("metrics", "edge_comparison", "summary"))
  expect_named(r$metrics, c("correlation", "rmse", "edge_diff"))
  expect_s3_class(r$edge_comparison, "data.frame")
  expect_true(all(c("from", "to", "weight1", "weight2", "diff", "abs_diff") %in%
                    names(r$edge_comparison)))
  expect_equal(nrow(r$edge_comparison), 9L)
  expect_type(r$summary, "character")
})

test_that("compare_networks: determinism (same seeds -> identical)", {
  a <- compare_networks(make_weight_model(1), make_weight_model(2))
  b <- compare_networks(make_weight_model(1), make_weight_model(2))
  expect_equal(a, b)
})

# ---- compare_centralities -------------------------------------------------

make_centrality_seqs <- function(seed) {
  set.seed(seed)
  as.data.frame(matrix(sample(c("A", "B", "C"), 20 * 8, replace = TRUE), 20, 8))
}

test_that("compare_centralities: smoke + structure (shared-alphabet tna models)", {
  skip_if_not_installed("tna")
  t1 <- tna::tna(make_centrality_seqs(1))
  t2 <- tna::tna(make_centrality_seqs(2))
  r <- compare_centralities(t1, t2)
  expect_s3_class(r, "centrality_comparison")
  expect_named(r, c("correlations", "centrality_comparison", "summary"))
  expect_type(r$correlations, "list")
  expect_s3_class(r$centrality_comparison, "data.frame")
  expect_type(r$summary, "character")
})

test_that("compare_centralities: determinism (same seeds -> identical)", {
  skip_if_not_installed("tna")
  a <- compare_centralities(tna::tna(make_centrality_seqs(1)),
                            tna::tna(make_centrality_seqs(2)))
  b <- compare_centralities(tna::tna(make_centrality_seqs(1)),
                            tna::tna(make_centrality_seqs(2)))
  expect_equal(a, b)
})

# ---- compare_edge_recovery ------------------------------------------------

test_that("compare_edge_recovery: smoke + structure (shared-state list models)", {
  r <- compare_edge_recovery(make_weight_model(1), make_weight_model(2))
  expect_s3_class(r, "edge_recovery")
  expect_named(r, c("true_positives", "false_positives", "false_negatives",
                    "true_negatives", "precision", "recall", "f1_score",
                    "accuracy", "jaccard"))
  expect_type(r$precision, "double")
  expect_type(r$recall, "double")
  expect_type(r$f1_score, "double")
})

test_that("compare_edge_recovery: determinism (same seeds -> identical)", {
  a <- compare_edge_recovery(make_weight_model(1), make_weight_model(2))
  b <- compare_edge_recovery(make_weight_model(1), make_weight_model(2))
  expect_equal(a, b)
})

# ---- calculate_edge_recovery (alias) --------------------------------------

test_that("calculate_edge_recovery: alias matches compare_edge_recovery", {
  alias  <- calculate_edge_recovery(make_weight_model(1), make_weight_model(2))
  canon  <- compare_edge_recovery(make_weight_model(1), make_weight_model(2))
  expect_s3_class(alias, "edge_recovery")
  expect_named(alias, names(canon))
  expect_equal(alias, canon)
})

# ---- compare_network_estimation ------------------------------------------

test_that("compare_network_estimation: smoke + structure", {
  seqs <- make_seqs()
  r <- compare_network_estimation(
    seqs, model_types = c("tna", "ftna"),
    sampling_percent = 0.3, iterations = 2, seed = 42, verbose = FALSE
  )
  expect_s3_class(r, "network_estimation")
  expect_named(r, c("individual", "aggregated", "ranking", "winner", "params"))
  expect_s3_class(r$individual, "data.frame")
  expect_true(all(c("category", "metric", "value", "iteration", "model_type") %in%
                    names(r$individual)))
  expect_true(all(c("tna", "ftna") %in% unique(r$individual$model_type)))
  expect_true(r$winner %in% c("tna", "ftna"))
  expect_setequal(r$ranking, c("tna", "ftna"))
})

test_that("compare_network_estimation: determinism (same seed -> identical)", {
  seqs <- make_seqs()
  a <- compare_network_estimation(seqs, iterations = 2, seed = 7, verbose = FALSE)
  b <- compare_network_estimation(seqs, iterations = 2, seed = 7, verbose = FALSE)
  expect_equal(a$individual, b$individual)
  expect_equal(a$winner, b$winner)
})

# ---- compare_tna_models (alias of compare_network_estimation) -------------

test_that("compare_tna_models: alias smoke + structure", {
  seqs <- make_seqs()
  r <- compare_tna_models(
    seqs, model_types = c("tna", "ftna"),
    iterations = 2, seed = 42, verbose = FALSE
  )
  expect_s3_class(r, "network_estimation")
  expect_named(r, c("individual", "aggregated", "ranking", "winner", "params"))
})

# ---- cross_validate_tna ---------------------------------------------------

test_that("cross_validate_tna: smoke + structure", {
  seqs <- make_seqs()
  r <- cross_validate_tna(
    seqs, model_types = c("relative", "frequency"),
    sampling_percent = 0.3, iterations = 2, seed = 42, verbose = FALSE
  )
  expect_type(r, "list")
  expect_named(r, c("relative", "frequency"))
  expect_true(all(c("model", "results") %in% names(r$relative)))
  expect_s3_class(r$relative$model, "tna")
})

test_that("cross_validate_tna: determinism (same seed -> identical results)", {
  seqs <- make_seqs()
  a <- cross_validate_tna(seqs, model_types = "relative", iterations = 2,
                          seed = 9, verbose = FALSE)
  b <- cross_validate_tna(seqs, model_types = "relative", iterations = 2,
                          seed = 9, verbose = FALSE)
  expect_equal(a$relative$results, b$relative$results)
})

# ---- compare_estimation (self-simulating) ---------------------------------

test_that("compare_estimation: smoke + structure", {
  r <- compare_estimation(
    models = c("tna", "ftna"), n_simulations = 2,
    n_sequences = 30, seq_length = 12, n_states = 4,
    scaling = "minmax", seed = 42, verbose = FALSE
  )
  expect_type(r, "list")
  expect_named(r, c("comparison", "summary", "raw_results", "ranking",
                    "winner", "params"))
  expect_s3_class(r$comparison, "data.frame")
  expect_true(all(c("category", "metric", "tna_mean", "tna_sd",
                    "ftna_mean", "ftna_sd") %in% names(r$comparison)))
  expect_true(r$winner %in% c("tna", "ftna"))
  expect_setequal(r$ranking, c("tna", "ftna"))
})

test_that("compare_estimation: determinism (same seed -> identical raw)", {
  args <- list(models = c("tna", "ftna"), n_simulations = 2,
               n_sequences = 30, seq_length = 12, n_states = 4,
               seed = 5, verbose = FALSE)
  a <- do.call(compare_estimation, args)
  b <- do.call(compare_estimation, args)
  expect_equal(a$raw_results, b$raw_results)
})

# ---- compare_reliability (self-simulating) --------------------------------

test_that("compare_reliability: smoke + structure", {
  r <- compare_reliability(
    n_simulations = 2, n_sequences = c(30, 40), seq_length = c(10, 12),
    n_states = 4, model_type = "tna", reliability_iter = 2,
    scaling = "none", seed = 42, verbose = FALSE
  )
  expect_s3_class(r, "tna_reliability_comparison")
  expect_named(r, c("raw_results", "summary", "params"))
  expect_s3_class(r$raw_results, "data.frame")
  expect_true(all(c("metric", "mean", "sd", "sim_id") %in% names(r$raw_results)))
})

test_that("compare_reliability: determinism (same seed -> identical raw)", {
  args <- list(n_simulations = 2, n_sequences = c(30, 40),
               seq_length = c(10, 12), n_states = 4, reliability_iter = 2,
               seed = 3, verbose = FALSE)
  a <- do.call(compare_reliability, args)
  b <- do.call(compare_reliability, args)
  expect_equal(a$raw_results, b$raw_results)
})
