# Characterization / regression tests for batch / grid / sampling orchestration
# functions. These pin the CURRENTLY-OBSERVED behavior (structure + key
# names/dims) so future refactors are caught. Everything is kept tiny and
# sequential (cores = 1, <= 2 runs, short sequences, 2 states) because these
# orchestrators are slow and may use future/parallel/progressr frameworks.
#
# Functions covered:
#   generate_param_grid, create_param_grid (alias)
#   batch_fit_models, batch_apply
#   run_bootstrap_iteration, run_bootstrap_simulation
#   run_grid_simulation, run_network_simulation, run_sampling_analysis

# ---------------------------------------------------------------------------
# Shared tiny fixtures
# ---------------------------------------------------------------------------

make_trans <- function() {
  tm <- matrix(c(0.7, 0.3, 0.4, 0.6), nrow = 2, byrow = TRUE)
  rownames(tm) <- colnames(tm) <- c("A", "B")
  tm
}
make_inits <- function() c(A = 0.5, B = 0.5)
make_stable <- function() list(c("A", "B"))
make_model_list <- function() list(weights = make_trans(), inits = make_inits())

make_seqs <- function(n = 30, len = 10, seed = 1) {
  simulate_sequences(
    trans_matrix = make_trans(),
    init_probs = make_inits(),
    seq_length = len,
    n_sequences = n,
    seed = seed
  )
}

# ===========================================================================
# generate_param_grid / create_param_grid
# ===========================================================================

test_that("generate_param_grid (random) smoke + structure", {
  ranges <- list(num_rows = c(50, 100), seq_length = c(10, 20),
                 stability = c(0.7, 1.0))
  g <- generate_param_grid(ranges, n = 4, method = "random")
  expect_s3_class(g, "data.frame")
  expect_equal(nrow(g), 4)
  expect_setequal(names(g), c("num_rows", "seq_length", "stability"))
  # integer ranges (both endpoints integer) get rounded to integers
  expect_true(all(g$num_rows == round(g$num_rows)))
  expect_true(all(g$seq_length == round(g$seq_length)))
  # continuous range stays inside [0.7, 1.0]
  expect_true(all(g$stability >= 0.7 & g$stability <= 1.0))
})

test_that("generate_param_grid default (no args) returns 10-row demo grid", {
  g <- generate_param_grid()
  expect_s3_class(g, "data.frame")
  expect_equal(nrow(g), 10)
  expect_setequal(names(g), c("n_sequences", "seq_length", "n_states"))
})

test_that("generate_param_grid grid method returns <= n rows with all params", {
  ranges <- list(num_rows = c(50, 100), seq_length = c(10, 20),
                 stability = c(0.7, 1.0))
  g <- generate_param_grid(ranges, n = 4, method = "grid")
  expect_s3_class(g, "data.frame")
  expect_lte(nrow(g), 4)
  expect_setequal(names(g), c("num_rows", "seq_length", "stability"))
})

test_that("generate_param_grid lhs method returns n rows", {
  skip_if_not_installed("lhs")
  ranges <- list(num_rows = c(50, 100), seq_length = c(10, 20),
                 stability = c(0.7, 1.0))
  g <- generate_param_grid(ranges, n = 4, method = "lhs")
  expect_s3_class(g, "data.frame")
  expect_equal(nrow(g), 4)
  expect_setequal(names(g), c("num_rows", "seq_length", "stability"))
})

test_that("generate_param_grid (random) is deterministic under set.seed", {
  ranges <- list(num_rows = c(50, 100), seq_length = c(10, 20))
  set.seed(1); a <- generate_param_grid(ranges, n = 4, method = "random")
  set.seed(1); b <- generate_param_grid(ranges, n = 4, method = "random")
  expect_identical(a, b)
})

test_that("create_param_grid is an alias for generate_param_grid", {
  ranges <- list(num_rows = c(50, 100), seq_length = c(10, 20))
  set.seed(2); a <- generate_param_grid(ranges, n = 4, method = "random")
  set.seed(2); b <- create_param_grid(ranges, n = 4, method = "random")
  expect_identical(a, b)
})

# ===========================================================================
# batch_fit_models
# ===========================================================================

test_that("batch_fit_models smoke + structure (sequential)", {
  skip_if_not_installed("tna")
  dl <- list(make_seqs(seed = 1), make_seqs(seed = 2))
  m <- suppressWarnings(batch_fit_models(dl, model_type = "tna", progress = FALSE))
  expect_type(m, "list")
  expect_length(m, 2)
  expect_s3_class(m[[1]], "tna")
  expect_s3_class(m[[2]], "tna")
})

test_that("batch_fit_models preserves names of a named list", {
  skip_if_not_installed("tna")
  dl <- list(first = make_seqs(seed = 1), second = make_seqs(seed = 2))
  m <- suppressWarnings(batch_fit_models(dl, model_type = "tna", progress = FALSE))
  expect_equal(names(m), c("first", "second"))
})

test_that("batch_fit_models returns empty list for empty input", {
  expect_identical(batch_fit_models(list(), progress = FALSE), list())
})

test_that("batch_fit_models errors when data_list is not a list", {
  expect_error(batch_fit_models(42), "data_list must be a list")
})

# ===========================================================================
# batch_apply
# ===========================================================================

test_that("batch_apply applies a function over a list and returns a list", {
  skip_if_not_installed("tna")
  dl <- list(make_seqs(seed = 1), make_seqs(seed = 2))
  m <- suppressWarnings(batch_fit_models(dl, model_type = "tna", progress = FALSE))
  r <- suppressWarnings(batch_apply(m, function(x) x$weights, progress = FALSE))
  expect_type(r, "list")
  expect_length(r, 2)
})

test_that("batch_apply simplify collapses length-1 atomic results to a vector", {
  skip_if_not_installed("tna")
  dl <- list(make_seqs(seed = 1), make_seqs(seed = 2))
  m <- suppressWarnings(batch_fit_models(dl, model_type = "tna", progress = FALSE))
  r <- suppressWarnings(
    batch_apply(m, function(x) nrow(x$weights), progress = FALSE, simplify = TRUE)
  )
  expect_true(is.atomic(r))
  expect_length(r, 2)
  expect_equal(unname(r), c(2L, 2L))
})

test_that("batch_apply returns empty list for empty input", {
  expect_identical(batch_apply(list(), function(x) x, progress = FALSE), list())
})

test_that("batch_apply errors on bad arguments", {
  expect_error(batch_apply(42, identity), "object_list must be a list")
  expect_error(batch_apply(list(1), "not a function"), "fun must be a function")
})

# ===========================================================================
# run_bootstrap_iteration
# ===========================================================================

test_that("run_bootstrap_iteration smoke + structure", {
  skip_if_not_installed("tna")
  res <- suppressWarnings(run_bootstrap_iteration(
    trans_matrix = make_trans(), init_probs = make_inits(),
    stable_transitions = make_stable(),
    seq_length = 12, n_sequences = 25, level = 0.05, include_na = FALSE
  ))
  expect_type(res, "list")
  expect_setequal(
    names(res),
    c("per_edge", "bootstrap_summary_raw", "TP_matrix", "TN_matrix",
      "FP_matrix", "FN_matrix")
  )
  expect_s3_class(res$per_edge, "data.frame")
  # 2 states -> 4 directed edges (rows), 14 metric columns
  expect_equal(nrow(res$per_edge), 4)
  expect_true(all(c("from", "to", "ground_truth_stable",
                    "bootstrap_significant_run", "TP", "TN", "FP", "FN",
                    "p_value", "weight") %in% names(res$per_edge)))
  expect_equal(dim(res$TP_matrix), c(2L, 2L))
  expect_equal(dim(res$FN_matrix), c(2L, 2L))
})

test_that("run_bootstrap_iteration requires trans_matrix and init_probs", {
  expect_error(
    run_bootstrap_iteration(init_probs = make_inits(),
                            stable_transitions = make_stable()),
    "trans_matrix is required"
  )
  expect_error(
    run_bootstrap_iteration(trans_matrix = make_trans(),
                            stable_transitions = make_stable()),
    "init_probs is required"
  )
})

test_that("run_bootstrap_iteration is deterministic under outer set.seed", {
  skip_if_not_installed("tna")
  set.seed(123)
  r1 <- suppressWarnings(run_bootstrap_iteration(
    trans_matrix = make_trans(), init_probs = make_inits(),
    stable_transitions = make_stable(),
    seq_length = 12, n_sequences = 25, include_na = FALSE
  ))
  set.seed(123)
  r2 <- suppressWarnings(run_bootstrap_iteration(
    trans_matrix = make_trans(), init_probs = make_inits(),
    stable_transitions = make_stable(),
    seq_length = 12, n_sequences = 25, include_na = FALSE
  ))
  expect_identical(r1$per_edge, r2$per_edge)
})

# ===========================================================================
# run_bootstrap_simulation
# ===========================================================================

test_that("run_bootstrap_simulation smoke + structure (cores = 1)", {
  skip_if_not_installed("tna")
  res <- suppressWarnings(run_bootstrap_simulation(
    Model = make_model_list(), stable_transitions = make_stable(),
    num_runs = 2, seq_length = 12, n_sequences = 25,
    include_na = FALSE, num_cores = 1
  ))
  expect_type(res, "list")
  expect_setequal(
    names(res),
    c("aggregated_summary", "individual_runs", "successful_runs")
  )
  expect_setequal(
    names(res$aggregated_summary),
    c("overall_performance", "edge_significance")
  )
  op <- res$aggregated_summary$overall_performance
  expect_s3_class(op, "data.frame")
  expect_setequal(names(op), c("Metric", "Value"))
  expect_setequal(
    op$Metric,
    c("Sensitivity (TPR)", "Specificity (TNR)", "FPR", "FNR")
  )
  es <- res$aggregated_summary$edge_significance
  expect_true(all(c("from", "to", "n_significant", "recovery_rate",
                    "ground_truth_stable") %in% names(es)))
  expect_lte(res$successful_runs, 2)
})

test_that("run_bootstrap_simulation validates the Model argument", {
  expect_error(
    suppressWarnings(run_bootstrap_simulation(
      Model = list(weights = make_trans()),  # missing inits
      stable_transitions = make_stable(), num_runs = 1, num_cores = 1
    ))
  )
})

# ===========================================================================
# run_grid_simulation
# ===========================================================================

test_that("run_grid_simulation smoke + structure (1 grid point, cores = 1)", {
  skip_if_not_installed("tna")
  res <- suppressWarnings(run_grid_simulation(
    Model = make_model_list(), stable_transitions = make_stable(),
    num_runs = 2, n_sequences_vec = c(25), seq_length_vec = c(12),
    na_range_list = list(list(min = 0, max = 0)),
    include_na = FALSE, num_cores = 1
  ))
  expect_type(res, "list")
  expect_length(res, 1)
  expect_equal(names(res), "nr25_sl12_na0-0")
  expect_setequal(
    names(res[[1]]),
    c("aggregated_summary", "individual_runs", "successful_runs", "parameters")
  )
  expect_lte(res[[1]]$successful_runs, 2)
})

test_that("run_grid_simulation produces one element per grid combination", {
  skip_if_not_installed("tna")
  res <- suppressWarnings(run_grid_simulation(
    Model = make_model_list(), stable_transitions = make_stable(),
    num_runs = 2, n_sequences_vec = c(25), seq_length_vec = c(12, 14),
    na_range_list = list(list(min = 0, max = 0)),
    include_na = FALSE, num_cores = 1
  ))
  # 1 n_seq x 2 seq_len x 1 na_range = 2 combinations
  expect_length(res, 2)
  expect_setequal(names(res), c("nr25_sl12_na0-0", "nr25_sl14_na0-0"))
})

# ===========================================================================
# run_network_simulation
# ===========================================================================

test_that("run_network_simulation smoke + structure (sequential)", {
  skip_if_not_installed("tna")
  orig <- make_seqs(n = 30, len = 10, seed = 4)
  res <- suppressWarnings(run_network_simulation(
    original_data_list = orig,
    sim_params = list(seq_length = 10, n_sequences = 30, min_na = 0, max_na = 0),
    models = c("tna"), comparisons = c("original"),
    num_runs = 1, parallel = FALSE
  ))
  expect_type(res, "list")
  expect_setequal(names(res), c("metrics", "summary_stats", "parameters"))
  expect_s3_class(res$metrics, "data.frame")
  expect_gt(nrow(res$metrics), 0)
  expect_true(all(c("model_type", "comparison_type", "metric_category",
                    "metric_name", "value") %in% names(res$metrics)))
  expect_s3_class(res$summary_stats, "data.frame")
  expect_true(all(c("mean_value", "sd_value", "median_value",
                    "n_obs") %in% names(res$summary_stats)))
  expect_setequal(
    names(res$parameters),
    c("models", "comparisons", "num_runs", "scaling")
  )
})

test_that("run_network_simulation errors on empty data list", {
  expect_error(
    suppressWarnings(run_network_simulation(original_data_list = list())),
    "non-empty list"
  )
})

# ===========================================================================
# run_sampling_analysis
# ===========================================================================

test_that("run_sampling_analysis smoke + structure", {
  skip_if_not_installed("tna")
  model <- tna::tna(make_seqs(n = 40, len = 10, seed = 3))
  res <- suppressWarnings(run_sampling_analysis(
    model, sampling_percent = 0.4, iterations = 2, seed = 5, verbose = FALSE
  ))
  expect_type(res, "list")
  expect_setequal(names(res), c("aggregated", "individual", "params"))
  expect_s3_class(res$aggregated, "data.frame")
  expect_true(all(c("category", "metric", "mean", "sd", "median",
                    "q25", "q75", "n_iterations") %in% names(res$aggregated)))
  expect_s3_class(res$individual, "data.frame")
  expect_true(all(c("category", "metric", "value", "iteration") %in%
                    names(res$individual)))
  expect_true(all(c("sampling_percent", "iterations", "successful") %in%
                    names(res$params)))
})

test_that("run_sampling_analysis is deterministic under seed", {
  skip_if_not_installed("tna")
  model <- tna::tna(make_seqs(n = 40, len = 10, seed = 3))
  a <- suppressWarnings(run_sampling_analysis(
    model, sampling_percent = 0.4, iterations = 2, seed = 99, verbose = FALSE
  ))
  b <- suppressWarnings(run_sampling_analysis(
    model, sampling_percent = 0.4, iterations = 2, seed = 99, verbose = FALSE
  ))
  expect_identical(a$individual, b$individual)
})

test_that("run_sampling_analysis validates inputs", {
  skip_if_not_installed("tna")
  model <- tna::tna(make_seqs(n = 40, len = 10, seed = 3))
  expect_error(
    run_sampling_analysis(model, sampling_percent = 1.5),
    "between 0 and 1"
  )
  expect_error(
    run_sampling_analysis(model, iterations = 0),
    "at least 1"
  )
  expect_error(
    run_sampling_analysis("not a model"),
    "must be a TNA object"
  )
})
