# Characterization / regression tests for summary, plotting, utility, and
# miscellaneous-simulator functions. Assertions match ACTUAL observed behavior
# (probed against devtools::load_all on 2026-06-20), not assumptions.
#
# PREVIOUSLY-BROKEN (now FIXED, characterized below):
#   summarize_grid_results() / analyze_grid_results() used to call the undefined
#   helper check_val_in_range() and errored on any valid run_grid_simulation()
#   output. The missing helpers (check_val_in_range(), safe_bind_rows()) are now
#   implemented, so both functions work and are characterized with real tests.

# ---------------------------------------------------------------------------
# Shared cheap fixtures (built once, reused across tests)
# ---------------------------------------------------------------------------

make_bootstrap_summary_df <- function() {
  mat <- matrix(c(0.6, 0.4, 0.3, 0.7), 2, 2, byrow = TRUE,
                dimnames = list(c("A", "B"), c("A", "B")))
  bs <- suppressWarnings(suppressMessages(run_bootstrap_simulation(
    Model = list(weights = mat, inits = c(A = 0.5, B = 0.5)),
    stable_transitions = list(c("A", "B")),
    num_runs = 4, seq_length = 12, n_sequences = 40, num_cores = 1
  )))
  do.call(rbind, bs$individual_runs$list_of_summaries)
}

make_sampling_individual <- function() {
  nets <- suppressWarnings(suppressMessages(simulate_tna_networks(
    n_networks = 1, n_states = 4, n_sequences = 40, seq_length = 12,
    seed = 42, verbose = FALSE
  )))
  sa <- suppressWarnings(suppressMessages(run_sampling_analysis(
    nets[[1]], sampling_percent = 0.3, iterations = 12, seed = 7, verbose = FALSE
  )))
  sa$individual
}


# ---------------------------------------------------------------------------
# GROUP_REGULATION_ACTIONS (exported constant)
# ---------------------------------------------------------------------------

test_that("GROUP_REGULATION_ACTIONS is the expected character vector", {
  expect_true(exists("GROUP_REGULATION_ACTIONS"))
  expect_type(GROUP_REGULATION_ACTIONS, "character")
  expect_length(GROUP_REGULATION_ACTIONS, 9)
  expect_setequal(
    GROUP_REGULATION_ACTIONS,
    c("adapt", "cohesion", "consensus", "coregulate", "discuss",
      "emotion", "monitor", "plan", "synthesis")
  )
  expect_false(any(duplicated(GROUP_REGULATION_ACTIONS)))
})


# ---------------------------------------------------------------------------
# validate_sim_params
# ---------------------------------------------------------------------------

test_that("validate_sim_params(list()) fills defaults including legacy names", {
  v <- validate_sim_params(list())
  expect_type(v, "list")
  expect_setequal(
    names(v),
    c("seq_length", "n_sequences", "min_na", "max_na",
      "max_seq_length", "num_rows")
  )
  expect_equal(v$seq_length, 30)
  expect_equal(v$n_sequences, 100)
  expect_equal(v$min_na, 0)
  # legacy aliases mirror the canonical values
  expect_equal(v$max_seq_length, v$seq_length)
  expect_equal(v$num_rows, v$n_sequences)
})

test_that("validate_sim_params honours supplied values", {
  v <- validate_sim_params(list(seq_length = 15, n_sequences = 42))
  expect_equal(v$seq_length, 15)
  expect_equal(v$n_sequences, 42)
  expect_equal(v$max_seq_length, 15)
  expect_equal(v$num_rows, 42)
})


# ---------------------------------------------------------------------------
# smart_select_states (alias of select_states)
# ---------------------------------------------------------------------------

test_that("smart_select_states returns n distinct states, reproducibly", {
  s <- smart_select_states(5, seed = 1)
  expect_type(s, "character")
  expect_length(s, 5)
  expect_false(any(duplicated(s)))

  s2 <- smart_select_states(5, seed = 1)
  expect_identical(s, s2)
})


# ---------------------------------------------------------------------------
# simulate_long_data
# ---------------------------------------------------------------------------

test_that("simulate_long_data returns expected long tibble and is reproducible", {
  ld <- simulate_long_data(n_groups = 2, n_actors = 3, n_courses = 1,
                           n_states = 4, seq_length_range = c(5, 8), seed = 42)
  expect_s3_class(ld, "data.frame")
  expect_setequal(names(ld),
                  c("Actor", "Achiever", "Group", "Course", "Time", "Action"))
  expect_gt(nrow(ld), 0)

  ld2 <- simulate_long_data(n_groups = 2, n_actors = 3, n_courses = 1,
                            n_states = 4, seq_length_range = c(5, 8), seed = 42)
  expect_equal(as.data.frame(ld), as.data.frame(ld2))
})


# ---------------------------------------------------------------------------
# simulate_onehot_data
# ---------------------------------------------------------------------------

test_that("simulate_onehot_data returns binary state columns, reproducibly", {
  oh <- simulate_onehot_data(n_groups = 2, n_actors = 3, n_courses = 1,
                             n_states = 4, seq_length_range = c(5, 8), seed = 42)
  expect_s3_class(oh, "data.frame")
  meta_cols <- c("Actor", "Achiever", "Group", "Course", "Time")
  expect_true(all(meta_cols %in% names(oh)))
  expect_gt(nrow(oh), 0)

  state_cols <- setdiff(names(oh), meta_cols)
  expect_gt(length(state_cols), 0)
  state_vals <- unlist(oh[state_cols], use.names = FALSE)
  expect_true(all(state_vals %in% c(0, 1)))

  oh2 <- simulate_onehot_data(n_groups = 2, n_actors = 3, n_courses = 1,
                              n_states = 4, seq_length_range = c(5, 8), seed = 42)
  expect_equal(as.data.frame(oh), as.data.frame(oh2))
})


# ---------------------------------------------------------------------------
# simulate_mtna / simulate_mlna (aliases of simulate_htna)
# ---------------------------------------------------------------------------

test_that("simulate_mtna returns a hierarchical matrix list, reproducibly", {
  mt <- simulate_mtna(seed = 42)
  expect_type(mt, "list")
  expect_setequal(names(mt),
                  c("matrix", "node_types", "type_names", "n_nodes_per_type"))
  expect_true(is.matrix(mt$matrix))

  mt2 <- simulate_mtna(seed = 42)
  expect_identical(mt, mt2)
})

test_that("simulate_mlna is an alias producing the same output as simulate_mtna", {
  expect_identical(simulate_mlna(seed = 42), simulate_mtna(seed = 42))
})


# ---------------------------------------------------------------------------
# summarize_networks
# ---------------------------------------------------------------------------

test_that("summarize_networks summarizes a list of TNA models", {
  nets <- suppressWarnings(suppressMessages(simulate_tna_networks(
    n_networks = 3, n_states = 4, n_sequences = 40, seq_length = 12,
    seed = 42, verbose = FALSE
  )))
  sn <- suppressWarnings(summarize_networks(nets))

  expect_s3_class(sn, "network_summary")
  expect_setequal(names(sn), c("summary_table", "aggregate", "n_networks"))
  expect_equal(sn$n_networks, 3L)
  expect_s3_class(sn$summary_table, "data.frame")
  expect_equal(nrow(sn$summary_table), 3L)
  expect_type(sn$aggregate, "list")
})


# ---------------------------------------------------------------------------
# summarize_simulation
# ---------------------------------------------------------------------------

test_that("summarize_simulation returns a one-row numeric summary", {
  df <- make_bootstrap_summary_df()
  ss <- summarize_simulation(df)
  expect_s3_class(ss, "data.frame")
  expect_equal(nrow(ss), 1L)
  expect_true("weight_mean" %in% names(ss))
  expect_true("weight_ci_lower" %in% names(ss))
})

test_that("summarize_simulation with 'by' groups the summary", {
  df <- make_bootstrap_summary_df()
  ss_by <- summarize_simulation(df, by = "from")
  expect_s3_class(ss_by, "data.frame")
  expect_true("from" %in% names(ss_by))
  # one row per distinct 'from' level
  expect_equal(nrow(ss_by), length(unique(df$from)))
})


# ---------------------------------------------------------------------------
# evaluate_bootstrap (alias of run_bootstrap_iteration)
# ---------------------------------------------------------------------------

test_that("evaluate_bootstrap returns per-edge + confusion-matrix components", {
  mat <- matrix(c(0.6, 0.4, 0.3, 0.7), 2, 2, byrow = TRUE,
                dimnames = list(c("A", "B"), c("A", "B")))
  eb <- suppressWarnings(suppressMessages(evaluate_bootstrap(
    trans_matrix = mat, init_probs = c(A = 0.5, B = 0.5),
    stable_transitions = list(c("A", "B")),
    seq_length = 12, n_sequences = 40
  )))
  expect_type(eb, "list")
  expect_true(all(c("per_edge", "TP_matrix", "TN_matrix",
                    "FP_matrix", "FN_matrix") %in% names(eb)))
  expect_s3_class(eb$per_edge, "data.frame")
})


# ---------------------------------------------------------------------------
# summarize_grid_results / analyze_grid_results  (FIXED: helpers implemented)
# ---------------------------------------------------------------------------

make_grid_results <- function() {
  trans <- matrix(c(.7, .3, .4, .6), 2, 2, byrow = TRUE,
                  dimnames = list(c("A", "B"), c("A", "B")))
  Model <- list(weights = trans, inits = c(A = .5, B = .5))
  suppressWarnings(suppressMessages(run_grid_simulation(
    Model = Model, stable_transitions = list(c("A", "B")),
    num_runs = 2, n_sequences_vec = 20, seq_length_vec = 10, num_cores = 1
  )))
}

test_that("summarize_grid_results returns the expected 4-element list", {
  skip_if_not_installed("tna")
  g <- make_grid_results()
  s <- suppressWarnings(summarize_grid_results(g, print_output = FALSE))

  expect_type(s, "list")
  expect_setequal(
    names(s),
    c("n_selected", "aggregated_summary",
      "selected_settings_summary_df", "compiled_individual_runs")
  )
  expect_equal(s$n_selected, 1L)
  expect_s3_class(s$selected_settings_summary_df, "data.frame")
  expect_equal(nrow(s$selected_settings_summary_df), 1L)
  expect_type(s$aggregated_summary, "list")
  expect_type(s$compiled_individual_runs, "list")
})

test_that("analyze_grid_results is an alias matching summarize_grid_results", {
  skip_if_not_installed("tna")
  g <- make_grid_results()
  s <- suppressWarnings(summarize_grid_results(g, print_output = FALSE))
  a <- suppressWarnings(analyze_grid_results(g, print_output = FALSE))

  expect_type(a, "list")
  expect_identical(class(a), class(s))
  expect_setequal(names(a), names(s))
  expect_equal(a$n_selected, s$n_selected)
})


# ---------------------------------------------------------------------------
# Plot functions (assert object type only; never inspect pixels)
# ---------------------------------------------------------------------------

test_that("plot_sampling_distribution returns a ggplot", {
  skip_if_not_installed("ggplot2")
  ind <- make_sampling_individual()
  cat_val <- unique(ind$category)[1]
  met_val <- unique(ind$metric)[1]
  p <- plot_sampling_distribution(ind, category = cat_val, metric = met_val)
  expect_s3_class(p, "ggplot")
})

test_that("plot_tna_comparison returns a ggplot", {
  skip_if_not_installed("ggplot2")
  ind <- make_sampling_individual()
  met_val <- unique(ind$metric)[1]
  p <- plot_tna_comparison(ind, plot_type = "histogram", metric_name = met_val)
  expect_s3_class(p, "ggplot")
})

test_that("plot_network_estimation returns a ggplot", {
  skip_if_not_installed("ggplot2")
  seqs <- suppressWarnings(suppressMessages(simulate_sequences(
    n_sequences = 40, seq_length = 12, n_states = 4, seed = 11
  )))
  ne <- suppressWarnings(suppressMessages(compare_network_estimation(
    seqs, model_types = c("tna", "ftna"), iterations = 10, seed = 5,
    verbose = FALSE
  )))
  p <- plot_network_estimation(ne, metric_name = "Pearson", plot_type = "histogram")
  expect_s3_class(p, "ggplot")
})
