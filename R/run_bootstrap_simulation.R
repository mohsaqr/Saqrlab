#' Run Bootstrap Simulation Across Multiple Runs
#'
#' @description
#' Simulate and analyze edge recovery across multiple bootstrap runs. Aggregates
#' results to compute overall performance metrics and edge-level significance.
#'
#' @param Model A TNA model object with `weights` (transition matrix) and
#'   `inits` (initial probabilities).
#' @param stable_transitions List of character vectors defining ground truth
#'   stable transitions. Each vector contains two state names (from, to).
#' @param num_runs Integer. Number of simulation runs to perform.
#' @param max_seq_length Integer. Maximum sequence length.
#' @param num_rows Integer. Number of sequences per run.
#' @param stability_prob Numeric in (0 to 1). Probability of following stable
#'   transitions. Default: 0.95.
#' @param unstable_mode Character. Mode for unstable transitions:
#'   "random_jump", "perturb_prob", or "unlikely_jump". Default: "random_jump".
#' @param unstable_random_transition_prob Numeric in (0 to 1). Probability of
#'   unstable action. Default: 0.5.
#' @param unstable_perturb_noise Numeric in (0 to 1). Noise for perturbation mode.
#'   Default: 0.5.
#' @param unlikely_prob_threshold Numeric in (0 to 1). Threshold for unlikely
#'   transitions. Default: 0.1.
#' @param min_na Integer. Minimum NAs per sequence. Default: 0.
#' @param max_na Integer. Maximum NAs per sequence. Default: 0.
#' @param include_na Logical. Whether to include NAs. Default: TRUE.
#' @param consistency_range Numeric vector of length 2. Bootstrap consistency
#'   range. Default: c(0.75, 1.25).
#' @param level Numeric in (0,1). Significance level. Default: 0.05.
#' @param num_cores Integer. Number of cores for parallel processing.
#'   Default: detectCores() - 1.
#'
#' @return A list containing:
#' \describe{
#'   \item{aggregated_summary}{List with:
#'     \itemize{
#'       \item overall_performance: Data frame with Sensitivity, Specificity, FPR, FNR.
#'       \item edge_significance: Data frame with per-edge recovery rates and statistics.
#'     }
#'   }
#'   \item{individual_runs}{List with:
#'     \itemize{
#'       \item list_of_summaries: Bootstrap summaries from each run.
#'       \item list_of_per_edge_performance: Per-edge results from each run.
#'     }
#'   }
#'   \item{successful_runs}{Number of successfully completed runs.}
#' }
#'
#' @details
#' For each run, the function calls `evaluate_bootstrap()` and accumulates
#' TP/TN/FP/FN counts. After all runs, it computes:
#'
#' - **Overall Performance**: Aggregated sensitivity (TPR), specificity (TNR),
#'   false positive rate (FPR), and false negative rate (FNR).
#'
#' - **Edge Significance**: For each edge, the number of times it was detected
#'   as significant, average p-value, average weight, and recovery rate.
#'
#' Parallel processing uses `parallel::mclapply()`.
#'
#' @examples
#' \dontrun{
#' # First create a TNA model from some data
#' trans_mat <- matrix(c(
#'   0.6, 0.3, 0.1,
#'   0.2, 0.6, 0.2,
#'   0.1, 0.2, 0.7
#' ), nrow = 3, byrow = TRUE)
#' rownames(trans_mat) <- colnames(trans_mat) <- c("A", "B", "C")
#'
#' # Create a mock Model object
#' Model <- list(
#'   weights = trans_mat,
#'   inits = c(A = 0.33, B = 0.34, C = 0.33)
#' )
#'
#' # Define stable transitions
#' stable <- list(c("A", "B"), c("B", "C"))
#'
#' # Run bootstrap simulation
#' results <- run_bootstrap_simulation(
#'   Model = Model,
#'   stable_transitions = stable,
#'   num_runs = 50,
#'   max_seq_length = 30,
#'   num_rows = 100,
#'   stability_prob = 0.95,
#'   unstable_mode = "random_jump",
#'   num_cores = 4
#' )
#'
#' # View overall performance
#' results$aggregated_summary$overall_performance
#'
#' # View edge-level results
#' results$aggregated_summary$edge_significance
#' }
#'
#' @import dplyr
#' @importFrom tidyr complete expand_grid
#' @importFrom parallel mclapply detectCores
#' @export
run_bootstrap_simulation <- function(Model,
                                      stable_transitions,
                                      num_runs,
                                      max_seq_length,
                                      num_rows,
                                      stability_prob = 0.95,
                                      unstable_mode = "random_jump",
                                      unstable_random_transition_prob = 0.5,
                                      unstable_perturb_noise = 0.5,
                                      unlikely_prob_threshold = 0.1,
                                      min_na = 0,
                                      max_na = 0,
                                      include_na = TRUE,
                                      consistency_range = c(0.75, 1.25),
                                      level = 0.05,
                                      num_cores = parallel::detectCores() - 1) {
  # --- Input Validation ---
  stopifnot(
    is.list(Model), all(c("weights", "inits") %in% names(Model)), is.matrix(Model$weights),
    is.numeric(Model$weights), is.numeric(Model$inits), nrow(Model$weights) == length(Model$inits),
    is.list(stable_transitions) || is.null(stable_transitions),
    is.numeric(num_runs), num_runs > 0,
    is.numeric(num_cores), num_cores > 0
  )

  transition_matrix <- Model$weights
  initial_probabilities <- Model$inits
  states_names <- rownames(transition_matrix)
  if (is.null(states_names)) stop("Model$weights must have rownames.")
  num_states <- length(states_names)

  # Initialize accumulators
  total_TP <- matrix(0L, num_states, num_states, dimnames = list(states_names, states_names))
  total_TN <- matrix(0L, num_states, num_states, dimnames = list(states_names, states_names))
  total_FP <- matrix(0L, num_states, num_states, dimnames = list(states_names, states_names))
  total_FN <- matrix(0L, num_states, num_states, dimnames = list(states_names, states_names))

  # Prepare parameter list
  run_params_list <- lapply(1:num_runs, function(X) {
    list(
      run_id = X,
      transition_matrix = transition_matrix, initial_probabilities = initial_probabilities,
      stable_transitions = stable_transitions, max_seq_length = max_seq_length,
      num_rows = num_rows, stability_prob = stability_prob,
      unstable_mode = unstable_mode, unstable_random_transition_prob = unstable_random_transition_prob,
      unstable_perturb_noise = unstable_perturb_noise, unlikely_prob_threshold = unlikely_prob_threshold,
      min_na = min_na, max_na = max_na, include_na = include_na,
      consistency_range = consistency_range, level = level
    )
  })

  # Parallel execution
  ListofRunResults <- parallel::mclapply(
    run_params_list, function(params) {
      run_result <- tryCatch(
        {
          evaluate_bootstrap(
            # Pass all params by name
            transition_matrix = params$transition_matrix, initial_probabilities = params$initial_probabilities,
            stable_transitions = params$stable_transitions, max_seq_length = params$max_seq_length,
            num_rows = params$num_rows, stability_prob = params$stability_prob,
            unstable_mode = params$unstable_mode, unstable_random_transition_prob = params$unstable_random_transition_prob,
            unstable_perturb_noise = params$unstable_perturb_noise, unlikely_prob_threshold = params$unlikely_prob_threshold,
            min_na = params$min_na, max_na = params$max_na, include_na = params$include_na,
            consistency_range = params$consistency_range, level = params$level
          )
        },
        error = function(e) {
          warning(sprintf("Run %d evaluation failed: %s", params$run_id, e$message))
          return(NULL)
        }
      )
      if (!is.null(run_result)) {
        run_result$run_id <- params$run_id
      }
      return(run_result)
    },
    mc.cores = num_cores
  )

  # --- Post-processing ---
  ListofRunResults <- Filter(Negate(is.null), ListofRunResults)
  successful_runs <- length(ListofRunResults)

  empty_overall_perf <- data.frame(Metric = c("Sensitivity (TPR)", "Specificity (TNR)", "FPR", "FNR"), Value = NA_real_)
  empty_edge_sig <- data.frame(
    from = character(), to = character(), ground_truth_stable = logical(),
    n_significant = integer(), avg_p_value = numeric(), avg_weight = numeric(),
    recovery_rate = numeric(), stringsAsFactors = FALSE
  )
  empty_return_list <- list(
    aggregated_summary = list(overall_performance = empty_overall_perf, edge_significance = empty_edge_sig),
    individual_runs = list(list_of_summaries = list(), list_of_per_edge_performance = list()),
    successful_runs = 0
  )

  if (successful_runs == 0) {
    warning("No successful simulation runs.")
    return(empty_return_list)
  }

  # Create Lists of Raw Data
  ListOfSummaries <- lapply(ListofRunResults, function(x) if (!is.null(x$bootstrap_summary_raw)) x$bootstrap_summary_raw else NULL)
  ListOfPerEdgePerformance <- lapply(ListofRunResults, function(x) if (!is.null(x$per_edge)) x$per_edge else NULL)
  ListOfSummaries <- Filter(Negate(is.null), ListOfSummaries)
  ListOfPerEdgePerformance <- Filter(Negate(is.null), ListOfPerEdgePerformance)
  run_ids <- sapply(ListofRunResults, `[[`, "run_id")
  if (length(ListOfSummaries) == successful_runs) names(ListOfSummaries) <- paste0("run_", run_ids)
  if (length(ListOfPerEdgePerformance) == successful_runs) names(ListOfPerEdgePerformance) <- paste0("run_", run_ids)

  # Accumulate TP/TN/FP/FN
  valid_matrices_found <- FALSE
  for (run_result in ListofRunResults) {
    if (!is.null(run_result$TP_matrix) && is.matrix(run_result$TP_matrix) && all(dim(run_result$TP_matrix) == c(num_states, num_states))) {
      total_TP <- total_TP + run_result$TP_matrix
      valid_matrices_found <- TRUE
    } else {
      if (!is.null(run_result$TP_matrix)) warning("Invalid TP_matrix.")
    }
    if (!is.null(run_result$TN_matrix) && is.matrix(run_result$TN_matrix) && all(dim(run_result$TN_matrix) == c(num_states, num_states))) {
      total_TN <- total_TN + run_result$TN_matrix
    } else {
      if (!is.null(run_result$TN_matrix)) warning("Invalid TN_matrix.")
    }
    if (!is.null(run_result$FP_matrix) && is.matrix(run_result$FP_matrix) && all(dim(run_result$FP_matrix) == c(num_states, num_states))) {
      total_FP <- total_FP + run_result$FP_matrix
    } else {
      if (!is.null(run_result$FP_matrix)) warning("Invalid FP_matrix.")
    }
    if (!is.null(run_result$FN_matrix) && is.matrix(run_result$FN_matrix) && all(dim(run_result$FN_matrix) == c(num_states, num_states))) {
      total_FN <- total_FN + run_result$FN_matrix
    } else {
      if (!is.null(run_result$FN_matrix)) warning("Invalid FN_matrix.")
    }
  }

  if (!valid_matrices_found) {
    warning("Could not accumulate valid TP/TN/FP/FN counts.")
    empty_return_list$successful_runs <- successful_runs
    return(empty_return_list)
  }

  # Create Aggregated Edge Significance Summary
  CombinedPerEdgePerformance <- suppressWarnings(tryCatch(
    dplyr::bind_rows(ListOfPerEdgePerformance),
    error = function(e) {
      warning("Failed to bind per_edge performance: ", e$message)
      data.frame()
    }
  ))

  # Helper function for ground truth DF (used in aggregation)
  create_ground_truth_df <- function(states, transitions) {
    gt_df <- tidyr::expand_grid(from = states, to = states) %>% dplyr::mutate(gt_stable = FALSE) # Use distinct name
    if (!is.null(transitions)) {
      for (pair in transitions) {
        if (length(pair) == 2 && all(pair %in% states)) {
          gt_df$gt_stable[gt_df$from == pair[1] & gt_df$to == pair[2]] <- TRUE
        }
      }
    }
    return(gt_df)
  }

  aggregated_edge_summary <- if (!is.null(CombinedPerEdgePerformance) && nrow(CombinedPerEdgePerformance) > 0) {
    gt_col_present <- "ground_truth_stable" %in% names(CombinedPerEdgePerformance)
    if (!gt_col_present) warning("Column 'ground_truth_stable' missing in combined performance data.")

    CombinedPerEdgePerformance %>%
      dplyr::mutate(
        p_value_num = suppressWarnings(as.numeric(p_value)),
        weight_num = suppressWarnings(as.numeric(weight))
      ) %>%
      dplyr::filter(!is.na(from), !is.na(to)) %>%
      dplyr::group_by(from, to) %>%
      dplyr::summarise(
        n_significant = sum(bootstrap_significant_run, na.rm = TRUE),
        avg_p_value = mean(p_value_num, na.rm = TRUE),
        avg_weight = mean(weight_num, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        recovery_rate = n_significant / successful_runs
      ) %>%
      # Ensure all edges present and definitive ground truth
      tidyr::complete(from = states_names, to = states_names, fill = list(n_significant = 0L, recovery_rate = 0.0)) %>%
      dplyr::left_join(create_ground_truth_df(states_names, stable_transitions), by = c("from", "to")) %>%
      # Use the definitive ground truth from create_ground_truth_df
      dplyr::mutate(ground_truth_stable = gt_stable) %>%
      dplyr::select(-dplyr::any_of("gt_stable"), -dplyr::any_of("ground_truth_stable_agg")) %>% # Clean up temp columns
      # Final coalesce for numeric summaries
      dplyr::mutate(
        n_significant = dplyr::coalesce(n_significant, 0L),
        recovery_rate = dplyr::coalesce(recovery_rate, 0.0)
      ) %>%
      dplyr::arrange(dplyr::desc(ground_truth_stable), dplyr::desc(recovery_rate)) # Sort stable first
  } else {
    empty_edge_sig
  }

  # Calculate OVERALL Performance Metrics
  overall_sensitivity <- if (sum(total_TP) + sum(total_FN) > 0) sum(total_TP) / (sum(total_TP) + sum(total_FN)) else 0
  overall_specificity <- if (sum(total_TN) + sum(total_FP) > 0) sum(total_TN) / (sum(total_TN) + sum(total_FP)) else 0
  overall_fpr <- if (sum(total_FP) + sum(total_TN) > 0) sum(total_FP) / (sum(total_FP) + sum(total_TN)) else 0
  overall_fnr <- if (sum(total_FN) + sum(total_TP) > 0) sum(total_FN) / (sum(total_FN) + sum(total_TP)) else 0

  overall_performance_summary <- data.frame(
    Metric = c("Sensitivity (TPR)", "Specificity (TNR)", "FPR", "FNR"),
    Value = c(overall_sensitivity, overall_specificity, overall_fpr, overall_fnr)
  )

  # Return Structured List
  return(list(
    aggregated_summary = list(
      overall_performance = overall_performance_summary,
      edge_significance = aggregated_edge_summary
    ),
    individual_runs = list(
      list_of_summaries = ListOfSummaries,
      list_of_per_edge_performance = ListOfPerEdgePerformance
    ),
    successful_runs = successful_runs
  ))
}
