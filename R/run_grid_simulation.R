#' Run Grid Search Over Simulation Parameters
#'
#' @description
#' Execute bootstrap simulations across a grid of parameter combinations.
#' Useful for studying how different settings affect edge recovery performance.
#'
#' @param Model A TNA model object with `weights` (transition matrix) and
#'   `inits` (initial probabilities).
#' @param stable_transitions List of character vectors defining ground truth
#'   stable transitions.
#' @param num_runs Integer. Number of bootstrap runs per parameter combination.
#' @param n_sequences_vec Numeric vector. Values of n_sequences to test.
#' @param seq_length_vec Numeric vector. Values of seq_length to test.
#' @param na_range_list List of lists. Each inner list has `min` and `max`
#'   elements defining an NA range to test.
#' @param stability_prob Numeric. Fixed stability probability. Default: 0.95.
#' @param unstable_mode Character. Fixed unstable mode. Default: "random_jump".
#' @param unstable_random_transition_prob Numeric. Fixed unstable probability.
#'   Default: 0.5.
#' @param unstable_perturb_noise Numeric. Fixed perturbation noise. Default: 0.5.
#' @param unlikely_prob_threshold Numeric. Fixed unlikely threshold. Default: 0.1.
#' @param include_na Logical. Whether to include NAs. Default: TRUE.
#' @param consistency_range Numeric vector of length 2. Bootstrap consistency
#'   range. Default: c(0.75, 1.25).
#' @param level Numeric. Significance level. Default: 0.05.
#' @param num_cores Integer. Number of cores for parallel processing.
#'   Default: detectCores() - 1.
#'
#' @param num_rows_vec Deprecated. Use `n_sequences_vec` instead.
#' @param max_seq_length_vec Deprecated. Use `seq_length_vec` instead.
#'
#' @return A named list where each element corresponds to a parameter combination.
#'   Names follow the pattern `nr<n_sequences>_sl<seq_length>_na<min>-<max>`.
#'   Each element contains:
#' \describe{
#'   \item{aggregated_summary}{Aggregated performance and edge significance.}
#'   \item{individual_runs}{Per-run details.}
#'   \item{successful_runs}{Number of successful runs.}
#'   \item{parameters}{The parameter values used for this combination.}
#' }
#'
#' @details
#' The function creates a full factorial grid from:
#' - `n_sequences_vec` x `seq_length_vec` x `na_range_list`
#'
#' For each combination, it runs `run_bootstrap_simulation()` and stores
#' the results along with the parameter values used.
#'
#' Progress messages are printed to track execution.
#'
#' @examples
#' \dontrun{
#' # Create a model
#' trans_mat <- matrix(c(
#'   0.6, 0.3, 0.1,
#'   0.2, 0.6, 0.2,
#'   0.1, 0.2, 0.7
#' ), nrow = 3, byrow = TRUE)
#' rownames(trans_mat) <- colnames(trans_mat) <- c("A", "B", "C")
#'
#' Model <- list(
#'   weights = trans_mat,
#'   inits = c(A = 0.33, B = 0.34, C = 0.33)
#' )
#'
#' stable <- list(c("A", "B"), c("B", "C"))
#'
#' # Run grid search
#' grid_results <- run_grid_simulation(
#'   Model = Model,
#'   stable_transitions = stable,
#'   num_runs = 20,
#'   n_sequences_vec = c(50, 100, 200),
#'   seq_length_vec = c(20, 30, 50),
#'   na_range_list = list(
#'     list(min = 0, max = 0),
#'     list(min = 0, max = 5),
#'     list(min = 5, max = 10)
#'   ),
#'   num_cores = 4
#' )
#'
#' # Analyze results
#' analyze_grid_results(grid_results)
#'
#' # Old parameter names still work
#' grid_results <- run_grid_simulation(
#'   Model = Model,
#'   stable_transitions = stable,
#'   num_runs = 20,
#'   num_rows_vec = c(50, 100, 200),
#'   max_seq_length_vec = c(20, 30, 50)
#' )
#' }
#'
#' @import dplyr
#' @importFrom parallel detectCores
#' @export
run_grid_simulation <- function(Model,
                                 stable_transitions,
                                 num_runs,
                                 n_sequences_vec = NULL,
                                 seq_length_vec = NULL,
                                 na_range_list = list(list(min = 0, max = 0)),
                                 stability_prob = 0.95,
                                 unstable_mode = "random_jump",
                                 unstable_random_transition_prob = 0.5,
                                 unstable_perturb_noise = 0.5,
                                 unlikely_prob_threshold = 0.1,
                                 include_na = TRUE,
                                 consistency_range = c(0.75, 1.25),
                                 level = 0.05,
                                 num_cores = parallel::detectCores() - 1,
                                 # Backward compatibility - old parameter names
                                 num_rows_vec = NULL,
                                 max_seq_length_vec = NULL) {
  # --- Backward compatibility: map old names to new names ---
  if (!is.null(num_rows_vec) && is.null(n_sequences_vec)) n_sequences_vec <- num_rows_vec
  if (!is.null(max_seq_length_vec) && is.null(seq_length_vec)) seq_length_vec <- max_seq_length_vec

  # Set defaults if still NULL
  if (is.null(n_sequences_vec)) n_sequences_vec <- c(50, 100, 200)
  if (is.null(seq_length_vec)) seq_length_vec <- c(20, 30)

  # --- Input Validation ---
  stopifnot(
    is.numeric(n_sequences_vec), all(n_sequences_vec > 0),
    is.numeric(seq_length_vec), all(seq_length_vec > 0),
    is.list(na_range_list)
  )

  # Create parameter grid
  na_params <- dplyr::bind_rows(lapply(na_range_list, function(x) data.frame(min_na = x$min, max_na = x$max)))
  param_grid <- expand.grid(n_sequences = n_sequences_vec, seq_length = seq_length_vec, stringsAsFactors = FALSE)
  param_grid <- dplyr::cross_join(param_grid, na_params)

  results_list <- vector("list", nrow(param_grid))
  total_grid_points <- nrow(param_grid)
  message(sprintf("Starting simulation grid with %d parameter combinations.", total_grid_points))

  # Loop through parameter grid
  for (i in 1:total_grid_points) {
    current_params_from_grid <- param_grid[i, ]
    message(sprintf(
      "Running grid point %d/%d: n_seq=%d, seq_len=%d, na=%d-%d",
      i, total_grid_points, current_params_from_grid$n_sequences, current_params_from_grid$seq_length,
      current_params_from_grid$min_na, current_params_from_grid$max_na
    ))
    current_include_na <- if (current_params_from_grid$max_na > 0) TRUE else include_na

    # Call the engine function (uses updated defaults if not overridden)
    sim_results <- run_bootstrap_simulation(
      Model = Model, stable_transitions = stable_transitions, num_runs = num_runs,
      seq_length = current_params_from_grid$seq_length, n_sequences = current_params_from_grid$n_sequences,
      na_range = c(current_params_from_grid$min_na, current_params_from_grid$max_na),
      include_na = current_include_na, stability_prob = stability_prob,
      unstable_mode = unstable_mode, unstable_random_transition_prob = unstable_random_transition_prob,
      unstable_perturb_noise = unstable_perturb_noise, unlikely_prob_threshold = unlikely_prob_threshold,
      consistency_range = consistency_range, level = level, num_cores = num_cores
    )

    # Add parameters to results (include both old and new names for compatibility)
    sim_results$parameters <- c(
      as.list(current_params_from_grid),
      list(
        # Old names for backward compatibility with analyze_grid_results
        num_rows = current_params_from_grid$n_sequences,
        max_seq_length = current_params_from_grid$seq_length,
        stability_prob = stability_prob, unstable_mode = unstable_mode,
        unstable_random_transition_prob = unstable_random_transition_prob,
        unstable_perturb_noise = unstable_perturb_noise,
        unlikely_prob_threshold = unlikely_prob_threshold,
        include_na = current_include_na,
        consistency_range = consistency_range, level = level,
        num_runs_attempted = num_runs, stable_transitions = stable_transitions
      )
    )

    results_list[[i]] <- sim_results
  }
  names(results_list) <- sprintf(
    "nr%d_sl%d_na%d-%d", param_grid$n_sequences, param_grid$seq_length, param_grid$min_na, param_grid$max_na
  )
  message("Simulation grid finished.")
  return(results_list)
}
