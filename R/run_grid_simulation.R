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
#' @param num_rows_vec Numeric vector. Values of num_rows to test.
#' @param max_seq_length_vec Numeric vector. Values of max_seq_length to test.
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
#' @return A named list where each element corresponds to a parameter combination.
#'   Names follow the pattern `nr<num_rows>_sl<seq_length>_na<min>-<max>`.
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
#' - `num_rows_vec` x `max_seq_length_vec` x `na_range_list`
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
#'   num_rows_vec = c(50, 100, 200),
#'   max_seq_length_vec = c(20, 30, 50),
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
#' }
#'
#' @import dplyr
#' @importFrom parallel detectCores
#' @export
run_grid_simulation <- function(Model,
                                 stable_transitions,
                                 num_runs,
                                 num_rows_vec,
                                 max_seq_length_vec,
                                 na_range_list,
                                 stability_prob = 0.95,
                                 unstable_mode = "random_jump",
                                 unstable_random_transition_prob = 0.5,
                                 unstable_perturb_noise = 0.5,
                                 unlikely_prob_threshold = 0.1,
                                 include_na = TRUE,
                                 consistency_range = c(0.75, 1.25),
                                 level = 0.05,
                                 num_cores = parallel::detectCores() - 1) {
  # --- Input Validation ---
  stopifnot(
    is.numeric(num_rows_vec), all(num_rows_vec > 0),
    is.numeric(max_seq_length_vec), all(max_seq_length_vec > 0),
    is.list(na_range_list)
  )

  # Create parameter grid
  na_params <- dplyr::bind_rows(lapply(na_range_list, function(x) data.frame(min_na = x$min, max_na = x$max)))
  param_grid <- expand.grid(num_rows = num_rows_vec, max_seq_length = max_seq_length_vec, stringsAsFactors = FALSE)
  param_grid <- dplyr::cross_join(param_grid, na_params)

  results_list <- vector("list", nrow(param_grid))
  total_grid_points <- nrow(param_grid)
  message(sprintf("Starting simulation grid with %d parameter combinations.", total_grid_points))

  # Loop through parameter grid
  for (i in 1:total_grid_points) {
    current_params_from_grid <- param_grid[i, ]
    message(sprintf(
      "Running grid point %d/%d: nr=%d, sl=%d, na=%d-%d",
      i, total_grid_points, current_params_from_grid$num_rows, current_params_from_grid$max_seq_length,
      current_params_from_grid$min_na, current_params_from_grid$max_na
    ))
    current_include_na <- if (current_params_from_grid$max_na > 0) TRUE else include_na

    # Call the engine function (uses updated defaults if not overridden)
    sim_results <- run_bootstrap_simulation(
      Model = Model, stable_transitions = stable_transitions, num_runs = num_runs,
      max_seq_length = current_params_from_grid$max_seq_length, num_rows = current_params_from_grid$num_rows,
      min_na = current_params_from_grid$min_na, max_na = current_params_from_grid$max_na,
      include_na = current_include_na, stability_prob = stability_prob,
      unstable_mode = unstable_mode, unstable_random_transition_prob = unstable_random_transition_prob,
      unstable_perturb_noise = unstable_perturb_noise, unlikely_prob_threshold = unlikely_prob_threshold,
      consistency_range = consistency_range, level = level, num_cores = num_cores
    )

    # Add parameters to results
    sim_results$parameters <- c(
      as.list(current_params_from_grid),
      list(
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
    "nr%d_sl%d_na%d-%d", param_grid$num_rows, param_grid$max_seq_length, param_grid$min_na, param_grid$max_na
  )
  message("Simulation grid finished.")
  return(results_list)
}
