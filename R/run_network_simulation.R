#' Run Network Analysis Simulations
#'
#' @description
#' Run TNA simulations comparing fitted models to original data. Supports
#' multiple model types, comparison modes, and parallel processing.
#'
#' @param original_data_list A list of data frames containing original sequence
#'   data, or a single data frame. Used as reference for comparisons.
#' @param sim_params Simulation parameters. Can be:
#' \itemize{
#'   \item NULL (uses defaults).
#'   \item A list of parameters.
#'   \item A data frame where each row is a parameter combination.
#' }
#' @param models Character vector of model types to simulate. Options:
#'   "tna", "ftna", "ctna", "atna". Default: c("tna").
#' @param comparisons Character vector of comparison types. Options:
#' \describe{
#'   \item{"original"}{Compare simulated models to original data.}
#'   \item{"across_models"}{Compare different model types against each other.}
#'   \item{"ftna_reference"}{Compare all models to their ftna counterparts.}
#' }
#' Default: c("original").
#' @param num_runs Integer. Number of simulation runs per parameter set.
#'   Default: 3.
#' @param parallel Logical. Whether to use parallel processing. Default: FALSE.
#' @param scaling Character. Scaling method for comparisons.
#'   Options: "none", "minmax", "zscore". Default: "minmax".
#'
#' @return A list containing:
#' \describe{
#'   \item{metrics}{Data frame of all comparison metrics from all runs.}
#'   \item{summary_stats}{Aggregated summary statistics by model and comparison type.}
#'   \item{parameters}{List of simulation parameters used.}
#' }
#'
#' @details
#' For each parameter combination and original dataset, the function:
#' \enumerate{
#'   \item Extracts transition probabilities from the original TNA model.
#'   \item Generates simulated sequences using those probabilities.
#'   \item Fits the specified model types to the simulated sequences.
#'   \item Computes comparison metrics against reference models.
#' }
#'
#' Parallel processing uses the future.apply package when enabled.
#'
#' @examples
#' \dontrun{
#' # Create some original sequence data
#' trans_mat <- matrix(c(
#'   0.7, 0.2, 0.1,
#'   0.3, 0.5, 0.2,
#'   0.2, 0.3, 0.5
#' ), nrow = 3, byrow = TRUE)
#' rownames(trans_mat) <- colnames(trans_mat) <- c("A", "B", "C")
#' init_probs <- c(A = 0.5, B = 0.3, C = 0.2)
#'
#' original_data <- simulate_sequences(
#'   trans_matrix = trans_mat,
#'   init_probs = init_probs,
#'   seq_length = 30,
#'   n_sequences = 100
#' )
#'
#' # Run simulations comparing tna and ftna models
#' results <- run_network_simulation(
#'   original_data_list = original_data,
#'   sim_params = list(seq_length = 30, n_sequences = 100),
#'   models = c("tna", "ftna"),
#'   comparisons = c("original"),
#'   num_runs = 5
#' )
#'
#' # View summary statistics
#' results$summary_stats
#' }
#'
#' @import dplyr
#' @import tna
#' @importFrom future plan
#' @importFrom future.apply future_lapply
#' @importFrom progressr with_progress progressor
#' @importFrom stats median
#' @export
run_network_simulation <- function(original_data_list,
                                    sim_params = NULL,
                                    models = c("tna"),
                                    comparisons = c("original"),
                                    num_runs = 3,
                                    parallel = FALSE,
                                    scaling = "minmax") {
  # --- Setup parallel processing if requested ---
  if (parallel && requireNamespace("future.apply", quietly = TRUE)) {
    future::plan(future::multisession)
    cat("Parallel processing enabled\n")
  }

  # --- Handle single data frame input ---
  if (is.data.frame(original_data_list)) {
    original_data_list <- list(original_data_list)
  }

  # --- Validate original data list ---
  if (!is.list(original_data_list) || length(original_data_list) == 0) {
    stop("original_data_list must be a non-empty list or a data frame.")
  }

  # --- Process simulation parameters ---
  if (is.null(sim_params)) {
    # Default parameter if none provided
    sim_params <- list(
      list(seq_length = 30, n_sequences = 100, min_na = 0, max_na = 5)
    )
  } else if (is.data.frame(sim_params)) {
    # Convert data frame to list of parameter sets
    sim_params <- lapply(1:nrow(sim_params), function(i) {
      as.list(sim_params[i, ])
    })
  } else if (!is.list(sim_params[[1]])) {
    # Single parameter set provided as a flat list
    sim_params <- list(sim_params)
  }

  # --- Define valid model types and comparisons ---
  valid_models <- c("tna", "ftna", "ctna", "atna")
  valid_comparisons <- c(
    "original", # Compare simulated models to original data
    "across_models", # Compare different model types against each other
    "ftna_reference" # Compare all models to their ftna counterparts
  )

  # Validate inputs
  models <- match.arg(models, valid_models, several.ok = TRUE)
  comparisons <- match.arg(comparisons, valid_comparisons, several.ok = TRUE)

  # --- Prepare results containers ---
  metrics_list <- list()


  # --- Define a single run function for potential parallelization ---
  run_single_param_set <- function(param_set_idx) {
    param_set <- sim_params[[param_set_idx]]

    # Ensure parameters are valid
    param_set <- validate_sim_params(param_set)
    cat(sprintf(
      "Running simulation %d/%d with parameters: %s\n",
      param_set_idx, length(sim_params),
      paste(names(param_set), param_set, sep = "=", collapse = ", ")
    ))

    # Container for metrics from this parameter set
    param_metrics_list <- list()
    # --- Iterate over original datasets ---
    for (data_idx in 1:length(original_data_list)) {
      original_data <- original_data_list[[data_idx]]

      # --- Pre-compute original models for this dataset ---
      original_tna_model <- tna::tna(original_data)
      original_ftna_model <- tna::ftna(original_data)

      # --- Extract transition probabilities once ---
      original_transition_probs <- original_tna_model$weights
      initial_probabilities <- original_tna_model$inits

      # Get number of states from the original data
      num_states <- length(initial_probabilities)


      # --- Process each model type ---
      for (model_type in models) {
        # Run multiple simulations for this model type
        for (run_idx in 1:num_runs) {
          # Generate sequences only once per run
          sim_sequences <- simulate_sequences(
            trans_matrix = original_transition_probs,
            init_probs = initial_probabilities,
            seq_length = param_set$seq_length,
            n_sequences = param_set$n_sequences,
            na_range = c(param_set$min_na, param_set$max_na)
          )

          # Check if sequences are valid
          if (all(is.na(sim_sequences))) {
            warning(paste("Invalid sequences generated for model:", model_type, "run:", run_idx, "param_set:", param_set_idx, "dataset:", data_idx))
            next
          }

          # Fit models - only compute what's needed
          sim_model <- fit_network_model(sim_sequences, model_type)

          # Check if model is valid
          if (is.null(sim_model) || any(is.na(sim_model$weights))) {
            warning(paste("Invalid model generated for model:", model_type, "run:", run_idx, "param_set:", param_set_idx, "dataset:", data_idx))
            next
          }

          # Only compute ftna if needed
          sim_ftna_model <- if ("ftna_reference" %in% comparisons && model_type != "ftna") {
            tna::ftna(sim_sequences)
          } else if (model_type == "ftna") {
            sim_model # If it's already an ftna model, use it directly
          } else {
            NULL
          }

          # Process requested comparisons
          if ("original" %in% comparisons) {
            # Compare to original models with appropriate scaling
            vs_orig_tna_metrics <- tryCatch(
              {
                vs_orig_tna <- tna::compare(sim_model, original_tna_model, scaling = scaling)
                if (!is.null(vs_orig_tna) && !is.null(vs_orig_tna$summary_metrics)) {
                  metrics_df <- vs_orig_tna$summary_metrics
                  metrics_df$model_type <- model_type
                  metrics_df$comparison_type <- "vs_original_tna"
                  metrics_df$run_idx <- run_idx
                  metrics_df$param_set_idx <- param_set_idx
                  metrics_df$data_idx <- data_idx # Add dataset index
                  metrics_df$num_states <- num_states # Add number of states

                  # Add parameters
                  for (param_name in names(param_set)) {
                    metrics_df[[param_name]] <- param_set[[param_name]]
                  }

                  # Calculate median sequence length and median NA per row
                  metrics_df$median_seq_length <- median(apply(sim_sequences, 1, function(row) sum(!is.na(row))))
                  metrics_df$median_na <- median(apply(sim_sequences, 1, function(row) sum(is.na(row))))

                  list(metrics_df) # Return as a list
                } else {
                  list() # Return empty list if no metrics
                }
              },
              error = function(e) {
                warning(paste("Comparison vs original_tna failed for model:", model_type, "run:", run_idx, "param_set:", param_set_idx, "dataset:", data_idx, "Error:", e$message))
                list() # Return empty list in case of error
              }
            )
            param_metrics_list <- c(param_metrics_list, vs_orig_tna_metrics)

            # If ftna comparison is also requested
            if ("ftna_reference" %in% comparisons && !is.null(sim_ftna_model)) {
              vs_orig_ftna_metrics <- tryCatch(
                {
                  vs_orig_ftna <- tna::compare(sim_ftna_model, original_ftna_model, scaling = scaling)

                  if (!is.null(vs_orig_ftna) && !is.null(vs_orig_ftna$summary_metrics)) {
                    metrics_df <- vs_orig_ftna$summary_metrics
                    metrics_df$model_type <- model_type
                    metrics_df$comparison_type <- "vs_original_ftna"
                    metrics_df$run_idx <- run_idx
                    metrics_df$param_set_idx <- param_set_idx
                    metrics_df$data_idx <- data_idx # Add dataset index
                    metrics_df$num_states <- num_states # Add number of states

                    # Add parameters
                    for (param_name in names(param_set)) {
                      metrics_df[[param_name]] <- param_set[[param_name]]
                    }

                    # Calculate median sequence length and median NA per row
                    metrics_df$median_seq_length <- median(apply(sim_sequences, 1, function(row) sum(!is.na(row))))
                    metrics_df$median_na <- median(apply(sim_sequences, 1, function(row) sum(is.na(row))))

                    list(metrics_df) # Return as a list
                  } else {
                    list() # Return empty list if no metrics
                  }
                },
                error = function(e) {
                  warning(paste("Comparison vs original_ftna failed for model:", model_type, "run:", run_idx, "param_set:", param_set_idx, "dataset:", data_idx, "Error:", e$message))
                  list() # Return empty list in case of error
                }
              )
              param_metrics_list <- c(param_metrics_list, vs_orig_ftna_metrics)
            }
          }

          # If model vs ftna comparison is requested
          if ("ftna_reference" %in% comparisons && model_type != "ftna" && !is.null(sim_ftna_model)) {
            vs_ftna_metrics <- tryCatch(
              {
                vs_ftna <- tna::compare(sim_model, sim_ftna_model, scaling = scaling)

                if (!is.null(vs_ftna) && !is.null(vs_ftna$summary_metrics)) {
                  metrics_df <- vs_ftna$summary_metrics
                  metrics_df$model_type <- model_type
                  metrics_df$comparison_type <- "vs_ftna"
                  metrics_df$run_idx <- run_idx
                  metrics_df$param_set_idx <- param_set_idx
                  metrics_df$data_idx <- data_idx # Add dataset index
                  metrics_df$num_states <- num_states # Add number of states

                  # Add parameters
                  for (param_name in names(param_set)) {
                    metrics_df[[param_name]] <- param_set[[param_name]]
                  }

                  # Calculate median sequence length and median NA per row
                  metrics_df$median_seq_length <- median(apply(sim_sequences, 1, function(row) sum(!is.na(row))))
                  metrics_df$median_na <- median(apply(sim_sequences, 1, function(row) sum(is.na(row))))

                  list(metrics_df) # Return as a list
                } else {
                  list() # Return empty list if no metrics
                }
              },
              error = function(e) {
                warning(paste("Comparison vs ftna failed for model:", model_type, "run:", run_idx, "param_set:", param_set_idx, "dataset:", data_idx, "Error:", e$message))
                list() # Return empty list in case of error
              }
            )
            param_metrics_list <- c(param_metrics_list, vs_ftna_metrics)
          }
        }
      }

      # Process across-model comparisons if requested
      if ("across_models" %in% comparisons && length(models) > 1) {
        # Generate sequences once for cross-model comparison, specific to this dataset
        cross_model_sequences <- simulate_sequences(
          trans_matrix = original_transition_probs,
          init_probs = initial_probabilities,
          seq_length = param_set$seq_length,
          n_sequences = param_set$n_sequences,
          na_range = c(param_set$min_na, param_set$max_na)
        )

        # Fit all models once for this dataset
        model_objects <- lapply(models, function(m) fit_network_model(cross_model_sequences, m))
        names(model_objects) <- models

        # Compare each pair of models
        if (length(models) >= 2) {
          for (i in 1:(length(models) - 1)) {
            for (j in (i + 1):length(models)) {
              model_i <- models[i]
              model_j <- models[j]

              cross_comp_metrics <- tryCatch(
                {
                  cross_comp <- tna::compare(model_objects[[model_i]], model_objects[[model_j]], scaling = scaling)

                  if (!is.null(cross_comp) && !is.null(cross_comp$summary_metrics)) {
                    metrics_df <- cross_comp$summary_metrics
                    metrics_df$model_type <- paste(model_i, "vs", model_j)
                    metrics_df$comparison_type <- "cross_model"
                    metrics_df$run_idx <- 1 # Only one run for cross-model
                    metrics_df$param_set_idx <- param_set_idx
                    metrics_df$data_idx <- data_idx # Add dataset index
                    metrics_df$num_states <- num_states # Add number of states


                    # Add parameters
                    for (param_name in names(param_set)) {
                      metrics_df[[param_name]] <- param_set[[param_name]]
                    }

                    # Calculate median sequence length and median NA per row
                    metrics_df$median_seq_length <- median(apply(cross_model_sequences, 1, function(row) sum(!is.na(row))))
                    metrics_df$median_na <- median(apply(cross_model_sequences, 1, function(row) sum(is.na(row))))

                    list(metrics_df) # Return as a list
                  } else {
                    list() # Return empty list if no metrics
                  }
                },
                error = function(e) {
                  warning(paste("Cross-model comparison failed for models:", model_i, "vs", model_j, "param_set:", param_set_idx, "dataset:", data_idx, "Error:", e$message))
                  list() # Return empty list in case of error
                }
              )
              param_metrics_list <- c(param_metrics_list, cross_comp_metrics)
            }
          }
        }
      }
    } # End of dataset loop
    return(param_metrics_list)
  }


  # --- Execute runs (in parallel if requested) ---
  if (parallel && requireNamespace("future.apply", quietly = TRUE)) {
    # With progress reporting
    progressr::with_progress({
      p <- progressr::progressor(steps = length(sim_params))

      all_metrics_lists <- future.apply::future_lapply(1:length(sim_params), function(i) {
        result <- run_single_param_set(i)
        p(sprintf("Parameter set %d complete", i))
        return(result)
      })
    })
  } else {
    # Sequential processing
    all_metrics_lists <- lapply(1:length(sim_params), run_single_param_set)
  }

  # Flatten list of lists
  metrics_list <- unlist(all_metrics_lists, recursive = FALSE)

  # Combine all metrics
  if (length(metrics_list) > 0) {
    metrics_df <- dplyr::bind_rows(metrics_list)

    # Clean up column names and structure
    metrics_df <- metrics_df %>%
      dplyr::rename(metric_category = category, metric_name = metric) %>%
      dplyr::select(
        model_type, comparison_type, metric_category, metric_name, value,
        data_idx, num_states, dplyr::everything()
      ) # Include data_idx and num_states
  } else {
    metrics_df <- data.frame()
  }

  # --- Generate summary statistics ---
  if (nrow(metrics_df) > 0) {
    summary_stats <- metrics_df %>%
      dplyr::group_by(model_type, comparison_type, metric_category, metric_name, data_idx, num_states) %>% # Group by dataset
      dplyr::summarise(
        mean_value = mean(value, na.rm = TRUE),
        sd_value = sd(value, na.rm = TRUE),
        median_value = median(value, na.rm = TRUE),
        min_value = min(value, na.rm = TRUE),
        max_value = max(value, na.rm = TRUE),
        n_obs = sum(!is.na(value)),
        .groups = "drop"
      )
  } else {
    summary_stats <- data.frame()
  }

  return(list(
    metrics = metrics_df,
    summary_stats = summary_stats,
    parameters = list(
      models = models,
      comparisons = comparisons,
      num_runs = num_runs,
      scaling = scaling
    )
  ))
}
