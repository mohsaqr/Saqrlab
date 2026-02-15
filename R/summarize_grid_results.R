#' Summarize Grid Simulation Results
#'
#' @description
#' Summarize results from `run_grid_simulation()`, filtering by parameter ranges
#' and computing aggregated performance metrics. Provides detailed summaries
#' at both the setting level and across all selected settings.
#'
#' @param grid_results_list List output from `run_grid_simulation()`.
#' @param n_sequences_range Numeric vector of length 2 (min, max) or NULL.
#'   Filter settings by n_sequences.
#' @param seq_length_range Numeric vector of length 2 (min, max) or NULL.
#'   Filter settings by seq_length.
#' @param min_na_range Numeric vector of length 2 (min, max) or NULL.
#'   Filter settings by min_na.
#' @param max_na_range Numeric vector of length 2 (min, max) or NULL.
#'   Filter settings by max_na.
#'
#' @param num_rows_range Deprecated. Use `n_sequences_range` instead.
#' @param max_seq_length_range Deprecated. Use `seq_length_range` instead.
#' @param level_context Numeric. Significance level used for calculations
#'   (for context in output). Default: 0.05.
#' @param print_output Logical. Master switch for console printing. Default: TRUE.
#' @param print_aggregated_overall Logical. Print averaged overall performance
#'   metrics across selected settings. Default: TRUE.
#' @param print_aggregated_edges Logical. Print aggregated edge significance
#'   summary. Default: TRUE.
#' @param print_settings_summary Logical. Print summary table for each selected
#'   setting. Default: TRUE.
#'
#' @return A list containing:
#' \describe{
#'   \item{n_selected}{Number of settings matching the filter criteria.}
#'   \item{aggregated_summary}{List with:
#'     \itemize{
#'       \item overall_performance: Averaged metrics across settings.
#'       \item edge_significance: Aggregated edge-level statistics.
#'     }
#'   }
#'   \item{selected_settings_summary_df}{Data frame with per-setting metrics
#'     calculated from total TP/TN/FP/FN counts.}
#'   \item{compiled_individual_runs}{List with detailed run-level data:
#'     \itemize{
#'       \item all_raw_summaries: Combined bootstrap summaries.
#'       \item all_per_edge_performance: Combined per-edge results.
#'       \item run_level_performance_metrics: Metrics per run.
#'       \item setting_level_summary_stats: Mean/Median/SD of run metrics.
#'     }
#'   }
#' }
#' Returns NULL (invisibly) if no settings match.
#'
#' @details
#' The function performs comprehensive analysis:
#'
#' **Filtering**: Selects settings where all parameters fall within specified ranges.
#'
#' **Aggregation from Input Summaries**: Extracts and averages overall performance
#' and edge significance from the original `aggregated_summary` in each setting.
#'
#' **Run-Level Metrics**: Recomputes TP/TN/FP/FN counts from per-edge data,
#' then calculates Sensitivity, Specificity, FPR, FNR, Accuracy, and MCC per run.
#'
#' **Setting-Level Metrics**: Two approaches:
#' 1. Mean/Median/SD of run-level metrics.
#' 2. Metrics computed from total counts across all runs (more robust).
#'
#' The `selected_settings_summary_df` uses approach #2 for the most accurate
#' overall metrics per setting.
#'
#' @examples
#' \dontrun{
#' # After running grid simulation
#' grid_results <- run_grid_simulation(...)
#'
#' # Analyze all results
#' analysis <- summarize_grid_results(grid_results)
#'
#' # Filter to specific parameter ranges
#' analysis_filtered <- summarize_grid_results(
#'   grid_results,
#'   n_sequences_range = c(100, 200),
#'   seq_length_range = c(30, 50),
#'   max_na_range = c(0, 5)
#' )
#'
#' # Access the detailed run-level metrics
#' run_metrics <- analysis$compiled_individual_runs$run_level_performance_metrics
#'
#' # Access setting-level summary
#' setting_summary <- analysis$selected_settings_summary_df
#'
#' # Old parameter names still work
#' analysis_filtered <- summarize_grid_results(
#'   grid_results,
#'   num_rows_range = c(100, 200),
#'   max_seq_length_range = c(30, 50)
#' )
#' }
#'
#' @import dplyr
#' @importFrom tidyr complete
#' @importFrom stats sd
#' @export
summarize_grid_results <- function(grid_results_list,
                                  n_sequences_range = NULL,
                                  seq_length_range = NULL,
                                  min_na_range = NULL,
                                  max_na_range = NULL,
                                  level_context = 0.05,
                                  print_output = TRUE,
                                  print_aggregated_overall = TRUE,
                                  print_aggregated_edges = TRUE,
                                  print_settings_summary = TRUE,
                                  # Backward compatibility - old parameter names
                                  num_rows_range = NULL,
                                  max_seq_length_range = NULL) {
  # --- Backward compatibility: map old names to new names ---
  if (!is.null(num_rows_range) && is.null(n_sequences_range)) n_sequences_range <- num_rows_range
  if (!is.null(max_seq_length_range) && is.null(seq_length_range)) seq_length_range <- max_seq_length_range

  # --- Input Validation ---
  if (!is.list(grid_results_list) || length(grid_results_list) == 0) {
    stop("`grid_results_list` must be a non-empty list.")
  }

  validate_range <- function(range_arg, arg_name) {
    if (!is.null(range_arg)) {
      if (!is.numeric(range_arg) || length(range_arg) != 2 || range_arg[1] > range_arg[2]) {
        stop(sprintf("'%s' must be NULL or a numeric vector of length 2 with min <= max", arg_name))
      }
    }
    TRUE
  }

  validate_range(n_sequences_range, "n_sequences_range")
  validate_range(seq_length_range, "seq_length_range")
  validate_range(min_na_range, "min_na_range")
  validate_range(max_na_range, "max_na_range")

  if (!is.numeric(level_context) || length(level_context) != 1 || level_context <= 0 || level_context >= 1) {
    stop("'level_context' must be a single numeric value between 0 and 1")
  }

  logical_params <- list(
    print_output = print_output,
    print_aggregated_overall = print_aggregated_overall,
    print_aggregated_edges = print_aggregated_edges,
    print_settings_summary = print_settings_summary
  )

  for (param_name in names(logical_params)) {
    param_value <- logical_params[[param_name]]
    if (!is.logical(param_value) || length(param_value) != 1 || is.na(param_value)) {
      stop(sprintf("'%s' must be a single TRUE or FALSE value", param_name))
    }
  }

  # --- Filter the results list ---
  selected_settings <- list()
  original_names <- names(grid_results_list)
  if (is.null(original_names)) original_names <- paste0("setting_", seq_along(grid_results_list))

  for (i in seq_along(grid_results_list)) {
    name <- original_names[i]
    res <- grid_results_list[[i]]

    if (is.null(res)) {
      warning(sprintf("Element %d ('%s') is NULL - skipping", i, name))
      next
    }
    if (!is.list(res)) {
      warning(sprintf("Element %d ('%s') is not a list - skipping", i, name))
      next
    }
    required_components <- c("parameters", "aggregated_summary", "individual_runs", "successful_runs")
    missing_components <- setdiff(required_components, names(res))
    if (length(missing_components) > 0) {
      warning(sprintf(
        "Element %d ('%s') is missing required components: %s - skipping",
        i, name, paste(missing_components, collapse = ", ")
      ))
      next
    }
    # Check for required parameters (support both old and new names)
    params <- res$parameters
    has_n_sequences <- !is.null(params$n_sequences) || !is.null(params$num_rows)
    has_seq_length <- !is.null(params$seq_length) || !is.null(params$max_seq_length)
    has_min_na <- !is.null(params$min_na)
    has_max_na <- !is.null(params$max_na)

    if (!has_n_sequences || !has_seq_length || !has_min_na || !has_max_na) {
      warning(sprintf(
        "Element %d ('%s') parameters missing required values - skipping",
        i, name
      ))
      next
    }
    if (!is.list(res$aggregated_summary) ||
      !is.list(res$individual_runs) ||
      !("edge_significance" %in% names(res$aggregated_summary)) ||
      !("list_of_summaries" %in% names(res$individual_runs)) ||
      !("list_of_per_edge_performance" %in% names(res$individual_runs))) {
      warning(sprintf("Element %d ('%s') has invalid structure in aggregated_summary or individual_runs - skipping", i, name))
      next
    }

    # Get parameter values (support both old and new names)
    param_n_sequences <- if (!is.null(params$n_sequences)) params$n_sequences else params$num_rows
    param_seq_length <- if (!is.null(params$seq_length)) params$seq_length else params$max_seq_length

    if (check_val_in_range(param_n_sequences, n_sequences_range) &&
      check_val_in_range(param_seq_length, seq_length_range) &&
      check_val_in_range(params$min_na, min_na_range) &&
      check_val_in_range(params$max_na, max_na_range)) {
      selected_settings[[name]] <- res
    }
  }

  n_selected <- length(selected_settings)
  if (n_selected == 0) {
    if (print_output) message("No simulation settings matched the specified criteria.")
    return(invisible(NULL))
  }
  if (print_output) message(sprintf("Found %d matching simulation settings. Processing...", n_selected))

  # --- Prepare Aggregated Summary Data (Setting Level) from INPUTS ---
  all_overall_perf_list <- lapply(names(selected_settings), function(name) {
    res <- selected_settings[[name]]
    perf_df <- res$aggregated_summary$overall_performance
    if (!is.data.frame(perf_df) || nrow(perf_df) == 0) {
      return(NULL)
    }
    perf_df$setting_name <- name
    perf_df
  })
  all_overall_perf <- safe_bind_rows(all_overall_perf_list, "input overall performance")

  aggregated_overall_metrics_from_input <- tryCatch(
    {
      if (nrow(all_overall_perf) > 0 && all(c("Metric", "Value", "setting_name") %in% names(all_overall_perf))) {
        all_overall_perf %>%
          dplyr::mutate(Value = suppressWarnings(as.numeric(Value))) %>%
          dplyr::filter(!is.na(Value)) %>%
          dplyr::group_by(Metric) %>%
          dplyr::summarise(
            Average_Value = mean(Value, na.rm = TRUE),
            SD_Value = sd(Value, na.rm = TRUE),
            Min_Value = min(Value, na.rm = TRUE),
            Max_Value = max(Value, na.rm = TRUE),
            N_Settings = dplyr::n_distinct(setting_name),
            .groups = "drop"
          )
      } else {
        data.frame()
      }
    },
    error = function(e) {
      warning(sprintf("Failed to aggregate overall metrics from input: %s", e$message))
      data.frame()
    }
  )

  # Extract edge significance from input summaries
  extract_edge_df <- function(res, setting_name) {
    edge_df <- res$aggregated_summary$edge_significance
    if (!is.data.frame(edge_df) || nrow(edge_df) == 0) {
      return(NULL)
    }
    edge_df$setting_name <- setting_name
    edge_df
  }
  all_edge_summaries_list <- lapply(names(selected_settings), function(name) {
    extract_edge_df(selected_settings[[name]], name)
  })
  all_edge_summaries <- safe_bind_rows(Filter(Negate(is.null), all_edge_summaries_list), "edge significance")

  aggregated_edge_summary <- tryCatch(
    {
      if (nrow(all_edge_summaries) > 0 && all(c("from", "to", "ground_truth_stable", "recovery_rate", "n_significant", "avg_p_value", "avg_weight", "setting_name") %in% names(all_edge_summaries))) {
        all_edge_summaries %>%
          dplyr::filter(!is.na(recovery_rate)) %>%
          dplyr::group_by(from, to, ground_truth_stable) %>%
          dplyr::summarise(
            avg_recovery_rate = mean(recovery_rate, na.rm = TRUE),
            sd_recovery_rate = sd(recovery_rate, na.rm = TRUE),
            avg_n_significant = mean(n_significant, na.rm = TRUE),
            overall_avg_p_value = mean(avg_p_value, na.rm = TRUE),
            overall_avg_weight = mean(avg_weight, na.rm = TRUE),
            n_settings_included = dplyr::n_distinct(setting_name),
            .groups = "drop"
          ) %>%
          dplyr::arrange(dplyr::desc(ground_truth_stable), dplyr::desc(avg_recovery_rate))
      } else {
        data.frame()
      }
    },
    error = function(e) {
      warning(sprintf("Failed to aggregate edge significance data: %s", e$message))
      data.frame()
    }
  )

  aggregated_summary_output <- list(
    overall_performance = aggregated_overall_metrics_from_input,
    edge_significance = aggregated_edge_summary
  )

  # --- Collate Individual Run Details Safely ---
  collated_summaries_list <- list()
  collated_perf_list <- list()
  setting_params_list <- list()

  for (setting_name in names(selected_settings)) {
    tryCatch(
      {
        res <- selected_settings[[setting_name]]
        params <- res$parameters
        individual_runs <- res$individual_runs

        # Store parameters for this setting
        setting_params_list[[setting_name]] <- data.frame(
          setting_name = setting_name,
          num_rows = as.numeric(params$num_rows),
          max_seq_length = as.numeric(params$max_seq_length),
          min_na = as.numeric(params$min_na),
          max_na = as.numeric(params$max_na),
          successful_runs = as.numeric(res$successful_runs),
          stringsAsFactors = FALSE
        )

        if (is.null(individual_runs) || !is.list(individual_runs)) {
          next
        }

        list_of_summaries <- individual_runs$list_of_summaries
        list_of_per_edge <- individual_runs$list_of_per_edge_performance

        params_df_row <- data.frame(setting_name = setting_name)
        params_df_row$num_rows <- as.numeric(params$num_rows)

        # Process raw summaries
        if (is.list(list_of_summaries) && length(list_of_summaries) > 0) {
          summaries_with_id <- lapply(seq_along(list_of_summaries), function(id) {
            df <- list_of_summaries[[id]]
            if (is.data.frame(df) && nrow(df) > 0) {
              df$run_id <- id
              return(df)
            } else {
              return(NULL)
            }
          })
          setting_summaries_df <- safe_bind_rows(Filter(Negate(is.null), summaries_with_id), sprintf("raw summaries '%s'", setting_name))
          if (nrow(setting_summaries_df) > 0) {
            params_repped <- params_df_row[rep(1, nrow(setting_summaries_df)), , drop = FALSE]
            cols_to_remove <- intersect(names(setting_summaries_df), names(params_repped))
            if (length(cols_to_remove) > 0) params_repped <- params_repped[, !names(params_repped) %in% cols_to_remove, drop = FALSE]
            combined_df <- cbind(setting_summaries_df, params_repped)
            collated_summaries_list[[setting_name]] <- combined_df
          }
        }

        # Process per-edge performance
        if (is.list(list_of_per_edge) && length(list_of_per_edge) > 0) {
          perf_with_id <- lapply(seq_along(list_of_per_edge), function(id) {
            df <- list_of_per_edge[[id]]
            required_cols <- c("from", "to", "p_value", "ground_truth_stable")
            if (is.data.frame(df) && nrow(df) > 0 && all(required_cols %in% names(df))) {
              df$run_id <- id
              return(df)
            } else {
              return(NULL)
            }
          })
          setting_perf_df <- safe_bind_rows(Filter(Negate(is.null), perf_with_id), sprintf("per-edge perf '%s'", setting_name))
          if (nrow(setting_perf_df) > 0) {
            params_repped <- params_df_row[rep(1, nrow(setting_perf_df)), , drop = FALSE]
            cols_to_remove <- intersect(names(setting_perf_df), names(params_repped))
            if (length(cols_to_remove) > 0) params_repped <- params_repped[, !names(params_repped) %in% cols_to_remove, drop = FALSE]
            combined_df <- cbind(setting_perf_df, params_repped)
            collated_perf_list[[setting_name]] <- combined_df
          }
        } else {
          warning(sprintf("Setting '%s': No 'list_of_per_edge_performance' found or empty. Cannot calculate run-level/setting metrics.", setting_name))
        }
      },
      error = function(e) {
        warning(sprintf("Failed processing runs for setting '%s': %s", setting_name, e$message))
      }
    )
  }

  # Combine collated results
  final_collated_summaries <- safe_bind_rows(collated_summaries_list, "all collated summaries")
  final_collated_performance <- safe_bind_rows(collated_perf_list, "all collated per-edge performance")
  all_setting_params <- safe_bind_rows(setting_params_list, "all setting parameters")

  # --- Calculate Run-Level Performance Metrics ---
  run_level_performance_metrics <- tryCatch(
    {
      if (nrow(final_collated_performance) > 0 &&
        all(c("run_id", "setting_name", "p_value", "ground_truth_stable") %in% names(final_collated_performance))) {
        perf_data <- final_collated_performance %>%
          dplyr::mutate(
            ground_truth_stable = as.logical(ground_truth_stable),
            is_significant = !is.na(p_value) & p_value < level_context
          ) %>%
          dplyr::filter(!is.na(ground_truth_stable))

        run_counts <- perf_data %>%
          dplyr::group_by(setting_name, run_id, dplyr::across(dplyr::any_of(c("num_rows")))) %>%
          dplyr::summarise(
            TP = sum(ground_truth_stable & is_significant, na.rm = TRUE),
            FN = sum(ground_truth_stable & !is_significant, na.rm = TRUE),
            TN = sum(!ground_truth_stable & !is_significant, na.rm = TRUE),
            FP = sum(!ground_truth_stable & is_significant, na.rm = TRUE),
            .groups = "drop"
          )

        run_metrics <- run_counts %>%
          dplyr::mutate(
            Sensitivity = TP / (TP + FN),
            Specificity = TN / (TN + FP),
            FPR = FP / (TN + FP),
            FNR = FN / (TP + FN),
            Accuracy = (TP + TN) / (TP + TN + FP + FN),
            mcc_denom_sq = as.numeric(TP + FP) * as.numeric(TP + FN) * as.numeric(TN + FP) * as.numeric(TN + FN),
            MCC = dplyr::if_else(mcc_denom_sq > 0,
              (as.numeric(TP) * as.numeric(TN) - as.numeric(FP) * as.numeric(FN)) / sqrt(mcc_denom_sq),
              0
            ),
            dplyr::across(c(Sensitivity, Specificity, FPR, FNR, Accuracy, MCC), ~ dplyr::if_else(is.nan(.x), NA_real_, .x))
          ) %>%
          dplyr::select(
            setting_name, run_id,
            dplyr::any_of(c("num_rows")),
            Sensitivity, Specificity, FPR, FNR, Accuracy, MCC,
            TP, FN, TN, FP
          )
        run_metrics
      } else {
        warning("Could not calculate run-level performance metrics. Input data missing/invalid.")
        data.frame(
          setting_name = character(), run_id = integer(), Sensitivity = numeric(), Specificity = numeric(),
          FPR = numeric(), FNR = numeric(), Accuracy = numeric(), MCC = numeric(),
          TP = integer(), FN = integer(), TN = integer(), FP = integer()
        )
      }
    },
    error = function(e) {
      warning(sprintf("Failed run-level metrics calculation: %s", e$message))
      data.frame(
        setting_name = character(), run_id = integer(), Sensitivity = numeric(), Specificity = numeric(),
        FPR = numeric(), FNR = numeric(), Accuracy = numeric(), MCC = numeric(),
        TP = integer(), FN = integer(), TN = integer(), FP = integer()
      )
    }
  )

  # --- Calculate Setting-Level Aggregated Stats (Mean/Median/SD of Run Metrics) ---
  setting_level_summary_stats <- tryCatch(
    {
      if (nrow(run_level_performance_metrics) > 0) {
        metrics_to_summarize <- c("Sensitivity", "Specificity", "FPR", "FNR", "Accuracy", "MCC", "TP", "FN", "TN", "FP")
        cols_present <- intersect(metrics_to_summarize, names(run_level_performance_metrics))
        if (length(cols_present) == 0) stop("No metric columns found in run_level_performance_metrics")

        summary_data <- run_level_performance_metrics %>%
          dplyr::select(setting_name, dplyr::all_of(cols_present)) %>%
          dplyr::group_by(setting_name)

        summary_median <- summary_data %>% dplyr::summarise(dplyr::across(dplyr::all_of(cols_present), safe_median, .names = "Median_{.col}"), .groups = "drop")
        summary_mean <- summary_data %>% dplyr::summarise(dplyr::across(dplyr::all_of(cols_present), safe_mean, .names = "Mean_{.col}"), .groups = "drop")
        summary_sd <- summary_data %>% dplyr::summarise(dplyr::across(dplyr::all_of(cols_present), safe_sd, .names = "SD_{.col}"), .groups = "drop")
        summary_n <- summary_data %>% dplyr::summarise(N_Runs = dplyr::n(), .groups = "drop")

        setting_summary <- summary_n %>%
          dplyr::left_join(summary_median, by = "setting_name") %>%
          dplyr::left_join(summary_mean, by = "setting_name") %>%
          dplyr::left_join(summary_sd, by = "setting_name")

        final_setting_summary <- all_setting_params %>%
          dplyr::select(setting_name, num_rows, max_seq_length, min_na, max_na) %>%
          dplyr::distinct() %>%
          dplyr::left_join(setting_summary, by = "setting_name")

        final_setting_summary
      } else {
        data.frame()
      }
    },
    error = function(e) {
      warning(sprintf("Failed setting-level aggregation of run metrics: %s", e$message))
      data.frame()
    }
  )

  # --- Calculate Setting-Level Overall Metrics from TOTAL Counts ---
  setting_total_counts_and_metrics <- tryCatch(
    {
      if (nrow(run_level_performance_metrics) > 0 && all(c("setting_name", "TP", "FN", "TN", "FP") %in% names(run_level_performance_metrics))) {
        total_counts <- run_level_performance_metrics %>%
          dplyr::group_by(setting_name) %>%
          dplyr::summarise(
            Total_TP = sum(TP, na.rm = TRUE),
            Total_FN = sum(FN, na.rm = TRUE),
            Total_TN = sum(TN, na.rm = TRUE),
            Total_FP = sum(FP, na.rm = TRUE),
            .groups = "drop"
          )

        metrics_from_totals <- total_counts %>%
          dplyr::mutate(
            Sensitivity = Total_TP / (Total_TP + Total_FN),
            Specificity = Total_TN / (Total_TN + Total_FP),
            FPR = Total_FP / (Total_TN + Total_FP),
            FNR = Total_FN / (Total_TP + Total_FN),
            Accuracy = (Total_TP + Total_TN) / (Total_TP + Total_TN + Total_FP + Total_FN),
            mcc_denom_sq = as.numeric(Total_TP + Total_FP) * as.numeric(Total_TP + Total_FN) * as.numeric(Total_TN + Total_FP) * as.numeric(Total_TN + Total_FN),
            MCC = dplyr::if_else(mcc_denom_sq > 0,
              (as.numeric(Total_TP) * as.numeric(Total_TN) - as.numeric(Total_FP) * as.numeric(Total_FN)) / sqrt(mcc_denom_sq),
              0
            ),
            dplyr::across(c(Sensitivity, Specificity, FPR, FNR, Accuracy, MCC), ~ dplyr::if_else(is.nan(.x), NA_real_, .x))
          ) %>%
          dplyr::select(setting_name, Sensitivity, Specificity, FPR, FNR, Accuracy, MCC)

        metrics_from_totals
      } else {
        warning("Cannot calculate setting metrics from total counts - run level metrics missing or incomplete.")
        data.frame(
          setting_name = character(), Sensitivity = numeric(), Specificity = numeric(),
          FPR = numeric(), FNR = numeric(), Accuracy = numeric(), MCC = numeric()
        )
      }
    },
    error = function(e) {
      warning(sprintf("Failed calculating setting metrics from total counts: %s", e$message))
      data.frame(
        setting_name = character(), Sensitivity = numeric(), Specificity = numeric(),
        FPR = numeric(), FNR = numeric(), Accuracy = numeric(), MCC = numeric()
      )
    }
  )

  # --- Build Selected Settings Summary Data Frame ---
  selected_settings_summary_df <- tryCatch(
    {
      if (nrow(all_setting_params) > 0 && nrow(setting_total_counts_and_metrics) > 0) {
        all_setting_params %>%
          dplyr::left_join(setting_total_counts_and_metrics, by = "setting_name") %>%
          dplyr::select(
            setting_name, num_rows, max_seq_length, min_na, max_na, successful_runs,
            Sensitivity, Specificity, FPR, FNR, Accuracy, MCC
          )
      } else {
        warning("Cannot create final settings summary table - parameter or calculated metrics data missing.")
        if (nrow(all_setting_params) > 0) {
          all_setting_params %>%
            dplyr::mutate(Sensitivity = NA_real_, Specificity = NA_real_, FPR = NA_real_, FNR = NA_real_, Accuracy = NA_real_, MCC = NA_real_) %>%
            dplyr::select(setting_name, num_rows, max_seq_length, min_na, max_na, successful_runs, Sensitivity, Specificity, FPR, FNR, Accuracy, MCC)
        } else {
          data.frame(
            setting_name = character(), num_rows = numeric(), max_seq_length = numeric(),
            min_na = numeric(), max_na = numeric(), successful_runs = numeric(),
            Sensitivity = numeric(), Specificity = numeric(), FPR = numeric(),
            FNR = numeric(), Accuracy = numeric(), MCC = numeric()
          )
        }
      }
    },
    error = function(e) {
      warning(sprintf("Failed creating final selected_settings_summary_df: %s", e$message))
      data.frame(
        setting_name = character(), num_rows = numeric(), max_seq_length = numeric(),
        min_na = numeric(), max_na = numeric(), successful_runs = numeric(),
        Sensitivity = numeric(), Specificity = numeric(), FPR = numeric(),
        FNR = numeric(), Accuracy = numeric(), MCC = numeric()
      )
    }
  )

  # --- Package Compiled Data ---
  compiled_data_output <- list(
    all_raw_summaries = final_collated_summaries,
    all_per_edge_performance = final_collated_performance,
    run_level_performance_metrics = run_level_performance_metrics,
    setting_level_summary_stats = setting_level_summary_stats
  )

  # --- Conditional Printing ---
  if (print_output) {
    cat("\n------------------------------------------------------------\n")
    cat(sprintf("Analysis Summary for %d Matching Simulation Settings\n", n_selected))
    cat("------------------------------------------------------------\n")
    cat("Applied Filters:\n")
    print_filter <- function(name, range_val) {
      if (!is.null(range_val)) cat(sprintf("  %s: [%s, %s]\n", name, range_val[1], range_val[2]))
    }
    print_filter("n_sequences", n_sequences_range)
    print_filter("seq_length", seq_length_range)
    print_filter("min_na", min_na_range)
    print_filter("max_na", max_na_range)
    cat(sprintf("  Significance Level (alpha) for Run Metrics: %.3f\n", level_context))
    cat("------------------------------------------------------------\n")

    if (print_aggregated_overall && nrow(aggregated_summary_output$overall_performance) > 0) {
      cat("\n--- Average Overall Performance Across Selected Settings (from ORIGINAL input summaries) ---\n")
      print(aggregated_summary_output$overall_performance, row.names = FALSE)
    } else if (print_aggregated_overall) {
      cat("\n--- Average Overall Performance (from ORIGINAL input summaries): Not Available or Empty ---\n")
    }

    if (print_aggregated_edges && nrow(aggregated_summary_output$edge_significance) > 0) {
      cat("\n--- Aggregated Edge Significance/Recovery Across Selected Settings (from input summaries) ---\n")
      print(aggregated_summary_output$edge_significance, row.names = FALSE)
    } else if (print_aggregated_edges) {
      cat("\n--- Aggregated Edge Significance/Recovery: Not Available ---\n")
    }

    if (print_settings_summary && nrow(selected_settings_summary_df) > 0) {
      cat("\n--- Overall Metrics for Each Selected Setting (Calculated from Run Totals) ---\n")
      print(selected_settings_summary_df, row.names = FALSE)
    } else if (print_settings_summary) {
      cat("\n--- Overall Metrics for Each Selected Setting: Not Available (Calculation Failed - check warnings) ---\n")
    }

    if (nrow(compiled_data_output$run_level_performance_metrics) > 0) {
      cat("\n--- Calculated Run-Level Performance Metrics (based on p < alpha) ---\n")
      cat(sprintf("  Data frame 'compiled_individual_runs$run_level_performance_metrics' created.\n"))
      cat(sprintf("  Dimensions: %d runs x %d columns\n", nrow(compiled_data_output$run_level_performance_metrics), ncol(compiled_data_output$run_level_performance_metrics)))
    } else {
      cat("\n--- Calculated Run-Level Performance Metrics: Not Available (check warnings) ---\n")
    }

    if (nrow(compiled_data_output$setting_level_summary_stats) > 0) {
      cat("\n--- Aggregated Setting-Level Statistics (Mean/Median/SD of run metrics) ---\n")
      cat(sprintf("  Data frame 'compiled_individual_runs$setting_level_summary_stats' created.\n"))
      cat(sprintf("  Dimensions: %d settings x %d columns\n", nrow(compiled_data_output$setting_level_summary_stats), ncol(compiled_data_output$setting_level_summary_stats)))
    } else {
      cat("\n--- Aggregated Setting-Level Statistics: Not Available (check warnings) ---\n")
    }
    cat("------------------------------------------------------------\n\n")
  }

  # --- Return Structured List ---
  return(invisible(list(
    n_selected = n_selected,
    aggregated_summary = aggregated_summary_output,
    selected_settings_summary_df = selected_settings_summary_df,
    compiled_individual_runs = compiled_data_output
  )))
}


# Backward compatibility alias (silent - no warning)

#' @rdname summarize_grid_results
#' @param ... Arguments passed to \code{summarize_grid_results}.
#' @export
analyze_grid_results <- function(...) summarize_grid_results(...)
