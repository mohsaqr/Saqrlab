#' @title Internal Helper Functions for Saqrlab
#' @name utils
#' @description Internal utility functions used by other Saqrlab functions.
#' @keywords internal
#' @importFrom utils tail head
NULL

# Global variable declarations to avoid R CMD check notes for dplyr NSE
# These are column names used in dplyr pipelines
utils::globalVariables(c(
  # Common column names
  ".", "from", "to", "value", "run_id", "setting_name",
  # Metrics
  "TP", "TN", "FP", "FN", "Sensitivity", "Specificity", "FPR", "FNR",
  "Accuracy", "MCC", "mcc_denom_sq", "Metric", "Value",
  # Totals
  "Total_TP", "Total_TN", "Total_FP", "Total_FN",
  # Bootstrap/simulation
  "p_value", "weight", "p_value_num", "weight_num",
  "bootstrap_significant_run", "is_significant",
  "n_significant", "recovery_rate", "avg_recovery_rate",
  "avg_p_value", "avg_weight",
  "ground_truth_stable", "gt_stable",
  # Parameters
  "num_rows", "max_seq_length", "min_na", "max_na", "num_states",
  "successful_runs",
  # Network simulation
  "category", "metric", "model_type", "comparison_type",
  "metric_category", "metric_name", "data_idx",
  # Data conversion
  "id", "Time", "Action", ".time_idx", ".time_name",
  # Summary functions
  ".n", ".sd", "network_id"
))

#' Safe Bind Rows
#'
#' Safely bind multiple data frames, handling errors and missing columns.
#'
#' @param df_list List of data frames to bind.
#' @param context Character string describing the context for warning messages.
#'
#' @return A single data frame with all rows bound together.
#'
#' @keywords internal
safe_bind_rows <- function(df_list, context = "") {
  if (length(df_list) == 0) return(data.frame())
  valid_dfs <- Filter(function(x) is.data.frame(x) && nrow(x) > 0 && ncol(x) > 0, df_list)
  if (length(valid_dfs) == 0) {
    warning(sprintf("No valid data frames to bind for %s", context))
    return(data.frame())
  }

  tryCatch({
    dplyr::bind_rows(valid_dfs)
  }, error = function(e) {
    warning(sprintf("Failed to bind rows for %s: %s", context, e$message))
    all_cols <- unique(unlist(lapply(valid_dfs, names)))
    aligned_dfs <- lapply(valid_dfs, function(df) {
      for (col in setdiff(all_cols, names(df))) df[[col]] <- NA
      df[, all_cols, drop = FALSE]
    })
    tryCatch({
      dplyr::bind_rows(aligned_dfs)
    }, error = function(e2) {
      warning(sprintf("Failed to bind rows for %s even after alignment attempt: %s", context, e2$message))
      return(data.frame())
    })
  })
}

#' Check Value in Range
#'
#' Check if a value falls within a specified range.
#'
#' @param value Numeric value to check.
#' @param range_val Numeric vector of length 2 with min and max, or NULL.
#'
#' @return Logical indicating whether value is in range.
#'
#' @keywords internal
check_val_in_range <- function(value, range_val) {
  if (is.null(range_val)) return(TRUE)
  if (is.na(value) || !is.numeric(value)) return(FALSE)
  return(value >= range_val[1] && value <= range_val[2])
}

#' Safe Median
#'
#' Calculate median with handling for empty vectors.
#'
#' @param x Numeric vector.
#'
#' @return Median value or NA if vector is empty.
#'
#' @keywords internal
safe_median <- function(x) {
  if (length(x) > 0) median(x, na.rm = TRUE) else NA_real_
}

#' Safe Mean
#'
#' Calculate mean with handling for empty vectors.
#'
#' @param x Numeric vector.
#'
#' @return Mean value or NA if vector is empty.
#'
#' @keywords internal
safe_mean <- function(x) {
  if (length(x) > 0) mean(x, na.rm = TRUE) else NA_real_
}

#' Safe Standard Deviation
#'
#' Calculate standard deviation with handling for single-value vectors.
#'
#' @param x Numeric vector.
#'
#' @return Standard deviation or NA if vector has fewer than 2 elements.
#'
#' @keywords internal
safe_sd <- function(x) {
  if (length(x) > 1) sd(x, na.rm = TRUE) else NA_real_
}
