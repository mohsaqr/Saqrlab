#' Summary Functions
#'
#' @description
#' Functions for summarizing simulation results and network statistics.
#'
#' @name summary_functions
#' @keywords internal
NULL

#' Summarize Simulation Results
#'
#' @description
#' Compute summary statistics for simulation results, optionally grouped
#' by parameter combinations or other factors.
#'
#' @param results A data frame of simulation results, typically from
#'   `run_network_simulation()` or `run_bootstrap_simulation()`.
#' @param by Character vector. Column names to group by when computing
#'   summaries. Default: NULL (compute overall summary).
#' @param metrics Character vector. Which summary metrics to compute.
#'   Options: "mean" (arithmetic mean), "sd" (standard deviation),
#'   "median" (median value), "ci" (95 percent confidence interval),
#'   "min" (minimum), "max" (maximum), "n" (count), or "all" to compute
#'   all metrics. Default: c("mean", "sd", "ci").
#' @param value_cols Character vector. Names of columns containing numeric
#'   values to summarize. If NULL, auto-detects numeric columns.
#'   Default: NULL.
#' @param na.rm Logical. Whether to remove NA values when computing
#'   statistics. Default: TRUE.
#'
#' @return A data frame with summary statistics. If `by` is specified,
#'   contains one row per unique combination of grouping variables.
#'
#' @details
#' This function provides flexible summarization of simulation results.
#' It automatically detects numeric columns and computes requested statistics.
#'
#' Confidence intervals are computed as mean +/- 1.96 * SE, where
#' SE = sd / sqrt(n).
#'
#' @examples
#' \dontrun{
#' # Run a simulation
#' results <- run_network_simulation(trans_mat, init_probs,
#'                                   num_rows_values = c(50, 100, 200),
#'                                   n_runs = 20)
#'
#' # Overall summary
#' summarize_simulation(results)
#'
#' # Summary by sample size
#' summarize_simulation(results, by = "num_rows")
#'
#' # Summary with all metrics
#' summarize_simulation(results, by = "num_rows", metrics = "all")
#' }
#'
#' @seealso [summarize_networks()] for summarizing network metrics,
#'   [run_network_simulation()] for running simulations.
#'
#' @importFrom stats qt sd median
#' @import dplyr
#' @export
summarize_simulation <- function(results,
                                 by = NULL,
                                 metrics = c("mean", "sd", "ci"),
                                 value_cols = NULL,
                                 na.rm = TRUE) {
  # Input validation
  if (!is.data.frame(results)) {
    stop("results must be a data frame.")
  }
  if (nrow(results) == 0) {
    warning("results data frame is empty.")
    return(data.frame())
  }

  if ("all" %in% metrics) {
    metrics <- c("mean", "sd", "median", "ci", "min", "max", "n")
  }

  # Auto-detect numeric columns
  if (is.null(value_cols)) {
    value_cols <- names(results)[sapply(results, is.numeric)]
    # Exclude common non-metric columns
    exclude_cols <- c("run_id", "setting_id", "run", "rep", "replicate", "seed")
    value_cols <- setdiff(value_cols, exclude_cols)
  }

  if (length(value_cols) == 0) {
    warning("No numeric columns found to summarize.")
    return(data.frame())
  }

  # Define summary functions
  summary_funs <- list()

  if ("mean" %in% metrics) {
    summary_funs$mean <- function(x) mean(x, na.rm = na.rm)
  }
  if ("sd" %in% metrics) {
    summary_funs$sd <- function(x) sd(x, na.rm = na.rm)
  }
  if ("median" %in% metrics) {
    summary_funs$median <- function(x) median(x, na.rm = na.rm)
  }
  if ("min" %in% metrics) {
    summary_funs$min <- function(x) min(x, na.rm = na.rm)
  }
  if ("max" %in% metrics) {
    summary_funs$max <- function(x) max(x, na.rm = na.rm)
  }
  if ("n" %in% metrics) {
    summary_funs$n <- function(x) sum(!is.na(x))
  }

  # Group if needed
  if (!is.null(by)) {
    # Check that grouping columns exist
    missing_by <- setdiff(by, names(results))
    if (length(missing_by) > 0) {
      stop("Grouping columns not found: ", paste(missing_by, collapse = ", "))
    }
    grouped_data <- dplyr::group_by(results, dplyr::across(dplyr::all_of(by)))
  } else {
    grouped_data <- results
  }

  # Compute summaries for each value column
  summary_list <- list()

  for (col in value_cols) {
    col_summaries <- dplyr::summarise(
      grouped_data,
      dplyr::across(
        dplyr::all_of(col),
        summary_funs,
        .names = "{.fn}"
      ),
      .groups = "drop"
    )

    # Add CI if requested
    if ("ci" %in% metrics && "mean" %in% names(col_summaries)) {
      col_data <- grouped_data[[col]]
      if (is.null(by)) {
        n_obs <- sum(!is.na(col_data))
        se <- sd(col_data, na.rm = na.rm) / sqrt(n_obs)
        col_summaries$ci_lower <- col_summaries$mean - 1.96 * se
        col_summaries$ci_upper <- col_summaries$mean + 1.96 * se
      } else {
        # Compute CI for each group
        ci_data <- dplyr::summarise(
          grouped_data,
          .n = sum(!is.na(.data[[col]])),
          .sd = sd(.data[[col]], na.rm = na.rm),
          .groups = "drop"
        )
        se <- ci_data$.sd / sqrt(ci_data$.n)
        col_summaries$ci_lower <- col_summaries$mean - 1.96 * se
        col_summaries$ci_upper <- col_summaries$mean + 1.96 * se
      }
    }

    # Rename columns to include original column name
    summary_cols <- setdiff(names(col_summaries), by)
    names(col_summaries)[names(col_summaries) %in% summary_cols] <-
      paste0(col, "_", summary_cols)

    summary_list[[col]] <- col_summaries
  }

  # Combine all summaries
  if (length(summary_list) == 1) {
    result <- summary_list[[1]]
  } else {
    result <- summary_list[[1]]
    for (i in 2:length(summary_list)) {
      if (!is.null(by)) {
        result <- dplyr::left_join(
          result,
          summary_list[[i]],
          by = by
        )
      } else {
        result <- cbind(result, summary_list[[i]])
      }
    }
  }

  # Convert to data frame
  result <- as.data.frame(result)
  return(result)
}


#' Summarize Multiple Networks
#'
#' @description
#' Compute aggregate statistics across multiple TNA network models.
#'
#' @param model_list A list of TNA model objects.
#' @param include Character vector. What to include in the summary.
#'   Options: "density" (network density), "centrality" (centrality measures),
#'   "edges" (edge weight statistics), or "all" for everything.
#'   Default: c("density", "centrality", "edges").
#' @param threshold Numeric. Minimum edge weight to consider present for
#'   density calculation. Default: 0.01.
#' @param centrality_measures Character vector. Which centrality measures
#'   to summarize. Default: c("OutStrength", "InStrength").
#'
#' @return A list containing:
#' \itemize{
#'   \item summary_table: Data frame with per-network statistics.
#'   \item aggregate: Named list of aggregate statistics across all networks.
#'   \item n_networks: Number of networks summarized.
#' }
#'
#' @details
#' This function is useful for summarizing results from simulation studies
#' where many networks are generated. It provides both per-network and
#' aggregate statistics.
#'
#' @examples
#' \dontrun{
#' # Generate and fit multiple networks
#' datasets <- lapply(1:20, function(i) {
#'   simulate_sequences(trans_mat, init_probs, 20, 100)
#' })
#' models <- batch_fit_models(datasets)
#'
#' # Summarize networks
#' network_summary <- summarize_networks(models)
#' print(network_summary$aggregate)
#' }
#'
#' @seealso [summarize_simulation()] for simulation result summaries,
#'   [batch_fit_models()] for fitting multiple models.
#'
#' @importFrom stats sd median
#' @export
summarize_networks <- function(model_list,
                               include = c("density", "centrality", "edges"),
                               threshold = 0.01,
                               centrality_measures = c("OutStrength", "InStrength")) {
  # Input validation
  if (!is.list(model_list)) {
    stop("model_list must be a list of TNA model objects.")
  }

  # Filter out NULL models
  valid_models <- Filter(Negate(is.null), model_list)
  n_models <- length(valid_models)

  if (n_models == 0) {
    warning("No valid models to summarize.")
    return(list(
      summary_table = data.frame(),
      aggregate = list(),
      n_networks = 0
    ))
  }

  if ("all" %in% include) {
    include <- c("density", "centrality", "edges")
  }

  # Initialize results storage
  per_network <- data.frame(
    network_id = seq_len(n_models),
    stringsAsFactors = FALSE
  )

  # Compute density for each network
  if ("density" %in% include) {
    densities <- sapply(valid_models, function(model) {
      tryCatch({
        weights <- extract_transition_matrix(model)
        n_edges <- nrow(weights) * ncol(weights)
        n_present <- sum(weights >= threshold)
        n_present / n_edges
      }, error = function(e) NA_real_)
    })
    per_network$density <- densities
  }

  # Compute edge statistics
  if ("edges" %in% include) {
    edge_stats <- t(sapply(valid_models, function(model) {
      tryCatch({
        weights <- extract_transition_matrix(model)
        c(
          edge_mean = mean(weights, na.rm = TRUE),
          edge_sd = sd(as.vector(weights), na.rm = TRUE),
          edge_max = max(weights, na.rm = TRUE),
          n_strong_edges = sum(weights >= 0.1)
        )
      }, error = function(e) c(
        edge_mean = NA_real_,
        edge_sd = NA_real_,
        edge_max = NA_real_,
        n_strong_edges = NA_integer_
      ))
    }))
    per_network <- cbind(per_network, as.data.frame(edge_stats))
  }

  # Compute centrality summaries
  if ("centrality" %in% include) {
    centrality_summaries <- lapply(valid_models, function(model) {
      tryCatch({
        cent <- tna::centralities(model)
        result <- list()
        for (measure in centrality_measures) {
          if (measure %in% colnames(cent)) {
            result[[paste0(measure, "_mean")]] <- mean(cent[, measure], na.rm = TRUE)
            result[[paste0(measure, "_sd")]] <- sd(cent[, measure], na.rm = TRUE)
          }
        }
        as.data.frame(result)
      }, error = function(e) {
        df <- data.frame(row.names = 1)
        for (measure in centrality_measures) {
          df[[paste0(measure, "_mean")]] <- NA_real_
          df[[paste0(measure, "_sd")]] <- NA_real_
        }
        df
      })
    })
    centrality_df <- do.call(rbind, centrality_summaries)
    per_network <- cbind(per_network, centrality_df)
  }

  # Compute aggregate statistics
  numeric_cols <- names(per_network)[sapply(per_network, is.numeric)]
  numeric_cols <- setdiff(numeric_cols, "network_id")

  aggregate_stats <- list()
  for (col in numeric_cols) {
    vals <- per_network[[col]]
    aggregate_stats[[col]] <- list(
      mean = mean(vals, na.rm = TRUE),
      sd = sd(vals, na.rm = TRUE),
      median = median(vals, na.rm = TRUE),
      min = min(vals, na.rm = TRUE),
      max = max(vals, na.rm = TRUE)
    )
  }

  # Create result
  result <- list(
    summary_table = per_network,
    aggregate = aggregate_stats,
    n_networks = n_models
  )
  class(result) <- c("network_summary", "list")

  return(result)
}
