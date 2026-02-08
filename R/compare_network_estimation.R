#' Compare Network Estimation Across TNA Model Types
#'
#' @description
#' Perform side-by-side comparison of multiple TNA model types (tna, atna, ctna, ftna)
#' using comprehensive sampling analysis. Uses a statistically correct approach:
#' comparing models built on independent data subsets (sample vs remaining).
#'
#' @param data A data frame containing sequence data.
#' @param model_types Character vector of model types to compare. Valid types are
#'   "tna", "atna", "ctna", "ftna". Default: c("tna", "ftna").
#' @param model_scaling A named list specifying scaling for each model type.
#'   Example: list(tna = "minmax", ftna = "skip"). If NULL, uses default scaling.
#' @param sampling_percent Numeric value between 0 and 1 indicating the
#'   proportion of data to sample. Default: 0.3.
#' @param iterations Integer specifying the number of sampling iterations
#'   for each model. Default: 100.
#' @param seed Integer seed for reproducible random sampling. Default: NULL.
#' @param verbose Logical. Print progress messages. Default: TRUE.
#'
#' @return A list containing:
#' \describe{
#'   \item{individual}{Data frame with all iteration results including model_type column.}
#'   \item{aggregated}{Data frame with aggregated metrics by model type.}
#'   \item{ranking}{Models ranked by mean Pearson correlation.}
#'   \item{winner}{Model with highest Pearson correlation.}
#'   \item{params}{Parameters used for the comparison.}
#' }
#'
#' @details
#' The function workflow:
#' 1. For each model type, builds a model from the data.
#' 2. Runs sampling analysis comparing sample vs remaining subsets.
#' 3. Combines all results for comparison.
#'
#' This approach is statistically correct because it compares models built on
#' independent data subsets, providing a true measure of estimation stability.
#'
#' @examples
#' \dontrun{
#' library(tna)
#' data(group_regulation)
#'
#' # Compare tna and ftna models
#' results <- compare_network_estimation(group_regulation)
#' results$ranking
#' results$winner
#'
#' # Compare all model types
#' results <- compare_network_estimation(
#'   group_regulation,
#'   model_types = c("tna", "atna", "ctna", "ftna"),
#'   iterations = 100
#' )
#'
#' # Plot results
#' plot(results)  # Default histogram
#' plot(results, metric_name = "Pearson")
#' }
#'
#' @import dplyr
#' @import tna
#' @export
compare_network_estimation <- function(data,
                                        model_types = c("tna", "ftna"),
                                        model_scaling = NULL,
                                        sampling_percent = 0.3,
                                        iterations = 100,
                                        seed = NULL,
                                        verbose = TRUE) {
  # Input validation
  if (!is.data.frame(data)) {
    stop("data must be a data frame")
  }

  valid_types <- c("tna", "atna", "ctna", "ftna")
  if (!all(model_types %in% valid_types)) {
    invalid <- setdiff(model_types, valid_types)
    stop(paste("Invalid model type(s):", paste(invalid, collapse = ", "),
               "\nValid types:", paste(valid_types, collapse = ", ")))
  }

  if (length(model_types) < 2) {
    stop("At least 2 model types required for comparison")
  }

  if (!is.null(seed)) set.seed(seed)

  if (verbose) {
    cat("=== Network Estimation Comparison ===\n")
    cat("Models:", paste(toupper(model_types), collapse = " vs "), "\n")
    cat("Iterations:", iterations, "per model\n")
    cat("Sample size:", paste0(sampling_percent * 100, "%"), "\n\n")
  }

  # Model type mapping
  model_map <- list(
    "tna" = "relative",
    "ftna" = "frequency",
    "ctna" = "co-occurrence"
  )

  all_individual <- list()
  all_aggregated <- list()

  for (model_type in model_types) {
    if (verbose) cat("Building", toupper(model_type), "model...\n")

    # Determine scaling for current model
    current_scaling <- if (!is.null(model_scaling) && model_type %in% names(model_scaling)) {
      model_scaling[[model_type]]
    } else {
      NULL
    }

    # Handle "none" string
    scaling_arg <- if (!is.null(current_scaling) && current_scaling == "none") {
      NULL
    } else {
      current_scaling
    }

    # Build the model
    tryCatch({
      if (model_type == "atna") {
        model <- tna::atna(data)
      } else {
        type_arg <- model_map[[model_type]]
        if (!is.null(scaling_arg)) {
          model <- tna::build_model(data, type = type_arg, scaling = scaling_arg)
        } else {
          model <- tna::build_model(data, type = type_arg)
        }
      }

      if (verbose) cat("  Running", iterations, "sampling iterations...\n")

      # Run sampling analysis
      results <- run_sampling_analysis(
        model = model,
        sampling_percent = sampling_percent,
        iterations = iterations,
        seed = NULL,  # Don't reset seed for each model
        model_scaling = current_scaling,
        verbose = FALSE
      )

      # Add model_type to results
      individual <- results$individual %>%
        dplyr::mutate(model_type = model_type)

      aggregated <- results$aggregated %>%
        dplyr::mutate(model_type = model_type)

      all_individual[[model_type]] <- individual
      all_aggregated[[model_type]] <- aggregated

      if (verbose) {
        pearson_mean <- aggregated$mean[aggregated$metric == "Pearson"]
        cat("  ", toupper(model_type), "complete. Pearson r =", round(pearson_mean, 3), "\n\n")
      }

    }, error = function(e) {
      warning("Error processing ", model_type, ": ", e$message)
    })
  }

  if (length(all_individual) == 0) {
    stop("All model types failed to process")
  }

  # Combine results
  combined_individual <- dplyr::bind_rows(all_individual)
  combined_aggregated <- dplyr::bind_rows(all_aggregated)

  # Create ranking based on Pearson correlation
  ranking_df <- combined_aggregated %>%
    dplyr::filter(metric == "Pearson") %>%
    dplyr::arrange(dplyr::desc(mean)) %>%
    dplyr::select(model_type, mean, sd)

  if (verbose) {
    cat("=== RESULTS ===\n")
    cat("Ranking by Pearson correlation:\n")
    for (i in seq_len(nrow(ranking_df))) {
      cat(sprintf("  %d. %s: r = %.4f (sd = %.4f)\n",
                  i,
                  toupper(ranking_df$model_type[i]),
                  ranking_df$mean[i],
                  ranking_df$sd[i]))
    }
    cat("\nWinner:", toupper(ranking_df$model_type[1]), "\n")
  }

  result <- list(
    individual = combined_individual,
    aggregated = combined_aggregated,
    ranking = ranking_df$model_type,
    winner = ranking_df$model_type[1],
    params = list(
      model_types = model_types,
      sampling_percent = sampling_percent,
      iterations = iterations,
      seed = seed,
      n_sequences = nrow(data)
    )
  )

  class(result) <- c("network_estimation", "list")
  result
}


#' Plot Network Estimation Comparison Results
#'
#' @description
#' Plot method for network_estimation objects. Default is histogram with
#' density overlay showing mean and median values.
#'
#' @param x A network_estimation object from compare_network_estimation().
#' @param metric_name Character. Metric to plot. Default: "Pearson".
#' @param plot_type Character. Plot type: "histogram" (default), "boxplot",
#'   "density", "bar", or "ridgeline".
#' @param ... Additional arguments (currently unused).
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' results <- compare_network_estimation(group_regulation)
#' plot(results)  # Histogram of Pearson correlation
#' plot(results, metric_name = "Spearman")
#' plot(results, plot_type = "boxplot")
#' }
#'
#' @import ggplot2
#' @import dplyr
#' @export
plot.network_estimation <- function(x, metric_name = "Pearson", plot_type = "histogram", ...) {
  plot_network_estimation(x, metric_name = metric_name, plot_type = plot_type)
}


#' Plot Network Estimation Results
#'
#' @description
#' Create publication-ready plots for network estimation comparison results.
#' Default is a histogram optimized for high iteration counts.
#'
#' @param results A network_estimation object or data frame with individual results.
#' @param metric_name Character. Metric to plot. Default: "Pearson".
#' @param plot_type Character. One of "histogram" (default), "boxplot", "density",
#'   "bar", or "ridgeline".
#' @param show_stats Logical. Show mean/median statistics. Default: TRUE.
#' @param bins Integer. Number of bins for histogram. Default: auto-calculated
#'   based on iterations (uses Sturges' rule).
#'
#' @return A ggplot object.
#'
#' @details
#' The histogram plot is optimized for high iteration counts (100+) and shows:
#' - Turquoise bars with density overlay
#' - Green vertical line for mean
#' - Orange dotted line for median
#' - Subtitle with exact mean and median values
#'
#' @examples
#' \dontrun{
#' results <- compare_network_estimation(group_regulation, iterations = 100)
#'
#' # Default histogram
#' plot_network_estimation(results)
#'
#' # Different metrics
#' plot_network_estimation(results, metric_name = "Spearman")
#' plot_network_estimation(results, metric_name = "Euclidean")
#'
#' # Different plot types
#' plot_network_estimation(results, plot_type = "boxplot")
#' plot_network_estimation(results, plot_type = "density")
#' }
#'
#' @import ggplot2
#' @import dplyr
#' @export
plot_network_estimation <- function(results,
                                     metric_name = "Pearson",
                                     plot_type = "histogram",
                                     show_stats = TRUE,
                                     bins = NULL) {

  # Extract data frame from network_estimation object
  if (inherits(results, "network_estimation")) {
    plot_data <- results$individual
    iterations <- results$params$iterations
  } else if (is.data.frame(results)) {
    plot_data <- results
    iterations <- length(unique(plot_data$iteration))
  } else {
    stop("results must be a network_estimation object or data frame")
  }

  # Filter for the specific metric
  plot_data <- plot_data %>%
    dplyr::filter(metric == metric_name)

  if (nrow(plot_data) == 0) {
    stop(paste("No data found for metric:", metric_name))
  }

  # Calculate bins based on iterations (Sturges' rule, capped)
  if (is.null(bins)) {
    bins <- min(50, max(20, ceiling(log2(iterations) + 1) * 3))
  }

  # Calculate statistics per model
  stats_df <- plot_data %>%
    dplyr::group_by(model_type) %>%
    dplyr::summarise(
      mean_val = mean(value, na.rm = TRUE),
      median_val = median(value, na.rm = TRUE),
      sd_val = sd(value, na.rm = TRUE),
      .groups = "drop"
    )

  # Get overall stats for subtitle
  overall_mean <- mean(plot_data$value, na.rm = TRUE)
  overall_median <- median(plot_data$value, na.rm = TRUE)

  # Create plot based on type
  if (plot_type == "histogram") {
    binwidth <- diff(range(plot_data$value, na.rm = TRUE)) / bins

    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = value)) +
      ggplot2::geom_histogram(
        bins = bins,
        fill = "turquoise3",
        color = "gray30",
        alpha = 0.8
      ) +
      ggplot2::geom_density(
        ggplot2::aes(y = ggplot2::after_stat(density) * binwidth * nrow(plot_data) / length(unique(plot_data$model_type))),
        color = "firebrick",
        linewidth = 1
      ) +
      ggplot2::geom_vline(
        data = stats_df,
        ggplot2::aes(xintercept = mean_val),
        color = "darkgreen",
        linewidth = 1.2
      ) +
      ggplot2::geom_vline(
        data = stats_df,
        ggplot2::aes(xintercept = median_val),
        color = "orange",
        linetype = "dashed",
        linewidth = 1
      ) +
      ggplot2::facet_wrap(~ toupper(model_type), ncol = 1, scales = "free_y") +
      ggplot2::labs(
        title = paste("Distribution of", metric_name, "Correlation"),
        subtitle = if (show_stats) {
          paste0("Green line = Mean | Orange dashed = Median")
        } else NULL,
        x = paste(metric_name, "Correlation"),
        y = "Frequency"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 14),
        strip.text = ggplot2::element_text(face = "bold", size = 12),
        axis.title = ggplot2::element_text(face = "bold")
      )

    # Add stat annotations
    if (show_stats) {
      p <- p +
        ggplot2::geom_text(
          data = stats_df,
          ggplot2::aes(
            x = mean_val,
            y = Inf,
            label = sprintf("Mean=%.4f", mean_val)
          ),
          vjust = 2,
          hjust = -0.1,
          size = 3.5,
          color = "darkgreen",
          fontface = "bold"
        )
    }

  } else if (plot_type == "boxplot") {
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = toupper(model_type), y = value, fill = model_type)) +
      ggplot2::geom_boxplot(alpha = 0.8, show.legend = FALSE) +
      ggplot2::stat_summary(
        fun = mean,
        geom = "point",
        shape = 18,
        size = 4,
        color = "darkgreen"
      ) +
      ggplot2::labs(
        title = paste("Distribution of", metric_name, "by Model Type"),
        subtitle = "Diamond = Mean",
        x = "Model Type",
        y = paste(metric_name, "Correlation")
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 14),
        axis.title = ggplot2::element_text(face = "bold")
      )

  } else if (plot_type == "density") {
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = value, fill = toupper(model_type))) +
      ggplot2::geom_density(alpha = 0.5) +
      ggplot2::geom_vline(
        data = stats_df,
        ggplot2::aes(xintercept = mean_val, color = toupper(model_type)),
        linewidth = 1,
        show.legend = FALSE
      ) +
      ggplot2::labs(
        title = paste("Density of", metric_name, "by Model Type"),
        subtitle = "Vertical lines = Means",
        x = paste(metric_name, "Correlation"),
        y = "Density",
        fill = "Model"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 14),
        axis.title = ggplot2::element_text(face = "bold"),
        legend.position = "bottom"
      )

  } else if (plot_type == "bar") {
    p <- ggplot2::ggplot(stats_df, ggplot2::aes(x = toupper(model_type), y = mean_val, fill = model_type)) +
      ggplot2::geom_col(alpha = 0.8, show.legend = FALSE) +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val),
        width = 0.2
      ) +
      ggplot2::geom_text(
        ggplot2::aes(label = sprintf("%.4f", mean_val)),
        vjust = -0.5,
        size = 4,
        fontface = "bold"
      ) +
      ggplot2::labs(
        title = paste("Mean", metric_name, "by Model Type"),
        subtitle = "Error bars = 1 SD",
        x = "Model Type",
        y = paste("Mean", metric_name)
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 14),
        axis.title = ggplot2::element_text(face = "bold")
      )

  } else if (plot_type == "ridgeline") {
    if (!requireNamespace("ggridges", quietly = TRUE)) {
      stop("Package 'ggridges' required. Install with: install.packages('ggridges')")
    }

    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = value, y = toupper(model_type), fill = model_type)) +
      ggridges::geom_density_ridges(alpha = 0.8, show.legend = FALSE) +
      ggplot2::labs(
        title = paste("Distribution of", metric_name, "by Model Type"),
        x = paste(metric_name, "Correlation"),
        y = "Model Type"
      ) +
      ggridges::theme_ridges() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 14)
      )

  } else {
    stop("Invalid plot_type. Choose: 'histogram', 'boxplot', 'density', 'bar', 'ridgeline'")
  }

  p
}


#' Print Network Estimation Results
#'
#' @param x A network_estimation object.
#' @param ... Additional arguments (unused).
#'
#' @export
print.network_estimation <- function(x, ...) {
  cat("Network Estimation Comparison\n")
  cat("=============================\n\n")

  cat("Models compared:", paste(toupper(x$params$model_types), collapse = " vs "), "\n")
  cat("Iterations:", x$params$iterations, "\n")
  cat("Sample size:", paste0(x$params$sampling_percent * 100, "%"), "\n")
  cat("Sequences:", x$params$n_sequences, "\n\n")

  cat("Ranking (by Pearson correlation):\n")
  ranking_df <- x$aggregated %>%
    dplyr::filter(metric == "Pearson") %>%
    dplyr::arrange(dplyr::desc(mean))

  for (i in seq_len(nrow(ranking_df))) {
    cat(sprintf("  %d. %s: r = %.4f (sd = %.4f)\n",
                i,
                toupper(ranking_df$model_type[i]),
                ranking_df$mean[i],
                ranking_df$sd[i]))
  }

  cat("\nWinner:", toupper(x$winner), "\n")

  invisible(x)
}


# Backward compatibility alias
#' @rdname compare_network_estimation
#' @export
compare_tna_models <- compare_network_estimation


#' Cross-Validate TNA Model Types
#'
#' @description
#' Perform cross-validation by testing different TNA model types on the same
#' data, allowing comparison of model performance.
#'
#' @param data A data frame containing sequence data.
#' @param model_types Character vector of model types to test.
#'   Default: c("relative", "frequency", "co-occurrence").
#' @param sampling_percent Proportion of data for sampling. Default: 0.3.
#' @param iterations Number of iterations per model type. Default: 50.
#' @param seed Random seed for reproducibility. Default: NULL.
#' @param verbose Logical. Print progress. Default: TRUE.
#'
#' @return A list containing cross-validation results for each model type.
#'
#' @examples
#' \dontrun{
#' data(group_regulation, package = "tna")
#' cv_results <- cross_validate_tna(
#'   group_regulation,
#'   model_types = c("relative", "frequency"),
#'   iterations = 30
#' )
#' }
#'
#' @import tna
#' @import dplyr
#' @export
cross_validate_tna <- function(data,
                                model_types = c("relative", "frequency", "co-occurrence"),
                                sampling_percent = 0.3,
                                iterations = 50,
                                seed = NULL,
                                verbose = TRUE) {

  if (!is.null(seed)) set.seed(seed)

  if (!is.data.frame(data)) {
    stop("data must be a data frame")
  }

  cv_results <- list()

  for (model_type in model_types) {
    if (verbose) cat("Cross-validating model type:", model_type, "\n")

    tryCatch({
      # Build the model
      model <- tna::build_model(data, type = model_type)

      # Perform sampling analysis
      results <- run_sampling_analysis(
        model = model,
        sampling_percent = sampling_percent,
        iterations = iterations,
        seed = NULL,
        verbose = FALSE
      )

      cv_results[[model_type]] <- list(
        model = model,
        results = results
      )

    }, error = function(e) {
      warning("Error processing model type ", model_type, ": ", e$message)
      cv_results[[model_type]] <- list(error = e$message)
    })
  }

  cv_results
}
