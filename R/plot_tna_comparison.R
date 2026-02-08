#' Plot TNA Model Comparison Results
#'
#' @description
#' Visualize the distribution of metrics from TNA model comparisons.
#' Supports multiple plot types including histogram (default), boxplot,
#' ridgeline, bar, and density plots.
#'
#' @param results A data frame containing comparison results from
#'   `compare_network_estimation()` (the individual element) or `run_sampling_analysis()`.
#' @param plot_type Character string specifying the plot type:
#'   "histogram" (default), "boxplot", "ridgeline", "bar", or "density".
#' @param metric_name Character string specifying a single metric to plot.
#'   Default: "Pearson". If NULL, all metrics will be plotted.
#' @param model_types Character vector specifying the models to include.
#'   If NULL (default), all models in the data will be used.
#' @param binwidth Numeric value for histogram bin width. If NULL (default),
#'   calculated automatically.
#'
#' @return A ggplot object.
#'
#' @details
#' Available plot types:
#' \describe{
#'   \item{boxplot}{Box plots showing distribution by model type.}
#'   \item{ridgeline}{Ridge line plots for density comparison.}
#'   \item{bar}{Bar plots with mean values and error bars.}
#'   \item{histogram}{Histogram with density overlay (requires metric_name).}
#'   \item{density}{Density plots with mean/median lines (requires metric_name).}
#' }
#'
#' @examples
#' \dontrun{
#' library(tna)
#' data(group_regulation)
#'
#' # Run comparison
#' results <- compare_tna_models(
#'   group_regulation,
#'   model_types = c("tna", "atna"),
#'   iterations = 30
#' )
#'
#' # Boxplot of Pearson correlation
#' plot_tna_comparison(results$individual, "boxplot", metric_name = "Pearson")
#'
#' # Histogram of single metric
#' plot_tna_comparison(results$individual, "histogram", metric_name = "Pearson")
#'
#' # All metrics as bar chart
#' plot_tna_comparison(results$individual, "bar")
#' }
#'
#' @import ggplot2
#' @import dplyr
#' @export
plot_tna_comparison <- function(results,
                                 plot_type = "histogram",
                                 metric_name = "Pearson",
                                 model_types = NULL,
                                 binwidth = NULL) {

  if (!is.data.frame(results)) {
    stop("results must be a data frame")
  }

  # Filter data based on user input
  plot_data <- results

  if (!is.null(metric_name)) {
    plot_data <- plot_data %>%
      dplyr::filter(metric == metric_name)
    if (nrow(plot_data) == 0) {
      stop(paste("No data found for metric:", metric_name))
    }
  }

  if (!is.null(model_types)) {
    plot_data <- plot_data %>%
      dplyr::filter(model_type %in% model_types)
    if (nrow(plot_data) == 0) {
      stop(paste("No data found for model type(s):", paste(model_types, collapse = ", ")))
    }
  }

  # Generate plot based on type
  p <- NULL

  if (plot_type == "boxplot") {
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = model_type, y = value, fill = model_type)) +
      ggplot2::geom_boxplot(alpha = 0.8) +
      ggplot2::labs(
        title = if (!is.null(metric_name)) paste("Stability of", metric_name) else "Comparison of All Metrics",
        x = "Model Type",
        y = "Metric Value"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "bottom")

    if (is.null(metric_name)) {
      p <- p +
        ggplot2::facet_wrap(~ paste(category, metric, sep = " - "), scales = "free_y") +
        ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())
    }

  } else if (plot_type == "ridgeline") {
    if (!requireNamespace("ggridges", quietly = TRUE)) {
      stop("Package 'ggridges' required for ridgeline plots. Install with: install.packages('ggridges')")
    }

    if (is.null(metric_name) && length(unique(plot_data$metric)) > 1) {
      p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = value, y = metric, fill = model_type)) +
        ggridges::geom_density_ridges(alpha = 0.8, rel_min_height = 0.01) +
        ggplot2::labs(
          title = "Distribution of All Metrics",
          x = "Metric Value",
          y = "Metric"
        ) +
        ggridges::theme_ridges() +
        ggplot2::theme(legend.position = "bottom")
    } else {
      p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = value, y = model_type, fill = model_type)) +
        ggridges::geom_density_ridges(alpha = 0.8, rel_min_height = 0.01) +
        ggplot2::labs(
          title = if (!is.null(metric_name)) paste("Distribution of", metric_name) else "Distribution of Metrics",
          x = "Metric Value",
          y = "Model Type"
        ) +
        ggridges::theme_ridges() +
        ggplot2::theme(legend.position = "none")
    }

  } else if (plot_type == "bar") {
    # Aggregate data for bar plot
    summary_data <- plot_data %>%
      dplyr::group_by(model_type, category, metric) %>%
      dplyr::summarise(
        mean_value = mean(value, na.rm = TRUE),
        sd_value = sd(value, na.rm = TRUE),
        .groups = "drop"
      )

    p <- ggplot2::ggplot(summary_data, ggplot2::aes(x = model_type, y = mean_value, fill = model_type)) +
      ggplot2::geom_col(alpha = 0.8) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value),
                             width = 0.2) +
      ggplot2::labs(
        title = if (!is.null(metric_name)) paste("Mean Value of", metric_name) else "Mean Value of All Metrics",
        subtitle = "With Standard Deviation Error Bars",
        x = "Model Type",
        y = "Mean Metric Value"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "bottom")

    if (is.null(metric_name)) {
      p <- p +
        ggplot2::facet_wrap(~ paste(category, metric, sep = " - "), scales = "free_y") +
        ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())
    }

  } else if (plot_type == "density") {
    if (is.null(metric_name)) {
      stop("Density plot requires metric_name to be specified")
    }

    mean_val <- mean(plot_data$value, na.rm = TRUE)
    median_val <- median(plot_data$value, na.rm = TRUE)

    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = value, color = model_type)) +
      ggplot2::geom_density(linewidth = 1.2) +
      ggplot2::geom_vline(xintercept = mean_val, color = "darkgreen", linetype = "solid", linewidth = 1.2) +
      ggplot2::geom_vline(xintercept = median_val, color = "orange", linetype = "dotted", linewidth = 1.2) +
      ggplot2::labs(
        title = paste("Distribution of", metric_name),
        subtitle = paste0("Mean: ", round(mean_val, 3), " | Median: ", round(median_val, 3)),
        x = "Metric Value",
        y = "Density"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        legend.position = "bottom",
        plot.title = ggplot2::element_text(face = "bold"),
        axis.title = ggplot2::element_text(face = "bold")
      )

    if ("model_type" %in% names(plot_data) && length(unique(plot_data$model_type)) > 1) {
      p <- p +
        ggplot2::facet_wrap(~ model_type, ncol = 1, scales = "free_y")
    }

  } else if (plot_type == "histogram") {
    if (is.null(metric_name)) {
      stop("Histogram requires metric_name to be specified")
    }

    if (is.null(binwidth)) {
      binwidth <- diff(range(plot_data$value, na.rm = TRUE)) / 30
    }

    mean_val <- mean(plot_data$value, na.rm = TRUE)
    median_val <- median(plot_data$value, na.rm = TRUE)

    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = value)) +
      ggplot2::geom_histogram(binwidth = binwidth, fill = "turquoise", color = "black", alpha = 0.7) +
      ggplot2::geom_density(ggplot2::aes(y = ggplot2::after_stat(density) * binwidth * nrow(plot_data)),
                            color = "firebrick", linetype = "solid", linewidth = 1) +
      ggplot2::geom_vline(xintercept = mean_val, color = "darkgreen", linetype = "solid", linewidth = 1) +
      ggplot2::geom_vline(xintercept = median_val, color = "orange", linetype = "dotted", linewidth = 1) +
      ggplot2::labs(
        title = paste("Distribution of", metric_name),
        subtitle = paste0("Mean: ", round(mean_val, 3), " | Median: ", round(median_val, 3)),
        x = "Metric Value",
        y = "Frequency"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        legend.position = "none",
        plot.title = ggplot2::element_text(face = "bold"),
        axis.title = ggplot2::element_text(face = "bold")
      )

    if ("model_type" %in% names(plot_data) && length(unique(plot_data$model_type)) > 1) {
      p <- p +
        ggplot2::facet_wrap(~ model_type, ncol = 1, scales = "free_y")
    }

  } else {
    stop("Invalid plot_type. Choose from: 'boxplot', 'ridgeline', 'bar', 'histogram', 'density'")
  }

  p
}


#' Plot Sampling Distribution for a Single Metric
#'
#' @description
#' Create a visualization showing the distribution of a single metric's
#' values across all sampling iterations.
#'
#' @param individual_results Data frame of raw metric values from iterations,
#'   as returned by `run_sampling_analysis()` or `compare_tna_models()`.
#' @param category The category of the metric to plot (e.g., "Correlations").
#' @param metric The specific metric name to plot (e.g., "Pearson").
#'
#' @return A ggplot object showing the distribution.
#'
#' @examples
#' \dontrun{
#' library(tna)
#' model <- tna(group_regulation)
#' results <- run_sampling_analysis(model, iterations = 100)
#'
#' # Plot Pearson correlation distribution
#' plot_sampling_distribution(
#'   results$individual,
#'   category = "Correlations",
#'   metric = "Pearson"
#' )
#' }
#'
#' @import ggplot2
#' @import dplyr
#' @export
plot_sampling_distribution <- function(individual_results, category, metric) {
  if (!is.data.frame(individual_results)) {
    stop("individual_results must be a data frame")
  }

  # Filter for the specific metric
  filtered_data <- individual_results %>%
    dplyr::filter(category == !!category, metric == !!metric)

  if (nrow(filtered_data) == 0) {
    stop(paste("No data found for metric:", metric, "in category:", category))
  }

  # Calculate statistics
  mean_val <- mean(filtered_data$value, na.rm = TRUE)
  median_val <- median(filtered_data$value, na.rm = TRUE)

  # Create the plot
  p <- ggplot2::ggplot(filtered_data, ggplot2::aes(x = value)) +
    ggplot2::geom_histogram(bins = 30, fill = "dodgerblue", color = "white", alpha = 0.7) +
    ggplot2::geom_density(
      ggplot2::aes(y = ggplot2::after_stat(density) * nrow(filtered_data) *
                     (max(filtered_data$value) - min(filtered_data$value)) / 30),
      color = "red",
      linetype = "dashed",
      linewidth = 1
    ) +
    ggplot2::geom_vline(xintercept = mean_val, color = "darkgreen", linetype = "solid", linewidth = 1.2) +
    ggplot2::labs(
      title = paste("Distribution of", metric, "across Iterations"),
      subtitle = paste("Mean:", round(mean_val, 3), " | Median:", round(median_val, 3)),
      x = paste(category, "-", metric, "Value"),
      y = "Frequency"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      axis.title = ggplot2::element_text(face = "bold")
    )

  p
}
