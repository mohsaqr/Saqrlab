#' Compare Reliability Across Data Conditions
#'
#' @description
#' Run multiple simulations assessing split-half reliability under varying
#' data conditions (number of sequences, sequence length, number of states).
#' Uses `tna::reliability()` for each simulation with randomly sampled
#' parameters from user-specified ranges.
#'
#' @param n_simulations Integer. Number of simulations to run. Default: 1000.
#' @param n_sequences Integer or integer vector of length 2. Number of
#'   sequences per simulation. If length 2, each simulation samples a random
#'   integer in that range. Default: c(50, 500).
#' @param seq_length Integer or integer vector of length 2. Sequence length.
#'   If length 2, each simulation samples a random integer in that range.
#'   Default: c(10, 50).
#' @param n_states Integer or integer vector of length 2. Number of states.
#'   If length 2, each simulation samples a random integer in that range.
#'   Default: 6.
#' @param alpha Numeric or numeric vector of length 2. Dirichlet concentration
#'   parameter for [generate_probabilities()]. If length 2, each simulation
#'   samples a random value from `runif(1, alpha[1], alpha[2])`.
#'   Small values (e.g., 0.1) produce sparse transition matrices; large values
#'   (e.g., 10) produce near-uniform matrices. Default: 1.
#' @param diag_c Numeric or numeric vector of length 2. Diagonal boost for
#'   [generate_probabilities()]. If length 2, each simulation samples a random
#'   value from `runif(1, diag_c[1], diag_c[2])`. Higher values create
#'   "sticky" states with strong self-transitions. Default: 0.
#' @param na_range Integer vector of length 2. Range of NAs per sequence.
#'   Default: c(0, 5).
#' @param model_type Character. Which model to fit. Default: "tna".
#' @param reliability_iter Integer. Number of iterations for
#'   `tna::reliability()`. Default: 100.
#' @param reliability_split Numeric. Split proportion for
#'   `tna::reliability()`. Default: 0.5.
#' @param scaling Character. Scaling for `tna::reliability()`.
#'   Default: "none".
#' @param seed Integer or NULL. Random seed. Default: NULL.
#' @param verbose Logical. Print progress. Default: TRUE.
#' @param parallel Logical. Use parallel processing. Default: FALSE.
#' @param cores Integer. Number of cores for parallel. Default:
#'   detectCores() - 1.
#'
#' @return A list of class `"tna_reliability_comparison"` containing:
#' \describe{
#'   \item{raw_results}{Data frame with all simulation results (one row per
#'     metric per simulation, with parameter columns).}
#'   \item{summary}{Data frame with mean/sd aggregated across simulations
#'     grouped by metric.}
#'   \item{params}{List of input parameters for reference.}
#' }
#'
#' @details
#' For each simulation:
#' \enumerate{
#'   \item Randomly sample `n_sequences`, `seq_length`, `n_states` from
#'     user-specified ranges.
#'   \item Generate ground truth via `generate_probabilities()`.
#'   \item Simulate sequences via `simulate_sequences()` with those
#'     parameters and NAs.
#'   \item Fit a TNA model (e.g., `tna::tna(sequences)`).
#'   \item Run `tna::reliability()` on the fitted model.
#'   \item Extract the `$summary` table (22 metrics with
#'     mean/sd/median/min/max/q25/q75).
#'   \item Append simulation parameters as columns.
#' }
#'
#' @examples
#' \dontrun{
#' # Quick test
#' res <- compare_reliability(
#'   n_simulations = 10,
#'   reliability_iter = 20,
#'   seed = 42
#' )
#' print(res)
#'
#' # Vary number of states too
#' res <- compare_reliability(
#'   n_simulations = 100,
#'   n_states = c(3, 8),
#'   reliability_iter = 50,
#'   seed = 42
#' )
#' plot(res)
#'
#' # Parallel execution
#' res <- compare_reliability(
#'   n_simulations = 500,
#'   parallel = TRUE,
#'   seed = 42
#' )
#' }
#'
#' @importFrom parallel detectCores mclapply
#' @import tna
#' @import dplyr
#' @export
compare_reliability <- function(n_simulations = 1000,
                                n_sequences = c(50, 500),
                                seq_length = c(10, 50),
                                n_states = 6,
                                alpha = 1,
                                diag_c = 0,
                                na_range = c(0, 5),
                                model_type = "tna",
                                reliability_iter = 100,
                                reliability_split = 0.5,
                                scaling = "none",
                                seed = NULL,
                                verbose = TRUE,
                                parallel = FALSE,
                                cores = parallel::detectCores() - 1) {
  # Validate model_type
  if (!exists(model_type, where = asNamespace("tna"), mode = "function")) {
    stop(sprintf(
      "Model '%s' not found in tna package. Available: tna, ftna, ctna, atna, sna, tsn",
      model_type
    ))
  }

  if (!is.null(seed)) set.seed(seed)

  # Normalize ranges
  if (length(na_range) == 1) na_range <- c(na_range, na_range)
  if (length(n_sequences) == 1) n_sequences <- c(n_sequences, n_sequences)
  if (length(seq_length) == 1) seq_length <- c(seq_length, seq_length)
  if (length(n_states) == 1) n_states <- c(n_states, n_states)
  if (length(alpha) == 1) alpha <- c(alpha, alpha)
  if (length(diag_c) == 1) diag_c <- c(diag_c, diag_c)

  # Helper to sample integer from range
  sample_param <- function(range) {
    if (range[1] == range[2]) range[1] else sample(range[1]:range[2], 1)
  }

  # Helper to sample continuous value from range
  sample_continuous <- function(range) {
    if (range[1] == range[2]) range[1] else runif(1, range[1], range[2])
  }

  # Function to run single simulation
  run_single_sim <- function(sim_id) {
    # Sample parameters for this simulation
    sim_n_seq <- sample_param(n_sequences)
    sim_seq_len <- sample_param(seq_length)
    sim_n_states <- sample_param(n_states)
    sim_alpha <- sample_continuous(alpha)
    sim_diag_c <- sample_continuous(diag_c)

    tryCatch({
      # Generate ground truth
      true_probs <- generate_probabilities(
        n_states = sim_n_states, alpha = sim_alpha, diag_c = sim_diag_c
      )

      # Simulate sequences
      sequences <- simulate_sequences(
        trans_matrix = true_probs$transition_probs,
        init_probs = true_probs$initial_probs,
        n_sequences = sim_n_seq,
        seq_length = sim_seq_len,
        na_range = na_range
      )

      # Fit model
      model_fn <- get(model_type, envir = asNamespace("tna"))
      fitted_model <- model_fn(sequences)

      # Run reliability
      rel <- tna::reliability(
        fitted_model,
        iter = reliability_iter,
        split = reliability_split,
        scaling = scaling
      )

      # Extract summary and append simulation params
      summary_tbl <- rel$summary
      summary_tbl$sim_id <- sim_id
      summary_tbl$sim_n_sequences <- sim_n_seq
      summary_tbl$sim_seq_length <- sim_seq_len
      summary_tbl$sim_n_states <- sim_n_states
      summary_tbl$sim_alpha <- sim_alpha
      summary_tbl$sim_diag_c <- sim_diag_c
      summary_tbl$na_min <- na_range[1]
      summary_tbl$na_max <- na_range[2]

      summary_tbl
    }, error = function(e) {
      NULL
    })
  }

  # Run simulations
  if (verbose) {
    cat(sprintf(
      "Reliability comparison: %d simulations (model=%s, iter=%d)\n",
      n_simulations, model_type, reliability_iter
    ))
    cat(sprintf(
      "  n_sequences: %d-%d, seq_length: %d-%d, n_states: %d-%d, NA: %d-%d\n",
      n_sequences[1], n_sequences[2],
      seq_length[1], seq_length[2],
      n_states[1], n_states[2],
      na_range[1], na_range[2]
    ))
    cat(sprintf(
      "  alpha: %.2f-%.2f, diag_c: %.2f-%.2f\n",
      alpha[1], alpha[2], diag_c[1], diag_c[2]
    ))
  }

  if (parallel && cores > 1) {
    results_list <- parallel::mclapply(seq_len(n_simulations), function(i) {
      if (verbose && i %% 100 == 0) {
        cat(sprintf("  Completed %d/%d\n", i, n_simulations))
      }
      run_single_sim(i)
    }, mc.cores = cores)
  } else {
    results_list <- lapply(seq_len(n_simulations), function(i) {
      if (verbose && i %% 100 == 0) {
        cat(sprintf("  Completed %d/%d\n", i, n_simulations))
      }
      run_single_sim(i)
    })
  }

  # Combine results
  results_list <- Filter(Negate(is.null), results_list)
  if (length(results_list) == 0) {
    stop("All simulations failed.")
  }

  raw_results <- dplyr::bind_rows(results_list)

  # Summarize across all simulations grouped by metric
  summary_df <- raw_results %>%
    dplyr::group_by(metric) %>%
    dplyr::summarise(
      overall_mean = mean(mean, na.rm = TRUE),
      overall_sd = sd(mean, na.rm = TRUE),
      overall_median = median(mean, na.rm = TRUE),
      n = dplyr::n(),
      .groups = "drop"
    )

  if (verbose) {
    cat(sprintf(
      "\nCompleted %d/%d simulations successfully.\n",
      length(results_list), n_simulations
    ))
    pearson_row <- summary_df[summary_df$metric == "Pearson", ]
    if (nrow(pearson_row) > 0) {
      cat(sprintf(
        "Mean Pearson reliability: %.3f (sd=%.3f)\n",
        pearson_row$overall_mean, pearson_row$overall_sd
      ))
    }
  }

  result <- list(
    raw_results = raw_results,
    summary = summary_df,
    params = list(
      n_simulations = n_simulations,
      n_sequences = n_sequences,
      seq_length = seq_length,
      n_states = n_states,
      alpha = alpha,
      diag_c = diag_c,
      na_range = na_range,
      model_type = model_type,
      reliability_iter = reliability_iter,
      reliability_split = reliability_split,
      scaling = scaling,
      successful = length(results_list)
    )
  )
  class(result) <- "tna_reliability_comparison"
  result
}

#' Print Method for TNA Reliability Comparison
#'
#' @param x A `tna_reliability_comparison` object.
#' @param ... Additional arguments (ignored).
#'
#' @export
print.tna_reliability_comparison <- function(x, ...) {
  p <- x$params
  cat("TNA Reliability Comparison\n")
  cat(sprintf("  Model: %s\n", p$model_type))
  cat(sprintf(
    "  Simulations: %d successful / %d total\n",
    p$successful, p$n_simulations
  ))
  cat(sprintf(
    "  Reliability iterations: %d (split=%.1f, scaling=%s)\n",
    p$reliability_iter, p$reliability_split, p$scaling
  ))
  cat("\nParameter ranges:\n")
  cat(sprintf("  n_sequences: %d - %d\n", p$n_sequences[1], p$n_sequences[2]))
  cat(sprintf("  seq_length:  %d - %d\n", p$seq_length[1], p$seq_length[2]))
  cat(sprintf("  n_states:    %d - %d\n", p$n_states[1], p$n_states[2]))
  cat(sprintf("  alpha:       %.2f - %.2f\n", p$alpha[1], p$alpha[2]))
  cat(sprintf("  diag_c:      %.2f - %.2f\n", p$diag_c[1], p$diag_c[2]))
  cat(sprintf("  na_range:    %d - %d\n", p$na_range[1], p$na_range[2]))

  cat("\nReliability summary (mean across simulations):\n")
  s <- x$summary
  key_metrics <- c("Pearson", "Spearman", "Kendall", "Euclidean", "Frobenius")
  key_rows <- s[s$metric %in% key_metrics, ]
  key_rows <- key_rows[match(key_metrics, key_rows$metric), ]
  key_rows <- key_rows[!is.na(key_rows$metric), ]
  for (i in seq_len(nrow(key_rows))) {
    cat(sprintf(
      "  %-12s  %.3f (sd=%.3f)\n",
      key_rows$metric[i], key_rows$overall_mean[i], key_rows$overall_sd[i]
    ))
  }
  invisible(x)
}

#' Plot Method for TNA Reliability Comparison
#'
#' @param x A `tna_reliability_comparison` object.
#' @param metric Character. Metric to plot. Default: "Pearson".
#' @param type Character. Plot type: "scatter", "boxplot", or "histogram".
#'   Default: "scatter".
#' @param ... Additional arguments (ignored).
#'
#' @import ggplot2
#' @export
plot.tna_reliability_comparison <- function(x,
                                            metric = "Pearson",
                                            type = "scatter",
                                            ...) {
  raw <- x$raw_results
  metric_data <- raw[raw$metric == metric, ]

  if (nrow(metric_data) == 0) {
    stop(sprintf("Metric '%s' not found in results.", metric))
  }

  if (type == "scatter") {
    # Build long-format data for faceting
    plot_data <- data.frame(
      value = rep(metric_data$mean, 5),
      param_value = c(
        metric_data$sim_n_sequences,
        metric_data$sim_seq_length,
        metric_data$sim_n_states,
        metric_data$sim_alpha,
        metric_data$sim_diag_c
      ),
      param = rep(
        c("n_sequences", "seq_length", "n_states", "alpha", "diag_c"),
        each = nrow(metric_data)
      )
    )

    p <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(x = .data$param_value, y = .data$value)
    ) +
      ggplot2::geom_point(alpha = 0.3, size = 1) +
      ggplot2::geom_smooth(method = "loess", se = TRUE, color = "steelblue") +
      ggplot2::facet_wrap(~param, scales = "free_x") +
      ggplot2::labs(
        title = sprintf("Reliability: %s vs Data Parameters", metric),
        x = "Parameter Value",
        y = sprintf("%s (mean per sim)", metric)
      ) +
      ggplot2::theme_minimal()
  } else if (type == "boxplot") {
    p <- ggplot2::ggplot(
      metric_data,
      ggplot2::aes(x = factor(.data$sim_n_states), y = .data$mean)
    ) +
      ggplot2::geom_boxplot(fill = "steelblue", alpha = 0.5) +
      ggplot2::labs(
        title = sprintf("Reliability: %s by Number of States", metric),
        x = "Number of States",
        y = sprintf("%s (mean per sim)", metric)
      ) +
      ggplot2::theme_minimal()
  } else if (type == "histogram") {
    p <- ggplot2::ggplot(
      metric_data,
      ggplot2::aes(x = .data$mean)
    ) +
      ggplot2::geom_histogram(
        bins = 30, fill = "steelblue", alpha = 0.7, color = "white"
      ) +
      ggplot2::labs(
        title = sprintf("Distribution of %s Reliability", metric),
        x = sprintf("%s (mean per sim)", metric),
        y = "Count"
      ) +
      ggplot2::theme_minimal()
  } else {
    stop("type must be 'scatter', 'boxplot', or 'histogram'.")
  }

  p
}
