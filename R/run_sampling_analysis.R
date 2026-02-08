#' Run Sampling Analysis on TNA Models
#'
#' @description
#' Perform comprehensive TNA model analysis through repeated sampling and
#' statistical comparison. Uses a statistically correct sampling procedure
#' by comparing models built on independent data subsets (sample vs remaining).
#'
#' @param model A TNA model object created with tna() or related functions.
#' @param sampling_percent Numeric value between 0 and 1 specifying the
#'   proportion of data to use for the sample model. Default: 0.3.
#' @param iterations Integer specifying the number of sampling iterations.
#'   Default: 100.
#' @param seed Integer seed for reproducible random sampling. Default: NULL.
#' @param model_scaling Character string to override the model's scaling:
#'   NULL (default) uses original model's scaling, "skip" forces no scaling,
#'   or valid TNA scaling options ("minmax", "max", "rank").
#' @param verbose Logical. Print progress messages. Default: TRUE.
#'
#' @return A list containing:
#' \describe{
#'   \item{aggregated}{Data frame with aggregated summary statistics by metric.}
#'   \item{individual}{Data frame with raw metric values from each iteration.}
#'   \item{params}{Parameters used for the analysis.}
#' }
#'
#' @details
#' This function implements a statistically correct sampling procedure:
#'
#' 1. **Data Splitting**: For each iteration, splits the original data into:
#'    - Sample set (specified percentage)
#'    - Remaining set (complement)
#'
#' 2. **Model Building**: Builds TNA models on both datasets using the
#'    original model's type and configurable scaling.
#'
#' 3. **Model Comparison**: Compares sample model vs remaining model using
#'    TNA's compare() function, providing metrics across categories:
#'    - Correlations (Pearson, Spearman, Kendall)
#'    - Dissimilarities (Euclidean, Manhattan, etc.)
#'    - Similarities (Cosine, Jaccard, etc.)
#'
#' 4. **Result Aggregation**: Collects metrics across all iterations
#'    and computes summary statistics (mean, sd, median, quartiles).
#'
#' This approach is statistically superior to comparing sample vs original
#' because it compares models built on independent data subsets, providing
#' a true measure of sampling variability and model stability.
#'
#' @examples
#' \dontrun{
#' library(tna)
#' model <- tna(group_regulation)
#'
#' # Basic analysis with defaults
#' results <- run_sampling_analysis(model, iterations = 50)
#' results$aggregated
#'
#' # With custom parameters
#' results <- run_sampling_analysis(
#'   model,
#'   sampling_percent = 0.4,
#'   iterations = 100,
#'   seed = 42
#' )
#' }
#'
#' @importFrom stats sd median quantile complete.cases
#' @import dplyr
#' @import tna
#' @export
run_sampling_analysis <- function(model,
                                   sampling_percent = 0.3,
                                   iterations = 100,
                                   seed = NULL,
                                   model_scaling = NULL,
                                   verbose = TRUE) {

  # Set seed if provided
  if (!is.null(seed)) set.seed(seed)

  # Validate inputs
  if (!inherits(model, "tna")) {
    stop("model must be a TNA object")
  }
  if (sampling_percent <= 0 || sampling_percent >= 1) {
    stop("sampling_percent must be between 0 and 1 (exclusive)")
  }
  if (iterations < 1) {
    stop("iterations must be at least 1")
  }

  # Get original model parameters
  model_type <- attr(model, "type")
  original_scaling <- attr(model, "scaling")

  # Handle scaling override
  if (!is.null(model_scaling)) {
    if (model_scaling == "skip") {
      final_scaling <- NULL
      if (verbose) cat("User override: skipping scaling\n")
    } else {
      final_scaling <- model_scaling
      if (verbose) cat("User override: using scaling =", model_scaling, "\n")
    }
  } else {
    final_scaling <- if (is.null(original_scaling) || length(original_scaling) == 0 || original_scaling == "") {
      NULL
    } else {
      original_scaling
    }
  }

  if (verbose) {
    cat("=== SAMPLING ANALYSIS ===\n")
    cat("Type:", ifelse(is.null(model_type), "default", model_type), "\n")
    cat("Scaling:", ifelse(is.null(final_scaling), "none", final_scaling), "\n")
  }

  # Get original data from model
  original_data <- model$data
  if (!is.data.frame(original_data)) {
    original_data <- as.data.frame(original_data)
  }

  if (verbose) {
    cat("Data dimensions:", nrow(original_data), "x", ncol(original_data), "\n")
    cat("Running", iterations, "iterations...\n")
  }

  # Collect all metrics

all_metrics <- list()

  for (i in 1:iterations) {
    # Split data into sample and remaining
    n_total <- nrow(original_data)
    n_sample <- round(n_total * sampling_percent)

    sample_indices <- sample(1:n_total, n_sample, replace = FALSE)
    remaining_indices <- setdiff(1:n_total, sample_indices)

    sample_data <- original_data[sample_indices, , drop = FALSE]
    remaining_data <- original_data[remaining_indices, , drop = FALSE]

    # Build models on both datasets
    tryCatch({
      if (!is.null(model_type)) {
        if (!is.null(final_scaling)) {
          sample_model <- tna::build_model(sample_data, type = model_type, scaling = final_scaling)
          remaining_model <- tna::build_model(remaining_data, type = model_type, scaling = final_scaling)
        } else {
          sample_model <- tna::build_model(sample_data, type = model_type)
          remaining_model <- tna::build_model(remaining_data, type = model_type)
        }
      } else {
        if (!is.null(final_scaling)) {
          sample_model <- tna::build_model(sample_data, scaling = final_scaling)
          remaining_model <- tna::build_model(remaining_data, scaling = final_scaling)
        } else {
          sample_model <- tna::build_model(sample_data)
          remaining_model <- tna::build_model(remaining_data)
        }
      }

      # Compare the two models (sample vs remaining)
      comparison <- tna::compare(sample_model, remaining_model)

      # Store the metrics
      if (!is.null(comparison$summary_metrics)) {
        metrics_df <- comparison$summary_metrics
        metrics_df$iteration <- i
        all_metrics[[i]] <- metrics_df
      }

    }, error = function(e) {
      if (verbose) warning("Iteration ", i, " failed: ", e$message)
    })

    if (verbose && i %% 10 == 0) cat("  Completed", i, "/", iterations, "\n")
  }

  # Filter out NULL results
  valid_metrics <- Filter(Negate(is.null), all_metrics)

  if (length(valid_metrics) == 0) {
    stop("All iterations failed. Check your model and data.")
  }

  if (verbose) {
    cat("Aggregating", length(valid_metrics), "valid results...\n")
  }

  # Combine all results
  individual_results <- dplyr::bind_rows(valid_metrics)

  # Aggregate by category and metric
  aggregated_results <- individual_results %>%
    dplyr::group_by(category, metric) %>%
    dplyr::summarise(
      mean = mean(value, na.rm = TRUE),
      sd = sd(value, na.rm = TRUE),
      median = median(value, na.rm = TRUE),
      min = min(value, na.rm = TRUE),
      max = max(value, na.rm = TRUE),
      q25 = quantile(value, 0.25, na.rm = TRUE),
      q75 = quantile(value, 0.75, na.rm = TRUE),
      n_iterations = dplyr::n(),
      .groups = "drop"
    )

  if (verbose) {
    cat("Analysis complete:", nrow(aggregated_results), "metrics aggregated\n")
  }

  list(
    aggregated = aggregated_results,
    individual = individual_results,
    params = list(
      sampling_percent = sampling_percent,
      iterations = iterations,
      model_type = model_type,
      scaling = final_scaling,
      successful = length(valid_metrics)
    )
  )
}


#' Sample and Re-estimate TNA Model
#'
#' @description
#' Take a TNA model, sample a percentage of its data, and re-estimate a new
#' TNA model from the sampled data using the same parameters.
#'
#' @param model A TNA model object or data frame containing sequence data.
#' @param sampling_percent Numeric value between 0 and 1 indicating the
#'   proportion of data to sample. Default: 0.3.
#' @param model_type Character string specifying the model type. If NULL,
#'   uses the type from the original model. Default: NULL.
#' @param scaling Character string specifying the scaling method.
#'   If NULL, uses the scaling from the original model. Default: NULL.
#'
#' @return A new TNA model object estimated from the sampled data.
#'
#' @examples
#' \dontrun{
#' library(tna)
#' model <- tna(group_regulation)
#' sampled <- sample_tna(model, 0.3)
#' }
#'
#' @import tna
#' @export
sample_tna <- function(model, sampling_percent = 0.3, model_type = NULL, scaling = NULL) {
  # Input validation
  if (sampling_percent <= 0 || sampling_percent > 1) {
    stop("sampling_percent must be between 0 and 1")
  }

  if (inherits(model, "tna")) {
    original_data <- as.data.frame(model$data)
    labels <- model$labels
    num_sequences <- nrow(original_data)

    if (is.null(model_type)) {
      model_type <- attr(model, "type")
    }

    # Randomly sample row indices
    num_samples <- round(num_sequences * sampling_percent)
    sampled_indices <- sample(1:num_sequences, size = num_samples, replace = FALSE)
    sampled_df <- original_data[sampled_indices, , drop = FALSE]

    # Convert labels if needed
    sampled_matrix <- as.matrix(sampled_df)
    character_matrix <- matrix(
      labels[as.vector(sampled_matrix)],
      nrow = nrow(sampled_matrix),
      ncol = ncol(sampled_matrix)
    )
    sampled_data <- as.data.frame(character_matrix)

  } else if (is.data.frame(model)) {
    original_data <- model
    num_sequences <- nrow(original_data)

    if (is.null(model_type)) {
      model_type <- "relative"
    }

    num_samples <- round(num_sequences * sampling_percent)
    sampled_indices <- sample(1:num_sequences, size = num_samples, replace = FALSE)
    sampled_data <- original_data[sampled_indices, , drop = FALSE]

  } else {
    stop("model must be a TNA object or a data frame")
  }

  # Remove rows with NA values
  sampled_data <- sampled_data[complete.cases(sampled_data), ]

  if (nrow(sampled_data) == 0) {
    stop("No complete cases remain after sampling. Try a larger sample percentage.")
  }

  # Re-estimate the model
  if (!is.null(model_type)) {
    if (!is.null(scaling) && scaling != "" && scaling != "none") {
      new_model <- tna::build_model(sampled_data, type = model_type, scaling = scaling)
    } else {
      new_model <- tna::build_model(sampled_data, type = model_type)
    }
  } else {
    if (!is.null(scaling) && scaling != "" && scaling != "none") {
      new_model <- tna::build_model(sampled_data, scaling = scaling)
    } else {
      new_model <- tna::build_model(sampled_data)
    }
  }

  new_model
}
