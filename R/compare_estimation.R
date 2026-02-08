#' Compare Model Estimation Across Simulations
#'
#' @description
#' Run multiple simulations comparing how well different TNA model types
#' (tna, ftna, ctna, atna) recover the true transition structure,
#' using `tna::compare()` for each comparison.
#'
#' @param models Character vector. Model types to compare.
#'   Options: "tna", "ftna", "ctna", "atna". Default: c("tna", "ftna").
#' @param n_simulations Integer. Number of simulations to run. Default: 1000.
#' @param n_sequences Integer. Number of sequences per simulation. Default: 200.
#' @param seq_length Integer. Maximum sequence length. Default: 25.
#' @param n_states Integer. Number of states. Default: 6.
#' @param na_range Integer vector of length 2. Range of NAs per sequence for
#'   varying lengths. Default: c(0, 5).
#' @param scaling Character. Scaling for tna::compare(). Default: "minmax".
#' @param seed Integer or NULL. Random seed. Default: NULL.
#' @param verbose Logical. Print progress. Default: TRUE.
#' @param parallel Logical. Use parallel processing. Default: FALSE.
#' @param cores Integer. Number of cores for parallel. Default: detectCores() - 1.
#'
#' @return A list containing:
#' \describe{
#'   \item{comparison}{Side-by-side comparison of key metrics across models.}
#'   \item{summary}{Data frame with mean/sd of all metrics by model type.}
#'   \item{raw_results}{Data frame with all simulation results.}
#'   \item{ranking}{Models ranked by Pearson correlation (best to worst).}
#'   \item{winner}{Model with highest Pearson correlation.}
#'   \item{params}{Parameters used for the simulation.}
#' }
#'
#' @details
#' For each simulation:
#' 1. Generate random transition probabilities (ground truth)
#' 2. Simulate sequences with optional NAs (varying lengths)
#' 3. Fit all specified model types
#' 4. Compare each to ground truth using `tna::compare()`
#' 5. Collect metrics (Pearson, Spearman, Kendall, Euclidean, etc.)
#'
#' Available model types:
#' \itemize{
#'   \item \code{tna}: Standard transition network analysis (probabilities)
#'   \item \code{ftna}: Frequency-based TNA (raw counts)
#'   \item \code{ctna}: Concurrent TNA
#'   \item \code{atna}: Absorbing TNA
#' }
#'
#' @examples
#' \dontrun{
#' # Compare 2 models (default: tna vs ftna)
#' results <- compare_estimation(n_simulations = 100, seed = 42)
#'
#' # Compare 3 models
#' results <- compare_estimation(
#'   models = c("tna", "ftna", "ctna"),
#'   n_simulations = 100,
#'   seed = 42
#' )
#'
#' # Compare all 4 models
#' results <- compare_estimation(
#'   models = c("tna", "ftna", "ctna", "atna"),
#'   n_simulations = 500,
#'   parallel = TRUE,
#'   seed = 42
#' )
#'
#' # View results
#' results$comparison  # Side-by-side metrics
#' results$ranking     # Best to worst
#' results$winner      # Top performer
#' }
#'
#' @importFrom parallel detectCores mclapply
#' @importFrom tidyr pivot_wider
#' @import tna
#' @import dplyr
#' @export
compare_estimation <- function(models = c("tna", "ftna"),
                                n_simulations = 1000,
                                n_sequences = 200,
                                seq_length = 25,
                                n_states = 6,
                                na_range = c(0, 5),
                                scaling = "minmax",
                                seed = NULL,
                                verbose = TRUE,
                                parallel = FALSE,
                                cores = parallel::detectCores() - 1) {

  # Validate models
 valid_models <- c("tna", "ftna", "ctna", "atna")
  models <- match.arg(models, valid_models, several.ok = TRUE)
  if (length(models) < 2) {
    stop("At least 2 models required for comparison.")
  }

  if (!is.null(seed)) set.seed(seed)

  # Normalize na_range
  if (length(na_range) == 1) na_range <- c(na_range, na_range)

  # Model fitting functions
  fit_model <- function(sequences, model_type) {
    tryCatch(
      switch(model_type,
        "tna" = tna::tna(sequences),
        "ftna" = tna::ftna(sequences),
        "ctna" = tna::ctna(sequences),
        "atna" = tna::atna(sequences)
      ),
      error = function(e) NULL
    )
  }

  # Function to run single simulation
  run_single_sim <- function(sim_id) {
    # Generate ground truth
    true_probs <- generate_probabilities(n_states = n_states)
    true_trans <- true_probs$transition_probs
    true_init <- true_probs$initial_probs

    # Simulate sequences with varying lengths (NAs)
    sequences <- simulate_sequences(
      trans_matrix = true_trans,
      init_probs = true_init,
      n_sequences = n_sequences,
      seq_length = seq_length,
      na_range = na_range
    )

    # Create ground truth "model" for comparison
    true_model <- list(weights = true_trans)
    class(true_model) <- "tna"

    # Fit and compare each model
    all_metrics <- list()

    for (model_type in models) {
      fitted_model <- fit_model(sequences, model_type)

      if (is.null(fitted_model)) next

      # Compare using tna::compare()
      comp <- tryCatch(
        tna::compare(fitted_model, true_model, scaling = scaling),
        error = function(e) NULL
      )

      if (!is.null(comp) && !is.null(comp$summary_metrics)) {
        df <- comp$summary_metrics
        df$model_type <- model_type
        df$sim_id <- sim_id
        all_metrics[[model_type]] <- df
      } else {
        # Fallback: calculate basic correlation
        model_weights <- fitted_model$weights
        if (model_type == "ftna") {
          model_weights <- model_weights / rowSums(model_weights)
          model_weights[is.nan(model_weights)] <- 0
        }
        cor_val <- tryCatch(
          cor(as.vector(true_trans), as.vector(model_weights)),
          error = function(e) NA
        )
        all_metrics[[model_type]] <- data.frame(
          sim_id = sim_id,
          model_type = model_type,
          category = "Correlations",
          metric = "Pearson",
          value = cor_val
        )
      }
    }

    if (length(all_metrics) == 0) return(NULL)
    dplyr::bind_rows(all_metrics)
  }

  # Run simulations
  if (verbose) {
    cat(sprintf("Comparing %d models: %s\n", length(models), paste(toupper(models), collapse = " vs ")))
    cat(sprintf("Running %d simulations (n=%d, len=%d, states=%d, NA=%d-%d)...\n",
                n_simulations, n_sequences, seq_length, n_states, na_range[1], na_range[2]))
  }

  if (parallel && cores > 1) {
    results_list <- parallel::mclapply(1:n_simulations, function(i) {
      if (verbose && i %% 100 == 0) cat(sprintf("  Completed %d/%d\n", i, n_simulations))
      run_single_sim(i)
    }, mc.cores = cores)
  } else {
    results_list <- lapply(1:n_simulations, function(i) {
      if (verbose && i %% 100 == 0) cat(sprintf("  Completed %d/%d\n", i, n_simulations))
      run_single_sim(i)
    })
  }

  # Combine results
  results_list <- Filter(Negate(is.null), results_list)
  if (length(results_list) == 0) {
    stop("All simulations failed.")
  }

  raw_results <- dplyr::bind_rows(results_list)

  # Summarize by model type
  summary_df <- raw_results %>%
    dplyr::group_by(model_type, category, metric) %>%
    dplyr::summarise(
      mean = mean(value, na.rm = TRUE),
      sd = sd(value, na.rm = TRUE),
      median = median(value, na.rm = TRUE),
      n = dplyr::n(),
      .groups = "drop"
    )

  # Create ranking based on Pearson correlation
  pearson_summary <- summary_df %>%
    dplyr::filter(metric == "Pearson") %>%
    dplyr::arrange(dplyr::desc(mean)) %>%
    dplyr::select(model_type, mean, sd)

  ranking <- pearson_summary$model_type
  winner <- ranking[1]

  # Create side-by-side comparison for key metrics
  key_metrics <- c("Pearson", "Spearman", "Kendall", "Euclidean", "Frobenius")

  comparison_df <- summary_df %>%
    dplyr::filter(metric %in% key_metrics) %>%
    tidyr::pivot_wider(
      id_cols = c(category, metric),
      names_from = model_type,
      values_from = c(mean, sd),
      names_glue = "{model_type}_{.value}"
    )

  # Reorder columns for readability
  model_cols <- unlist(lapply(models, function(m) c(paste0(m, "_mean"), paste0(m, "_sd"))))
  existing_cols <- intersect(model_cols, names(comparison_df))
  comparison_df <- comparison_df %>%
    dplyr::select(category, metric, dplyr::all_of(existing_cols))

  if (verbose) {
    cat(sprintf("\nCompleted %d/%d simulations successfully.\n",
                length(results_list), n_simulations))
    cat("\nPearson correlation by model:\n")
    for (i in seq_along(ranking)) {
      m <- ranking[i]
      vals <- pearson_summary[pearson_summary$model_type == m, ]
      cat(sprintf("  %d. %s: %.3f (sd=%.3f)\n", i, toupper(m), vals$mean, vals$sd))
    }
    cat(sprintf("\nWinner: %s\n", toupper(winner)))
  }

  list(
    comparison = comparison_df,
    summary = summary_df,
    raw_results = raw_results,
    ranking = ranking,
    winner = winner,
    params = list(
      models = models,
      n_simulations = n_simulations,
      n_sequences = n_sequences,
      seq_length = seq_length,
      n_states = n_states,
      na_range = na_range,
      scaling = scaling,
      successful = length(results_list)
    )
  )
}
