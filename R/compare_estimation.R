#' Compare TNA vs FTNA Estimation Across Simulations
#'
#' @description
#' Run multiple simulations comparing how well TNA and FTNA recover the true
#' transition structure, using `tna::compare()` for each comparison.
#'
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
#'   \item{summary}{Data frame with mean/sd of metrics by model type.}
#'   \item{raw_results}{Data frame with all simulation results.}
#'   \item{winner}{Which model performed better on average.}
#'   \item{params}{Parameters used for the simulation.}
#' }
#'
#' @details
#' For each simulation:
#' 1. Generate random transition probabilities (ground truth)
#' 2. Simulate sequences with optional NAs (varying lengths)
#' 3. Fit both TNA and FTNA models
#' 4. Compare each to ground truth using `tna::compare()`
#' 5. Collect metrics (correlation, RMSE, cosine similarity, etc.)
#'
#' @examples
#' \dontrun{
#' # Quick test with 100 simulations
#' results <- compare_estimation(n_simulations = 100, seed = 42)
#' results$summary
#' results$winner
#'
#' # Full study with 1000 simulations
#' results <- compare_estimation(
#'   n_simulations = 1000,
#'   n_sequences = 300,
#'   seq_length = 30,
#'   na_range = c(0, 10),
#'   parallel = TRUE,
#'   seed = 42
#' )
#' }
#'
#' @importFrom parallel detectCores mclapply
#' @importFrom tidyr pivot_wider
#' @import tna
#' @import dplyr
#' @export
compare_estimation <- function(n_simulations = 1000,
                                n_sequences = 200,
                                seq_length = 25,
                                n_states = 6,
                                na_range = c(0, 5),
                                scaling = "minmax",
                                seed = NULL,
                                verbose = TRUE,
                                parallel = FALSE,
                                cores = parallel::detectCores() - 1) {

  if (!is.null(seed)) set.seed(seed)

  # Normalize na_range
 if (length(na_range) == 1) na_range <- c(na_range, na_range)

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

    # Fit models
    model_tna <- tryCatch(tna::tna(sequences), error = function(e) NULL)
    model_ftna <- tryCatch(tna::ftna(sequences), error = function(e) NULL)

    if (is.null(model_tna) || is.null(model_ftna)) {
      return(NULL)
    }

    # Create ground truth "model" for comparison
    true_model <- list(weights = true_trans, class = "tna")
    class(true_model) <- "tna"

    # Compare using tna::compare()
    tna_comp <- tryCatch(
      tna::compare(model_tna, true_model, scaling = scaling),
      error = function(e) NULL
    )

    ftna_comp <- tryCatch(
      tna::compare(model_ftna, true_model, scaling = scaling),
      error = function(e) NULL
    )

    if (is.null(tna_comp) || is.null(ftna_comp)) {
      return(NULL)
    }

    # Extract metrics from tna::compare() output
    extract_metrics <- function(comp, model_type) {
      if (!is.null(comp$summary_metrics)) {
        df <- comp$summary_metrics
        df$model_type <- model_type
        df$sim_id <- sim_id
        return(df)
      }
      return(NULL)
    }

    tna_metrics <- extract_metrics(tna_comp, "tna")
    ftna_metrics <- extract_metrics(ftna_comp, "ftna")

    if (is.null(tna_metrics) && is.null(ftna_metrics)) {
      # Fallback: calculate basic metrics manually
      tna_cor <- cor(as.vector(true_trans), as.vector(model_tna$weights))
      ftna_norm <- model_ftna$weights / rowSums(model_ftna$weights)
      ftna_norm[is.nan(ftna_norm)] <- 0
      ftna_cor <- cor(as.vector(true_trans), as.vector(ftna_norm))

      return(data.frame(
        sim_id = sim_id,
        model_type = c("tna", "ftna"),
        correlation = c(tna_cor, ftna_cor)
      ))
    }

    rbind(tna_metrics, ftna_metrics)
  }

  # Run simulations
  if (verbose) {
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
  if ("value" %in% names(raw_results)) {
    # tna::compare() format
    summary_df <- raw_results %>%
      dplyr::group_by(model_type, category, metric) %>%
      dplyr::summarise(
        mean = mean(value, na.rm = TRUE),
        sd = sd(value, na.rm = TRUE),
        median = median(value, na.rm = TRUE),
        n = dplyr::n(),
        .groups = "drop"
      )
  } else {
    # Fallback format
    summary_df <- raw_results %>%
      dplyr::group_by(model_type) %>%
      dplyr::summarise(
        correlation_mean = mean(correlation, na.rm = TRUE),
        correlation_sd = sd(correlation, na.rm = TRUE),
        n = dplyr::n(),
        .groups = "drop"
      )
  }

  # Create comparison summary (TNA vs FTNA side by side)
  if ("metric" %in% names(summary_df)) {
    comparison <- summary_df %>%
      dplyr::filter(metric == "Pearson") %>%
      dplyr::select(model_type, mean, sd)

    tna_cor <- comparison$mean[comparison$model_type == "tna"]
    ftna_cor <- comparison$mean[comparison$model_type == "ftna"]
  } else if ("correlation_mean" %in% names(summary_df)) {
    tna_cor <- summary_df$correlation_mean[summary_df$model_type == "tna"]
    ftna_cor <- summary_df$correlation_mean[summary_df$model_type == "ftna"]
  } else {
    tna_cor <- ftna_cor <- NA
 }

  winner <- if (length(tna_cor) > 0 && length(ftna_cor) > 0 && !is.na(tna_cor) && !is.na(ftna_cor)) {
    if (abs(tna_cor - ftna_cor) < 0.001) "TIE" else if (tna_cor > ftna_cor) "TNA" else "FTNA"
  } else {
    "TIE"
  }

  # Create side-by-side comparison for key metrics
  if ("metric" %in% names(summary_df)) {
    key_metrics <- c("Pearson", "Spearman", "Kendall", "Euclidean", "Frobenius")
    comparison_df <- summary_df %>%
      dplyr::filter(metric %in% key_metrics) %>%
      tidyr::pivot_wider(
        id_cols = c(category, metric),
        names_from = model_type,
        values_from = c(mean, sd)
      ) %>%
      dplyr::select(category, metric, mean_tna, sd_tna, mean_ftna, sd_ftna) %>%
      dplyr::rename(
        TNA_mean = mean_tna, TNA_sd = sd_tna,
        FTNA_mean = mean_ftna, FTNA_sd = sd_ftna
      )
  } else {
    comparison_df <- summary_df
  }

  if (verbose) {
    cat(sprintf("\nCompleted %d/%d simulations successfully.\n",
                length(results_list), n_simulations))
    cat(sprintf("Pearson correlation: TNA=%.3f, FTNA=%.3f\n",
                ifelse(length(tna_cor) > 0, tna_cor, NA),
                ifelse(length(ftna_cor) > 0, ftna_cor, NA)))
    cat(sprintf("Winner: %s\n", winner))
  }

  list(
    comparison = comparison_df,
    summary = summary_df,
    raw_results = raw_results,
    winner = winner,
    params = list(
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
