#' Evaluate Bootstrap Performance for a Single Simulation Run
#'
#' @description
#' Performs a single bootstrap evaluation run: generates sequences, fits a TNA
#' model, runs bootstrap analysis, and computes classification metrics for edge
#' detection compared to ground truth stable transitions.
#'
#' @param transition_matrix Square numeric matrix of transition probabilities.
#'   Rows must sum to 1. Row names define state names.
#' @param initial_probabilities Named numeric vector of initial state
#'   probabilities. Must sum to 1.
#' @param stable_transitions List of character vectors. Each vector contains
#'   two state names defining a stable (ground truth) transition pair.
#' @param max_seq_length Integer. Maximum length of each sequence.
#' @param num_rows Integer. Number of sequences to generate.
#' @param stability_prob Numeric in (0 to 1). Probability of following stable transitions.
#' @param unstable_mode Character. Mode for unstable transitions:
#'   "random_jump", "perturb_prob", or "unlikely_jump".
#' @param unstable_random_transition_prob Numeric in (0 to 1). Probability of
#'   unstable action.
#' @param unstable_perturb_noise Numeric in (0 to 1). Noise factor for perturbation mode.
#' @param unlikely_prob_threshold Numeric in (0 to 1). Threshold for unlikely transitions.
#' @param min_na Integer. Minimum NAs per sequence.
#' @param max_na Integer. Maximum NAs per sequence.
#' @param include_na Logical. Whether to include NAs.
#' @param consistency_range Numeric vector of length 2. Range for bootstrap
#'   consistency analysis (e.g., c(0.75, 1.25)).
#' @param level Numeric in (0,1). Significance level for p-value threshold
#'   (e.g., 0.05).
#'
#' @return A list containing:
#' \describe{
#'   \item{per_edge}{Data frame with per-edge results including ground truth,
#'     bootstrap significance, TP/TN/FP/FN, p-values, and metrics.}
#'   \item{bootstrap_summary_raw}{Raw bootstrap summary from tna::bootstrap.}
#'   \item{TP_matrix}{Matrix of true positives by edge.}
#'   \item{TN_matrix}{Matrix of true negatives by edge.}
#'   \item{FP_matrix}{Matrix of false positives by edge.}
#'   \item{FN_matrix}{Matrix of false negatives by edge.}
#' }
#' Returns NULL if bootstrap fails.
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Generates sequences using `simulate_sequences_advanced()`.
#'   \item Fits a TNA model and runs bootstrap analysis.
#'   \item Creates a ground truth matrix from stable_transitions.
#'   \item Compares bootstrap-significant edges to ground truth.
#'   \item Computes TP, TN, FP, FN and derived metrics per edge.
#' }
#'
#' An edge is considered "significant" if its p-value < level.
#'
#' @examples
#' \dontrun{
#' # Create transition matrix
#' trans_mat <- matrix(c(
#'   0.6, 0.3, 0.1,
#'   0.2, 0.6, 0.2,
#'   0.1, 0.2, 0.7
#' ), nrow = 3, byrow = TRUE)
#' rownames(trans_mat) <- colnames(trans_mat) <- c("A", "B", "C")
#' init_probs <- c(A = 0.33, B = 0.34, C = 0.33)
#'
#' # Define ground truth stable edges
#' stable <- list(c("A", "B"), c("B", "C"))
#'
#' # Run single bootstrap evaluation
#' result <- evaluate_bootstrap(
#'   transition_matrix = trans_mat,
#'   initial_probabilities = init_probs,
#'   stable_transitions = stable,
#'   max_seq_length = 30,
#'   num_rows = 100,
#'   stability_prob = 0.95,
#'   unstable_mode = "random_jump",
#'   unstable_random_transition_prob = 0.5,
#'   unstable_perturb_noise = 0.5,
#'   unlikely_prob_threshold = 0.1,
#'   min_na = 0,
#'   max_na = 5,
#'   include_na = TRUE,
#'   consistency_range = c(0.75, 1.25),
#'   level = 0.05
#' )
#' }
#'
#' @import dplyr
#' @import tna
#' @export
evaluate_bootstrap <- function(transition_matrix,
                                initial_probabilities,
                                stable_transitions,
                                max_seq_length,
                                num_rows,
                                stability_prob,
                                unstable_mode,
                                unstable_random_transition_prob,
                                unstable_perturb_noise,
                                unlikely_prob_threshold,
                                min_na,
                                max_na,
                                include_na,
                                consistency_range,
                                level) {
  # --- Input Validation ---
  stopifnot(
    is.numeric(consistency_range), length(consistency_range) == 2,
    consistency_range[1] <= consistency_range[2], all(consistency_range > 0),
    is.numeric(level), length(level) == 1, level > 0, level < 1
  )

  states_names <- rownames(transition_matrix)
  if (is.null(states_names)) stop("transition_matrix requires rownames.")
  num_states <- length(states_names)

  # 1. Generate Sequences
  sequences_df <- simulate_sequences_advanced(
    transition_matrix = transition_matrix,
    initial_probabilities = initial_probabilities,
    max_seq_length = max_seq_length,
    num_rows = num_rows,
    stable_transitions = stable_transitions,
    stability_prob = stability_prob,
    unstable_mode = unstable_mode,
    unstable_random_transition_prob = unstable_random_transition_prob,
    unstable_perturb_noise = unstable_perturb_noise,
    unlikely_prob_threshold = unlikely_prob_threshold,
    min_na = min_na,
    max_na = max_na,
    include_na = include_na
  )

  # 2. Bootstrap the network
  bootstrap_result <- NULL
  tna_obj <- NULL
  tryCatch(
    {
      tna_obj <- tna::tna(sequences_df)
      if (!inherits(tna_obj, "tna") || is.null(tna_obj$data) || nrow(tna_obj$data) < 2 || ncol(tna_obj$data) < 2 || is.null(tna_obj$weights)) {
        stop(sprintf(
          "Invalid tna object created from sequences (rows=%d, cols=%d).",
          if (is.null(tna_obj$data)) 0 else nrow(tna_obj$data),
          if (is.null(tna_obj$data)) 0 else ncol(tna_obj$data)
        ))
      }
      bootstrap_result <- tna_obj |> tna::bootstrap(
        consistency_range = consistency_range, level = level, method = "stability", iter = 10000,
      )
    },
    error = function(e) {
      warning(paste("TNA object creation or Bootstrap failed:", e$message))
    }
  )

  if (is.null(bootstrap_result) || !inherits(bootstrap_result, "tna_bootstrap") || is.null(bootstrap_result$summary)) {
    return(NULL)
  }

  # 3. Create Ground Truth Matrix
  ground_truth_matrix <- matrix(FALSE, nrow = num_states, ncol = num_states, dimnames = list(states_names, states_names))
  if (!is.null(stable_transitions)) {
    for (stable_pair in stable_transitions) {
      if (length(stable_pair) == 2 && all(stable_pair %in% states_names)) {
        ground_truth_matrix[stable_pair[1], stable_pair[2]] <- TRUE # matrix[FROM, TO]
      }
    }
  }

  # 4. Extract Significant Edges (p_value based)
  p_value_matrix <- matrix(1.0, nrow = num_states, ncol = num_states, dimnames = list(states_names, states_names))
  bootstrap_summary <- bootstrap_result$summary
  significant_edges <- matrix(FALSE, nrow = num_states, ncol = num_states) # Default

  if (inherits(bootstrap_summary, "data.frame") && nrow(bootstrap_summary) > 0) {
    req_cols <- c("from", "to", "p_value")
    if (!all(req_cols %in% names(bootstrap_summary))) {
      warning("Bootstrap summary missing required columns (from, to, p_value).")
    } else {
      for (i in 1:nrow(bootstrap_summary)) {
        from_st <- bootstrap_summary$from[i]
        to_st <- bootstrap_summary$to[i]
        p_val <- bootstrap_summary$p_value[i]
        if (!is.na(from_st) && !is.na(to_st) && from_st %in% states_names && to_st %in% states_names && !is.na(p_val) && is.numeric(p_val)) {
          p_value_matrix[from_st, to_st] <- p_val
        }
      }
      significant_edges <- p_value_matrix < level
    }
  }

  # 5. Calculate TP, TN, FP, FN per edge
  TP <- (significant_edges == TRUE) & (ground_truth_matrix == TRUE)
  TN <- (significant_edges == FALSE) & (ground_truth_matrix == FALSE)
  FP <- (significant_edges == TRUE) & (ground_truth_matrix == FALSE)
  FN <- (significant_edges == FALSE) & (ground_truth_matrix == TRUE)

  # 6. Create base per_edge data frame with CORRECT ground truth vectorization
  per_edge_base <- data.frame(
    from = rep(states_names, each = num_states),
    to = rep(states_names, times = num_states),
    ground_truth_stable = as.vector(t(ground_truth_matrix)), # Row-wise vectorization
    bootstrap_significant_run = as.vector(t(significant_edges)), # Match vectorization
    TP = as.vector(t(TP)), TN = as.vector(t(TN)),
    FP = as.vector(t(FP)), FN = as.vector(t(FN)),
    stringsAsFactors = FALSE
  )

  # 7. Safely merge p_value and weight
  bootstrap_summary_sub <- if (inherits(bootstrap_summary, "data.frame") && nrow(bootstrap_summary) > 0) {
    cols_to_keep <- intersect(c("from", "to", "p_value", "weight"), names(bootstrap_summary))
    if (length(cols_to_keep) < 2) data.frame(from = character(), to = character()) else bootstrap_summary[, cols_to_keep, drop = FALSE]
  } else {
    data.frame(from = character(), to = character(), p_value = numeric(), weight = numeric())
  }

  if (!all(c("from", "to") %in% names(bootstrap_summary_sub))) {
    bootstrap_summary_sub <- data.frame(from = character(), to = character(), p_value = numeric(), weight = numeric())
  }

  # 8. Join and Calculate Metrics
  per_edge_results <- per_edge_base %>%
    dplyr::left_join(bootstrap_summary_sub, by = c("from", "to")) %>%
    dplyr::mutate(
      p_value = if ("p_value" %in% names(.)) suppressWarnings(as.numeric(p_value)) else NA_real_,
      weight = if ("weight" %in% names(.)) suppressWarnings(as.numeric(weight)) else NA_real_,
      sensitivity_run = dplyr::if_else(TP + FN > 0, TP / (TP + FN), 0),
      specificity_run = dplyr::if_else(TN + FP > 0, TN / (TN + FP), 0),
      fpr_run = dplyr::if_else(FP + TN > 0, FP / (FP + TN), 0),
      fnr_run = dplyr::if_else(TP + FN > 0, FN / (TP + FN), 0)
    )

  # Return list
  return(list(
    per_edge = per_edge_results,
    bootstrap_summary_raw = bootstrap_summary,
    TP_matrix = TP, TN_matrix = TN, FP_matrix = FP, FN_matrix = FN
  ))
}
