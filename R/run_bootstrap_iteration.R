#' Run Bootstrap Iteration for a Single Simulation Run
#'
#' @description
#' Performs a single bootstrap evaluation run: generates sequences, fits a TNA
#' model, runs bootstrap analysis, and computes classification metrics for edge
#' detection compared to ground truth stable transitions.
#'
#' @param trans_matrix Square numeric matrix of transition probabilities.
#'   Rows must sum to 1. Row names define state names.
#' @param init_probs Named numeric vector of initial state
#'   probabilities. Must sum to 1.
#' @param stable_transitions List of character vectors. Each vector contains
#'   two state names defining a stable (ground truth) transition pair.
#' @param seq_length Integer. Maximum length of each sequence.
#' @param n_sequences Integer. Number of sequences to generate.
#' @param stability_prob Numeric in (0 to 1). Probability of following stable transitions.
#' @param unstable_mode Character. Mode for unstable transitions:
#'   "random_jump", "perturb_prob", or "unlikely_jump".
#' @param unstable_random_transition_prob Numeric in (0 to 1). Probability of
#'   unstable action.
#' @param unstable_perturb_noise Numeric in (0 to 1). Noise factor for perturbation mode.
#' @param unlikely_prob_threshold Numeric in (0 to 1). Threshold for unlikely transitions.
#' @param na_range Integer vector of length 2 (min, max) or single integer (min=max).
#'   Range of NA values per sequence. Default: c(0, 0).
#' @param include_na Logical. Whether to include NAs.
#' @param consistency_range Numeric vector of length 2. Range for bootstrap
#'   consistency analysis (e.g., c(0.75, 1.25)).
#' @param level Numeric in (0,1). Significance level for p-value threshold
#'   (e.g., 0.05).
#'
#' @param transition_matrix Deprecated. Use `trans_matrix` instead.
#' @param initial_probabilities Deprecated. Use `init_probs` instead.
#' @param max_seq_length Deprecated. Use `seq_length` instead.
#' @param num_rows Deprecated. Use `n_sequences` instead.
#' @param min_na Deprecated. Use `na_range` instead.
#' @param max_na Deprecated. Use `na_range` instead.
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
#' result <- run_bootstrap_iteration(
#'   trans_matrix = trans_mat,
#'   init_probs = init_probs,
#'   stable_transitions = stable,
#'   seq_length = 30,
#'   n_sequences = 100,
#'   stability_prob = 0.95,
#'   unstable_mode = "random_jump",
#'   unstable_random_transition_prob = 0.5,
#'   na_range = c(0, 5),
#'   include_na = TRUE,
#'   consistency_range = c(0.75, 1.25),
#'   level = 0.05
#' )
#'
#' # Old parameter names still work
#' result <- run_bootstrap_iteration(
#'   transition_matrix = trans_mat,
#'   initial_probabilities = init_probs,
#'   stable_transitions = stable,
#'   max_seq_length = 30,
#'   num_rows = 100
#' )
#' }
#'
#' @import dplyr
#' @import tna
#' @export
run_bootstrap_iteration <- function(trans_matrix = NULL,
                                init_probs = NULL,
                                stable_transitions,
                                seq_length = 20,
                                n_sequences = 100,
                                stability_prob = 0.95,
                                unstable_mode = "random_jump",
                                unstable_random_transition_prob = 0.5,
                                unstable_perturb_noise = 0.5,
                                unlikely_prob_threshold = 0.1,
                                na_range = c(0, 0),
                                include_na = TRUE,
                                consistency_range = c(0.75, 1.25),
                                level = 0.05,
                                # Backward compatibility - old parameter names
                                transition_matrix = NULL,
                                initial_probabilities = NULL,
                                max_seq_length = NULL,
                                num_rows = NULL,
                                min_na = NULL,
                                max_na = NULL) {
  # --- Backward compatibility: map old names to new names ---
  if (!is.null(transition_matrix)) trans_matrix <- transition_matrix
  if (!is.null(initial_probabilities)) init_probs <- initial_probabilities
  if (!is.null(max_seq_length)) seq_length <- max_seq_length
  if (!is.null(num_rows)) n_sequences <- num_rows

  # Handle na_range from min_na/max_na
  if (!is.null(min_na) || !is.null(max_na)) {
    min_na_val <- if (!is.null(min_na)) min_na else 0
    max_na_val <- if (!is.null(max_na)) max_na else min_na_val
    na_range <- c(min_na_val, max_na_val)
  }

  # Normalize na_range to length 2
  if (length(na_range) == 1) {
    na_range <- c(na_range, na_range)
  }

  # --- Input Validation ---
  stopifnot(
    is.numeric(consistency_range), length(consistency_range) == 2,
    consistency_range[1] <= consistency_range[2], all(consistency_range > 0),
    is.numeric(level), length(level) == 1, level > 0, level < 1
  )

  if (is.null(trans_matrix)) stop("trans_matrix is required.")
  if (is.null(init_probs)) stop("init_probs is required.")

  states_names <- rownames(trans_matrix)
  if (is.null(states_names)) stop("trans_matrix requires rownames.")
  num_states <- length(states_names)

  # 1. Generate Sequences
  sequences_df <- simulate_sequences_advanced(
    trans_matrix = trans_matrix,
    init_probs = init_probs,
    seq_length = seq_length,
    n_sequences = n_sequences,
    stable_transitions = stable_transitions,
    stability_prob = stability_prob,
    unstable_mode = unstable_mode,
    unstable_random_transition_prob = unstable_random_transition_prob,
    unstable_perturb_noise = unstable_perturb_noise,
    unlikely_prob_threshold = unlikely_prob_threshold,
    na_range = na_range,
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


# Backward compatibility alias (silent - no warning)

#' @rdname run_bootstrap_iteration
#' @export
evaluate_bootstrap <- function(...) run_bootstrap_iteration(...)
