#' Fit a Temporal Network Analysis Model
#'
#' @description
#' Fit a TNA (Temporal Network Analysis) model to sequence data using one of
#' several available model types.
#'
#' @param sequences Data frame of sequences. Each row is a sequence, each
#'   column is a time point.
#' @param model_type Character. Type of model to fit. One of:
#' \describe{
#'   \item{"tna"}{Standard Temporal Network Analysis model.}
#'   \item{"ftna"}{Filtered TNA model.}
#'   \item{"ctna"}{Conditional TNA model.}
#'   \item{"atna"}{Aggregated TNA model.}
#' }
#'
#' @return A fitted TNA model object of the appropriate class.
#'
#' @details
#' This function is a wrapper around the tna package's model fitting functions.
#' It provides a unified interface for fitting different types of TNA models.
#'
#' @examples
#' \dontrun{
#' # Generate some sequence data
#' trans_mat <- matrix(c(
#'   0.7, 0.2, 0.1,
#'   0.3, 0.5, 0.2,
#'   0.2, 0.3, 0.5
#' ), nrow = 3, byrow = TRUE)
#' rownames(trans_mat) <- colnames(trans_mat) <- c("A", "B", "C")
#' init_probs <- c(A = 0.5, B = 0.3, C = 0.2)
#'
#' sequences <- simulate_sequences(
#'   transition_matrix = trans_mat,
#'   initial_probabilities = init_probs,
#'   max_seq_length = 20,
#'   num_rows = 100
#' )
#'
#' # Fit different model types
#' model_tna <- fit_network_model(sequences, "tna")
#' model_ftna <- fit_network_model(sequences, "ftna")
#' }
#'
#' @import tna
#' @export
fit_network_model <- function(sequences, model_type) {
  switch(model_type,
    "tna" = tna::tna(sequences),
    "ftna" = tna::ftna(sequences),
    "ctna" = tna::ctna(sequences),
    "atna" = tna::atna(sequences),
    stop("Invalid model type: ", model_type)
  )
}
