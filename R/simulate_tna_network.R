#' Simulate a Single TNA Network
#'
#' @description
#' Generate a single fitted TNA network model with learning state names.
#' This is the simplest way to get a ready-to-use tna model object.
#'
#' @param n_states Integer. Number of states (8-10 recommended). Default: 9.
#' @param n_sequences Integer or NULL. Number of sequences to simulate.
#'   If NULL (default), randomly selects between 600-1200.
#' @param seq_length Integer. Length of each sequence. Default: 25.
#' @param categories Character vector. Learning state categories to use.
#'   Options: "metacognitive", "cognitive", "behavioral", "social",
#'   "motivational", "affective", "group_regulation", "lms", or "all".
#'   Default: c("metacognitive", "cognitive").
#' @param seed Integer or NULL. Random seed for reproducibility. Default: NULL.
#'
#' @return A tna model object (class "tna") ready for use with all tna functions:
#'   \code{plot()}, \code{centralities()}, \code{communities()}, etc.
#'
#' @details
#' This function provides the simplest interface for generating a TNA network:
#' \enumerate{
#'   \item Selects learning state names from specified categories
#'   \item Generates random transition probabilities
#'   \item Simulates Markov chain sequences
#'   \item Fits and returns a tna model
#' }
#'
#' The returned model is 100% compatible with all tna package functions.
#'
#' @examples
#' library(Saqrlab)
#'
#' # Generate a network (simplest usage)
#' model <- simulate_tna_network(seed = 42)
#'
#' # Use with tna functions
#' \dontrun{
#' library(tna)
#' plot(model)
#' centralities(model)
#' communities(model)
#' }
#'
#' # Custom configuration
#' model <- simulate_tna_network(
#'   n_states = 8,
#'   n_sequences = 300,
#'   categories = "group_regulation",
#'   seed = 123
#' )
#'
#' # Access model components
#' model$weights
#' model$initial_probs
#'
#' @seealso
#' \code{\link{simulate_tna_networks}} for generating multiple networks,
#' \code{\link{simulate_sequences}} for generating sequences without fitting,
#' \code{\link{fit_network_model}} for fitting models to existing data.
#'
#' @importFrom seqHMM simulate_initial_probs simulate_transition_probs
#' @import tna
#' @export
simulate_tna_network <- function(n_states = 9,
                                  n_sequences = NULL,
                                  seq_length = 25,
                                  categories = c("metacognitive", "cognitive"),
                                  seed = NULL) {
  # Set seed
  if (!is.null(seed)) set.seed(seed)

 # Random n_sequences if not specified
  if (is.null(n_sequences)) {
    n_sequences <- sample(600:1200, 1)
  }


  # Get learning state names
states <- get_learning_states(categories = categories, n = n_states)

  # Generate random probabilities
  init_probs <- seqHMM::simulate_initial_probs(n_states)
  names(init_probs) <- states

  trans_probs <- seqHMM::simulate_transition_probs(n_states)
  dimnames(trans_probs) <- list(states, states)

  # Simulate sequences
  sequences <- simulate_sequences(
    trans_matrix = trans_probs,
    init_probs = init_probs,
    n_sequences = n_sequences,
    seq_length = seq_length
  )

  # Fit and return tna model
  model <- tna::tna(sequences)

  return(model)
}
