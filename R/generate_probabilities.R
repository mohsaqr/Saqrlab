#' Generate Transition and Initial Probabilities
#'
#' @description
#' Generate random transition probabilities and initial state probabilities
#' for a Markov chain using the seqHMM package.
#'
#' @param n_states Integer. Number of states in the Markov chain. Default: 5.
#' @param states Character vector. Names for the states.
#'   Must have at least `n_states` elements. If NULL, uses letters A, B, C, ...
#'   Default: NULL.
#' @param seed Integer or NULL. Random seed for reproducibility. Default: NULL.
#'
#' @param possible_state_names Deprecated. Use `states` instead.
#'
#' @return A list containing:
#' \describe{
#'   \item{initial_probs}{Named numeric vector of initial state probabilities.}
#'   \item{transition_probs}{Square matrix of transition probabilities with
#'     row and column names set to state names.}
#'   \item{state_names}{Character vector of state names used.}
#' }
#'
#' @details
#' The transition probabilities are generated using
#' `seqHMM::simulate_transition_probs()` and initial probabilities using
#' `seqHMM::simulate_initial_probs()`.
#'
#' @examples
#' \dontrun{
#' # Generate probabilities for 4 states with default names
#' probs <- generate_probabilities(n_states = 4, seed = 42)
#'
#' # Generate probabilities with custom state names
#' probs <- generate_probabilities(
#'   n_states = 4,
#'   states = c("A", "B", "C", "D", "E"),
#'   seed = 123
#' )
#'
#' # View initial probabilities
#' probs$initial_probs
#'
#' # View transition matrix
#' probs$transition_probs
#'
#' # Old parameter name still works
#' probs <- generate_probabilities(
#'   n_states = 3,
#'   possible_state_names = c("X", "Y", "Z")
#' )
#' }
#'
#' @importFrom seqHMM simulate_initial_probs simulate_transition_probs
#' @export
generate_probabilities <- function(n_states = 8,
                                    states = NULL,
                                    seed = NULL,
                                    # Backward compatibility
                                    possible_state_names = NULL) {
  # Backward compatibility
  if (!is.null(possible_state_names)) states <- possible_state_names

  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Determine state names
  if (is.null(states)) {
    if (n_states <= 26) {
      states <- LETTERS[1:n_states]
    } else {
      states <- paste0("S", 1:n_states)
    }
  } else {
    if (length(states) < n_states) {
      stop("states must have at least n_states elements.")
    }
    states <- states[1:n_states]
  }

  initial_probs <- seqHMM::simulate_initial_probs(n_states)
  names(initial_probs) <- states
  transition_probs <- seqHMM::simulate_transition_probs(n_states)
  dimnames(transition_probs) <- list(states, states)

  return(list(
    initial_probs = initial_probs,
    transition_probs = transition_probs,
    state_names = states
  ))
}
