#' Generate Transition and Initial Probabilities
#'
#' @description
#' Generate random transition probabilities and initial state probabilities
#' for a Markov chain using the seqHMM package.
#'
#' @param n_states Integer. Number of states in the Markov chain.
#' @param possible_state_names Character vector. Names for the states.
#'   Must have at least `n_states` elements.
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
#' This function uses a fixed seed (123) to ensure reproducible probability
#' generation. The transition probabilities are generated using
#' `seqHMM::simulate_transition_probs()` and initial probabilities using
#' `seqHMM::simulate_initial_probs()`.
#'
#' @examples
#' \dontrun{
#' # Generate probabilities for 4 states
#' probs <- generate_probabilities(
#'   n_states = 4,
#'   possible_state_names = c("A", "B", "C", "D", "E")
#' )
#'
#' # View initial probabilities
#' probs$initial_probs
#'
#' # View transition matrix
#' probs$transition_probs
#' }
#'
#' @importFrom seqHMM simulate_initial_probs simulate_transition_probs
#' @export
generate_probabilities <- function(n_states, possible_state_names) {
  set.seed(123) # Consistent seed
  state_names <- possible_state_names[1:n_states]
  initial_probs <- seqHMM::simulate_initial_probs(n_states)
  names(initial_probs) <- state_names
  transition_probs <- seqHMM::simulate_transition_probs(n_states)
  dimnames(transition_probs) <- list(state_names, state_names)
  return(list(
    initial_probs = initial_probs,
    transition_probs = transition_probs,
    state_names = state_names
  ))
}
