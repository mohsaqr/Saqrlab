# ===========================================================================
# Hidden Markov model simulation -- latent state paths + observed emissions
# Returns saqr_sim with $data and $params
# ===========================================================================


#' Simulate Hidden Markov Model Sequences
#'
#' @description Generate observed symbol sequences from a discrete hidden
#'   Markov model (HMM). For each sequence a latent state path is drawn from
#'   the Markov chain defined by \code{init} (initial distribution) and
#'   \code{trans} (state transition matrix); each hidden state then emits an
#'   observed symbol according to the \code{emission} matrix. The true latent
#'   paths are returned in \code{$params$hidden_paths} for parameter-recovery
#'   testing.
#'
#' @param n_sequences Integer. Number of sequences to simulate. Default: 50.
#' @param seq_length Integer. Length (number of time points) of each sequence.
#'   Default: 30.
#' @param n_states Integer. Number of hidden states. Default: 2.
#' @param n_symbols Integer. Number of observable symbols. Default: 3.
#' @param trans Numeric matrix (\code{n_states x n_states}) or NULL. Row-
#'   stochastic state transition matrix. When NULL, a diagonally dominant
#'   matrix (states persist) is auto-generated.
#' @param emission Numeric matrix (\code{n_states x n_symbols}) or NULL. Row-
#'   stochastic emission matrix. When NULL, an auto-generated matrix where each
#'   state favours a distinct symbol is used.
#' @param init Numeric vector (length \code{n_states}) or NULL. Initial state
#'   distribution. When NULL, a uniform distribution is used.
#' @param seed Integer or NULL. Random seed.
#'
#' @return A \code{\link{saqr_sim}} object with:
#'   \describe{
#'     \item{\code{$data}}{long-format data.frame with columns
#'       \code{sequence_id} (integer), \code{time} (integer, 1...seq_length),
#'       and \code{symbol} (integer observed symbol in 1...n_symbols).}
#'     \item{\code{$params}}{list with \code{trans}, \code{emission},
#'       \code{init}, \code{n_states}, \code{n_symbols}, and
#'       \code{hidden_paths} (an \code{n_sequences x seq_length} integer matrix
#'       of the TRUE latent state sequences).}
#'   }
#'
#' @details The within-sequence latent recursion is inherently sequential and
#'   is implemented with \code{Reduce(..., accumulate = TRUE)} over the time
#'   index; the work across sequences is vectorised with \code{Map}/\code{vapply}.
#'
#' @examples
#' r <- simulate_hmm(n_sequences = 100, seq_length = 40, seed = 1)
#' print(r)
#' r$params$trans
#' r$params$emission
#' head(as.data.frame(r))
#'
#' # Explicit three-state model
#' tr <- matrix(c(0.8, 0.1, 0.1,
#'                0.1, 0.8, 0.1,
#'                0.1, 0.1, 0.8), nrow = 3, byrow = TRUE)
#' em <- matrix(c(0.7, 0.2, 0.1,
#'                0.1, 0.7, 0.2,
#'                0.2, 0.1, 0.7), nrow = 3, byrow = TRUE)
#' r2 <- simulate_hmm(n_sequences = 80, seq_length = 50, n_states = 3,
#'                    n_symbols = 3, trans = tr, emission = em, seed = 42)
#' summary(r2)
#'
#' @export
simulate_hmm <- function(n_sequences = 50, seq_length = 30, n_states = 2,
                         n_symbols = 3, trans = NULL, emission = NULL,
                         init = NULL, seed = NULL) {
  stopifnot(
    is.numeric(n_sequences), length(n_sequences) == 1L, n_sequences >= 1L,
    is.numeric(seq_length), length(seq_length) == 1L, seq_length >= 1L,
    is.numeric(n_states), length(n_states) == 1L, n_states >= 1L,
    is.numeric(n_symbols), length(n_symbols) == 1L, n_symbols >= 1L
  )

  n_sequences <- as.integer(n_sequences)
  seq_length  <- as.integer(seq_length)
  n_states    <- as.integer(n_states)
  n_symbols   <- as.integer(n_symbols)

  # Auto-generate row-stochastic transition matrix (diagonally dominant).
  if (is.null(trans)) {
    off <- 0.2 / max(n_states - 1L, 1L)
    trans <- matrix(off, nrow = n_states, ncol = n_states)
    diag(trans) <- 0.8
    if (n_states == 1L) trans <- matrix(1, 1L, 1L)
  }
  stopifnot(
    is.matrix(trans), is.numeric(trans),
    nrow(trans) == n_states, ncol(trans) == n_states,
    all(trans >= 0), all(abs(rowSums(trans) - 1) < 1e-8)
  )

  # Auto-generate emission matrix: each state favours a distinct symbol.
  if (is.null(emission)) {
    low <- 0.2 / max(n_symbols - 1L, 1L)
    emission <- matrix(low, nrow = n_states, ncol = n_symbols)
    favoured <- ((seq_len(n_states) - 1L) %% n_symbols) + 1L
    emission[cbind(seq_len(n_states), favoured)] <- 0.8
    # Renormalise rows (handles n_symbols == 1 and rounding).
    emission <- emission / rowSums(emission)
  }
  stopifnot(
    is.matrix(emission), is.numeric(emission),
    nrow(emission) == n_states, ncol(emission) == n_symbols,
    all(emission >= 0), all(abs(rowSums(emission) - 1) < 1e-8)
  )

  # Auto-generate initial distribution (uniform).
  if (is.null(init)) init <- rep(1 / n_states, n_states)
  stopifnot(
    is.numeric(init), length(init) == n_states,
    all(init >= 0), abs(sum(init) - 1) < 1e-8
  )

  if (!is.null(seed)) set.seed(as.integer(seed))

  states <- seq_len(n_states)
  symbols <- seq_len(n_symbols)

  # --- Latent state paths -------------------------------------------------
  # One path per sequence. Within a sequence the recursion is sequential, so
  # we use Reduce(accumulate = TRUE): each step samples the next state from
  # the transition row of the current state. The "innovations" (uniforms) are
  # pre-drawn so the recursion stays a pure deterministic fold.
  simulate_path <- function(seq_idx) {
    first_state <- sample(states, size = 1L, prob = init)
    if (seq_length == 1L) return(first_state)
    # Pre-draw one uniform per transition step; fold over them.
    innovations <- stats::runif(seq_length - 1L)
    next_state <- function(current, u) {
      cum <- cumsum(trans[current, ])
      cum[length(cum)] <- 1  # guard FP drift so the last bucket always catches u
      states[which(u <= cum)[1L]]
    }
    Reduce(next_state, innovations, accumulate = TRUE, init = first_state)
  }

  path_list <- Map(simulate_path, seq_len(n_sequences))
  hidden_paths <- matrix(
    unlist(path_list, use.names = FALSE),
    nrow = n_sequences, ncol = seq_length, byrow = TRUE
  )

  # --- Emissions ----------------------------------------------------------
  # Given the full latent matrix, emit a symbol for every (sequence, time)
  # cell. Vectorised across all cells via vapply over the flattened states.
  flat_states <- as.integer(hidden_paths)
  emit_symbol <- function(state) {
    sample(symbols, size = 1L, prob = emission[state, ])
  }
  flat_symbols <- vapply(flat_states, emit_symbol, integer(1L))
  symbol_mat <- matrix(flat_symbols, nrow = n_sequences, ncol = seq_length)

  # --- Long-format data ---------------------------------------------------
  df <- data.frame(
    sequence_id = rep(seq_len(n_sequences), times = seq_length),
    time        = rep(seq_len(seq_length), each = n_sequences),
    symbol      = as.integer(symbol_mat)
  )
  # Order by sequence then time for readability.
  ord <- order(df$sequence_id, df$time)
  df <- df[ord, , drop = FALSE]
  rownames(df) <- NULL

  saqr_sim(
    data   = df,
    params = list(
      trans        = trans,
      emission     = emission,
      init         = init,
      n_states     = n_states,
      n_symbols    = n_symbols,
      hidden_paths = hidden_paths
    ),
    type = "hmm",
    seed = seed,
    call = match.call()
  )
}
