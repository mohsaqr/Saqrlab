# ===========================================================================
# Item Response Theory (IRT) simulation
# Returns saqr_sim with $data (person-by-item responses) and $params
# (ground-truth item and person parameters)
# ===========================================================================


#' Simulate Item Response Theory (IRT) Data
#'
#' @description Generate item-response data under a known IRT model with
#'   ground-truth item and person parameters. Supports the 1PL (Rasch), 2PL,
#'   3PL dichotomous models and the graded response model (GRM) for ordered
#'   polytomous items. Designed so that, at large \code{n_persons}, item
#'   proportion-correct recovers item difficulty (negatively) and person total
#'   score recovers true ability (positively).
#'
#'   Response probabilities are computed for the whole person-by-item grid in
#'   one vectorised pass (\code{outer()} / matrix ops) and responses are drawn
#'   with a single comparison against a uniform random matrix -- no loops.
#'
#' @param model Character. One of \code{"2PL"}, \code{"1PL"}, \code{"3PL"},
#'   \code{"GRM"}. Matched via \code{\link{match.arg}}.
#' @param n_persons Integer. Number of persons (rows). Default: 500.
#' @param n_items Integer. Number of items (columns). Default: 20.
#' @param n_categories Integer. Number of response categories. Must be 2 for
#'   the dichotomous models (1PL/2PL/3PL) and \code{> 2} for GRM. Default: 2.
#' @param a Numeric vector of length \code{n_items} or NULL. Item
#'   discriminations. NULL auto-generates \code{rlnorm(meanlog = 0,
#'   sdlog = 0.3)} (values near 1). For 1PL, discrimination is fixed at 1
#'   regardless of \code{a}.
#' @param b Numeric vector of length \code{n_items} or NULL. Item difficulties.
#'   NULL auto-generates \code{N(0, 1)}. For GRM, \code{b} is the item location
#'   around which ordered category thresholds are centred.
#' @param c Numeric vector of length \code{n_items} or NULL. Lower-asymptote
#'   (guessing) parameters. NULL auto-generates \code{0.2} for 3PL and
#'   \code{0} otherwise. Ignored (forced to 0) for 1PL/2PL/GRM.
#' @param theta Numeric vector of length \code{n_persons} or NULL. True person
#'   abilities. NULL auto-generates \code{N(0, 1)}.
#' @param seed Integer or NULL. Random seed. \code{set.seed()} is called only
#'   when not NULL.
#'
#' @return A \code{\link{saqr_sim}} object with:
#'   \describe{
#'     \item{\code{$data}}{data.frame with \code{n_persons} rows and columns
#'       \code{item1}...\code{itemJ}. Values are \code{0/1} for dichotomous
#'       models and \code{0}...\code{n_categories - 1} for GRM.}
#'     \item{\code{$params}}{list with \code{model}, \code{a}, \code{b},
#'       \code{c} (as used), \code{theta} (true abilities), \code{thresholds}
#'       (GRM category thresholds matrix, NULL otherwise), and
#'       \code{n_categories}.}
#'   }
#'
#' @examples
#' # 2PL with default auto-generated parameters
#' fit <- simulate_irt(n_persons = 500, n_items = 20, model = "2PL", seed = 1)
#' print(fit)
#' fit$params$b
#'
#' # Recover difficulty: easier items are answered correctly more often
#' prop_correct <- colMeans(as.data.frame(fit))
#' cor(prop_correct, fit$params$b)        # strongly negative
#'
#' # Recover ability: higher total score tracks higher true theta
#' total <- rowSums(as.data.frame(fit))
#' cor(total, fit$params$theta)           # strongly positive
#'
#' # Rasch (1PL): discrimination fixed at 1
#' simulate_irt(n_persons = 300, n_items = 10, model = "1PL", seed = 2)
#'
#' # Graded response model with 4 ordered categories
#' simulate_irt(n_persons = 400, n_items = 8, model = "GRM",
#'              n_categories = 4, seed = 3)
#'
#' @export
simulate_irt <- function(n_persons = 500, n_items = 20,
                         model = c("2PL", "1PL", "3PL", "GRM"),
                         n_categories = 2,
                         a = NULL, b = NULL, c = NULL, theta = NULL,
                         seed = NULL) {
  model <- match.arg(model)

  stopifnot(
    is.numeric(n_persons), length(n_persons) == 1L, n_persons >= 2L,
    is.numeric(n_items), length(n_items) == 1L, n_items >= 1L,
    is.numeric(n_categories), length(n_categories) == 1L, n_categories >= 2L
  )
  n_persons <- as.integer(n_persons)
  n_items   <- as.integer(n_items)
  n_categories <- as.integer(n_categories)

  if (model == "GRM") {
    stopifnot(n_categories > 2L)
  } else {
    stopifnot(n_categories == 2L)
  }

  if (!is.null(a)) stopifnot(is.numeric(a), length(a) == n_items, all(a > 0))
  if (!is.null(b)) stopifnot(is.numeric(b), length(b) == n_items)
  if (!is.null(c)) {
    stopifnot(is.numeric(c), length(c) == n_items, all(c >= 0), all(c < 1))
  }
  if (!is.null(theta)) {
    stopifnot(is.numeric(theta), length(theta) == n_persons)
  }

  if (!is.null(seed)) set.seed(as.integer(seed))

  # --- Person abilities -----------------------------------------------------
  if (is.null(theta)) theta <- stats::rnorm(n_persons, mean = 0, sd = 1)

  # --- Item parameters ------------------------------------------------------
  if (is.null(a)) a <- stats::rlnorm(n_items, meanlog = 0, sdlog = 0.3)
  if (is.null(b)) b <- stats::rnorm(n_items, mean = 0, sd = 1)
  if (is.null(c)) c <- if (model == "3PL") rep(0.2, n_items) else rep(0, n_items)

  # 1PL fixes discrimination at 1; 2PL/3PL keep a; GRM keeps a; no model uses c
  # except 3PL.
  if (model == "1PL") a <- rep(1, n_items)
  if (model != "3PL") c <- rep(0, n_items)

  item_names <- paste0("item", seq_len(n_items))

  if (model == "GRM") {
    out <- simulate_irt_grm(theta, a, b, n_categories)
    responses <- out$responses
    thresholds <- out$thresholds
  } else {
    # z[i, j] = a_j * (theta_i - b_j) for the whole grid, vectorised.
    z <- sweep(outer(theta, b, FUN = "-"), 2L, a, FUN = "*")
    p_correct <- sweep(stats::plogis(z), 2L,
                       (1 - c), FUN = "*")          # (1 - c) * plogis(z)
    p_correct <- sweep(p_correct, 2L, c, FUN = "+") # + c
    u <- matrix(stats::runif(n_persons * n_items),
                nrow = n_persons, ncol = n_items)
    responses <- (u < p_correct) * 1L
    thresholds <- NULL
  }

  colnames(responses) <- item_names
  df <- as.data.frame(responses)

  saqr_sim(
    data   = df,
    params = list(
      model        = model,
      a            = stats::setNames(a, item_names),
      b            = stats::setNames(b, item_names),
      c            = stats::setNames(c, item_names),
      theta        = theta,
      thresholds   = thresholds,
      n_categories = n_categories
    ),
    type = "irt",
    seed = seed,
    call = match.call()
  )
}


# ---------------------------------------------------------------------------
# Internal: graded response model draw (vectorised, no loops)
# ---------------------------------------------------------------------------

#' Draw Graded Response Model responses
#'
#' @description Internal helper. Builds per-item ordered category thresholds
#'   centred on each item location \code{b}, computes cumulative-logit
#'   boundary probabilities for the whole person-by-item grid, derives category
#'   probabilities, and draws an ordered category for every cell via inverse-CDF
#'   against one uniform random matrix.
#'
#' @param theta Numeric vector of person abilities.
#' @param a Numeric vector of item discriminations.
#' @param b Numeric vector of item locations.
#' @param n_categories Integer number of ordered categories (\code{> 2}).
#'
#' @return List with \code{responses} (n_persons x n_items integer matrix,
#'   values \code{0}...\code{n_categories - 1}) and \code{thresholds}
#'   (n_items x (n_categories - 1) matrix of ordered category thresholds).
#'
#' @keywords internal
#' @noRd
simulate_irt_grm <- function(theta, a, b, n_categories) {
  n_persons <- length(theta)
  n_items   <- length(b)
  n_thr     <- n_categories - 1L

  # Ordered thresholds per item: evenly spaced offsets around b, so each item
  # has increasing category boundaries b + offset (offsets centred at 0).
  offsets <- seq(-1.5, 1.5, length.out = n_thr)
  thresholds <- outer(b, offsets, FUN = "+")          # n_items x n_thr

  # Cumulative-logit boundary probabilities P(Y >= k) for every (person, item,
  # threshold). Compute one threshold-slice at a time across the full grid and
  # store the cell-wise category index by inverse-CDF.
  u <- matrix(stats::runif(n_persons * n_items),
              nrow = n_persons, ncol = n_items)

  # boundary_prob[[k]] = P(Y >= k) as a person x item matrix, for k = 1..n_thr.
  z_threshold <- function(k) {
    # a_j * (theta_i - thresholds_jk) over the grid
    sweep(outer(theta, thresholds[, k], FUN = "-"), 2L, a, FUN = "*")
  }
  boundary <- lapply(seq_len(n_thr), function(k) stats::plogis(z_threshold(k)))

  # P(Y >= k) is non-increasing in k. Response category = number of crossed
  # boundaries = sum over k of (u < P(Y >= k)). Vectorised via Reduce.
  crossings <- Reduce(`+`, lapply(boundary, function(p) (u < p) * 1L))
  responses <- crossings                       # values 0..n_thr == 0..K-1

  list(responses = responses, thresholds = thresholds)
}
