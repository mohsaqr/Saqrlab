# ===========================================================================
# Unified simulation dispatcher -- one discoverable entry point over all the
# explicit-parameter simulators. Each target returns a saqr_sim object.
# ===========================================================================

# Internal: the canonical type -> simulator-function map. A named list of the
# functions themselves (not a switch), so dispatch is a single do.call and the
# list also drives list_simulators().
.simulator_map <- function() {
  list(
    ttest        = simulate_ttest,
    anova        = simulate_anova,
    correlation  = simulate_correlation,
    clusters     = simulate_clusters,
    prediction   = simulate_prediction,
    regression   = simulate_regression,
    lpa          = simulate_lpa,
    lca          = simulate_lca,
    fa           = simulate_fa,
    seq_clusters = simulate_seq_clusters,
    longitudinal = simulate_longitudinal,
    mlm          = simulate_mlm,
    growth       = simulate_growth,
    irt          = simulate_irt,
    survival     = simulate_survival,
    hmm          = simulate_hmm
  )
}

# Internal: one-line descriptions, keyed by type (same keys as .simulator_map).
.simulator_descriptions <- function() {
  c(
    ttest        = "Two-group comparison with known means/SDs (t-test).",
    anova        = "Multi-group comparison with known group means (one-way ANOVA).",
    correlation  = "Multivariate normal data from an explicit correlation/covariance matrix.",
    clusters     = "Gaussian mixture with known cluster centres and SDs.",
    prediction   = "Regression data with continuous + categorical predictors and known coefficients.",
    regression   = "Linear regression data with known coefficients and residual SD.",
    lpa          = "Latent profile analysis: continuous indicators from known profiles.",
    lca          = "Latent class analysis: binary indicators from known class probabilities.",
    fa           = "Factor analysis: indicators from known loadings and factor structure.",
    seq_clusters = "Categorical sequences drawn from known cluster-specific generators.",
    longitudinal = "Longitudinal panel data with a known within-person VAR(1) process.",
    mlm          = "Multilevel (mixed-effects) data with known fixed and random effects.",
    growth       = "Latent growth-curve data with known intercept/slope means and variances.",
    irt          = "Item response data from a known IRT model (difficulty/discrimination).",
    survival     = "Survival/time-to-event data with known covariate effects (Cox-style).",
    hmm          = "Hidden Markov sequences from known transition/emission matrices."
  )
}


#' Simulate Data via a Unified Dispatcher
#'
#' @description A single discoverable entry point over every explicit-parameter
#'   simulator in the package. You pass the simulation \code{type} plus the
#'   named arguments that simulator expects, and \code{simulate()} forwards them
#'   straight through, returning that simulator's \code{\link{saqr_sim}} object
#'   unchanged. Use \code{\link{list_simulators}} to see the available types.
#'
#'   For random, seed-driven generation (where the seed picks both the data and
#'   the structural parameters) and a bare \code{data.frame} return, use the
#'   separate \code{\link{simulate_data}} entry point instead.
#'
#' @param type Character scalar. The simulation type. One of the values in the
#'   \code{type} column of \code{\link{list_simulators}}: \code{"ttest"},
#'   \code{"anova"}, \code{"correlation"}, \code{"clusters"},
#'   \code{"prediction"}, \code{"regression"}, \code{"lpa"}, \code{"lca"},
#'   \code{"fa"}, \code{"seq_clusters"}, \code{"longitudinal"}, \code{"mlm"},
#'   \code{"growth"}, \code{"irt"}, \code{"survival"}, \code{"hmm"}.
#' @param ... Named arguments passed straight through to the matching
#'   simulator. See the help page of the underlying simulator (e.g.
#'   \code{\link{simulate_ttest}}) for the arguments it accepts.
#' @param seed Integer or NULL. Random seed, forwarded to the target simulator.
#'
#' @return The \code{\link{saqr_sim}} object returned by the dispatched
#'   simulator (with \code{$data}, \code{$params}, \code{$type}, \code{$seed}).
#'
#' @seealso \code{\link{list_simulators}} for the catalogue of types,
#'   \code{\link{simulate_data}} for random seed-driven generation.
#'
#' @examples
#' simulate("ttest", n_a = 50, n_b = 50, mean_a = 0, mean_b = 0.5, seed = 1)
#'
#' simulate("anova", n = 30, means = c(10, 12, 15), seed = 1)
#'
#' @export
simulate <- function(type, ..., seed = NULL) {
  stopifnot(is.character(type), length(type) == 1L)

  map <- .simulator_map()
  valid <- names(map)

  if (!type %in% valid) {
    stop(sprintf(
      "Unknown simulation type \"%s\". Valid types are: %s.",
      type, paste(valid, collapse = ", ")
    ), call. = FALSE)
  }

  fn <- map[[type]]
  args <- list(...)
  if ("seed" %in% names(formals(fn))) args$seed <- seed

  do.call(fn, args)
}


#' List Available Simulators
#'
#' @description Return the catalogue of simulation types reachable through
#'   \code{\link{simulate}}, one row per type, with the underlying function name
#'   and a one-line description. Print it directly.
#'
#'   Note: this lists the explicit-parameter simulators (those returning a
#'   \code{\link{saqr_sim}}). The separate \code{\link{simulate_data}} entry
#'   point provides random, seed-driven generation and returns a bare
#'   \code{data.frame}; it is not part of this catalogue.
#'
#' @return A base \code{data.frame} with one row per simulator and columns:
#'   \describe{
#'     \item{\code{type}}{The \code{type} string passed to \code{\link{simulate}}.}
#'     \item{\code{function}}{The name of the underlying simulator function.}
#'     \item{\code{description}}{A short one-line description.}
#'   }
#'
#' @seealso \code{\link{simulate}}
#'
#' @examples
#' list_simulators()
#'
#' @export
list_simulators <- function() {
  map <- .simulator_map()
  descriptions <- .simulator_descriptions()
  types <- names(map)

  data.frame(
    type        = types,
    `function`  = paste0("simulate_", types),
    description = unname(descriptions[types]),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}
