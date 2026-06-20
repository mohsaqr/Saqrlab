# ===========================================================================
# Scenario presets + tidy results / export
#
# A "scenario" is a named, hard-coded recipe: an ordered list of argument-lists,
# each ready to do.call() against a target simulate_* function. Scenarios let a
# user reproduce a whole study design (power sweep, ICC sweep, ...) with one
# call, and flatten the results into a single tidy data.frame for export.
# ===========================================================================


# ---------------------------------------------------------------------------
# Internal registry of scenario recipes
#
# Each entry holds:
#   $description : one-line human description
#   $simulator   : name of the target simulate_* function
#   $cases       : a list of named argument-lists (one per case)
# Every argument-list uses ONLY arguments that exist on the target simulator.
# ---------------------------------------------------------------------------

.saqr_scenarios <- function() {
  list(

    # Power analysis for a two-group t-test: effect size x sample size grid.
    power_ttest = list(
      description = "Two-group t-test power grid: Cohen's d x sample size.",
      simulator   = "simulate_ttest",
      cases = list(
        list(n_a = 30,  n_b = 30,  mean_a = 0, mean_b = 0.2),
        list(n_a = 30,  n_b = 30,  mean_a = 0, mean_b = 0.5),
        list(n_a = 30,  n_b = 30,  mean_a = 0, mean_b = 0.8),
        list(n_a = 100, n_b = 100, mean_a = 0, mean_b = 0.2),
        list(n_a = 100, n_b = 100, mean_a = 0, mean_b = 0.5),
        list(n_a = 100, n_b = 100, mean_a = 0, mean_b = 0.8)
      )
    ),

    # Multilevel model with the intraclass correlation swept across values.
    mlm_icc_sweep = list(
      description = "Multilevel model sweeping the intraclass correlation.",
      simulator   = "simulate_mlm",
      cases = list(
        list(n_clusters = 30, cluster_size = 20, icc = 0.05),
        list(n_clusters = 30, cluster_size = 20, icc = 0.10),
        list(n_clusters = 30, cluster_size = 20, icc = 0.20),
        list(n_clusters = 30, cluster_size = 20, icc = 0.30)
      )
    ),

    # IRT response data under each common dichotomous model.
    irt_models = list(
      description = "IRT response data across 1PL, 2PL and 3PL models.",
      simulator   = "simulate_irt",
      cases = list(
        list(n_persons = 300, n_items = 15, model = "1PL"),
        list(n_persons = 300, n_items = 15, model = "2PL"),
        list(n_persons = 300, n_items = 15, model = "3PL")
      )
    ),

    # Survival data with the censoring rate swept across values.
    survival_censoring = list(
      description = "Survival data sweeping the censoring rate.",
      simulator   = "simulate_survival",
      cases = list(
        list(n = 200, n_covariates = 2, censoring_rate = 0.1),
        list(n = 200, n_covariates = 2, censoring_rate = 0.3),
        list(n = 200, n_covariates = 2, censoring_rate = 0.5)
      )
    )
  )
}


# ---------------------------------------------------------------------------
# Internal: fetch one scenario by name, with a helpful error
# ---------------------------------------------------------------------------

.get_scenario_entry <- function(scenario) {
  stopifnot(is.character(scenario), length(scenario) == 1L)
  registry <- .saqr_scenarios()
  if (!scenario %in% names(registry)) {
    stop(sprintf(
      "Unknown scenario '%s'. Valid scenarios: %s",
      scenario, paste(names(registry), collapse = ", ")
    ), call. = FALSE)
  }
  registry[[scenario]]
}


# ---------------------------------------------------------------------------
# list_scenarios
# ---------------------------------------------------------------------------

#' List Available Scenario Presets
#'
#' @description Return a tidy catalogue of the built-in scenario presets. Each
#'   scenario is a hard-coded recipe of argument-lists targeting one
#'   \code{simulate_*} function. Print the result directly.
#'
#' @return A base \code{data.frame} with one row per scenario and columns
#'   \code{scenario}, \code{description}, \code{simulator} (the target
#'   \code{simulate_*} function) and \code{n_cases} (number of argument-lists
#'   in the recipe).
#'
#' @examples
#' list_scenarios()
#'
#' @export
list_scenarios <- function() {
  registry <- .saqr_scenarios()
  data.frame(
    scenario    = names(registry),
    description = vapply(registry, function(s) s$description, character(1)),
    simulator   = vapply(registry, function(s) s$simulator, character(1)),
    n_cases     = vapply(registry, function(s) length(s$cases), integer(1)),
    row.names   = NULL,
    stringsAsFactors = FALSE
  )
}


# ---------------------------------------------------------------------------
# get_scenario
# ---------------------------------------------------------------------------

#' Get a Scenario Recipe
#'
#' @description Retrieve the argument-lists for a named scenario. With
#'   \code{case = NULL} the full list of argument-lists is returned; with an
#'   integer \code{case} the single argument-list for that case is returned,
#'   ready to \code{do.call()} against the scenario's simulator.
#'
#' @param scenario Character scalar. The scenario name (see
#'   \code{\link{list_scenarios}}).
#' @param case Integer scalar or \code{NULL}. When \code{NULL} (default) the
#'   full list of argument-lists is returned. When an integer, the single
#'   argument-list at that index is returned.
#'
#' @return When \code{case = NULL}, a list of named argument-lists. When
#'   \code{case} is an integer, a single named argument-list.
#'
#' @examples
#' get_scenario("power_ttest")
#' get_scenario("power_ttest", case = 1)
#'
#' @export
get_scenario <- function(scenario, case = NULL) {
  entry <- .get_scenario_entry(scenario)
  cases <- entry$cases
  if (is.null(case)) {
    return(cases)
  }
  stopifnot(is.numeric(case), length(case) == 1L)
  case <- as.integer(case)
  if (case < 1L || case > length(cases)) {
    stop(sprintf(
      "case must be between 1 and %d for scenario '%s'.",
      length(cases), scenario
    ), call. = FALSE)
  }
  cases[[case]]
}


# ---------------------------------------------------------------------------
# run_scenario
# ---------------------------------------------------------------------------

#' Run a Whole Scenario
#'
#' @description Execute every case of a scenario by \code{do.call()}-ing its
#'   argument-list against the scenario's simulator. Returns a named list of
#'   the resulting \code{\link{saqr_sim}} objects, labelled \code{case1},
#'   \code{case2}, ...
#'
#' @param scenario Character scalar. The scenario name (see
#'   \code{\link{list_scenarios}}).
#' @param seed Integer scalar or \code{NULL}. When supplied, each case is given
#'   a deterministic seed of \code{seed + case index} so the whole run is
#'   reproducible while every case still differs. When \code{NULL} the
#'   simulators use their own defaults.
#'
#' @return A named list of \code{\link{saqr_sim}} objects, one per case, named
#'   \code{case1}, \code{case2}, ...
#'
#' @examples
#' s <- run_scenario("power_ttest", seed = 1)
#' s$case1
#'
#' @export
run_scenario <- function(scenario, seed = NULL) {
  entry <- .get_scenario_entry(scenario)
  cases <- entry$cases
  simulator <- match.fun(entry$simulator)

  if (!is.null(seed)) {
    stopifnot(is.numeric(seed), length(seed) == 1L)
    seed <- as.integer(seed)
  }

  indices <- seq_along(cases)
  results <- Map(function(args, idx) {
    if (!is.null(seed)) args$seed <- seed + idx
    do.call(simulator, args)
  }, cases, indices)

  names(results) <- paste0("case", indices)
  results
}


# ---------------------------------------------------------------------------
# tidy_simulation_results
# ---------------------------------------------------------------------------

# Internal: scalar (length-1, atomic) params of one saqr_sim, as a 1-row df.
.scalar_params_row <- function(params) {
  if (length(params) == 0L) return(data.frame(row.names = NULL))
  is_scalar <- vapply(params, function(p) {
    is.atomic(p) && length(p) == 1L && !is.null(p)
  }, logical(1))
  scalars <- params[is_scalar]
  if (length(scalars) == 0L) return(data.frame(row.names = NULL))
  as.data.frame(scalars, stringsAsFactors = FALSE, row.names = NULL)
}

# Internal: union-column row-bind that tolerates differing column sets.
.safe_rbind <- function(frames) {
  all_cols <- unique(unlist(lapply(frames, names)))
  filled <- lapply(frames, function(df) {
    missing <- setdiff(all_cols, names(df))
    df[missing] <- if (length(missing) > 0L) NA else NULL
    df[all_cols]
  })
  do.call(rbind, filled)
}

#' Flatten Simulation Results into One Tidy Data Frame
#'
#' @description Stack the \code{$data} of several \code{\link{saqr_sim}} objects
#'   into a single tidy \code{data.frame}. Each sim contributes its rows plus a
#'   \code{.case} column (the list name) and one column per scalar (length-1)
#'   ground-truth parameter in \code{$params}. Non-scalar parameters (matrices,
#'   vectors) are skipped. Columns are unioned across cases, so cases with
#'   different parameter sets still bind safely.
#'
#' @param x A named list of \code{\link{saqr_sim}} objects (for example the
#'   output of \code{\link{run_scenario}}), or a single \code{saqr_sim} object.
#'
#' @return A single base \code{data.frame} with all data rows stacked, a
#'   \code{.case} column, and one column per scalar parameter.
#'
#' @examples
#' s <- run_scenario("power_ttest", seed = 1)
#' tidy_simulation_results(s)
#'
#' @export
tidy_simulation_results <- function(x) {
  if (inherits(x, "saqr_sim")) x <- list(case1 = x)
  stopifnot(is.list(x), length(x) > 0L)
  stopifnot(all(vapply(x, inherits, logical(1), what = "saqr_sim")))

  case_names <- names(x)
  if (is.null(case_names)) case_names <- paste0("case", seq_along(x))

  frames <- Map(function(sim, nm) {
    d <- as.data.frame(sim$data, stringsAsFactors = FALSE)
    param_row <- .scalar_params_row(sim$params)
    base <- data.frame(.case = rep(nm, nrow(d)), stringsAsFactors = FALSE)
    if (ncol(param_row) > 0L) {
      base <- cbind(base, param_row[rep(1L, nrow(d)), , drop = FALSE],
                    row.names = NULL)
    }
    cbind(base, d, row.names = NULL)
  }, x, case_names)

  out <- .safe_rbind(frames)
  rownames(out) <- NULL
  out
}


# ---------------------------------------------------------------------------
# export_simulation
# ---------------------------------------------------------------------------

#' Export Simulation Results to CSV
#'
#' @description Write a tidy simulation data.frame to CSV. A list of
#'   \code{\link{saqr_sim}} objects (or a single one) is first flattened with
#'   \code{\link{tidy_simulation_results}}; a \code{data.frame} is written
#'   as-is.
#'
#' @param x A \code{data.frame}, a single \code{\link{saqr_sim}} object, or a
#'   list of \code{saqr_sim} objects.
#' @param file Character scalar. Path of the CSV file to write.
#'
#' @return The \code{file} path, invisibly.
#'
#' @examples
#' s <- run_scenario("power_ttest", seed = 1)
#' path <- tempfile(fileext = ".csv")
#' export_simulation(s, path)
#'
#' @export
export_simulation <- function(x, file) {
  stopifnot(is.character(file), length(file) == 1L)
  df <- if (is.data.frame(x)) x else tidy_simulation_results(x)
  utils::write.csv(df, file = file, row.names = FALSE)
  invisible(file)
}
