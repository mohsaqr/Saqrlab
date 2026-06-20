# ===========================================================================
# validate_recovery -- parameter-recovery scoring
# "Did my method recover the truth I simulated?"
# ===========================================================================


# ---------------------------------------------------------------------------
# Internal: flatten the numeric entries of a $params list into a flat named
# numeric vector. A scalar `error_sd = 1` becomes "error_sd"; a named vector
# `betas = c(x1 = .5, x2 = -.3)` becomes "betas.x1", "betas.x2"; an unnamed
# vector `v = c(1, 2)` becomes "v.1", "v.2". Non-numeric entries are skipped.
# ---------------------------------------------------------------------------

flatten_params <- function(params) {
  stopifnot(is.list(params))

  pieces <- Map(
    function(value, name) {
      if (!is.numeric(value)) {
        return(numeric(0))
      }
      if (length(value) == 1L) {
        return(stats::setNames(as.numeric(value), name))
      }
      sub <- names(value)
      if (is.null(sub) || any(!nzchar(sub))) {
        sub <- as.character(seq_along(value))
      }
      stats::setNames(as.numeric(value), paste(name, sub, sep = "."))
    },
    params,
    names(params)
  )

  # Inner names are already fully qualified; drop the list-level names so
  # unlist() does not prepend them a second time.
  out <- unlist(unname(pieces), use.names = TRUE)
  if (is.null(out)) out <- numeric(0)
  out
}


# ---------------------------------------------------------------------------
# Internal: coerce the `sim` argument into a flat named numeric truth vector.
# Accepts a saqr_sim (uses $params), or a plain named list / named vector.
# ---------------------------------------------------------------------------

as_truth_vector <- function(sim) {
  if (inherits(sim, "saqr_sim")) {
    return(flatten_params(sim$params))
  }
  if (is.list(sim)) {
    return(flatten_params(sim))
  }
  stopifnot(is.numeric(sim), !is.null(names(sim)))
  stats::setNames(as.numeric(sim), names(sim))
}


#' Score Parameter Recovery Against Simulated Ground Truth
#'
#' @description Compare a method's point estimates to the known true values that
#'   generated a simulation, matched \strong{by name}. This is the core
#'   "did my method recover the truth I simulated?" check of the simulation
#'   laboratory: feed in a \code{\link{saqr_sim}} object (or a plain named list /
#'   vector of true values) and a named vector of estimates, and get back a tidy
#'   one-row-per-parameter table of errors and tolerance flags.
#'
#' @param sim A \code{\link{saqr_sim}} object (its \code{$params} supplies the
#'   ground truth) \emph{or} a plain named list / named numeric vector of true
#'   values.
#' @param estimates A named numeric vector (or named list of scalars) of point
#'   estimates produced by a fitted method. Names are matched against the
#'   (flattened) names of the truth.
#' @param params Optional character vector restricting which parameter names to
#'   compare. Default \code{NULL} uses the intersection of names available in
#'   both the truth and the estimates.
#' @param tolerance Positive numeric. Threshold a parameter must fall within to
#'   count as recovered. Interpreted as a relative tolerance when
#'   \code{relative = TRUE}, otherwise an absolute tolerance. Default \code{0.1}.
#' @param relative Logical. If \code{TRUE} (default), \code{within_tol} compares
#'   \code{rel_error <= tolerance}; if \code{FALSE}, compares
#'   \code{abs_error <= tolerance}.
#'
#' @return An S3 object of class \code{c("recovery_result", "data.frame")} with
#'   one row per compared parameter and columns:
#'   \describe{
#'     \item{\code{parameter}}{Parameter name (flattened, e.g. \code{coefs.x1}).}
#'     \item{\code{true}}{True value from the simulation.}
#'     \item{\code{estimate}}{Estimated value.}
#'     \item{\code{abs_error}}{\code{abs(estimate - true)}.}
#'     \item{\code{rel_error}}{\code{abs_error / abs(true)}, \code{NA} when
#'       \code{true == 0}.}
#'     \item{\code{within_tol}}{Logical recovery flag (see \code{relative}).}
#'   }
#'   Use \code{\link{summary.recovery_result}} for a one-row scorecard.
#'
#' @examples
#' # t-test: recover the true Cohen's d
#' sim <- simulate_ttest(n_a = 200, n_b = 200, mean_a = 0, mean_b = 0.5,
#'                       sd_a = 1, sd_b = 1, seed = 1)
#' fit <- t.test(score ~ group, data = sim$data)
#' d_hat <- unname(fit$estimate[2] - fit$estimate[1])   # mean_b - mean_a
#' validate_recovery(sim, estimates = c(cohens_d = d_hat))
#'
#' # Regression: recover known coefficients from lm()
#' rsim <- simulate_regression(
#'   coefs = c("(Intercept)" = 1, x1 = 2, x2 = -1.5),
#'   predictor_sds = c(x1 = 1, x2 = 1), error_sd = 1, n = 2000, seed = 7
#' )
#' fit_lm <- lm(y ~ x1 + x2, data = rsim$data)
#' est <- stats::setNames(coef(fit_lm),
#'                        paste0("coefs.", names(coef(fit_lm))))
#' validate_recovery(rsim, estimates = est)
#'
#' @export
validate_recovery <- function(sim, estimates, params = NULL,
                              tolerance = 0.1, relative = TRUE) {
  stopifnot(
    inherits(sim, "saqr_sim") || is.list(sim) ||
      (is.numeric(sim) && !is.null(names(sim))),
    is.numeric(tolerance), length(tolerance) == 1L, tolerance > 0,
    is.logical(relative), length(relative) == 1L
  )

  truth <- as_truth_vector(sim)

  if (is.list(estimates)) estimates <- unlist(estimates, use.names = TRUE)
  stopifnot(is.numeric(estimates), !is.null(names(estimates)))

  shared <- intersect(names(truth), names(estimates))
  if (!is.null(params)) {
    stopifnot(is.character(params))
    shared <- intersect(shared, params)
  }
  if (length(shared) == 0L) {
    stop("No parameter names are shared between truth and estimates",
         if (!is.null(params)) " (after applying `params`)" else "", ".")
  }

  true_vals <- unname(truth[shared])
  est_vals  <- unname(estimates[shared])

  abs_error <- abs(est_vals - true_vals)
  rel_error <- ifelse(true_vals == 0, NA_real_, abs_error / abs(true_vals))

  within_tol <- if (isTRUE(relative)) {
    rel_error <= tolerance
  } else {
    abs_error <= tolerance
  }

  out <- data.frame(
    parameter  = shared,
    true       = true_vals,
    estimate   = est_vals,
    abs_error  = abs_error,
    rel_error  = rel_error,
    within_tol = within_tol,
    stringsAsFactors = FALSE
  )

  attr(out, "tolerance") <- tolerance
  attr(out, "relative")  <- relative
  class(out) <- c("recovery_result", "data.frame")
  out
}


# ---------------------------------------------------------------------------
# print
# ---------------------------------------------------------------------------

#' @rdname validate_recovery
#' @param x A \code{recovery_result} object.
#' @param ... Further arguments (ignored).
#' @export
print.recovery_result <- function(x, ...) {
  n_params  <- nrow(x)
  n_within  <- sum(x$within_tol, na.rm = TRUE)
  pct       <- if (n_params > 0L) 100 * n_within / n_params else NA_real_
  tol       <- attr(x, "tolerance")
  rel       <- isTRUE(attr(x, "relative"))

  cat(sprintf(
    "Parameter recovery  (%s tolerance = %g)\n",
    if (rel) "relative" else "absolute", tol
  ))
  cat(sprintf(
    "  %d parameters | %d within tolerance (%.1f%%)\n",
    n_params, n_within, pct
  ))
  cat(sprintf(
    "  mean abs error = %.4g | mean rel error = %.4g\n\n",
    mean(x$abs_error), mean(x$rel_error, na.rm = TRUE)
  ))

  print(as.data.frame(x), row.names = FALSE)
  invisible(x)
}


# ---------------------------------------------------------------------------
# summary -- one-row scorecard
# ---------------------------------------------------------------------------

#' @rdname validate_recovery
#' @param object A \code{recovery_result} object.
#' @export
summary.recovery_result <- function(object, ...) {
  n_params <- nrow(object)
  n_within <- sum(object$within_tol, na.rm = TRUE)

  data.frame(
    n_params       = n_params,
    n_within_tol   = n_within,
    pct_within_tol = if (n_params > 0L) 100 * n_within / n_params else NA_real_,
    mean_abs_error = mean(object$abs_error),
    mean_rel_error = mean(object$rel_error, na.rm = TRUE),
    max_abs_error  = max(object$abs_error),
    stringsAsFactors = FALSE
  )
}


# ---------------------------------------------------------------------------
# Ensure S3 dispatch works even before NAMESPACE is regenerated. This is a
# no-op once roxygen2 writes the S3method() entries; it keeps print/summary
# dispatching correctly under devtools::load_all() in the meantime.
# ---------------------------------------------------------------------------

registerS3method("print", "recovery_result", print.recovery_result)
registerS3method("summary", "recovery_result", summary.recovery_result)
