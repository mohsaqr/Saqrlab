# ===========================================================================
# saqr_sim S3 class -- unified wrapper for simulation results
# ===========================================================================

#' Create a saqr_sim Object
#'
#' @description Wraps simulation output in a standardised S3 class. Every
#'   simulation function that returns numerical data uses this class, giving
#'   a consistent interface: \code{$data} (the data.frame), \code{$params}
#'   (ground-truth generating parameters), \code{$type}, \code{$seed}.
#'
#'   Backward-compatible: \code{$data} and \code{$params} still work exactly
#'   as before. Additionally, \code{[} and \code{as.data.frame()} delegate to
#'   the underlying data, so existing code that subscripts or coerces the
#'   result keeps working.
#'
#' @param data A data.frame (or matrix) of simulated data.
#' @param params A list of ground-truth generating parameters.
#' @param type Character label identifying the simulation type
#'   (e.g. \code{"lpa"}, \code{"regression"}, \code{"longitudinal"}).
#' @param seed The random seed actually used (integer or NULL).
#' @param call The matched call that created this object (optional).
#' @param extras Named list of additional fields to attach (optional).
#'
#' @return An object of class \code{saqr_sim} (inherits from \code{list}).
#'
#' @examples
#' sim <- saqr_sim(
#'   data = data.frame(x = rnorm(10), y = rnorm(10)),
#'   params = list(mean_x = 0, mean_y = 0),
#'   type = "demo", seed = 1
#' )
#' sim
#' names(sim)
#' head(sim)
#' dim(sim)
#' sim$params
#'
#' @keywords internal
#' @export
saqr_sim <- function(data, params, type, seed = NULL, call = NULL,
                     extras = list()) {
  stopifnot(is.data.frame(data) || is.matrix(data))
  stopifnot(is.list(params))
  stopifnot(is.character(type), length(type) == 1L)

  obj <- c(
    list(data = data, params = params, type = type, seed = seed, call = call),
    extras
  )
  class(obj) <- c("saqr_sim", "list")
  obj
}


# ---------------------------------------------------------------------------
# print
# ---------------------------------------------------------------------------

#' @export
print.saqr_sim <- function(x, ...) {
  type_label <- x$type
  d <- x$data
  n_row <- if (is.data.frame(d)) nrow(d) else NROW(d)
  n_col <- if (is.data.frame(d)) ncol(d) else NCOL(d)

  cat(sprintf("saqr_sim [%s]  %d x %d", type_label, n_row, n_col))
  if (!is.null(x$seed)) cat(sprintf("  (seed=%d)", x$seed))
  cat("\n")

  # Show param names

  pnames <- names(x$params)
  if (length(pnames) > 0L) {
    cat("  params:", paste(pnames, collapse = ", "), "\n")
  }

  # Show column names (truncated)
  cnames <- if (is.data.frame(d)) names(d) else colnames(d)
  if (length(cnames) > 8L) {
    cat("  cols:  ", paste(c(cnames[1:6], "...",
                             cnames[length(cnames)]), collapse = ", "), "\n")
  } else if (length(cnames) > 0L) {
    cat("  cols:  ", paste(cnames, collapse = ", "), "\n")
  }

  invisible(x)
}


# ---------------------------------------------------------------------------
# summary
# ---------------------------------------------------------------------------

#' @export
summary.saqr_sim <- function(object, ...) {
  cat(sprintf("Simulation type: %s\n", object$type))
  if (!is.null(object$seed)) cat(sprintf("Seed: %d\n", object$seed))
  cat("\n--- Data ---\n")
  print(summary(object$data))
  cat("\n--- Parameters ---\n")
  str(object$params, max.level = 1L)
  invisible(object)
}


# ---------------------------------------------------------------------------
# as.data.frame -- extract $data
# ---------------------------------------------------------------------------

#' @export
as.data.frame.saqr_sim <- function(x, ...) {
  x$data
}


# ---------------------------------------------------------------------------
# [ -- delegate to $data so result[1:5, ] works
# ---------------------------------------------------------------------------

#' @export
`[.saqr_sim` <- function(x, ...) {
  x$data[...]
}


# ---------------------------------------------------------------------------
# dim / nrow / ncol -- delegate to $data
# ---------------------------------------------------------------------------

#' @export
dim.saqr_sim <- function(x) {
  dim(x$data)
}

#' @export
names.saqr_sim <- function(x) {
  # Return list element names (data, params, type, seed, call, ...)
  # This preserves backward compat with $data, $params access
  NextMethod()
}

#' @export
head.saqr_sim <- function(x, n = 6L, ...) {
  head(x$data, n = n, ...)
}

#' @export
tail.saqr_sim <- function(x, n = 6L, ...) {
  tail(x$data, n = n, ...)
}

#' @export
str.saqr_sim <- function(object, ...) {
  cat(sprintf("saqr_sim [%s]", object$type))
  if (!is.null(object$seed)) cat(sprintf("  seed=%d", object$seed))
  cat("\n")
  cat(" $ data  :")
  str(object$data, indent.str = "   ", ...)
  cat(" $ params:")
  str(object$params, max.level = 1L, indent.str = "   ", ...)
  invisible(NULL)
}
