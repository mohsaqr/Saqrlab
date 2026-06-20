# ===========================================================================
# Missing-data-mechanism injection (MCAR / MAR / MNAR)
# A transformer: takes a data.frame and returns it with NAs injected, plus
# an attached "missing_info" attribute describing what was done.
# ===========================================================================


#' Inject Missing Values Under a Known Mechanism (MCAR / MAR / MNAR)
#'
#' @description Take an existing complete data.frame and inject missing values
#'   (\code{NA}) into target columns under a chosen, well-defined missingness
#'   mechanism. This is a *transformer*, not a generator: the returned object
#'   is the same data.frame with some cells replaced by \code{NA}, carrying an
#'   attached \code{"missing_info"} attribute that records exactly which cells
#'   were made missing and how. Useful for testing imputation methods and
#'   missing-data-aware estimators against a known ground truth.
#'
#'   The three mechanisms follow the standard Rubin taxonomy:
#'   \describe{
#'     \item{\code{"MCAR"}}{Missing Completely At Random. Each cell of the
#'       target columns is independently set missing with probability
#'       \code{prop}.}
#'     \item{\code{"MAR"}}{Missing At Random. The probability that a target
#'       cell is missing depends on an *observed* \code{predictor} column,
#'       never on the to-be-missing value itself. Probabilities are calibrated
#'       (via a rank/logistic weighting on the predictor) so the average
#'       missing fraction is approximately \code{prop}.}
#'     \item{\code{"MNAR"}}{Missing Not At Random. The probability that a
#'       target cell is missing depends on that cell's *own* (numeric) value:
#'       larger values are more likely to go missing. Calibrated so the
#'       average missing fraction per column is approximately \code{prop}.}
#'   }
#'
#' @param data A data.frame, or a \code{\link{saqr_sim}} object (its
#'   \code{$data} is used). The complete data into which to inject missingness.
#' @param mechanism Character. One of \code{"MCAR"}, \code{"MAR"},
#'   \code{"MNAR"} (matched via \code{\link[base]{match.arg}}).
#' @param prop Numeric scalar in \code{[0, 1]}. Target average fraction of
#'   cells in the target columns to set missing. Default: \code{0.1}.
#' @param cols Character vector of target column names, or \code{NULL} (the
#'   default) meaning all columns. For \code{"MAR"}, the \code{predictor}
#'   column is automatically excluded from the targets.
#' @param predictor Character scalar naming the observed column that drives
#'   missingness under \code{"MAR"}. Required for \code{"MAR"}, ignored
#'   otherwise.
#' @param seed Integer or \code{NULL}. If non-\code{NULL}, \code{set.seed()} is
#'   called for reproducibility.
#'
#' @return The input data.frame with \code{NA}s injected into the target
#'   columns. The returned object additionally carries an attribute
#'   \code{"missing_info"}, a list with elements:
#'   \describe{
#'     \item{\code{mechanism}}{The resolved mechanism string.}
#'     \item{\code{prop}}{The requested target proportion.}
#'     \item{\code{cols}}{The target column names actually used.}
#'     \item{\code{n_missing}}{Total number of cells set missing.}
#'     \item{\code{realized_prop}}{Realized fraction of target cells set
#'       missing (\code{n_missing / (nrow * length(cols))}).}
#'     \item{\code{indicator}}{A logical matrix (\code{nrow(data)} x
#'       \code{length(cols)}) marking the cells that were injected as
#'       \code{TRUE}, with column names matching \code{cols}.}
#'   }
#'
#' @examples
#' df <- simulate_prediction(
#'   n = 500,
#'   coefs = c("(Intercept)" = 1, x1 = 2, x2 = -1),
#'   seed = 1
#' )$data
#'
#' # MCAR: 20% of every column missing at random
#' mcar <- inject_missingness(df, mechanism = "MCAR", prop = 0.2, seed = 7)
#' attr(mcar, "missing_info")$realized_prop
#' str(attr(mcar, "missing_info"))
#'
#' # MAR: missingness in y driven by observed x1
#' mar <- inject_missingness(df, mechanism = "MAR", prop = 0.2,
#'                           cols = "y", predictor = "x1", seed = 7)
#' colSums(is.na(mar))
#'
#' # MNAR: large y-values more likely missing
#' mnar <- inject_missingness(df, mechanism = "MNAR", prop = 0.2,
#'                            cols = "y", seed = 7)
#' attr(mnar, "missing_info")$n_missing
#'
#' @export
inject_missingness <- function(data,
                               mechanism = c("MCAR", "MAR", "MNAR"),
                               prop = 0.1,
                               cols = NULL,
                               predictor = NULL,
                               seed = NULL) {
  # Accept saqr_sim by operating on its underlying data.frame
  if (inherits(data, "saqr_sim")) data <- data$data

  mechanism <- match.arg(mechanism)

  stopifnot(
    is.data.frame(data),
    is.numeric(prop), length(prop) == 1L, prop >= 0, prop <= 1,
    is.null(cols) || is.character(cols),
    is.null(predictor) || (is.character(predictor) &&
                             length(predictor) == 1L)
  )

  all_cols <- names(data)

  # Resolve target columns
  if (is.null(cols)) {
    target_cols <- all_cols
  } else {
    stopifnot(all(cols %in% all_cols))
    target_cols <- cols
  }

  # MAR requires an observed predictor; it is never itself a target
  if (mechanism == "MAR") {
    stopifnot(
      !is.null(predictor),
      predictor %in% all_cols
    )
    target_cols <- setdiff(target_cols, predictor)
    stopifnot(length(target_cols) >= 1L)
  }

  stopifnot(length(target_cols) >= 1L)

  if (!is.null(seed)) set.seed(as.integer(seed))

  n <- nrow(data)

  # Build a per-column logical "make missing" indicator using the chosen
  # mechanism. lapply over columns (no for-loops).
  ind_list <- switch(
    mechanism,
    MCAR = lapply(target_cols, function(nm) {
      stats::runif(n) < prop
    }),
    MAR = {
      pred_vals <- data[[predictor]]
      pred_w <- .missing_weight_from(pred_vals)
      lapply(target_cols, function(nm) {
        .bernoulli_calibrated(pred_w, prop)
      })
    },
    MNAR = lapply(target_cols, function(nm) {
      col_w <- .missing_weight_from(data[[nm]])
      .bernoulli_calibrated(col_w, prop)
    })
  )

  indicator <- matrix(unlist(ind_list, use.names = FALSE),
                      nrow = n, ncol = length(target_cols))
  colnames(indicator) <- target_cols

  # Apply the NAs, one column at a time (lapply, no loops)
  data[target_cols] <- Map(function(col, flag) {
    col[flag] <- NA
    col
  }, data[target_cols], split(indicator, col(indicator)))

  n_missing <- sum(indicator)

  attr(data, "missing_info") <- list(
    mechanism     = mechanism,
    prop          = prop,
    cols          = target_cols,
    n_missing     = n_missing,
    realized_prop = n_missing / (n * length(target_cols)),
    indicator     = indicator
  )

  data
}


# ---------------------------------------------------------------------------
# Internal helpers (live inside the package, hidden from the public surface)
# ---------------------------------------------------------------------------

# Map an observed variable to a [0, 1) "propensity weight" in which larger
# values yield larger weights. Uses ranks so it is robust to scale and works
# for numeric, ordered, and (coerced) categorical inputs. Ties and NAs are
# handled gracefully.
.missing_weight_from <- function(x) {
  v <- if (is.numeric(x)) x else as.numeric(as.factor(x))
  if (all(is.na(v)) || stats::sd(v, na.rm = TRUE) == 0) {
    return(rep(0.5, length(v)))
  }
  r <- rank(v, na.last = "keep", ties.method = "average")
  w <- (r - 0.5) / sum(!is.na(r))
  w[is.na(w)] <- 0.5
  w
}

# Given a vector of propensity weights in [0, 1] and a target average
# probability `prop`, return a logical "make missing" vector whose expected
# fraction is `prop`. Calibration is CLAMP-AWARE: it finds the multiplier `m`
# such that mean(pmin(m * weights, 1)) == prop, so the expected missingness
# stays on target even when some probabilities hit the [0, 1] ceiling (a plain
# prop/mean scale undershoots once clamping kicks in, e.g. for prop > 0.5).
.bernoulli_calibrated <- function(weights, prop) {
  if (prop <= 0) return(rep(FALSE, length(weights)))
  if (prop >= 1) return(rep(TRUE, length(weights)))
  w <- pmax(weights, 0)
  mw <- mean(w)
  if (mw <= 0) {
    p <- rep(prop, length(w))
  } else {
    f  <- function(m) mean(pmin(m * w, 1)) - prop
    hi <- prop / mw
    while (f(hi) < 0 && hi < 1e9) hi <- hi * 2
    m  <- if (f(hi) >= 0) stats::uniroot(f, c(0, hi))$root else hi
    p  <- pmin(m * w, 1)
  }
  stats::runif(length(p)) < p
}
