#' Partial Correlation Network via EBICglasso
#'
#' @description
#' Estimates a sparse partial correlation network from frequency data or a
#' correlation/covariance matrix using the graphical lasso with EBIC model
#' selection. The graphical lasso is fit across a path of penalty values
#' (lambda), and the model with the lowest Extended Bayesian Information
#' Criterion (EBIC) is selected.
#'
#' @param data A data frame of per-sequence action frequencies (e.g., output of
#'   \code{\link{convert_sequence_format}} with \code{format = "frequency"}),
#'   or a square symmetric correlation/covariance matrix.
#' @param id_col Character vector. Name(s) of ID column(s) to exclude when
#'   \code{data} is a data frame. If NULL, columns named \code{"rid"} and
#'   non-numeric columns are automatically excluded. Default: NULL.
#' @param n Integer. Sample size. Required when \code{data} is a matrix.
#'   Ignored when \code{data} is a data frame (n is taken from \code{nrow()}).
#'   Default: NULL.
#' @param gamma Numeric. EBIC tuning parameter controlling sparsity preference.
#'   Higher values (e.g., 0.5) favour sparser networks; 0 reduces to BIC.
#'   Default: 0.5.
#' @param nlambda Integer. Number of lambda values in the regularisation path.
#'   Default: 100.
#' @param lambda.min.ratio Numeric. Ratio of the smallest to the largest lambda
#'   value. Default: 0.01.
#' @param penalize.diagonal Logical. Whether to penalise the diagonal of the
#'   precision matrix. Default: FALSE.
#' @param threshold Numeric. Absolute values in the partial correlation matrix
#'   below this are set to zero. Default: 1e-4.
#' @param cor_method Character. Correlation method when computing from a data
#'   frame: \code{"pearson"}, \code{"spearman"}, or \code{"kendall"}.
#'   Default: \code{"pearson"}.
#' @param input_type Character. How to interpret a matrix input:
#'   \code{"auto"} (detect from diagonal), \code{"cor"}, or \code{"cov"}.
#'   Default: \code{"auto"}.
#'
#' @return An object of class \code{"pcor_network"}, a list containing:
#' \describe{
#'   \item{pcor_matrix}{Sparse partial correlation matrix (zero-diagonal).}
#'   \item{precision_matrix}{Estimated precision (inverse covariance) matrix.}
#'   \item{cor_matrix}{Correlation matrix used as input.}
#'   \item{edges}{Data frame of non-zero edges (from, to, weight).}
#'   \item{lambda_selected}{The lambda value of the selected model.}
#'   \item{ebic_selected}{The EBIC value of the selected model.}
#'   \item{lambda_path}{Numeric vector of all lambda values evaluated.}
#'   \item{ebic_path}{Numeric vector of EBIC values for each lambda.}
#'   \item{n}{Sample size.}
#'   \item{p}{Number of variables.}
#'   \item{gamma}{EBIC gamma used.}
#'   \item{n_edges}{Number of non-zero edges.}
#' }
#'
#' @details
#' When \code{data} is a data frame, the function automatically excludes
#' non-numeric columns, ID columns (\code{id_col} and \code{"rid"}), columns
#' with non-syntactic names (e.g., \code{"\%"}), zero-variance columns, and
#' columns that are entirely \code{NA}. Rows with any \code{NA} in the
#' remaining columns are dropped before computing correlations.
#'
#' The estimation proceeds as follows:
#' \enumerate{
#'   \item Compute the correlation matrix from frequency data (or use supplied
#'     matrix directly).
#'   \item Generate a logarithmically spaced path of lambda values.
#'   \item For each lambda, fit the graphical lasso via \code{glasso::glasso()}
#'     using warm starts.
#'   \item Compute the EBIC for each fit and select the model with the lowest
#'     EBIC.
#'   \item Convert the selected precision matrix to partial correlations.
#' }
#'
#' @examples
#' \dontrun{
#' library(tna)
#'
#' # From wide sequence data
#' freq <- convert_sequence_format(group_regulation, format = "frequency")
#' net <- pcor_network(freq)
#' print(net)
#'
#' # From a correlation matrix
#' cor_mat <- cor(freq[, -c(1, 2)])
#' net2 <- pcor_network(cor_mat, n = nrow(freq))
#' }
#'
#' @seealso \code{\link{convert_sequence_format}} for producing frequency data.
#'
#' @importFrom stats cor cov complete.cases
#' @export
pcor_network <- function(data,
                         id_col = NULL,
                         n = NULL,
                         gamma = 0.5,
                         nlambda = 100L,
                         lambda.min.ratio = 0.01,
                         penalize.diagonal = FALSE,
                         threshold = 1e-4,
                         cor_method = c("pearson", "spearman", "kendall"),
                         input_type = c("auto", "cor", "cov")) {
  cor_method <- match.arg(cor_method)
  input_type <- match.arg(input_type)
  stopifnot(is.numeric(gamma), length(gamma) == 1, gamma >= 0)
  stopifnot(is.numeric(nlambda), length(nlambda) == 1, nlambda >= 2)
  stopifnot(is.numeric(lambda.min.ratio), lambda.min.ratio > 0,
            lambda.min.ratio < 1)
  stopifnot(is.logical(penalize.diagonal), length(penalize.diagonal) == 1)
  stopifnot(is.numeric(threshold), threshold >= 0)

  # Prepare input: get correlation matrix S and sample size n
  prepared <- .prepare_pcor_input(
    data, id_col = id_col, n = n,
    cor_method = cor_method, input_type = input_type
  )
  S <- prepared$S
  n <- prepared$n
  p <- ncol(S)

  # Compute lambda path
  lambda_path <- .compute_lambda_path(S, nlambda, lambda.min.ratio)

  # Select best model via EBIC
  selected <- .select_ebic(S, lambda_path, n, gamma, penalize.diagonal)

  # Convert precision matrix to partial correlations
  pcor <- .precision_to_pcor(selected$wi, threshold)
  colnames(pcor) <- rownames(pcor) <- colnames(S)

  # Build edge data frame
  edges <- .pcor_to_edges(pcor)

  structure(
    list(
      pcor_matrix      = pcor,
      precision_matrix = selected$wi,
      cor_matrix       = S,
      edges            = edges,
      lambda_selected  = selected$lambda,
      ebic_selected    = selected$ebic,
      lambda_path      = lambda_path,
      ebic_path        = selected$ebic_path,
      n                = n,
      p                = p,
      gamma            = gamma,
      n_edges          = nrow(edges)
    ),
    class = "pcor_network"
  )
}


# ---- Input preparation ----

#' Validate and prepare input for pcor_network
#' @noRd
.prepare_pcor_input <- function(data, id_col, n, cor_method, input_type) {
  if (is.data.frame(data)) {
    # Exclude id columns, "rid", and non-numeric columns
    exclude <- c(id_col, "rid")
    numeric_cols <- vapply(data, is.numeric, logical(1))
    keep <- setdiff(names(data)[numeric_cols], exclude)

    # Drop columns with non-syntactic names (e.g. "%", "*", "NA")
    syntactic <- make.names(keep) == keep
    if (any(!syntactic)) {
      dropped <- keep[!syntactic]
      message("Dropping non-syntactic columns: ",
              paste(dropped, collapse = ", "))
      keep <- keep[syntactic]
    }

    if (length(keep) < 2) {
      stop("At least 2 numeric columns are required after cleaning.")
    }

    mat <- as.matrix(data[, keep, drop = FALSE])

    # Drop all-NA columns
    all_na <- apply(mat, 2, function(x) all(is.na(x)))
    if (any(all_na)) {
      message("Dropping all-NA columns: ",
              paste(keep[all_na], collapse = ", "))
      mat <- mat[, !all_na, drop = FALSE]
    }

    # Drop rows with any NA
    complete <- complete.cases(mat)
    if (!all(complete)) {
      n_dropped <- sum(!complete)
      message("Dropping ", n_dropped, " rows with NA values.")
      mat <- mat[complete, , drop = FALSE]
    }

    if (nrow(mat) < 3) {
      stop("Fewer than 3 complete rows remain after removing NAs.")
    }

    # Drop zero-variance columns silently
    col_vars <- apply(mat, 2, stats::var)
    zero_var <- colnames(mat)[col_vars == 0]
    if (length(zero_var) > 0) {
      message("Dropping zero-variance columns: ",
              paste(zero_var, collapse = ", "))
      mat <- mat[, col_vars > 0, drop = FALSE]
    }

    if (ncol(mat) < 2) {
      stop("At least 2 variable columns are required after cleaning.")
    }

    n <- nrow(mat)
    S <- cor(mat, method = cor_method)

  } else if (is.matrix(data)) {
    stopifnot(nrow(data) == ncol(data))
    if (!isSymmetric(unname(data), tol = 1e-8)) {
      stop("Matrix input must be symmetric.")
    }
    if (is.null(n)) {
      stop("Sample size 'n' is required when data is a matrix.")
    }
    stopifnot(is.numeric(n), length(n) == 1, n > 0)

    if (input_type == "auto") {
      diag_vals <- diag(data)
      input_type <- if (all(abs(diag_vals - 1) < 1e-8)) "cor" else "cov"
    }

    if (input_type == "cov") {
      d <- sqrt(diag(data))
      S <- data / outer(d, d)
    } else {
      S <- data
    }

    if (is.null(colnames(S))) {
      colnames(S) <- rownames(S) <- paste0("V", seq_len(ncol(S)))
    }
  } else {
    stop("data must be a data frame or a square symmetric matrix.")
  }

  list(S = S, n = n)
}


# ---- Lambda path ----

#' Compute log-spaced lambda path
#' @noRd
.compute_lambda_path <- function(S, nlambda, lambda.min.ratio) {
  lambda_max <- max(abs(S[upper.tri(S)]))
  if (lambda_max <= 0) {
    stop("All off-diagonal correlations are zero; nothing to regularise.")
  }
  lambda_min <- lambda_max * lambda.min.ratio
  exp(seq(log(lambda_max), log(lambda_min), length.out = nlambda))
}


# ---- EBIC model selection ----

#' Select best lambda via EBIC using glasso fits with warm starts
#' @noRd
.select_ebic <- function(S, lambda_path, n, gamma, penalize_diagonal) {
  p <- ncol(S)
  n_lambda <- length(lambda_path)
  ebic_vals <- numeric(n_lambda)

  # Initial fit (cold start)
  w_prev <- NULL
  wi_prev <- NULL
  best_idx <- 1L
  best_ebic <- Inf
  best_wi <- NULL

  for (i in seq_along(lambda_path)) {
    lam <- lambda_path[i]

    fit <- tryCatch(
      glasso::glasso(
        s = S,
        rho = lam,
        penalize.diagonal = penalize_diagonal,
        start = if (is.null(w_prev)) "cold" else "warm",
        w.init = w_prev,
        wi.init = wi_prev,
        trace = FALSE
      ),
      error = function(e) NULL
    )

    if (is.null(fit)) {
      ebic_vals[i] <- Inf
      next
    }

    # Warm start for next iteration
    w_prev <- fit$w
    wi_prev <- fit$wi

    # EBIC calculation
    log_det <- determinant(fit$wi, logarithm = TRUE)
    if (log_det$sign <= 0) {
      ebic_vals[i] <- Inf
      next
    }
    log_det_val <- as.numeric(log_det$modulus)

    loglik <- (n / 2) * (log_det_val - sum(diag(S %*% fit$wi)))
    npar <- sum(abs(fit$wi[upper.tri(fit$wi)]) > 1e-10)
    ebic_vals[i] <- -2 * loglik + npar * log(n) +
      4 * npar * gamma * log(p)

    if (ebic_vals[i] < best_ebic) {
      best_ebic <- ebic_vals[i]
      best_idx <- i
      best_wi <- fit$wi
    }
  }

  if (is.null(best_wi)) {
    stop("All glasso fits failed. Check your input data.")
  }

  colnames(best_wi) <- rownames(best_wi) <- colnames(S)

  list(
    wi        = best_wi,
    lambda    = lambda_path[best_idx],
    ebic      = best_ebic,
    ebic_path = ebic_vals
  )
}


# ---- Partial correlation conversion ----

#' Convert precision matrix to partial correlations
#' @noRd
.precision_to_pcor <- function(Wi, threshold) {
  d <- sqrt(diag(Wi))
  pcor <- -Wi / outer(d, d)
  diag(pcor) <- 0
  pcor[abs(pcor) < threshold] <- 0
  pcor
}


# ---- Edge extraction ----

#' Extract non-zero edges from partial correlation matrix
#' @noRd
.pcor_to_edges <- function(pcor) {
  idx <- which(upper.tri(pcor) & pcor != 0, arr.ind = TRUE)
  if (nrow(idx) == 0) {
    return(data.frame(
      from = character(0), to = character(0),
      weight = numeric(0), stringsAsFactors = FALSE
    ))
  }
  nms <- colnames(pcor)
  data.frame(
    from   = nms[idx[, 1]],
    to     = nms[idx[, 2]],
    weight = pcor[idx],
    stringsAsFactors = FALSE
  )
}


# ---- S3 methods ----

#' Print Method for Partial Correlation Network
#'
#' @param x A \code{pcor_network} object.
#' @param ... Additional arguments (ignored).
#'
#' @export
print.pcor_network <- function(x, ...) {
  total_possible <- x$p * (x$p - 1) / 2
  pct <- if (total_possible > 0) {
    sprintf("%.1f%%", 100 * x$n_edges / total_possible)
  } else {
    "0%"
  }
  cat("Partial Correlation Network (EBICglasso)\n")
  cat(sprintf("  Variables: %d  |  Sample size: %d\n", x$p, x$n))
  cat(sprintf(
    "  Non-zero edges: %d / %d (%s)\n",
    x$n_edges, total_possible, pct
  ))
  cat(sprintf(
    "  Gamma: %.2f  |  Lambda: %.4f\n",
    x$gamma, x$lambda_selected
  ))
  invisible(x)
}


#' Plot Method for Partial Correlation Network
#'
#' @description
#' Plots the partial correlation network using \code{cograph::tplot()}.
#' Requires the \pkg{cograph} package to be installed.
#'
#' @param x A \code{pcor_network} object.
#' @param ... Additional arguments passed to \code{cograph::tplot()}.
#'
#' @export
plot.pcor_network <- function(x, ...) {
  if (!requireNamespace("cograph", quietly = TRUE)) {
    stop(
      "Package 'cograph' is required for plotting. ",
      "Install it with: install.packages('cograph')"
    )
  }
  cograph::tplot(x$pcor_matrix, directed = FALSE, ...)
}
