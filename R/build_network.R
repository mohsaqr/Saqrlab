#' Build a Psychological Network
#'
#' @description
#' Estimates a network from frequency data or a correlation/covariance matrix.
#' Supports three methods:
#' \itemize{
#'   \item \code{"glasso"} (aliases: \code{"ebicglasso"}, \code{"regularized"}):
#'     Sparse partial correlations via graphical lasso with EBIC model selection.
#'   \item \code{"pcor"} (alias: \code{"partial"}):
#'     Unregularised partial correlations (inverts the correlation matrix
#'     directly; requires p < n).
#'   \item \code{"cor"} (alias: \code{"correlation"}):
#'     Pairwise correlation network.
#' }
#'
#' @param data A data frame of per-sequence action frequencies (e.g., output of
#'   \code{\link{convert_sequence_format}} with \code{format = "frequency"}),
#'   or a square symmetric correlation/covariance matrix.
#' @param method Character. Network estimation method. One of \code{"glasso"},
#'   \code{"ebicglasso"}, \code{"regularized"} (all equivalent),
#'   \code{"pcor"}, \code{"partial"} (equivalent),
#'   \code{"cor"}, \code{"correlation"} (equivalent).
#'   Default: \code{"glasso"}.
#' @param id_col Character vector. Name(s) of ID column(s) to exclude when
#'   \code{data} is a data frame. If NULL, columns named \code{"rid"} and
#'   non-numeric columns are automatically excluded. Default: NULL.
#' @param level Character or NULL. Decomposition level for multilevel networks
#'   when data has repeated measures per person. One of:
#'   \itemize{
#'     \item \code{NULL} (default): no decomposition, current behavior.
#'     \item \code{"between"}: aggregate to person means, estimate network
#'       of trait-level associations.
#'     \item \code{"within"}: person-mean center each variable, estimate
#'       network of state-level (occasion-to-occasion) associations.
#'     \item \code{"both"}: return both between and within networks in a
#'       \code{psych_network_ml} object.
#'   }
#'   Requires \code{id_col} and a data frame input.
#' @param n Integer. Sample size. Required when \code{data} is a matrix.
#'   Ignored when \code{data} is a data frame (n is taken from \code{nrow()}).
#'   Default: NULL.
#' @param gamma Numeric. EBIC tuning parameter controlling sparsity preference
#'   (only used when \code{method = "glasso"}). Higher values favour sparser
#'   networks; 0 reduces to BIC. Default: 0.5.
#' @param nlambda Integer. Number of lambda values in the regularisation path
#'   (only used when \code{method = "glasso"}). Default: 100.
#' @param lambda.min.ratio Numeric. Ratio of the smallest to the largest lambda
#'   value (only used when \code{method = "glasso"}). Default: 0.01.
#' @param penalize.diagonal Logical. Whether to penalise the diagonal of the
#'   precision matrix (only used when \code{method = "glasso"}).
#'   Default: FALSE.
#' @param threshold Numeric. Absolute values in the result matrix below this
#'   are set to zero. For \code{"glasso"} and \code{"pcor"} this applies to
#'   partial correlations; for \code{"cor"} it applies to correlations.
#'   Default: 1e-4.
#' @param cor_method Character. Correlation method when computing from a data
#'   frame: \code{"pearson"}, \code{"spearman"}, or \code{"kendall"}.
#'   Default: \code{"pearson"}.
#' @param input_type Character. How to interpret a matrix input:
#'   \code{"auto"} (detect from diagonal), \code{"cor"}, or \code{"cov"}.
#'   Default: \code{"auto"}.
#'
#' @return An object of class \code{"psych_network"} (or
#'   \code{"psych_network_ml"} when \code{level = "both"}), a list containing:
#' \describe{
#'   \item{network_matrix}{The estimated network weight matrix (partial
#'     correlations for \code{"glasso"} and \code{"pcor"}; correlations for
#'     \code{"cor"}). Diagonal is zero.}
#'   \item{cor_matrix}{Correlation matrix used as input.}
#'   \item{edges}{Data frame of non-zero edges (from, to, weight).}
#'   \item{n}{Sample size.}
#'   \item{p}{Number of variables.}
#'   \item{method}{The resolved method name (\code{"glasso"}, \code{"pcor"},
#'     or \code{"cor"}).}
#'   \item{n_edges}{Number of non-zero edges.}
#'   \item{level}{The decomposition level used (\code{NULL}, \code{"between"},
#'     or \code{"within"}).}
#' }
#'
#' Additional fields for \code{method = "glasso"}:
#' \describe{
#'   \item{precision_matrix}{Estimated precision (inverse covariance) matrix.}
#'   \item{lambda_selected}{The lambda value of the selected model.}
#'   \item{ebic_selected}{The EBIC value of the selected model.}
#'   \item{lambda_path}{Numeric vector of all lambda values evaluated.}
#'   \item{ebic_path}{Numeric vector of EBIC values for each lambda.}
#'   \item{gamma}{EBIC gamma used.}
#' }
#'
#' Additional fields for \code{method = "pcor"}:
#' \describe{
#'   \item{precision_matrix}{Precision (inverse correlation) matrix.}
#' }
#'
#' When \code{level = "both"}, a \code{"psych_network_ml"} object is returned
#' containing:
#' \describe{
#'   \item{between}{A \code{psych_network} object for the between-person level.}
#'   \item{within}{A \code{psych_network} object for the within-person level.}
#'   \item{method}{The resolved method name.}
#' }
#'
#' @details
#' When \code{data} is a data frame, the function automatically excludes
#' non-numeric columns, ID columns (\code{id_col} and \code{"rid"}), columns
#' with non-syntactic names (e.g., \code{"\%"}), zero-variance columns, and
#' columns that are entirely \code{NA}. Rows with any \code{NA} in the
#' remaining columns are dropped before computing correlations.
#'
#' \strong{Method details:}
#' \describe{
#'   \item{glasso}{Fits the graphical lasso across a log-spaced lambda path
#'     using \code{glasso::glasso()} with warm starts, selects the model with
#'     the lowest EBIC, and converts the precision matrix to partial
#'     correlations.}
#'   \item{pcor}{Directly inverts the correlation matrix to obtain the
#'     precision matrix, then converts to partial correlations. Requires
#'     \code{p < n} (more observations than variables).}
#'   \item{cor}{Returns the pairwise correlation matrix as the network.
#'     Values below \code{threshold} are zeroed out.}
#' }
#'
#' \strong{Multilevel decomposition:}
#' When \code{level} is specified, the data is decomposed before estimation:
#' \describe{
#'   \item{between}{Rows are aggregated to person means (using the first
#'     element of \code{id_col} as the grouping variable), then the network
#'     is estimated from the aggregated data.}
#'   \item{within}{Each variable is person-mean centered (group-mean centering).
#'     Persons with only one observation are dropped (no within-person
#'     variance). The network is estimated from the pooled centered residuals.}
#'   \item{both}{Both between and within networks are estimated and returned
#'     in a \code{psych_network_ml} object.}
#' }
#'
#' @examples
#' \dontrun{
#' library(tna)
#'
#' # From wide sequence data — regularised (default)
#' freq <- convert_sequence_format(group_regulation, format = "frequency")
#' net <- build_network(freq)
#' print(net)
#'
#' # Unregularised partial correlations
#' net_pcor <- build_network(freq, method = "pcor")
#'
#' # Correlation network
#' net_cor <- build_network(freq, method = "cor", threshold = 0.1)
#'
#' # From a correlation matrix
#' S <- cor(freq[, -c(1, 2)])
#' net2 <- build_network(S, n = nrow(freq), method = "glasso")
#'
#' # Multilevel: between-person network
#' freq_long <- convert_sequence_format(group_regulation_long,
#'   action = "Action", id_col = "Actor", format = "frequency")
#' net_between <- build_network(freq_long, id_col = "Actor",
#'   level = "between")
#'
#' # Multilevel: both levels
#' net_ml <- build_network(freq_long, id_col = "Actor", level = "both")
#' print(net_ml)
#' }
#'
#' @seealso \code{\link{convert_sequence_format}} for producing frequency data.
#'
#' @importFrom stats aggregate ave cor cov complete.cases var
#' @export
build_network <- function(data,
                          method = c("glasso", "ebicglasso", "regularized",
                                     "pcor", "partial",
                                     "cor", "correlation"),
                          id_col = NULL,
                          level = NULL,
                          n = NULL,
                          gamma = 0.5,
                          nlambda = 100L,
                          lambda.min.ratio = 0.01,
                          penalize.diagonal = FALSE,
                          threshold = 1e-4,
                          cor_method = c("pearson", "spearman", "kendall"),
                          input_type = c("auto", "cor", "cov")) {
  method <- match.arg(method)
  cor_method <- match.arg(cor_method)
  input_type <- match.arg(input_type)

  # Normalise method aliases
  method <- switch(method,
    ebicglasso  = ,
    regularized = "glasso",
    partial     = "pcor",
    correlation = "cor",
    method
  )

  stopifnot(is.numeric(threshold), threshold >= 0)

  # Validate level parameter
  if (!is.null(level)) {
    level <- match.arg(level, c("between", "within", "both"))
    if (is.null(id_col)) {
      stop("'id_col' is required when 'level' is specified.", call. = FALSE)
    }
    if (!is.data.frame(data)) {
      stop("'data' must be a data frame when 'level' is specified.",
           call. = FALSE)
    }
  }

  # level = "both": recursive dispatch

  if (identical(level, "both")) {
    between <- build_network(
      data, method = method, id_col = id_col, level = "between", n = n,
      gamma = gamma, nlambda = nlambda, lambda.min.ratio = lambda.min.ratio,
      penalize.diagonal = penalize.diagonal, threshold = threshold,
      cor_method = cor_method, input_type = input_type
    )
    within <- build_network(
      data, method = method, id_col = id_col, level = "within", n = n,
      gamma = gamma, nlambda = nlambda, lambda.min.ratio = lambda.min.ratio,
      penalize.diagonal = penalize.diagonal, threshold = threshold,
      cor_method = cor_method, input_type = input_type
    )
    result <- list(between = between, within = within, method = method)
    class(result) <- "psych_network_ml"
    return(result)
  }

  # Prepare input: get correlation matrix S and sample size n
  prepared <- .prepare_network_input(
    data, id_col = id_col, n = n,
    cor_method = cor_method, input_type = input_type,
    level = level
  )
  S <- prepared$S
  n <- prepared$n
  p <- ncol(S)

  # Dispatch to method-specific estimation
  result <- switch(method,
    glasso = .estimate_glasso(S, n, p, gamma, nlambda, lambda.min.ratio,
                              penalize.diagonal, threshold),
    pcor   = .estimate_pcor(S, p, threshold),
    cor    = .estimate_cor(S, threshold)
  )

  result$cor_matrix <- S
  result$n <- n
  result$p <- p
  result$method <- method
  result$level <- level

  structure(result, class = "psych_network")
}


# ---- Method: glasso (EBICglasso) ----

#' Regularised partial correlations via graphical lasso + EBIC
#' @noRd
.estimate_glasso <- function(S, n, p, gamma, nlambda, lambda.min.ratio,
                             penalize.diagonal, threshold) {
  stopifnot(is.numeric(gamma), length(gamma) == 1, gamma >= 0)
  stopifnot(is.numeric(nlambda), length(nlambda) == 1, nlambda >= 2)
  stopifnot(is.numeric(lambda.min.ratio), lambda.min.ratio > 0,
            lambda.min.ratio < 1)
  stopifnot(is.logical(penalize.diagonal), length(penalize.diagonal) == 1)

  lambda_path <- .compute_lambda_path(S, nlambda, lambda.min.ratio)
  selected <- .select_ebic(S, lambda_path, n, gamma, penalize.diagonal)

  pcor <- .precision_to_pcor(selected$wi, threshold)
  colnames(pcor) <- rownames(pcor) <- colnames(S)
  edges <- .network_to_edges(pcor)

  list(
    network_matrix   = pcor,
    precision_matrix = selected$wi,
    edges            = edges,
    lambda_selected  = selected$lambda,
    ebic_selected    = selected$ebic,
    lambda_path      = lambda_path,
    ebic_path        = selected$ebic_path,
    gamma            = gamma,
    n_edges          = nrow(edges)
  )
}


# ---- Method: pcor (unregularised) ----

#' Unregularised partial correlations via matrix inversion
#' @noRd
.estimate_pcor <- function(S, p, threshold) {
  Wi <- tryCatch(
    solve(S),
    error = function(e) {
      stop(
        "Correlation matrix is singular (p >= n or collinear variables). ",
        "Use method = 'glasso' for regularised estimation.",
        call. = FALSE
      )
    }
  )
  colnames(Wi) <- rownames(Wi) <- colnames(S)

  pcor <- .precision_to_pcor(Wi, threshold)
  colnames(pcor) <- rownames(pcor) <- colnames(S)
  edges <- .network_to_edges(pcor)

  list(
    network_matrix   = pcor,
    precision_matrix = Wi,
    edges            = edges,
    n_edges          = nrow(edges)
  )
}


# ---- Method: cor (correlation network) ----

#' Correlation network (threshold and return)
#' @noRd
.estimate_cor <- function(S, threshold) {
  net <- S
  diag(net) <- 0
  net[abs(net) < threshold] <- 0
  colnames(net) <- rownames(net) <- colnames(S)
  edges <- .network_to_edges(net)

  list(
    network_matrix = net,
    edges          = edges,
    n_edges        = nrow(edges)
  )
}


# ---- Input preparation ----

#' Validate and prepare input for build_network
#' @noRd
.prepare_network_input <- function(data, id_col, n, cor_method, input_type,
                                   level = NULL) {
  if (is.data.frame(data)) {
    # Extract id values before subsetting (needed for level decomposition)
    id_vals_raw <- if (!is.null(level) && !is.null(id_col)) {
      data[[id_col[1]]]
    } else {
      NULL
    }

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

    # Track which rows survive cleaning (for aligning id_vals)
    row_mask <- rep(TRUE, nrow(mat))

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
      row_mask[!complete] <- FALSE
    }

    if (nrow(mat) < 3) {
      stop("Fewer than 3 complete rows remain after removing NAs.")
    }

    # Drop zero-variance columns
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

    # Multilevel preprocessing
    if (!is.null(level) && !is.null(id_vals_raw)) {
      id_vals <- id_vals_raw[row_mask]

      if (level == "between") {
        # Aggregate to person means
        mat_df <- as.data.frame(mat)
        mat_df$.id <- id_vals
        agg <- aggregate(. ~ .id, data = mat_df, FUN = mean)
        mat <- as.matrix(agg[, names(agg) != ".id", drop = FALSE])
      } else if (level == "within") {
        # Drop persons with < 2 observations
        tab <- table(id_vals)
        multi <- names(tab[tab >= 2])
        keep_rows <- id_vals %in% multi
        if (any(!keep_rows)) {
          n_single <- sum(!keep_rows)
          message("Dropping ", n_single,
                  " single-observation rows (within-person centering).")
          mat <- mat[keep_rows, , drop = FALSE]
          id_vals <- id_vals[keep_rows]
        }

        if (nrow(mat) < 3) {
          stop("Fewer than 3 rows remain after dropping ",
               "single-observation persons.")
        }

        # Person-mean center each variable
        for (j in seq_len(ncol(mat))) {
          mat[, j] <- mat[, j] - ave(mat[, j], id_vals, FUN = mean)
        }

        # Re-check zero-variance columns after centering
        col_vars2 <- apply(mat, 2, stats::var)
        zero_var2 <- colnames(mat)[col_vars2 == 0]
        if (length(zero_var2) > 0) {
          message("Dropping zero-variance columns after centering: ",
                  paste(zero_var2, collapse = ", "))
          mat <- mat[, col_vars2 > 0, drop = FALSE]
        }

        if (ncol(mat) < 2) {
          stop("At least 2 variable columns required after within centering.")
        }
      }
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

    w_prev <- fit$w
    wi_prev <- fit$wi

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

#' Extract non-zero edges from a symmetric network matrix
#' @noRd
.network_to_edges <- function(net) {
  idx <- which(upper.tri(net) & net != 0, arr.ind = TRUE)
  if (nrow(idx) == 0) {
    return(data.frame(
      from = character(0), to = character(0),
      weight = numeric(0), stringsAsFactors = FALSE
    ))
  }
  nms <- colnames(net)
  data.frame(
    from   = nms[idx[, 1]],
    to     = nms[idx[, 2]],
    weight = net[idx],
    stringsAsFactors = FALSE
  )
}


# ---- S3 methods ----

#' Print Method for Psychological Network
#'
#' @param x A \code{psych_network} object.
#' @param ... Additional arguments (ignored).
#'
#' @export
print.psych_network <- function(x, ...) {
  total_possible <- x$p * (x$p - 1) / 2
  pct <- if (total_possible > 0) {
    sprintf("%.1f%%", 100 * x$n_edges / total_possible)
  } else {
    "0%"
  }

  label <- switch(x$method,
    glasso = "Partial Correlation Network (EBICglasso)",
    pcor   = "Partial Correlation Network (unregularised)",
    cor    = "Correlation Network"
  )

  level_label <- if (!is.null(x$level)) {
    sprintf(" [%s-person]", x$level)
  } else {
    ""
  }

  n_label <- if (identical(x$level, "between")) {
    sprintf("  Variables: %d  |  Sample size: %d (unique persons)\n",
            x$p, x$n)
  } else if (identical(x$level, "within")) {
    sprintf("  Variables: %d  |  Sample size: %d (observations)\n",
            x$p, x$n)
  } else {
    sprintf("  Variables: %d  |  Sample size: %d\n", x$p, x$n)
  }

  cat(label, level_label, "\n", sep = "")
  cat(n_label)
  cat(sprintf("  Non-zero edges: %d / %d (%s)\n",
              x$n_edges, total_possible, pct))

  if (x$method == "glasso") {
    cat(sprintf("  Gamma: %.2f  |  Lambda: %.4f\n",
                x$gamma, x$lambda_selected))
  }

  r2 <- predictability.psych_network(x)
  cat(sprintf("  Avg. predictability: %.3f  |  Range: [%.3f, %.3f]\n",
              mean(r2), min(r2), max(r2)))

  invisible(x)
}


#' Print Method for Multilevel Psychological Network
#'
#' @param x A \code{psych_network_ml} object.
#' @param ... Additional arguments (ignored).
#'
#' @export
print.psych_network_ml <- function(x, ...) {
  label <- switch(x$method,
    glasso = "glasso",
    pcor   = "pcor",
    cor    = "cor"
  )
  cat(sprintf("Multilevel Psychological Network (%s)\n", label))
  cat("-- Between-person --\n")
  b <- x$between
  total_b <- b$p * (b$p - 1) / 2
  pct_b <- if (total_b > 0) {
    sprintf("%.1f%%", 100 * b$n_edges / total_b)
  } else {
    "0%"
  }
  cat(sprintf("  Variables: %d  |  Sample size: %d (unique persons)\n",
              b$p, b$n))
  cat(sprintf("  Non-zero edges: %d / %d (%s)\n",
              b$n_edges, total_b, pct_b))
  if (b$method == "glasso") {
    cat(sprintf("  Gamma: %.2f  |  Lambda: %.4f\n",
                b$gamma, b$lambda_selected))
  }
  r2_b <- predictability.psych_network(b)
  cat(sprintf("  Avg. predictability: %.3f  |  Range: [%.3f, %.3f]\n",
              mean(r2_b), min(r2_b), max(r2_b)))
  cat("-- Within-person --\n")
  w <- x$within
  total_w <- w$p * (w$p - 1) / 2
  pct_w <- if (total_w > 0) {
    sprintf("%.1f%%", 100 * w$n_edges / total_w)
  } else {
    "0%"
  }
  cat(sprintf("  Variables: %d  |  Sample size: %d (observations)\n",
              w$p, w$n))
  cat(sprintf("  Non-zero edges: %d / %d (%s)\n",
              w$n_edges, total_w, pct_w))
  if (w$method == "glasso") {
    cat(sprintf("  Gamma: %.2f  |  Lambda: %.4f\n",
                w$gamma, w$lambda_selected))
  }
  r2_w <- predictability.psych_network(w)
  cat(sprintf("  Avg. predictability: %.3f  |  Range: [%.3f, %.3f]\n",
              mean(r2_w), min(r2_w), max(r2_w)))
  invisible(x)
}


#' Plot Method for Multilevel Psychological Network
#'
#' @description
#' Plots the between-person and within-person networks side by side using
#' \code{cograph::tplot()}.
#' Requires the \pkg{cograph} package to be installed.
#'
#' @param x A \code{psych_network_ml} object.
#' @param predictability Logical. If \code{TRUE}, display node predictability
#'   as pie chart rings (default: \code{FALSE}).
#' @param pie_color Character. Color for the predictability pie segments.
#'   Default: \code{"#377EB8"}.
#' @param ... Additional arguments passed to \code{cograph::tplot()}.
#'
#' @importFrom graphics par title
#' @export
plot.psych_network_ml <- function(x, predictability = FALSE,
                                  pie_color = "#377EB8", ...) {
  if (!requireNamespace("cograph", quietly = TRUE)) {
    stop(
      "Package 'cograph' is required for plotting. ",
      "Install it with: install.packages('cograph')"
    )
  }
  old_par <- graphics::par(mfrow = c(1, 2))
  on.exit(graphics::par(old_par))

  dots <- list(...)
  if (predictability) {
    r2 <- predictability.psych_network_ml(x)
    dots_b <- c(dots, list(pie = r2$between,
                           pieColor = rep(pie_color, length(r2$between))))
    dots_w <- c(dots, list(pie = r2$within,
                           pieColor = rep(pie_color, length(r2$within))))
  } else {
    dots_b <- dots
    dots_w <- dots
  }

  do.call(cograph::tplot,
          c(list(x$between$network_matrix, directed = FALSE), dots_b))
  graphics::title("Between-person")
  do.call(cograph::tplot,
          c(list(x$within$network_matrix, directed = FALSE), dots_w))
  graphics::title("Within-person")
}


#' Plot Method for Psychological Network
#'
#' @description
#' Plots the network using \code{cograph::tplot()}.
#' Requires the \pkg{cograph} package to be installed.
#' When \code{predictability = TRUE}, node predictability (R\eqn{^2}) is shown
#' as pie rings around each node.
#'
#' @param x A \code{psych_network} object.
#' @param predictability Logical. If \code{TRUE}, display node predictability
#'   as pie chart rings (default: \code{FALSE}).
#' @param pie_color Character. Color for the predictability pie segments.
#'   Default: \code{"#377EB8"} (blue, following mgm convention).
#' @param ... Additional arguments passed to \code{cograph::tplot()}.
#'
#' @export
plot.psych_network <- function(x, predictability = FALSE,
                               pie_color = "#377EB8", ...) {
  if (!requireNamespace("cograph", quietly = TRUE)) {
    stop(
      "Package 'cograph' is required for plotting. ",
      "Install it with: install.packages('cograph')"
    )
  }
  dots <- list(...)
  if (predictability) {
    r2 <- predictability.psych_network(x)
    dots$pie <- r2
    dots$pieColor <- rep(pie_color, length(r2))
  }
  do.call(cograph::tplot,
          c(list(x$network_matrix, directed = FALSE), dots))
}


# ---- Predictability ----

#' Compute Node Predictability
#'
#' @description
#' Computes the proportion of variance explained (R\eqn{^2}) for each node in
#' the network, following Haslbeck & Waldorp (2018).
#'
#' For \code{method = "glasso"} or \code{"pcor"}, predictability is computed
#' analytically from the precision matrix:
#' \deqn{R^2_j = 1 - 1 / \Omega_{jj}}
#' where \eqn{\Omega} is the precision (inverse correlation) matrix.
#'
#' For \code{method = "cor"}, predictability is the multiple R\eqn{^2} from
#' regressing each node on its network neighbors (nodes with non-zero edges).
#'
#' @param object A \code{psych_network} or \code{psych_network_ml} object.
#' @param ... Additional arguments (ignored).
#'
#' @return For \code{psych_network}: a named numeric vector of R\eqn{^2} values
#'   (one per node, between 0 and 1).
#'
#'   For \code{psych_network_ml}: a list with elements \code{$between} and
#'   \code{$within}, each a named numeric vector.
#'
#' @references
#' Haslbeck, J. M. B., & Waldorp, L. J. (2018). How well do network models
#' predict observations? On the importance of predictability in network models.
#' \emph{Behavior Research Methods}, 50(2), 853--861.
#' \doi{10.3758/s13428-017-0910-x}
#'
#' @examples
#' \dontrun{
#' freq <- convert_sequence_format(group_regulation, format = "frequency")
#' net <- build_network(freq)
#' predictability(net)
#'
#' # Plot with predictability rings
#' plot(net, predictability = TRUE)
#' }
#'
#' @export
predictability <- function(object, ...) {
  UseMethod("predictability")
}


#' @rdname predictability
#' @export
predictability.psych_network <- function(object, ...) {
  if (object$method %in% c("glasso", "pcor")) {
    # From precision matrix: R²_j = 1 - 1/Omega_jj
    omega_diag <- diag(object$precision_matrix)
    r2 <- 1 - 1 / omega_diag
  } else {
    # cor method: multiple R² from correlation matrix
    S <- object$cor_matrix
    net <- object$network_matrix
    p <- ncol(net)
    r2 <- vapply(seq_len(p), function(j) {
      neighbors <- which(net[j, ] != 0)
      if (length(neighbors) == 0L) return(0)
      if (length(neighbors) == 1L) return(S[neighbors, j]^2)
      r_vec <- S[neighbors, j]
      R_nn <- S[neighbors, neighbors]
      tryCatch(
        as.numeric(crossprod(r_vec, solve(R_nn, r_vec))),
        error = function(e) 0
      )
    }, numeric(1))
  }
  r2 <- pmin(pmax(r2, 0), 1)
  names(r2) <- colnames(object$network_matrix)
  r2
}


#' @rdname predictability
#' @export
predictability.psych_network_ml <- function(object, ...) {
  list(
    between = predictability(object$between),
    within  = predictability(object$within)
  )
}
