#' Estimate a Network
#'
#' @description
#' Universal network estimation function that supports both transition
#' networks (relative, frequency, co-occurrence) and association networks
#' (correlation, partial correlation, graphical lasso). Uses the global
#' estimator registry, so custom estimators can also be used.
#'
#' @param data Data frame (sequences or per-observation frequencies) or a
#'   square symmetric matrix (correlation or covariance).
#' @param method Character. Name of a registered estimator. Built-in methods:
#'   \code{"relative"} (default), \code{"frequency"}, \code{"co_occurrence"},
#'   \code{"cor"}, \code{"pcor"}, \code{"glasso"}.
#'   Aliases: \code{"ebicglasso"} and \code{"regularized"} map to
#'   \code{"glasso"}; \code{"partial"} maps to \code{"pcor"};
#'   \code{"correlation"} maps to \code{"cor"};
#'   \code{"transition"} maps to \code{"relative"};
#'   \code{"counts"} maps to \code{"frequency"}.
#' @param params Named list. Method-specific parameters passed to the estimator
#'   function (e.g. \code{list(gamma = 0.5)} for glasso, or
#'   \code{list(format = "wide")} for transition methods). This is the key
#'   composability feature: downstream functions like bootstrap or grid search
#'   can store and replay the full params list without knowing method internals.
#' @param scaling Character vector or NULL. Post-estimation scaling to apply
#'   (in order). Options: \code{"minmax"}, \code{"max"}, \code{"rank"},
#'   \code{"normalize"}. Can combine: \code{c("rank", "minmax")}.
#'   Default: \code{NULL} (no scaling).
#' @param threshold Numeric. Absolute values below this are set to zero in the
#'   result matrix. Default: 0 (no thresholding).
#' @param level Character or NULL. Multilevel decomposition for association
#'   methods. One of \code{NULL}, \code{"between"}, \code{"within"},
#'   \code{"both"}. Requires \code{id_col}. Default: \code{NULL}.
#' @param id_col Character. Name of the ID column for multilevel decomposition.
#'   Default: \code{NULL}.
#'
#' @return An object of class \code{"saqr_network"} containing:
#' \describe{
#'   \item{matrix}{The estimated network weight matrix.}
#'   \item{nodes}{Character vector of node names.}
#'   \item{directed}{Logical. Whether the network is directed.}
#'   \item{method}{The resolved method name.}
#'   \item{params}{The params list used (for reproducibility).}
#'   \item{scaling}{The scaling applied (or NULL).}
#'   \item{threshold}{The threshold applied.}
#'   \item{n_nodes}{Number of nodes.}
#'   \item{n_edges}{Number of non-zero edges.}
#'   \item{edges}{Data frame of non-zero edges (from, to, weight).}
#'   \item{level}{Decomposition level used (or NULL).}
#' }
#' Method-specific extras (e.g. \code{precision_matrix}, \code{cor_matrix},
#' \code{frequency_matrix}, \code{lambda_selected}, etc.) are preserved
#' from the estimator output.
#'
#' When \code{level = "both"}, returns an object of class
#' \code{"saqr_network_ml"} with \code{$between} and \code{$within}
#' sub-networks and a \code{$method} field.
#'
#' @details
#' The function works as follows:
#' \enumerate{
#'   \item Resolves method aliases to canonical names.
#'   \item Retrieves the estimator function from the global registry.
#'   \item For association methods with \code{level} specified, decomposes
#'     the data (between-person means or within-person centering).
#'   \item Calls the estimator: \code{do.call(fn, c(list(data = data), params))}.
#'   \item Applies scaling and thresholding to the result matrix.
#'   \item Extracts edges and constructs the \code{saqr_network} object.
#' }
#'
#' @examples
#' \dontrun{
#' library(tna)
#'
#' # Transition network (relative probabilities)
#' net <- estimate_network(group_regulation, method = "relative")
#' print(net)
#'
#' # Transition network (raw frequencies)
#' net_freq <- estimate_network(group_regulation, method = "frequency")
#'
#' # Association network (glasso)
#' freq_data <- convert_sequence_format(group_regulation, format = "frequency")
#' net_glasso <- estimate_network(freq_data, method = "glasso",
#'                                 params = list(gamma = 0.5, nlambda = 50))
#'
#' # With scaling
#' net_scaled <- estimate_network(group_regulation, method = "relative",
#'                                 scaling = c("rank", "minmax"))
#'
#' # Composable: replay config on new data
#' config <- net_glasso$params
#' net2 <- estimate_network(new_data, method = net_glasso$method,
#'                           params = config)
#'
#' # Custom estimator
#' register_estimator("my_method", my_fn, "My method", directed = FALSE)
#' net_custom <- estimate_network(data, method = "my_method")
#' }
#'
#' @seealso \code{\link{register_estimator}}, \code{\link{list_estimators}},
#'   \code{\link{build_network}}
#'
#' @importFrom stats aggregate ave cor complete.cases var
#' @export
estimate_network <- function(data,
                             method = "relative",
                             params = list(),
                             scaling = NULL,
                             threshold = 0,
                             level = NULL,
                             id_col = NULL) {
  stopifnot(is.character(method), length(method) == 1)
  stopifnot(is.list(params))
  stopifnot(is.numeric(threshold), length(threshold) == 1, threshold >= 0)

  # Resolve method aliases
  method <- .resolve_method_alias(method)

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

  # Validate scaling
  if (!is.null(scaling)) {
    valid_scaling <- c("minmax", "max", "rank", "normalize")
    bad <- setdiff(scaling, valid_scaling)
    if (length(bad) > 0) {
      stop("Unknown scaling method(s): ", paste(bad, collapse = ", "),
           ". Options: ", paste(valid_scaling, collapse = ", "),
           call. = FALSE)
    }
  }

  # Get estimator from registry
  estimator <- get_estimator(method)

  # level = "both": recursive dispatch
  if (identical(level, "both")) {
    between <- estimate_network(
      data, method = method, params = params, scaling = scaling,
      threshold = threshold, level = "between", id_col = id_col
    )
    within_net <- estimate_network(
      data, method = method, params = params, scaling = scaling,
      threshold = threshold, level = "within", id_col = id_col
    )
    result <- list(between = between, within = within_net, method = method)
    class(result) <- "saqr_network_ml"
    return(result)
  }

  # Multilevel decomposition for association methods
  if (!is.null(level) && !estimator$directed) {
    data <- .decompose_multilevel(data, id_col = id_col, level = level)
  }

  # Call estimator
  est_result <- do.call(estimator$fn, c(list(data = data), params))

  # Validate estimator output
  if (!is.list(est_result) ||
      is.null(est_result$matrix) ||
      is.null(est_result$nodes) ||
      is.null(est_result$directed)) {
    stop("Estimator '", method,
         "' must return a list with 'matrix', 'nodes', and 'directed'.",
         call. = FALSE)
  }

  net_matrix <- est_result$matrix
  nodes <- est_result$nodes
  directed <- est_result$directed

  # Apply scaling
  if (!is.null(scaling)) {
    net_matrix <- .apply_scaling(net_matrix, scaling)
  }

  # Apply threshold
  if (threshold > 0) {
    net_matrix[abs(net_matrix) < threshold] <- 0
  }

  # Extract edges
  edges <- .extract_edges_from_matrix(net_matrix, directed = directed)

  # Build saqr_network object
  result <- list(
    matrix = net_matrix,
    nodes = nodes,
    directed = directed,
    method = method,
    params = params,
    scaling = scaling,
    threshold = threshold,
    n_nodes = length(nodes),
    n_edges = nrow(edges),
    edges = edges,
    level = level
  )

  # Carry over method-specific extras
  extras <- setdiff(names(est_result), c("matrix", "nodes", "directed"))
  for (key in extras) {
    result[[key]] <- est_result[[key]]
  }

  structure(result, class = "saqr_network")
}


# ---- Method alias resolution ----

#' Resolve method aliases to canonical names
#' @noRd
.resolve_method_alias <- function(method) {
  aliases <- c(
    ebicglasso  = "glasso",
    regularized = "glasso",
    partial     = "pcor",
    correlation = "cor",
    transition  = "relative",
    counts      = "frequency"
  )
  if (method %in% names(aliases)) {
    aliases[[method]]
  } else {
    method
  }
}


# ---- Post-estimation scaling ----

#' Apply scaling transformations to a network matrix
#' @noRd
.apply_scaling <- function(mat, scaling) {
  for (s in scaling) {
    mat <- switch(s,
      minmax = {
        vals <- mat[mat != 0]
        if (length(vals) == 0) {
          mat
        } else {
          rng <- range(vals)
          if (rng[1] == rng[2]) mat
          else {
            mat[mat != 0] <- (mat[mat != 0] - rng[1]) / (rng[2] - rng[1])
            mat
          }
        }
      },
      max = {
        max_abs <- max(abs(mat))
        if (max_abs > 0) mat / max_abs else mat
      },
      rank = {
        nz <- mat != 0
        if (any(nz)) {
          mat[nz] <- rank(mat[nz])
          mat
        } else {
          mat
        }
      },
      normalize = {
        rs <- rowSums(abs(mat))
        nonzero_rows <- rs > 0
        mat[nonzero_rows, ] <- mat[nonzero_rows, ] / rs[nonzero_rows]
        mat
      },
      mat  # default: no change
    )
  }
  mat
}


# ---- Edge extraction ----

#' Extract non-zero edges from a network matrix
#'
#' For undirected networks, uses upper triangle only.
#' For directed networks, uses all non-diagonal non-zero entries.
#'
#' @noRd
.extract_edges_from_matrix <- function(mat, directed = FALSE) {
  nms <- rownames(mat)
  if (is.null(nms)) nms <- paste0("V", seq_len(nrow(mat)))

  if (directed) {
    idx <- which(mat != 0 & row(mat) != col(mat), arr.ind = TRUE)
  } else {
    idx <- which(upper.tri(mat) & mat != 0, arr.ind = TRUE)
  }

  if (nrow(idx) == 0) {
    return(data.frame(
      from = character(0), to = character(0),
      weight = numeric(0), stringsAsFactors = FALSE
    ))
  }

  data.frame(
    from   = nms[idx[, 1]],
    to     = nms[idx[, 2]],
    weight = mat[idx],
    stringsAsFactors = FALSE
  )
}


# ---- Multilevel decomposition ----

#' Decompose data for multilevel analysis
#'
#' @param data Data frame.
#' @param id_col Character. Grouping variable name.
#' @param level Character: "between" or "within".
#'
#' @return Transformed data frame.
#' @noRd
.decompose_multilevel <- function(data, id_col, level) {
  stopifnot(is.data.frame(data))
  grp_var <- id_col[1]

  if (!grp_var %in% names(data)) {
    stop("id_col '", grp_var, "' not found in data.", call. = FALSE)
  }

  # Get numeric columns (exclude id columns and "rid")
  exclude <- c(id_col, "rid")
  numeric_cols <- vapply(data, is.numeric, logical(1))
  keep <- setdiff(names(data)[numeric_cols], exclude)

  if (length(keep) < 2) {
    stop("At least 2 numeric columns are required for multilevel decomposition.")
  }

  mat <- data[, keep, drop = FALSE]
  id_vals <- data[[grp_var]]

  if (level == "between") {
    # Aggregate to person means
    mat$.id <- id_vals
    agg <- aggregate(. ~ .id, data = mat, FUN = mean)
    result <- agg[, names(agg) != ".id", drop = FALSE]
    return(as.data.frame(result))

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
    mat_m <- as.matrix(mat)
    for (j in seq_len(ncol(mat_m))) {
      mat_m[, j] <- mat_m[, j] - ave(mat_m[, j], id_vals, FUN = mean)
    }

    return(as.data.frame(mat_m))
  }

  data
}


# ---- S3 methods ----

#' Print Method for saqr_network
#'
#' @param x A \code{saqr_network} object.
#' @param ... Additional arguments (ignored).
#'
#' @export
print.saqr_network <- function(x, ...) {
  method_labels <- c(
    relative      = "Transition Network (relative probabilities)",
    frequency     = "Transition Network (frequency counts)",
    co_occurrence = "Co-occurrence Network",
    glasso        = "Partial Correlation Network (EBICglasso)",
    pcor          = "Partial Correlation Network (unregularised)",
    cor           = "Correlation Network"
  )

  label <- if (x$method %in% names(method_labels)) {
    method_labels[[x$method]]
  } else {
    sprintf("Network (method: %s)", x$method)
  }

  dir_label <- if (x$directed) " [directed]" else " [undirected]"

  level_label <- if (!is.null(x$level)) {
    sprintf(" [%s-person]", x$level)
  } else {
    ""
  }

  cat(label, dir_label, level_label, "\n", sep = "")
  cat(sprintf("  Nodes: %d  |  Edges: %d\n", x$n_nodes, x$n_edges))

  if (!is.null(x$n)) {
    cat(sprintf("  Sample size: %d\n", x$n))
  }

  if (x$method == "glasso" && !is.null(x$gamma)) {
    cat(sprintf("  Gamma: %.2f  |  Lambda: %.4f\n",
                x$gamma, x$lambda_selected))
  }

  if (!is.null(x$scaling)) {
    cat(sprintf("  Scaling: %s\n", paste(x$scaling, collapse = " -> ")))
  }

  if (x$threshold > 0) {
    cat(sprintf("  Threshold: %g\n", x$threshold))
  }

  invisible(x)
}


#' Print Method for saqr_network_ml
#'
#' @param x A \code{saqr_network_ml} object.
#' @param ... Additional arguments (ignored).
#'
#' @export
print.saqr_network_ml <- function(x, ...) {
  cat(sprintf("Multilevel Network (method: %s)\n", x$method))
  cat("-- Between-person --\n")
  b <- x$between
  cat(sprintf("  Nodes: %d  |  Edges: %d\n", b$n_nodes, b$n_edges))
  if (!is.null(b$n)) {
    cat(sprintf("  Sample size: %d (unique persons)\n", b$n))
  }
  cat("-- Within-person --\n")
  w <- x$within
  cat(sprintf("  Nodes: %d  |  Edges: %d\n", w$n_nodes, w$n_edges))
  if (!is.null(w$n)) {
    cat(sprintf("  Sample size: %d (observations)\n", w$n))
  }
  invisible(x)
}


#' Plot Method for saqr_network
#'
#' @description
#' Plots the network using \code{cograph::splot()}.
#' Requires the \pkg{cograph} package to be installed.
#'
#' @param x A \code{saqr_network} object.
#' @param ... Additional arguments passed to \code{cograph::splot()}.
#'
#' @export
plot.saqr_network <- function(x, ...) {
  if (!requireNamespace("cograph", quietly = TRUE)) {
    stop(
      "Package 'cograph' is required for plotting. ",
      "Install it with: install.packages('cograph')"
    )
  }

  node_cols <- .node_colors(x$n_nodes)

  dots <- list(
    x = x$matrix,
    directed = x$directed,
    node_fill = node_cols,
    edge_labels = TRUE,
    edge_label_size = 0.65,
    node_size = 8,
    theme = "colorblind",
    ...
  )

  do.call(cograph::splot, dots)
}


#' Plot Method for saqr_network_ml
#'
#' @description
#' Plots between and within networks side by side using
#' \code{cograph::splot()}.
#'
#' @param x A \code{saqr_network_ml} object.
#' @param ... Additional arguments passed to \code{cograph::splot()}.
#'
#' @importFrom graphics par
#' @export
plot.saqr_network_ml <- function(x, ...) {
  if (!requireNamespace("cograph", quietly = TRUE)) {
    stop(
      "Package 'cograph' is required for plotting. ",
      "Install it with: install.packages('cograph')"
    )
  }
  old_par <- graphics::par(mfrow = c(1, 2))
  on.exit(graphics::par(old_par))

  node_cols <- .node_colors(x$between$n_nodes)

  dots <- list(
    directed = x$between$directed,
    node_fill = node_cols,
    edge_labels = TRUE,
    edge_label_size = 0.65,
    node_size = 8,
    theme = "colorblind",
    ...
  )

  dots_b <- c(list(x = x$between$matrix, title = "Between-person"), dots)
  dots_w <- c(list(x = x$within$matrix, title = "Within-person"), dots)

  do.call(cograph::splot, dots_b)
  do.call(cograph::splot, dots_w)
}
