# ---- Built-in Estimator Implementations ----

# ---- Shared helpers ----

#' Select state columns from a data frame
#'
#' Resolves which columns contain state data based on explicit \code{cols},
#' exclusion of \code{id} columns, or all columns as fallback.
#'
#' @param data Data frame.
#' @param id Character vector or NULL. ID column(s) to exclude.
#' @param cols Character vector or NULL. Explicit state columns.
#' @return Character vector of state column names.
#' @noRd
.select_state_cols <- function(data, id = NULL, cols = NULL) {
  if (!is.null(cols)) {
    cols
  } else if (!is.null(id)) {
    setdiff(names(data), id)
  } else {
    names(data)
  }
}


# ---- Core transition counting engine ----

#' Count transitions from sequence data
#'
#' Dispatches to wide or long format counting. Returns a square integer matrix
#' of transition frequencies.
#'
#' @param data Data frame of sequence data.
#' @param format Character: \code{"auto"}, \code{"wide"}, or \code{"long"}.
#' @param action Character. Action column name (long format).
#' @param id Character vector. ID column(s).
#' @param time Character. Time column name (long format).
#' @param cols Character vector. State columns (wide format).
#'
#' @return Square integer matrix with row/column names = sorted unique states.
#' @noRd
.count_transitions <- function(data,
                               format = "auto",
                               action = "Action",
                               id = NULL,
                               time = "Time",
                               cols = NULL) {
  stopifnot(is.data.frame(data))

  if (format == "auto") {
    format <- if (action %in% names(data)) "long" else "wide"
  }

  if (format == "wide") {
    .count_transitions_wide(data, id = id, cols = cols)
  } else {
    .count_transitions_long(data, action = action, id = id, time = time)
  }
}


#' Count transitions from wide format (vectorized base R)
#'
#' Each row is a sequence, columns are consecutive time points.
#' Uses matrix slicing + tabulate for speed.
#'
#' @noRd
.count_transitions_wide <- function(data, id = NULL, cols = NULL) {
  state_cols <- .select_state_cols(data, id, cols)

  missing_cols <- setdiff(state_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("Columns not found: ", paste(missing_cols, collapse = ", "))
  }
  if (length(state_cols) < 2L) {
    stop("At least 2 state columns are required for wide format.")
  }

  mat <- as.matrix(data[, state_cols, drop = FALSE])
  nc <- ncol(mat)

  # from = all columns except last, to = all columns except first
  from_vec <- as.vector(mat[, -nc, drop = FALSE])
  to_vec <- as.vector(mat[, -1L, drop = FALSE])

  # Remove pairs where either is NA
  valid <- !is.na(from_vec) & !is.na(to_vec)
  from_vec <- from_vec[valid]
  to_vec <- to_vec[valid]

  # Integer encode + tabulate
  all_states <- sort(unique(c(from_vec, to_vec)))
  n_states <- length(all_states)

  if (n_states == 0L) {
    return(matrix(0L, 0, 0))
  }

  from_int <- match(from_vec, all_states)
  to_int <- match(to_vec, all_states)
  pair_idx <- (from_int - 1L) * n_states + to_int

  counts <- tabulate(pair_idx, nbins = n_states * n_states)

  matrix(
    as.integer(counts),
    nrow = n_states,
    ncol = n_states,
    byrow = TRUE,
    dimnames = list(all_states, all_states)
  )
}


#' Count transitions from long format (data.table)
#'
#' Uses data.table for fast grouping and lag computation.
#'
#' @importFrom data.table setDT setorderv
#' @noRd
.count_transitions_long <- function(data, action = "Action", id = NULL,
                                    time = "Time") {
  if (!action %in% names(data)) {
    stop("Action column '", action, "' not found in data.")
  }
  if (!is.null(id)) {
    missing_ids <- setdiff(id, names(data))
    if (length(missing_ids) > 0) {
      stop("ID column(s) not found: ", paste(missing_ids, collapse = ", "))
    }
  }

  dt <- data.table::as.data.table(data)

  # Order by id + time
  order_cols <- c(id, if (time %in% names(dt)) time)
  if (length(order_cols) > 0) {
    data.table::setorderv(dt, order_cols)
  }

  action_col <- action

  # Build group key for sequences
  if (is.null(id)) {
    # Single sequence: consecutive pairs from all rows
    actions <- dt[[action_col]]
    n <- length(actions)
    if (n < 2L) {
      all_vals <- unique(actions[!is.na(actions)])
      all_states <- sort(all_vals)
      n_states <- length(all_states)
      return(matrix(0L, nrow = n_states, ncol = n_states,
                    dimnames = list(all_states, all_states)))
    }
    from_vec <- actions[-n]
    to_vec <- actions[-1L]
    # Filter out pairs where either side is NA
    valid <- !is.na(from_vec) & !is.na(to_vec)
    from_vec <- from_vec[valid]
    to_vec <- to_vec[valid]
  } else {
    # Group by ID columns, extract consecutive pairs
    # Use data.table's fast grouping
    if (length(id) == 1L) {
      grp_col <- id
    } else {
      # Create composite group key
      dt[, .grp_key := do.call(paste, c(.SD, sep = "\x1f")),
         .SDcols = id]
      grp_col <- ".grp_key"
    }

    # Extract from/to pairs per group using data.table
    # NAs are kept in position; pairs with NA on either side are filtered out
    pairs <- dt[, {
      a <- get(action_col)
      n <- length(a)
      if (n < 2L) {
        list(from = character(0), to = character(0))
      } else {
        f <- a[-n]
        t <- a[-1L]
        ok <- !is.na(f) & !is.na(t)
        list(from = f[ok], to = t[ok])
      }
    }, by = grp_col]

    from_vec <- pairs$from
    to_vec <- pairs$to
  }

  if (length(from_vec) == 0L) {
    # Collect all states to set matrix dimensions
    all_vals <- unique(dt[[action_col]])
    all_vals <- all_vals[!is.na(all_vals)]
    all_states <- sort(all_vals)
    n_states <- length(all_states)
    return(matrix(0L, nrow = n_states, ncol = n_states,
                  dimnames = list(all_states, all_states)))
  }

  # Integer encode + tabulate
  all_states <- sort(unique(c(from_vec, to_vec)))
  n_states <- length(all_states)

  from_int <- match(from_vec, all_states)
  to_int <- match(to_vec, all_states)
  pair_idx <- (from_int - 1L) * n_states + to_int

  counts <- tabulate(pair_idx, nbins = n_states * n_states)

  matrix(
    as.integer(counts),
    nrow = n_states,
    ncol = n_states,
    byrow = TRUE,
    dimnames = list(all_states, all_states)
  )
}


# ---- Transition estimators ----

#' Frequency estimator: raw transition counts
#' @noRd
.estimator_frequency <- function(data,
                                 format = "auto",
                                 action = "Action",
                                 id = NULL,
                                 time = "Time",
                                 cols = NULL,
                                 ...) {
  freq_mat <- .count_transitions(
    data, format = format, action = action, id = id, time = time, cols = cols
  )
  states <- rownames(freq_mat)
  list(
    matrix = freq_mat,
    nodes = states,
    directed = TRUE,
    cleaned_data = data,
    frequency_matrix = freq_mat
  )
}


#' Relative estimator: row-normalized transition probabilities
#' @noRd
.estimator_relative <- function(data,
                                format = "auto",
                                action = "Action",
                                id = NULL,
                                time = "Time",
                                cols = NULL,
                                ...) {
  freq_mat <- .count_transitions(
    data, format = format, action = action, id = id, time = time, cols = cols
  )
  states <- rownames(freq_mat)

  # Row-normalize
  row_sums <- rowSums(freq_mat)
  rel_mat <- freq_mat
  storage.mode(rel_mat) <- "double"
  nonzero <- row_sums > 0
  rel_mat[nonzero, ] <- rel_mat[nonzero, ] / row_sums[nonzero]

  list(
    matrix = rel_mat,
    nodes = states,
    directed = TRUE,
    cleaned_data = data,
    frequency_matrix = freq_mat
  )
}


#' Co-occurrence estimator: positional co-occurrence within sequences
#'
#' Counts all positional column pairs (i, j) where i < j. For each pair,
#' if both positions have non-NA states, the co-occurrence count is incremented
#' for both (from, to) and (to, from). This matches tna::ctna() semantics.
#'
#' @noRd
.estimator_co_occurrence <- function(data,
                                     format = "auto",
                                     action = "Action",
                                     id = NULL,
                                     time = "Time",
                                     cols = NULL,
                                     ...) {
  stopifnot(is.data.frame(data))

  if (format == "auto") {
    format <- if (action %in% names(data)) "long" else "wide"
  }

  if (format == "wide") {
    cooc_mat <- .count_cooccurrence_wide(data, id = id, cols = cols)
  } else {
    cooc_mat <- .count_cooccurrence_long(
      data, action = action, id = id, time = time
    )
  }

  list(
    matrix = cooc_mat,
    nodes = rownames(cooc_mat),
    directed = FALSE,
    cleaned_data = data
  )
}


#' Count positional co-occurrences from wide format
#'
#' For each pair of column positions (i, j) where i < j, counts how many
#' sequences have non-NA values at both positions. Symmetric: both (from, to)
#' and (to, from) are incremented.
#'
#' Strategy: integer-encode the matrix once, then iterate over column pairs
#' with lightweight tabulate accumulation. Avoids materializing huge
#' n_rows x n_pairs matrices.
#' @noRd
.count_cooccurrence_wide <- function(data, id = NULL, cols = NULL) {
  state_cols <- .select_state_cols(data, id, cols)

  mat <- as.matrix(data[, state_cols, drop = FALSE])
  nc <- ncol(mat)
  all_states <- sort(unique(as.vector(mat[!is.na(mat)])))
  n_states <- length(all_states)

  if (n_states == 0L || nc < 2L) {
    cooc <- matrix(0, n_states, n_states,
                   dimnames = list(all_states, all_states))
    return(cooc)
  }

  nbins <- n_states * n_states

  # Integer-encode the entire matrix once (NA stays NA)
  int_mat <- matrix(match(mat, all_states), nrow = nrow(mat), ncol = nc)

  # Accumulate counts across column pairs
  counts <- integer(nbins)

  for (i in seq_len(nc - 1L)) {
    col_i <- int_mat[, i]
    for (j in seq(i + 1L, nc)) {
      col_j <- int_mat[, j]
      valid <- !is.na(col_i) & !is.na(col_j)
      fi <- col_i[valid]
      tj <- col_j[valid]
      # Forward direction: i -> j
      idx_fwd <- (fi - 1L) * n_states + tj
      # Reverse direction: j -> i
      idx_rev <- (tj - 1L) * n_states + fi
      counts <- counts + tabulate(idx_fwd, nbins) + tabulate(idx_rev, nbins)
    }
  }

  cooc <- matrix(
    as.numeric(counts),
    nrow = n_states,
    ncol = n_states,
    byrow = TRUE,
    dimnames = list(all_states, all_states)
  )

  # Self-pairs (A,A) are double-counted by the bidirectional approach
  diag(cooc) <- diag(cooc) / 2

  cooc
}


#' Count positional co-occurrences from long format
#'
#' Converts to wide-like structure per group, then counts column-pair
#' co-occurrences.
#' @noRd
.count_cooccurrence_long <- function(data, action = "Action", id = NULL,
                                     time = "Time") {
  if (!action %in% names(data)) {
    stop("Action column '", action, "' not found in data.")
  }

  dt <- data.table::as.data.table(data)

  # Order by id + time
  order_cols <- c(id, if (time %in% names(dt)) time)
  if (length(order_cols) > 0) {
    data.table::setorderv(dt, order_cols)
  }

  # Build group key
  if (is.null(id)) {
    dt[, .seq_grp := 1L]
    grp_col <- ".seq_grp"
  } else if (length(id) == 1L) {
    grp_col <- id
  } else {
    dt[, .grp_key := do.call(paste, c(.SD, sep = "\x1f")),
       .SDcols = id]
    grp_col <- ".grp_key"
  }

  action_col <- action

  # For each group, create all position pairs and collect (from, to)
  pairs <- dt[!is.na(get(action_col)), {
    a <- get(action_col)
    n <- length(a)
    if (n < 2L) {
      list(from = character(0), to = character(0))
    } else {
      cp <- utils::combn(n, 2)
      f <- a[cp[1, ]]
      t <- a[cp[2, ]]
      # Both directions
      list(from = c(f, t), to = c(t, f))
    }
  }, by = grp_col]

  from_vec <- pairs$from
  to_vec <- pairs$to

  # Collect all states
  all_vals <- unique(dt[[action_col]])
  all_vals <- all_vals[!is.na(all_vals)]
  all_states <- sort(all_vals)
  n_states <- length(all_states)

  if (length(from_vec) == 0L || n_states == 0L) {
    return(matrix(0, nrow = n_states, ncol = n_states,
                  dimnames = list(all_states, all_states)))
  }

  from_int <- match(from_vec, all_states)
  to_int <- match(to_vec, all_states)
  pair_idx <- (from_int - 1L) * n_states + to_int

  counts <- tabulate(pair_idx, nbins = n_states * n_states)

  cooc <- matrix(
    as.numeric(counts),
    nrow = n_states,
    ncol = n_states,
    byrow = TRUE,
    dimnames = list(all_states, all_states)
  )

  # Self-pairs (A,A) are double-counted by the bidirectional approach
  diag(cooc) <- diag(cooc) / 2

  cooc
}


# ---- Association estimators ----

#' Prepare association input: clean data frame or validate matrix
#'
#' Handles data frame cleaning (drop NA, zero-variance, non-syntactic cols)
#' and matrix input (symmetric check, cor/cov detection).
#'
#' @return List with \code{S} (correlation matrix), \code{n} (sample size).
#' @noRd
.prepare_association_input <- function(data, id_col = NULL, n = NULL,
                                       cor_method = "pearson",
                                       input_type = "auto") {
  if (is.data.frame(data)) {
    # Exclude id columns, "rid", and non-numeric columns
    exclude <- c(id_col, "rid")
    numeric_cols <- vapply(data, is.numeric, logical(1))
    keep <- setdiff(names(data)[numeric_cols], exclude)

    # Drop columns with non-syntactic names
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
              paste(colnames(mat)[all_na], collapse = ", "))
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
    mat <- NULL
  } else {
    stop("data must be a data frame or a square symmetric matrix.")
  }

  list(S = S, n = n, mat = mat)
}


#' Correlation estimator
#' @noRd
.estimator_cor <- function(data,
                           id_col = NULL,
                           n = NULL,
                           cor_method = "pearson",
                           input_type = "auto",
                           threshold = 0,
                           ...) {
  prepared <- .prepare_association_input(
    data, id_col = id_col, n = n,
    cor_method = cor_method, input_type = input_type
  )
  S <- prepared$S
  n_obs <- prepared$n

  net <- S
  diag(net) <- 0
  net[abs(net) < threshold] <- 0

  nodes <- colnames(net)
  list(
    matrix = net,
    nodes = nodes,
    directed = FALSE,
    cleaned_data = prepared$mat,
    cor_matrix = S,
    n = n_obs,
    p = ncol(S)
  )
}


#' Partial correlation estimator (unregularized)
#' @noRd
.estimator_pcor <- function(data,
                            id_col = NULL,
                            n = NULL,
                            cor_method = "pearson",
                            input_type = "auto",
                            threshold = 0,
                            ...) {
  prepared <- .prepare_association_input(
    data, id_col = id_col, n = n,
    cor_method = cor_method, input_type = input_type
  )
  S <- prepared$S
  n_obs <- prepared$n

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

  nodes <- colnames(pcor)
  list(
    matrix = pcor,
    nodes = nodes,
    directed = FALSE,
    cleaned_data = prepared$mat,
    precision_matrix = Wi,
    cor_matrix = S,
    n = n_obs,
    p = ncol(S)
  )
}


# ---- Shared association helpers ----

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


#' Convert precision matrix to partial correlations
#' @noRd
.precision_to_pcor <- function(Wi, threshold) {
  d <- sqrt(diag(Wi))
  pcor <- -Wi / outer(d, d)
  diag(pcor) <- 0
  pcor[abs(pcor) < threshold] <- 0
  pcor
}


#' EBICglasso estimator
#' @noRd
.estimator_glasso <- function(data,
                              id_col = NULL,
                              n = NULL,
                              gamma = 0.5,
                              nlambda = 100L,
                              lambda.min.ratio = 0.01,
                              penalize.diagonal = FALSE,
                              cor_method = "pearson",
                              input_type = "auto",
                              threshold = 0,
                              ...) {
  prepared <- .prepare_association_input(
    data, id_col = id_col, n = n,
    cor_method = cor_method, input_type = input_type
  )
  S <- prepared$S
  n_obs <- prepared$n
  p <- ncol(S)

  stopifnot(is.numeric(gamma), length(gamma) == 1, gamma >= 0)
  stopifnot(is.numeric(nlambda), length(nlambda) == 1, nlambda >= 2)
  stopifnot(is.numeric(lambda.min.ratio), lambda.min.ratio > 0,
            lambda.min.ratio < 1)
  stopifnot(is.logical(penalize.diagonal), length(penalize.diagonal) == 1)

  lambda_path <- .compute_lambda_path(S, nlambda, lambda.min.ratio)
  selected <- .select_ebic(S, lambda_path, n_obs, gamma, penalize.diagonal)

  pcor <- .precision_to_pcor(selected$wi, threshold)
  colnames(pcor) <- rownames(pcor) <- colnames(S)

  nodes <- colnames(pcor)
  list(
    matrix = pcor,
    nodes = nodes,
    directed = FALSE,
    cleaned_data = prepared$mat,
    precision_matrix = selected$wi,
    cor_matrix = S,
    lambda_selected = selected$lambda,
    ebic_selected = selected$ebic,
    lambda_path = lambda_path,
    ebic_path = selected$ebic_path,
    gamma = gamma,
    n = n_obs,
    p = p
  )
}
