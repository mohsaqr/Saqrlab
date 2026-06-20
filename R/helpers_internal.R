# Internal helper functions.
#
# These were referenced by several exported functions but never defined, which
# left those functions erroring on basic usage. They are implemented here from
# their call-site contracts. All are internal (@noRd) -- adding them does not
# change the package's public NAMESPACE.

#' Extract a transition / weight matrix from a model object
#'
#' Accepts a fitted `tna` model, a plain list carrying a `weights` element, or a
#' bare matrix, and returns the weight matrix (preserving state row/column
#' names). Used by the network-comparison functions.
#'
#' @param model A `tna` object, a list with a `weights` element, or a matrix.
#' @return A numeric matrix of edge weights.
#' @noRd
extract_transition_matrix <- function(model) {
  if (is.matrix(model)) {
    return(model)
  }
  w <- if (inherits(model, "tna")) {
    model$weights
  } else if (is.list(model) && !is.null(model[["weights"]])) {
    model[["weights"]]
  } else {
    stop(
      "Cannot extract a transition matrix: expected a 'tna' model, a list with ",
      "a 'weights' element, or a matrix.",
      call. = FALSE
    )
  }
  if (is.null(w)) {
    stop("Model has no 'weights' matrix to extract.", call. = FALSE)
  }
  as.matrix(w)
}

#' Test whether a value falls within an optional range filter
#'
#' Returns `TRUE` when `range` is `NULL` (interpreted as "no filter applied"),
#' otherwise tests `range[1] <= val <= range[2]`. A length-1 `range` is treated
#' as an exact-match target. Missing/empty `val` against an active range yields
#' `FALSE`. Used by the grid-result filtering in `summarize_grid_results()`.
#'
#' @param val A scalar numeric value (or `NULL`).
#' @param range `NULL`, a single value, or a `c(min, max)` pair.
#' @return A single logical.
#' @noRd
check_val_in_range <- function(val, range) {
  if (is.null(range)) {
    return(TRUE)
  }
  if (is.null(val) || length(val) == 0L || all(is.na(val))) {
    return(FALSE)
  }
  if (length(range) == 1L) {
    return(isTRUE(val == range))
  }
  isTRUE(val >= range[1L] && val <= range[2L])
}

#' Row-bind a list of data frames, tolerant of NULLs and differing columns
#'
#' Drops `NULL`/zero-row elements, takes the union of all columns (filling
#' absent ones with `NA`), and `rbind`s the result. Returns an empty
#' `data.frame` when nothing is bindable. `label` is used only to make any
#' failure message informative. Used throughout `summarize_grid_results()`.
#'
#' @param df_list A list of data frames (may contain `NULL`s).
#' @param label A short description of the contents, for error messages.
#' @return A single base `data.frame`.
#' @noRd
safe_bind_rows <- function(df_list, label = "data") {
  df_list <- Filter(function(x) !is.null(x) && NROW(x) > 0L, df_list)
  if (length(df_list) == 0L) {
    return(data.frame())
  }
  df_list <- lapply(df_list, function(d) {
    d <- as.data.frame(d, stringsAsFactors = FALSE)
    factor_cols <- vapply(d, is.factor, logical(1))
    if (any(factor_cols)) d[factor_cols] <- lapply(d[factor_cols], as.character)
    d
  })
  all_cols <- unique(unlist(lapply(df_list, names), use.names = FALSE))
  aligned <- lapply(df_list, function(d) {
    missing <- setdiff(all_cols, names(d))
    if (length(missing) > 0L) {
      d[missing] <- NA
    }
    d[, all_cols, drop = FALSE]
  })
  out <- tryCatch(
    do.call(rbind, aligned),
    error = function(e) {
      stop(
        sprintf("Failed to combine %s: %s", label, conditionMessage(e)),
        call. = FALSE
      )
    }
  )
  rownames(out) <- NULL
  out
}

#' Reshape long-format sequence data to wide (one row per id)
#'
#' Converts long event data (one row per id-time-state) into a wide sequence
#' table: one row per `id_col`, with the `state_col` values ordered by
#' `time_col` spread across columns `T1`, `T2`, ... Sequences shorter than the
#' longest are right-padded with `NA`. Used by `simulate_group_tna_networks()`
#' to feed `fit_network_model(..., "group_tna")`.
#'
#' @param long_data A data frame / tibble in long format.
#' @param id_col Name of the identifier column (default "Actor").
#' @param time_col Name of the ordering column (default "Time").
#' @param state_col Name of the state/action column (default "Action").
#' @return A base data frame with the id column followed by `T1..Tk` sequence
#'   columns.
#' @noRd
long_to_wide <- function(long_data, id_col = "Actor",
                         time_col = "Time", state_col = "Action") {
  long_data <- as.data.frame(long_data, stringsAsFactors = FALSE)
  stopifnot(all(c(id_col, time_col, state_col) %in% names(long_data)))

  ids <- unique(long_data[[id_col]])
  ids <- ids[!is.na(ids)]
  if (length(ids) == 0L) {
    empty <- data.frame(stringsAsFactors = FALSE)
    empty[[id_col]] <- ids
    return(empty)
  }
  seqs <- lapply(ids, function(id) {
    rows <- long_data[long_data[[id_col]] == id, , drop = FALSE]
    rows <- rows[order(rows[[time_col]]), , drop = FALSE]
    as.character(rows[[state_col]])
  })

  max_len <- max(vapply(seqs, length, integer(1L)))
  padded <- lapply(seqs, function(s) {
    length(s) <- max_len # extends with NA
    s
  })

  wide <- as.data.frame(
    do.call(rbind, padded),
    stringsAsFactors = FALSE
  )
  names(wide) <- paste0("T", seq_len(max_len))
  wide <- cbind(stats::setNames(data.frame(ids, stringsAsFactors = FALSE), id_col), wide)
  rownames(wide) <- NULL
  wide
}
