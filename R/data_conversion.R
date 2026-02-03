#' Data Format Conversion Functions
#'
#' @description
#' Functions for converting between wide and long data formats commonly used
#' in Temporal Network Analysis.
#'
#' @name data_conversion
#' @keywords internal
NULL

#' Convert Wide Sequences to Long Format
#'
#' @description
#' Convert sequence data from wide format (one row per sequence, columns as
#' time points) to long format (one row per action).
#'
#' @param data Data frame in wide format with sequences in rows.
#' @param id_col Character. Name of the ID column, or NULL to auto-generate IDs.
#'   Default: NULL.
#' @param time_prefix Character. Prefix for time point columns (e.g., "V" for V1,
#'   V2, ...). Default: "V".
#' @param action_col Character. Name of the action column in output.
#'   Default: "Action".
#' @param time_col Character. Name of the time column in output.
#'   Default: "Time".
#' @param drop_na Logical. Whether to drop NA values. Default: TRUE.
#'
#' @return A data frame in long format with columns:
#' \describe{
#'   \item{id}{Sequence identifier (integer).}
#'   \item{Time}{Time point within the sequence (integer).}
#'   \item{Action}{The action/state at that time point (character).}
#' }
#' Any additional columns from the original data are preserved.
#'
#' @details
#' This function converts data from the format produced by `simulate_sequences()`
#' to the long format used by many TNA functions and analyses.
#'
#' @examples
#' \dontrun{
#' # Generate wide format sequences
#' trans_mat <- matrix(c(0.7, 0.2, 0.1, 0.3, 0.5, 0.2, 0.2, 0.3, 0.5),
#'                     nrow = 3, byrow = TRUE)
#' rownames(trans_mat) <- colnames(trans_mat) <- c("A", "B", "C")
#' init_probs <- c(A = 0.5, B = 0.3, C = 0.2)
#'
#' wide_data <- simulate_sequences(
#'   transition_matrix = trans_mat,
#'   initial_probabilities = init_probs,
#'   max_seq_length = 10,
#'   num_rows = 50
#' )
#'
#' # Convert to long format
#' long_data <- wide_to_long(wide_data)
#' head(long_data)
#' }
#'
#' @seealso \code{\link{long_to_wide}} for the reverse conversion,
#'   \code{\link{prepare_for_tna}} for preparing data for TNA analysis.
#'
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr mutate arrange
#' @export
wide_to_long <- function(data,
                         id_col = NULL,
                         time_prefix = "V",
                         action_col = "Action",
                         time_col = "Time",
                         drop_na = TRUE) {
  # Input validation
  stopifnot(is.data.frame(data))
  stopifnot(is.character(time_prefix), length(time_prefix) == 1)
  stopifnot(is.character(action_col), length(action_col) == 1)
  stopifnot(is.character(time_col), length(time_col) == 1)
  stopifnot(is.logical(drop_na), length(drop_na) == 1)

  # Identify time columns
  time_cols <- grep(paste0("^", time_prefix, "[0-9]+$"), names(data), value = TRUE)
  if (length(time_cols) == 0) {
    stop("No columns matching time prefix '", time_prefix, "' found.")
  }

  # Add ID column if needed
  if (is.null(id_col)) {
    data$id <- seq_len(nrow(data))
    id_col <- "id"
  } else if (!id_col %in% names(data)) {
    stop("ID column '", id_col, "' not found in data.")
  }

  # Identify other columns to preserve
  other_cols <- setdiff(names(data), c(time_cols, id_col))

  # Create the long format data
  result <- tidyr::pivot_longer(
    data,
    cols = dplyr::all_of(time_cols),
    names_to = time_col,
    values_to = action_col,
    names_prefix = time_prefix,
    names_transform = list(Time = as.integer)
  )

  # Drop NA values if requested
  if (drop_na) {
    result <- result[!is.na(result[[action_col]]), ]
  }

  # Arrange by ID and time
  result <- dplyr::arrange(result, .data[[id_col]], .data[[time_col]])

  # Reset row names

  rownames(result) <- NULL

  return(result)
}


#' Convert Long Format to Wide Sequences
#'
#' @description
#' Convert sequence data from long format (one row per action) to wide format
#' (one row per sequence, columns as time points).
#'
#' @param data Data frame in long format.
#' @param id_col Character. Name of the column identifying sequences.
#'   Default: "Actor".
#' @param time_col Character. Name of the column identifying time points.
#'   Default: "Time".
#' @param action_col Character. Name of the column containing actions/states.
#'   Default: "Action".
#' @param time_prefix Character. Prefix for time point columns in output.
#'   Default: "V".
#' @param fill_na Logical. Whether to fill missing time points with NA.
#'   Default: TRUE.
#'
#' @return A data frame in wide format where each row is a sequence and
#'   columns V1, V2, ... contain the actions at each time point.
#'
#' @details
#' This function converts long format data (like that from `simulate_long_data()`)
#' to the wide format expected by `tna::tna()` and related functions.
#'
#' If `time_col` contains non-integer values (e.g., timestamps), the function
#' will use the ordering within each sequence to create time indices.
#'
#' @examples
#' \dontrun{
#' # Generate long format data
#' long_data <- simulate_long_data(
#'   n_groups = 5,
#'   actors_per_group = 10,
#'   seed = 42
#' )
#'
#' # Convert to wide format
#' wide_data <- long_to_wide(long_data, id_col = "Actor")
#' head(wide_data)
#' }
#'
#' @seealso \code{\link{wide_to_long}} for the reverse conversion,
#'   \code{\link{prepare_for_tna}} for preparing data for TNA analysis.
#'
#' @importFrom tidyr pivot_wider
#' @importFrom dplyr mutate group_by row_number ungroup
#' @export
long_to_wide <- function(data,
                         id_col = "Actor",
                         time_col = "Time",
                         action_col = "Action",
                         time_prefix = "V",
                         fill_na = TRUE) {
  # Input validation
  stopifnot(is.data.frame(data))
  stopifnot(is.character(id_col), length(id_col) == 1)
  stopifnot(is.character(time_col), length(time_col) == 1)
  stopifnot(is.character(action_col), length(action_col) == 1)
  stopifnot(is.character(time_prefix), length(time_prefix) == 1)

  # Check required columns exist
  required_cols <- c(id_col, action_col)
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # Create time index within each sequence
  # This handles both integer time and timestamp columns
  data <- dplyr::group_by(data, .data[[id_col]])

  if (time_col %in% names(data)) {
    # Sort by existing time column first
    data <- dplyr::arrange(data, .data[[time_col]])
  }

  data <- dplyr::mutate(data, .time_idx = dplyr::row_number())
  data <- dplyr::ungroup(data)

  # Create column names
  data[[".time_name"]] <- paste0(time_prefix, data[[".time_idx"]])

  # Pivot to wide format
  result <- tidyr::pivot_wider(
    data,
    id_cols = dplyr::all_of(id_col),
    names_from = ".time_name",
    values_from = dplyr::all_of(action_col),
    values_fill = if (fill_na) list(NA_character_) else NULL
  )

  # Reorder columns to ensure V1, V2, V3... order
  time_cols <- grep(paste0("^", time_prefix, "[0-9]+$"), names(result), value = TRUE)
  time_nums <- as.integer(gsub(time_prefix, "", time_cols))
  time_cols <- time_cols[order(time_nums)]
  other_cols <- setdiff(names(result), time_cols)
  result <- result[, c(other_cols, time_cols), drop = FALSE]

  # Convert to data frame (not tibble)
  result <- as.data.frame(result)
  rownames(result) <- NULL

  return(result)
}


#' Prepare Data for TNA Analysis
#'
#' @description
#' Prepare simulated or real data for use with `tna::tna()` and related
#' functions. Handles various input formats and ensures the output is
#' compatible with TNA models.
#'
#' @param data Data frame containing sequence data.
#' @param type Character. Type of input data:
#' \describe{
#'   \item{"sequences"}{Wide format with one row per sequence (default).}
#'   \item{"long"}{Long format with one row per action.}
#'   \item{"auto"}{Automatically detect format based on column names.}
#' }
#' @param state_names Character vector. Expected state names, or NULL to
#'   extract from data. Default: NULL.
#' @param id_col Character. Name of ID column for long format data.
#'   Default: "Actor".
#' @param time_col Character. Name of time column for long format data.
#'   Default: "Time".
#' @param action_col Character. Name of action column for long format data.
#'   Default: "Action".
#' @param validate Logical. Whether to validate that all actions are in
#'   state_names. Default: TRUE.
#'
#' @return A data frame ready for use with TNA functions. For "sequences" type,
#'   returns a data frame where each row is a sequence and columns are time
#'   points (V1, V2, ...). For "long" type, converts to wide format first.
#'
#' @details
#' This function performs several preparations:
#' \enumerate{
#'   \item Converts long format to wide format if needed.
#'   \item Validates that all actions/states are recognized.
#'   \item Removes any non-sequence columns (e.g., id, metadata).
#'   \item Converts factors to characters.
#'   \item Ensures consistent column naming (V1, V2, ...).
#' }
#'
#' @examples
#' \dontrun{
#' # From wide format sequences
#' sequences <- simulate_sequences(
#'   transition_matrix = trans_mat,
#'   initial_probabilities = init_probs,
#'   max_seq_length = 15,
#'   num_rows = 100
#' )
#' tna_data <- prepare_for_tna(sequences, type = "sequences")
#' model <- tna::tna(tna_data)
#'
#' # From long format
#' long_data <- simulate_long_data(n_groups = 10, seed = 42)
#' tna_data <- prepare_for_tna(long_data, type = "long")
#' model <- tna::tna(tna_data)
#' }
#'
#' @seealso \code{\link{wide_to_long}}, \code{\link{long_to_wide}} for
#'   format conversions, \code{\link{simulate_sequences}} and
#'   \code{\link{simulate_long_data}} for generating simulated data.
#'
#' @export
prepare_for_tna <- function(data,
                            type = c("sequences", "long", "auto"),
                            state_names = NULL,
                            id_col = "Actor",
                            time_col = "Time",
                            action_col = "Action",
                            validate = TRUE) {
  # Input validation
  stopifnot(is.data.frame(data))
  type <- match.arg(type)

  # Auto-detect format
  if (type == "auto") {
    has_time_cols <- any(grepl("^V[0-9]+$", names(data)))
    has_action_col <- action_col %in% names(data)

    if (has_time_cols && !has_action_col) {
      type <- "sequences"
    } else if (has_action_col && !has_time_cols) {
      type <- "long"
    } else if (has_time_cols && has_action_col) {
      # If both, use long if nrow >> number of unique IDs
      if (id_col %in% names(data)) {
        n_ids <- length(unique(data[[id_col]]))
        if (nrow(data) > n_ids * 2) {
          type <- "long"
        } else {
          type <- "sequences"
        }
      } else {
        type <- "sequences"
      }
    } else {
      stop("Cannot auto-detect data format. Please specify 'type' explicitly.")
    }
  }

  # Convert long to wide if needed
  if (type == "long") {
    if (!action_col %in% names(data)) {
      stop("Action column '", action_col, "' not found in data.")
    }
    if (!id_col %in% names(data)) {
      stop("ID column '", id_col, "' not found in data.")
    }

    data <- long_to_wide(
      data,
      id_col = id_col,
      time_col = time_col,
      action_col = action_col
    )
  }

  # Identify sequence columns (V1, V2, ...)
  seq_cols <- grep("^V[0-9]+$", names(data), value = TRUE)
  if (length(seq_cols) == 0) {
    stop("No sequence columns (V1, V2, ...) found in data.")
  }

  # Extract just the sequence columns
  result <- data[, seq_cols, drop = FALSE]

  # Convert factors to characters
  for (col in names(result)) {
    if (is.factor(result[[col]])) {
      result[[col]] <- as.character(result[[col]])
    }
  }

  # Validate state names if requested
  if (validate && !is.null(state_names)) {
    all_states <- unique(unlist(result, use.names = FALSE))
    all_states <- all_states[!is.na(all_states)]
    unknown_states <- setdiff(all_states, state_names)
    if (length(unknown_states) > 0) {
      warning("Unknown states found: ", paste(unknown_states, collapse = ", "))
    }
  }

  # Ensure proper column ordering
  col_nums <- as.integer(gsub("V", "", seq_cols))
  result <- result[, order(col_nums), drop = FALSE]

  # Convert to data frame
  result <- as.data.frame(result)
  rownames(result) <- NULL

  return(result)
}


#' Convert Action Column to One-Hot Encoding
#'
#' @description
#' Convert a categorical Action column to one-hot (binary indicator) columns.
#'
#' @param data Data frame containing an action column.
#' @param action_col Character. Name of the action column. Default: "Action".
#' @param states Character vector or NULL. States to include as columns.
#'   If NULL, uses all unique values. Default: NULL.
#' @param drop_action Logical. Remove the original action column. Default: TRUE.
#' @param sort_states Logical. Sort state columns alphabetically. Default: FALSE.
#' @param prefix Character. Prefix for state column names. Default: "".
#'
#' @return Data frame with one-hot encoded columns (0/1 integers).
#'
#' @examples
#' \dontrun{
#' # Generate long format data
#' long_data <- simulate_long_data(n_groups = 3, seed = 42)
#'
#' # Convert to one-hot encoding
#' onehot_data <- action_to_onehot(long_data)
#' head(onehot_data)
#'
#' # With custom prefix
#' onehot_data <- action_to_onehot(long_data, prefix = "state_")
#' names(onehot_data)
#' }
#'
#' @seealso \code{\link{simulate_onehot_data}} for directly generating one-hot data,
#'   \code{\link{simulate_long_data}} for generating long format data.
#'
#' @export
action_to_onehot <- function(data, action_col = "Action", states = NULL,
                             drop_action = TRUE, sort_states = FALSE,
                             prefix = "") {
  stopifnot(is.data.frame(data))
  stopifnot(action_col %in% names(data))

  if (is.null(states)) {
    states <- unique(data[[action_col]])
    states <- states[!is.na(states)]
  }

  if (sort_states) {
    states <- sort(states)
  }

  for (state in states) {
    col_name <- paste0(prefix, state)
    data[[col_name]] <- as.integer(data[[action_col]] == state)
  }

  if (drop_action) {
    data[[action_col]] <- NULL
  }

  return(data)
}
