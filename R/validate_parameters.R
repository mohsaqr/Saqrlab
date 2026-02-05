#' Validate and Standardize Simulation Parameters
#'
#' @description
#' Validate simulation parameters and set defaults for missing values.
#' Ensures parameter relationships are valid (e.g., min_na <= max_na).
#'
#' @param params List of simulation parameters. May include:
#' \describe{
#'   \item{seq_length}{Maximum sequence length (default: 30).}
#'   \item{n_sequences}{Number of sequences to generate (default: 100).}
#'   \item{min_na}{Minimum number of NA values per sequence (default: 0).}
#'   \item{max_na}{Maximum number of NA values per sequence (default: 5).}
#' }
#'
#' @return A list of validated parameters with defaults applied.
#'
#' @details
#' The function performs the following validations:
#' \itemize{
#'   \item Sets default values for missing parameters.
#'   \item Ensures `max_na <= seq_length`.
#'   \item Ensures `min_na <= max_na`.
#'   \item Ensures `min_na < seq_length`.
#' }
#'
#' Both new (`seq_length`, `n_sequences`) and old (`max_seq_length`, `num_rows`)
#' parameter names are supported for backward compatibility.
#'
#' @examples
#' # Validate with defaults
#' params <- validate_sim_params(list())
#'
#' # Validate custom parameters (new names)
#' params <- validate_sim_params(list(
#'   seq_length = 50,
#'   n_sequences = 200,
#'   min_na = 2,
#'   max_na = 10
#' ))
#'
#' # Old parameter names still work
#' params <- validate_sim_params(list(
#'   max_seq_length = 50,
#'   num_rows = 200
#' ))
#'
#' @export
validate_sim_params <- function(params) {
  # --- Backward compatibility: map old names to new names ---
  if (!is.null(params$max_seq_length) && is.null(params$seq_length)) {
    params$seq_length <- params$max_seq_length
  }
  if (!is.null(params$num_rows) && is.null(params$n_sequences)) {
    params$n_sequences <- params$num_rows
  }

  # Set defaults for missing parameters (use new standardized names)
  defaults <- list(
    seq_length = 30,
    n_sequences = 100,
    min_na = 0,
    max_na = 5
  )

  # Apply defaults for any missing parameters
  for (param_name in names(defaults)) {
    if (is.null(params[[param_name]])) {
      params[[param_name]] <- defaults[[param_name]]
    }
  }

  # Validate parameter relationships
  if (params$max_na > params$seq_length) {
    params$max_na <- params$seq_length
  }

  if (params$min_na > params$max_na) {
    params$min_na <- params$max_na
  }

  # Ensure min_na is less than seq_length
  if (params$min_na >= params$seq_length) {
    params$min_na <- max(0, params$seq_length - 1) # Ensure it's at least 0
  }

  # Ensure max_na is not greater than seq_length (redundant but for clarity)
  params$max_na <- min(params$max_na, params$seq_length)

  # Also set old names for backward compatibility with callers expecting them

  params$max_seq_length <- params$seq_length
  params$num_rows <- params$n_sequences

  return(params)
}
