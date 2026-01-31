#' Validate and Standardize Simulation Parameters
#'
#' @description
#' Validate simulation parameters and set defaults for missing values.
#' Ensures parameter relationships are valid (e.g., min_na <= max_na).
#'
#' @param params List of simulation parameters. May include:
#' \describe{
#'   \item{max_seq_length}{Maximum sequence length (default: 30).}
#'   \item{num_rows}{Number of sequences to generate (default: 100).}
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
#'   \item Ensures `max_na <= max_seq_length`.
#'   \item Ensures `min_na <= max_na`.
#'   \item Ensures `min_na < max_seq_length`.
#' }
#'
#' @examples
#' # Validate with defaults
#' params <- validate_sim_params(list())
#'
#' # Validate custom parameters
#' params <- validate_sim_params(list(
#'   max_seq_length = 50,
#'   num_rows = 200,
#'   min_na = 2,
#'   max_na = 10
#' ))
#'
#' @export
validate_sim_params <- function(params) {
  # Set defaults for missing parameters
  defaults <- list(
    max_seq_length = 30,
    num_rows = 100,
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
  if (params$max_na > params$max_seq_length) {
    params$max_na <- params$max_seq_length
  }

  if (params$min_na > params$max_na) {
    params$min_na <- params$max_na
  }

  # Ensure min_na is less than max_seq_length
  if (params$min_na >= params$max_seq_length) {
    params$min_na <- max(0, params$max_seq_length - 1) # Ensure it's at least 0
  }

  # Ensure max_na is not greater than max_seq_length (redundant but for clarity)
  params$max_na <- min(params$max_na, params$max_seq_length)

  return(params)
}
