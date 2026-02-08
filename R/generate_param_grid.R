#' Generate Parameter Grid for Simulations
#'
#' @description
#' Generate a grid of parameter combinations for simulation studies using
#' various sampling methods.
#'
#' @param param_ranges Named list of parameter ranges. Each element can be:
#' \itemize{
#'   \item A numeric vector of length 2 (min, max) for continuous/integer parameters.
#'   \item A vector of values for categorical parameters.
#' }
#' @param n Integer. Number of parameter combinations to generate. Default: 10.
#' @param method Character. Sampling method. One of:
#' \describe{
#'   \item{"random"}{Random uniform sampling within ranges.}
#'   \item{"grid"}{Regular grid sampling (may exceed n, then subsampled).}
#'   \item{"lhs"}{Latin Hypercube Sampling for better coverage (requires lhs package).}
#' }
#' Default: "random".
#'
#' @return A data frame with `n` rows and columns for each parameter.
#'
#' @details
#' **Method Details:**
#'
#' - "random": Draws uniform random values within each range. Integer parameters
#'   (detected when min and max are both integers) are rounded.
#'
#' - "grid": Creates a regular grid with approximately `n^(1/d)` points per
#'   dimension (where d is the number of parameters). If the resulting grid
#'   exceeds n points, it is randomly subsampled.
#'
#' - "lhs": Uses Latin Hypercube Sampling for space-filling designs that
#'   provide better coverage of the parameter space than random sampling.
#'
#' Categorical parameters are sampled uniformly with replacement for all methods.
#'
#' @examples
#' # Define parameter ranges
#' ranges <- list(
#'   num_rows = c(50, 500),       # Integer parameter
#'   max_seq_length = c(10, 100), # Integer parameter
#'   stability_prob = c(0.7, 1.0) # Continuous parameter
#' )
#'
#' # Random sampling
#' grid_random <- generate_param_grid(ranges, n = 20, method = "random")
#'
#' # Latin Hypercube Sampling
#' grid_lhs <- generate_param_grid(ranges, n = 20, method = "lhs")
#'
#' # Grid sampling
#' grid_regular <- generate_param_grid(ranges, n = 20, method = "grid")
#'
#' @importFrom lhs randomLHS
#' @importFrom stats runif
#' @export
generate_param_grid <- function(param_ranges, n = 10, method = "random") {
  if (method == "lhs" && requireNamespace("lhs", quietly = TRUE)) {
    # Latin Hypercube Sampling for better coverage of parameter space
    numeric_params <- names(param_ranges)[sapply(param_ranges, function(x)
      is.numeric(x) && length(x) == 2)]

    if (length(numeric_params) > 0) {
      # Generate Latin Hypercube sample
      lhs_design <- lhs::randomLHS(n, length(numeric_params))

      # Scale to parameter ranges
      param_grid <- data.frame(matrix(nrow = n, ncol = length(numeric_params)))
      names(param_grid) <- numeric_params

      for (i in 1:length(numeric_params)) {
        param_name <- numeric_params[i]
        range <- param_ranges[[param_name]]
        param_grid[[param_name]] <- range[1] + lhs_design[, i] * (range[2] - range[1])

        # Round to integers if needed
        if (all(range == round(range))) {
          param_grid[[param_name]] <- round(param_grid[[param_name]])
        }
      }

      # Add non-numeric parameters
      for (param_name in setdiff(names(param_ranges), numeric_params)) {
        values <- param_ranges[[param_name]]
        param_grid[[param_name]] <- sample(values, n, replace = TRUE)
      }

      return(param_grid)
    }
  } else if (method == "grid") {
    # Create a regular grid with approximately n points
    numeric_params <- names(param_ranges)[sapply(param_ranges, function(x)
      is.numeric(x) && length(x) == 2)]

    if (length(numeric_params) > 0) {
      # Determine number of points per dimension
      points_per_dim <- ceiling(n^(1 / length(numeric_params)))

      # Generate grid for each parameter
      grid_lists <- list()
      for (param_name in numeric_params) {
        range <- param_ranges[[param_name]]
        if (all(range == round(range))) {
          # Integer parameter
          step <- max(1, ceiling((range[2] - range[1]) / (points_per_dim - 1)))
          grid_lists[[param_name]] <- seq(range[1], range[2], by = step)
        } else {
          # Continuous parameter
          grid_lists[[param_name]] <- seq(range[1], range[2], length.out = points_per_dim)
        }
      }

      # Create all combinations
      param_grid <- expand.grid(grid_lists)

      # Subsample if too many combinations
      if (nrow(param_grid) > n) {
        param_grid <- param_grid[sample(1:nrow(param_grid), n), ]
      }

      # Add non-numeric parameters
      for (param_name in setdiff(names(param_ranges), numeric_params)) {
        values <- param_ranges[[param_name]]
        param_grid[[param_name]] <- sample(values, nrow(param_grid), replace = TRUE)
      }

      return(param_grid)
    }
  }

  # Default to random sampling if other methods not applicable
  param_grid <- data.frame(matrix(nrow = n, ncol = length(param_ranges)))
  names(param_grid) <- names(param_ranges)

  for (param_name in names(param_ranges)) {
    range <- param_ranges[[param_name]]
    if (is.numeric(range) && length(range) == 2) {
      # Numeric range
      values <- runif(n, range[1], range[2])
      if (all(range == round(range))) {
        values <- round(values)
      }
      param_grid[[param_name]] <- values
    } else {
      # Categorical or specific set of values
      param_grid[[param_name]] <- sample(range, n, replace = TRUE)
    }
  }

  return(param_grid)
}
