#' Batch Processing Functions
#'
#' @description
#' Functions for fitting models and applying operations to multiple datasets.
#'
#' @name batch_processing
#' @keywords internal
NULL

#' Fit Models to Multiple Datasets
#'
#' @description
#' Fit TNA models to multiple datasets at once, with optional parallel
#' processing.
#'
#' @param data_list A list of data frames, each containing sequence data.
#' @param model_type Character. Type of model to fit: "tna", "ftna", "ctna",
#'   or "atna". Default: "tna".
#' @param parallel Logical. Whether to use parallel processing. Default: FALSE.
#' @param cores Integer or NULL. Number of cores to use for parallel processing.
#'   If NULL, uses `parallel::detectCores() - 1`. Ignored if `parallel = FALSE`.
#'   Default: NULL.
#' @param progress Logical. Whether to show progress messages. Default: TRUE.
#' @param ... Additional arguments passed to the model fitting function.
#'
#' @return A list of fitted model objects, with the same names/indices as
#'   `data_list`. Failed fits are returned as NULL with a warning.
#'
#' @details
#' This function provides a convenient way to fit TNA models to many datasets,
#' such as when running simulations or analyzing multiple groups.
#'
#' For parallel processing on Windows, set up a parallel backend first using
#' the `future` package.
#'
#' @examples
#' \dontrun{
#' # Generate multiple datasets
#' datasets <- lapply(1:10, function(i) {
#'   simulate_sequences(trans_mat, init_probs, max_seq_length = 20, num_rows = 100)
#' })
#'
#' # Fit models in sequence
#' models <- batch_fit_models(datasets, model_type = "tna")
#'
#' # Fit models in parallel
#' models <- batch_fit_models(datasets, model_type = "tna",
#'                            parallel = TRUE, cores = 4)
#'
#' # Check results
#' sapply(models, function(m) if (!is.null(m)) "OK" else "Failed")
#' }
#'
#' @seealso \code{\link{fit_network_model}} for fitting a single model,
#'   \code{\link{batch_apply}} for applying functions to model lists.
#'
#' @importFrom parallel detectCores mclapply
#' @importFrom future.apply future_lapply
#' @importFrom future plan
#' @import tna
#' @export
batch_fit_models <- function(data_list,
                             model_type = c("tna", "ftna", "ctna", "atna"),
                             parallel = FALSE,
                             cores = NULL,
                             progress = TRUE,
                             ...) {
  # Input validation
  if (!is.list(data_list)) {
    stop("data_list must be a list of data frames.")
  }
  model_type <- match.arg(model_type)

  n_datasets <- length(data_list)
  if (n_datasets == 0) {
    return(list())
  }

  # Set up cores
  if (parallel && is.null(cores)) {
    cores <- max(1, parallel::detectCores() - 1)
  }

  if (progress) {
    message(sprintf("Fitting %d %s models%s...",
                    n_datasets, model_type,
                    if (parallel) sprintf(" using %d cores", cores) else ""))
  }

  # Define the fitting function
  fit_one <- function(data, idx) {
    tryCatch({
      fit_network_model(data, model_type)
    }, error = function(e) {
      warning(sprintf("Failed to fit model for dataset %d: %s", idx, e$message))
      NULL
    })
  }

  # Run fitting
  if (parallel && cores > 1) {
    # Use future.apply for cross-platform parallel processing
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    future::plan("multisession", workers = cores)

    results <- future.apply::future_lapply(
      seq_along(data_list),
      function(i) fit_one(data_list[[i]], i),
      future.seed = TRUE
    )
  } else {
    # Sequential processing
    results <- lapply(seq_along(data_list), function(i) {
      if (progress && i %% 10 == 0) {
        message(sprintf("  Fitted %d/%d models", i, n_datasets))
      }
      fit_one(data_list[[i]], i)
    })
  }

  # Preserve names if data_list was named
  if (!is.null(names(data_list))) {
    names(results) <- names(data_list)
  }

  # Count successes
  n_success <- sum(sapply(results, function(x) !is.null(x)))
  if (progress) {
    message(sprintf("Successfully fitted %d/%d models", n_success, n_datasets))
  }

  return(results)
}


#' Apply Function to Multiple Models or Datasets
#'
#' @description
#' Apply a function to a list of TNA models or datasets, with optional
#' parallel processing.
#'
#' @param object_list A list of TNA model objects or data frames.
#' @param fun A function to apply to each element.
#' @param parallel Logical. Whether to use parallel processing. Default: FALSE.
#' @param cores Integer or NULL. Number of cores for parallel processing.
#'   Default: NULL.
#' @param progress Logical. Whether to show progress messages. Default: TRUE.
#' @param simplify Logical. Whether to simplify results if possible.
#'   Default: FALSE.
#' @param ... Additional arguments passed to `fun`.
#'
#' @return A list of results from applying `fun` to each element.
#'   If `simplify = TRUE` and results are atomic, returns a simplified vector.
#'
#' @details
#' This is a general-purpose function for batch operations on model lists.
#' It handles errors gracefully, returning NULL for failed operations.
#'
#' @examples
#' \dontrun{
#' # Extract transition matrices from multiple models
#' trans_mats <- batch_apply(models, extract_transition_matrix)
#'
#' # Compare each model to a reference
#' ref_model <- models[[1]]
#' comparisons <- batch_apply(
#'   models[-1],
#'   function(m) compare_networks(ref_model, m)$metrics$correlation
#' )
#'
#' # Extract centralities with simplification
#' centralities <- batch_apply(
#'   models,
#'   function(m) tna::centralities(m),
#'   simplify = FALSE
#' )
#' }
#'
#' @seealso \code{\link{batch_fit_models}} for fitting multiple models.
#'
#' @importFrom parallel detectCores
#' @importFrom future.apply future_lapply
#' @importFrom future plan
#' @export
batch_apply <- function(object_list,
                        fun,
                        parallel = FALSE,
                        cores = NULL,
                        progress = TRUE,
                        simplify = FALSE,
                        ...) {
  # Input validation
  if (!is.list(object_list)) {
    stop("object_list must be a list.")
  }
  if (!is.function(fun)) {
    stop("fun must be a function.")
  }

  n_objects <- length(object_list)
  if (n_objects == 0) {
    return(list())
  }

  # Set up cores
  if (parallel && is.null(cores)) {
    cores <- max(1, parallel::detectCores() - 1)
  }

  if (progress) {
    message(sprintf("Applying function to %d objects%s...",
                    n_objects,
                    if (parallel) sprintf(" using %d cores", cores) else ""))
  }

  # Define wrapper function
  apply_one <- function(obj, idx) {
    tryCatch({
      fun(obj, ...)
    }, error = function(e) {
      warning(sprintf("Function failed for object %d: %s", idx, e$message))
      NULL
    })
  }

  # Run application
  if (parallel && cores > 1) {
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    future::plan("multisession", workers = cores)

    results <- future.apply::future_lapply(
      seq_along(object_list),
      function(i) apply_one(object_list[[i]], i),
      future.seed = TRUE
    )
  } else {
    results <- lapply(seq_along(object_list), function(i) {
      if (progress && i %% 10 == 0) {
        message(sprintf("  Processed %d/%d objects", i, n_objects))
      }
      apply_one(object_list[[i]], i)
    })
  }

  # Preserve names
  if (!is.null(names(object_list))) {
    names(results) <- names(object_list)
  }

  # Simplify if requested and possible
  if (simplify) {
    is_atomic <- all(sapply(results, function(x) {
      is.atomic(x) && length(x) == 1
    }))
    if (is_atomic) {
      results <- unlist(results)
    }
  }

  return(results)
}
