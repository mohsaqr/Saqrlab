#' @keywords internal
"_PACKAGE"

# Global variables to avoid R CMD check notes (ggplot2/dplyr NSE)
utils::globalVariables(c(
  "density", "mean_val", "mean_value", "median_val", "sd_val", "sd_value",
  "value", "metric", "category", "model_type", "iteration", "mean", "sd"
))
