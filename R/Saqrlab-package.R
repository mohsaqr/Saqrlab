#' @keywords internal
"_PACKAGE"

# Global variables to avoid R CMD check notes (ggplot2/dplyr NSE)
utils::globalVariables(c(
  "density", "mean_val", "mean_value", "median_val", "sd_val", "sd_value",
  "value", "metric", "category", "model_type", "iteration", "mean", "sd",
  "sim_alpha", "sim_diag_c",
  ".grp_key", ".", ".SD", ":="
))


#' @noRd
.onLoad <- function(libname, pkgname) {
  .register_builtin_estimators()
}
