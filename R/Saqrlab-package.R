#' @keywords internal
"_PACKAGE"

# Global variables to avoid R CMD check notes (ggplot2/dplyr NSE)
utils::globalVariables(c(
  "density", "mean_val", "mean_value", "median_val", "sd_val", "sd_value",
  "value", "metric", "category", "model_type", "iteration", "mean", "sd",
  "sim_alpha", "sim_diag_c",
  ".grp_key", ".", ".SD", ":=",
  "cr_lower", "cr_upper", "weight_x", "weight_y", ".seq_grp"
))


#' @noRd
.onLoad <- function(libname, pkgname) {
  .register_builtin_estimators()
}
