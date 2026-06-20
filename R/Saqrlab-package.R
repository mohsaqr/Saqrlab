#' Saqrlab: A Modern Laboratory for Data Simulation
#'
#' @description
#' Saqrlab generates synthetic datasets with known ground-truth parameters so
#' you can test and validate statistical methods. It spans classic designs
#' (t-test, ANOVA, regression, factor analysis), modern data-generating
#' processes (IRT, survival, multilevel, growth curves, hidden Markov,
#' missing-data mechanisms), and a full Temporal Network Analysis toolkit.
#'
#' @details
#' Two complementary tiers:
#' \itemize{
#'   \item \strong{Explicit-parameter} simulators (e.g. [simulate_ttest()],
#'     [simulate_irt()], [simulate_mlm()], and the [simulate()] dispatcher)
#'     return a [saqr_sim] object with `$data` and the true `$params`, for
#'     parameter-recovery testing.
#'   \item \strong{Random-parameter} generation via [simulate_data()] returns a
#'     bare data frame whose structure is invented from the seed, for
#'     stress/robustness testing.
#' }
#'
#' Start with [list_simulators()] for the catalogue, [simulate()] to generate,
#' and [validate_recovery()] to score whether a fitted method recovered the
#' truth. Browse the per-category HTML reference in the package's
#' `reference/` directory.
#'
#' @keywords internal
#' @importFrom utils head tail str
#' @importFrom stats plogis rlnorm
"_PACKAGE"

# Global variables to avoid R CMD check notes (ggplot2/dplyr NSE)
utils::globalVariables(c(
  "density", "mean_val", "mean_value", "median_val", "sd_val", "sd_value",
  "value", "metric", "category", "model_type", "iteration", "mean", "sd",
  "sim_alpha", "sim_diag_c",
  ".grp_key", ".", ".SD", ":=",
  "cr_lower", "cr_upper", "weight_x", "weight_y", ".seq_grp",
  # dplyr/NSE variables in run_bootstrap_*/run_network_simulation/summarize_grid_results
  "FN", "FNR", "FP", "FPR", "MCC", "Metric", "TN", "TP", "Value",
  "from", "to", "max_na", "min_na", "run_id", "p_value", "weight",
  "Accuracy", "Sensitivity", "Specificity", "Total_FN", "Total_FP",
  "Total_TN", "Total_TP", "avg_p_value", "avg_recovery_rate", "avg_weight",
  "bootstrap_significant_run", "comparison_type", "data_idx",
  "ground_truth_stable", "gt_stable", "is_significant", "max_seq_length",
  "mcc_denom_sq", "metric_category", "metric_name", "n_significant",
  "num_rows", "num_states", "p_value_num", "recovery_rate",
  "safe_mean", "safe_median", "safe_sd", "successful_runs", "weight_num"
))


#' @noRd
.onLoad <- function(libname, pkgname) invisible(NULL)
