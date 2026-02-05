#' @keywords internal
"_PACKAGE"

#' Saqrlab: Simulation and Analysis Tools for Temporal Network Analysis
#'
#' @description
#' A comprehensive toolkit for simulating, analyzing, and comparing Temporal
#' Network Analysis (TNA) models. Designed for educational research with
#' built-in support for realistic learning state names and hierarchical data
#' structures.
#'
#' @section Data Simulation:
#' Functions for generating synthetic TNA data:
#' \itemize{
#'   \item \code{\link{simulate_matrix}} - Simple transition matrices
#'   \item \code{\link{simulate_htna}} - Multi-type matrices (HTNA/MLNA/MTNA)
#'   \item \code{\link{simulate_sequences}} - Markov chain sequences
#'   \item \code{\link{simulate_sequences_advanced}} - Sequences with stability modes
#'   \item \code{\link{simulate_long_data}} - Hierarchical group data
#'   \item \code{\link{simulate_onehot_data}} - One-hot encoded sequences
#'   \item \code{\link{simulate_edge_list}} - Social network edge lists
#' }
#'
#' @section Network Generation:
#' Functions for creating complete TNA datasets and networks:
#' \itemize{
#'   \item \code{\link{generate_tna_datasets}} - Complete datasets with parameters
#'   \item \code{\link{generate_tna_networks}} - Fitted TNA models
#'   \item \code{\link{generate_group_tna_networks}} - Group TNA models
#'   \item \code{\link{generate_tna_matrix}} - HTNA/MLNA transition matrices
#'   \item \code{\link{generate_probabilities}} - Random transition probabilities
#' }
#'
#' @section Model Fitting & Extraction:
#' Functions for fitting models and extracting components:
#' \itemize{
#'   \item \code{\link{fit_network_model}} - Fit TNA/fTNA/cTNA/aTNA models
#'   \item \code{\link{extract_transition_matrix}} - Get transition matrix
#'   \item \code{\link{extract_initial_probs}} - Get initial probabilities
#'   \item \code{\link{extract_edges}} - Get edge list
#' }
#'
#' @section Network Comparison:
#' Functions for comparing and evaluating networks:
#' \itemize{
#'   \item \code{\link{compare_networks}} - Compare two networks
#'   \item \code{\link{compare_centralities}} - Compare centrality profiles
#'   \item \code{\link{calculate_edge_recovery}} - Edge recovery metrics
#' }
#'
#' @section Data Conversion:
#' Functions for transforming data formats:
#' \itemize{
#'   \item \code{\link{wide_to_long}} - Wide to long format
#'   \item \code{\link{long_to_wide}} - Long to wide format
#'   \item \code{\link{prepare_for_tna}} - Prepare for tna package
#'   \item \code{\link{action_to_onehot}} - Convert to one-hot encoding
#' }
#'
#' @section Batch Processing:
#' Functions for processing multiple datasets:
#' \itemize{
#'   \item \code{\link{batch_fit_models}} - Fit models to multiple datasets
#'   \item \code{\link{batch_apply}} - Apply function to list of objects
#' }
#'
#' @section Bootstrap & Simulation:
#' Functions for simulation studies:
#' \itemize{
#'   \item \code{\link{run_bootstrap_simulation}} - Bootstrap analysis
#'   \item \code{\link{run_grid_simulation}} - Parameter grid search
#'   \item \code{\link{run_network_simulation}} - Model comparison studies
#'   \item \code{\link{evaluate_bootstrap}} - Evaluate bootstrap results
#'   \item \code{\link{analyze_grid_results}} - Analyze grid output
#' }
#'
#' @section Learning States:
#' Datasets and functions for educational research:
#' \itemize{
#'   \item \code{\link{LEARNING_STATES}} - 180+ learning verbs in 8 categories
#'   \item \code{\link{get_learning_states}} - Get verbs by category
#'   \item \code{\link{list_learning_categories}} - Show available categories
#'   \item \code{\link{smart_select_states}} - Intelligent state selection
#'   \item \code{\link{GLOBAL_NAMES}} - 300 diverse names for simulation
#'   \item \code{\link{get_global_names}} - Retrieve global names
#' }
#'
#' @section Utilities:
#' Helper functions:
#' \itemize{
#'   \item \code{\link{create_param_grid}} - Create parameter grids
#'   \item \code{\link{validate_sim_params}} - Validate parameters
#'   \item \code{\link{summarize_simulation}} - Summary statistics
#'   \item \code{\link{summarize_networks}} - Network summaries
#' }
#'
#' @section Getting Started:
#' Start with a simple workflow:
#' \preformatted{
#' library(Saqrlab)
#' library(tna)
#'
#' # Generate sequences with learning state names
#' sequences <- simulate_sequences(
#'   n_sequences = 100,
#'   seq_length = 20,
#'   n_states = 5,
#'   seed = 42
#' )
#'
#' # Fit a TNA model
#' model <- fit_network_model(sequences, "tna")
#'
#' # Visualize
#' plot(model)
#' }
#'
#' @seealso
#' \itemize{
#'   \item \code{\link[tna:tna]{tna}} package for Temporal Network Analysis
#'   \item \code{vignette("introduction", package = "Saqrlab")} for getting started
#'   \item \code{vignette("simulation-guide", package = "Saqrlab")} for simulation details
#' }
#'
#' @author Mohammed Saqr \email{saqr@@saqr.me}
#'
#' @docType package
#' @name Saqrlab-package
#' @aliases Saqrlab
NULL
