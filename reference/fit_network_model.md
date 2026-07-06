# Fit a Temporal Network Analysis Model

Fit a TNA (Temporal Network Analysis) model to sequence data using one
of several available model types.

## Usage

``` r
fit_network_model(sequences, model_type, group = NULL)
```

## Arguments

- sequences:

  Data frame of sequences. Each row is a sequence, each column is a time
  point. For group models, should include a grouping column.

- model_type:

  Character. Type of model to fit. One of:

  "tna"

  :   Standard Temporal Network Analysis model.

  "ftna"

  :   Filtered TNA model.

  "ctna"

  :   Conditional TNA model.

  "atna"

  :   Aggregated TNA model.

  "group_tna"

  :   Group TNA model for multi-group analysis.

- group:

  Character or NULL. Name of the grouping column for group_tna models.
  Required when `model_type = "group_tna"`. Default: NULL.

## Value

A fitted TNA model object of the appropriate class.

## Details

This function is a wrapper around the tna package's model fitting
functions. It provides a unified interface for fitting different types
of TNA models.

For group_tna models, the `sequences` data frame must contain a column
identifying the group membership of each sequence, specified by the
`group` parameter.

## Examples

``` r
if (FALSE) { # \dontrun{
# Generate some sequence data
trans_mat <- matrix(c(
  0.7, 0.2, 0.1,
  0.3, 0.5, 0.2,
  0.2, 0.3, 0.5
), nrow = 3, byrow = TRUE)
rownames(trans_mat) <- colnames(trans_mat) <- c("A", "B", "C")
init_probs <- c(A = 0.5, B = 0.3, C = 0.2)

sequences <- simulate_sequences(
  transition_matrix = trans_mat,
  initial_probabilities = init_probs,
  max_seq_length = 20,
  num_rows = 100
)

# Fit different model types
model_tna <- fit_network_model(sequences, "tna")
model_ftna <- fit_network_model(sequences, "ftna")

# Fit group model
long_data <- simulate_long_data(n_groups = 5, actors_per_group = 10)
wide_data <- long_to_wide(long_data, id_col = "Actor")
# Add group column back
actor_groups <- unique(long_data[, c("Actor", "Group")])
wide_data <- merge(wide_data, actor_groups, by = "Actor")
model_group <- fit_network_model(wide_data, "group_tna", group = "Group")
} # }
```
