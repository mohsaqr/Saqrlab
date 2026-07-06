# Sample and Re-estimate TNA Model

Take a TNA model, sample a percentage of its data, and re-estimate a new
TNA model from the sampled data using the same parameters.

## Usage

``` r
sample_tna(model, sampling_percent = 0.3, model_type = NULL, scaling = NULL)
```

## Arguments

- model:

  A TNA model object or data frame containing sequence data.

- sampling_percent:

  Numeric value between 0 and 1 indicating the proportion of data to
  sample. Default: 0.3.

- model_type:

  Character string specifying the model type. If NULL, uses the type
  from the original model. Default: NULL.

- scaling:

  Character string specifying the scaling method. If NULL, uses the
  scaling from the original model. Default: NULL.

## Value

A new TNA model object estimated from the sampled data.

## Examples

``` r
if (FALSE) { # \dontrun{
library(tna)
model <- tna(group_regulation)
sampled <- sample_tna(model, 0.3)
} # }
```
