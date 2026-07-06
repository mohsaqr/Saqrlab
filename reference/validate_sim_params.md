# Validate and Standardize Simulation Parameters

Validate simulation parameters and set defaults for missing values.
Ensures parameter relationships are valid (e.g., min_na \<= max_na).

## Usage

``` r
validate_sim_params(params)
```

## Arguments

- params:

  List of simulation parameters. May include:

  seq_length

  :   Maximum sequence length (default: 30).

  n_sequences

  :   Number of sequences to generate (default: 100).

  min_na

  :   Minimum number of NA values per sequence (default: 0).

  max_na

  :   Maximum number of NA values per sequence (default: 5).

## Value

A list of validated parameters with defaults applied.

## Details

The function performs the following validations:

- Sets default values for missing parameters.

- Ensures `max_na <= seq_length`.

- Ensures `min_na <= max_na`.

- Ensures `min_na < seq_length`.

Both new (`seq_length`, `n_sequences`) and old (`max_seq_length`,
`num_rows`) parameter names are supported for backward compatibility.

## Examples

``` r
# Validate with defaults
params <- validate_sim_params(list())

# Validate custom parameters (new names)
params <- validate_sim_params(list(
  seq_length = 50,
  n_sequences = 200,
  min_na = 2,
  max_na = 10
))

# Old parameter names still work
params <- validate_sim_params(list(
  max_seq_length = 50,
  num_rows = 200
))
```
