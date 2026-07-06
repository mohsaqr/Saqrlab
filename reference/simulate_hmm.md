# Simulate Hidden Markov Model Sequences

Generate observed symbol sequences from a discrete hidden Markov model
(HMM). For each sequence a latent state path is drawn from the Markov
chain defined by `init` (initial distribution) and `trans` (state
transition matrix); each hidden state then emits an observed symbol
according to the `emission` matrix. The true latent paths are returned
in `$params$hidden_paths` for parameter-recovery testing.

## Usage

``` r
simulate_hmm(
  n_sequences = 50,
  seq_length = 30,
  n_states = 2,
  n_symbols = 3,
  trans = NULL,
  emission = NULL,
  init = NULL,
  seed = NULL
)
```

## Arguments

- n_sequences:

  Integer. Number of sequences to simulate. Default: 50.

- seq_length:

  Integer. Length (number of time points) of each sequence. Default: 30.

- n_states:

  Integer. Number of hidden states. Default: 2.

- n_symbols:

  Integer. Number of observable symbols. Default: 3.

- trans:

  Numeric matrix (`n_states x n_states`) or NULL. Row- stochastic state
  transition matrix. When NULL, a diagonally dominant matrix (states
  persist) is auto-generated.

- emission:

  Numeric matrix (`n_states x n_symbols`) or NULL. Row- stochastic
  emission matrix. When NULL, an auto-generated matrix where each state
  favours a distinct symbol is used.

- init:

  Numeric vector (length `n_states`) or NULL. Initial state
  distribution. When NULL, a uniform distribution is used.

- seed:

  Integer or NULL. Random seed.

## Value

A [`saqr_sim`](https://pak.dynasite.org/Saqrlab/reference/saqr_sim.md)
object with:

- `$data`:

  long-format data.frame with columns `sequence_id` (integer), `time`
  (integer, 1...seq_length), and `symbol` (integer observed symbol in
  1...n_symbols).

- `$params`:

  list with `trans`, `emission`, `init`, `n_states`, `n_symbols`, and
  `hidden_paths` (an `n_sequences x seq_length` integer matrix of the
  TRUE latent state sequences).

## Details

The within-sequence latent recursion is inherently sequential and is
implemented with `Reduce(..., accumulate = TRUE)` over the time index;
the work across sequences is vectorised with `Map`/`vapply`.

## Examples

``` r
r <- simulate_hmm(n_sequences = 100, seq_length = 40, seed = 1)
print(r)
#> saqr_sim [hmm]  4000 x 3  (seed=1)
#>   params: trans, emission, init, n_states, n_symbols, hidden_paths 
#>   cols:   sequence_id, time, symbol 
r$params$trans
#>      [,1] [,2]
#> [1,]  0.8  0.2
#> [2,]  0.2  0.8
r$params$emission
#>      [,1] [,2] [,3]
#> [1,]  0.8  0.1  0.1
#> [2,]  0.1  0.8  0.1
head(as.data.frame(r))
#>   sequence_id time symbol
#> 1           1    1      2
#> 2           1    2      1
#> 3           1    3      2
#> 4           1    4      2
#> 5           1    5      3
#> 6           1    6      3

# Explicit three-state model
tr <- matrix(c(0.8, 0.1, 0.1,
               0.1, 0.8, 0.1,
               0.1, 0.1, 0.8), nrow = 3, byrow = TRUE)
em <- matrix(c(0.7, 0.2, 0.1,
               0.1, 0.7, 0.2,
               0.2, 0.1, 0.7), nrow = 3, byrow = TRUE)
r2 <- simulate_hmm(n_sequences = 80, seq_length = 50, n_states = 3,
                   n_symbols = 3, trans = tr, emission = em, seed = 42)
summary(r2)
#> Simulation type: hmm
#> Seed: 42
#> 
#> --- Data ---
#>   sequence_id         time          symbol     
#>  Min.   : 1.00   Min.   : 1.0   Min.   :1.000  
#>  1st Qu.:20.75   1st Qu.:13.0   1st Qu.:1.000  
#>  Median :40.50   Median :25.5   Median :2.000  
#>  Mean   :40.50   Mean   :25.5   Mean   :2.002  
#>  3rd Qu.:60.25   3rd Qu.:38.0   3rd Qu.:3.000  
#>  Max.   :80.00   Max.   :50.0   Max.   :3.000  
#> 
#> --- Parameters ---
#> List of 6
#>  $ trans       : num [1:3, 1:3] 0.8 0.1 0.1 0.1 0.8 0.1 0.1 0.1 0.8
#>  $ emission    : num [1:3, 1:3] 0.7 0.1 0.2 0.2 0.7 0.1 0.1 0.2 0.7
#>  $ init        : num [1:3] 0.333 0.333 0.333
#>  $ n_states    : int 3
#>  $ n_symbols   : int 3
#>  $ hidden_paths: int [1:80, 1:50] 1 3 3 1 1 3 3 3 2 3 ...
```
