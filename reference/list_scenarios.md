# List Available Scenario Presets

Return a tidy catalogue of the built-in scenario presets. Each scenario
is a hard-coded recipe of argument-lists targeting one `simulate_*`
function. Print the result directly.

## Usage

``` r
list_scenarios()
```

## Value

A base `data.frame` with one row per scenario and columns `scenario`,
`description`, `simulator` (the target `simulate_*` function) and
`n_cases` (number of argument-lists in the recipe).

## Examples

``` r
list_scenarios()
#>             scenario                                           description
#> 1        power_ttest Two-group t-test power grid: Cohen's d x sample size.
#> 2      mlm_icc_sweep Multilevel model sweeping the intraclass correlation.
#> 3         irt_models     IRT response data across 1PL, 2PL and 3PL models.
#> 4 survival_censoring            Survival data sweeping the censoring rate.
#>           simulator n_cases
#> 1    simulate_ttest       6
#> 2      simulate_mlm       4
#> 3      simulate_irt       3
#> 4 simulate_survival       3
```
