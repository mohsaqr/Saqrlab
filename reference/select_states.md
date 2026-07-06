# Select States for Networks

Intelligently select learning states based on the number of nodes
needed, optionally biasing toward specific categories.

## Usage

``` r
select_states(
  n_states,
  primary_categories = NULL,
  secondary_categories = NULL,
  primary_ratio = 0.6,
  seed = NULL
)

smart_select_states(...)
```

## Arguments

- n_states:

  Integer. Number of states needed.

- primary_categories:

  Character vector. Categories to prioritize. Default: NULL (balanced
  selection).

- secondary_categories:

  Character vector. Categories to supplement with. Default: NULL (use
  all remaining).

- primary_ratio:

  Numeric (0 to 1). Proportion of states from primary categories.
  Default: 0.6.

- seed:

  Integer or NULL. Random seed. Default: NULL.

- ...:

  Arguments passed to `select_states`.

## Value

Character vector of selected learning states.

## Details

Selection logic:

- If n_states \<= 5: Single category or balanced small set

- If n_states 6-10: 1-2 categories prioritized

- If n_states 11-20: Multiple categories with primary focus

- If n_states \> 20: All categories combined

## Examples

``` r
# 5 states focused on self-regulation
select_states(5, primary_categories = "metacognitive")
#> [1] "Relief"   "Adapt"    "Excite"   "Reflect"  "Diagnose"

# 10 states: mostly cognitive, some behavioral
select_states(10,
  primary_categories = "cognitive",
  secondary_categories = "behavioral"
)
#>  [1] "Diagram"    "Conclude"   "List"       "Write"      "Encode"    
#>  [6] "Study"      "Generalize" "Attempt"    "Synthesize" "Retrieve"  

# 15 states: balanced across categories
select_states(15, seed = 42)
#>  [1] "Anxious"    "Develop"    "Improve"    "Analyze"    "Reread"    
#>  [6] "Link"       "Maintain"   "Forecast"   "Classify"   "Strive"    
#> [11] "Module"     "Brainstorm" "Celebrate"  "Create"     "Reconsider"
```
