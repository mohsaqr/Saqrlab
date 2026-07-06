# Get Learning State Verbs

Retrieve learning action verbs by category, with options for random
selection and combining multiple categories.

## Usage

``` r
get_learning_states(categories = "all", n = NULL, seed = NULL)
```

## Arguments

- categories:

  Character vector specifying which categories to include. Options:
  "metacognitive", "cognitive", "behavioral", "social", "motivational",
  "affective", or "all" for all categories. Default: "all".

- n:

  Integer or NULL. If specified, randomly sample n verbs from the
  selected categories. If NULL, return all verbs. Default: NULL.

- seed:

  Integer or NULL. Random seed for reproducible sampling. Default: NULL.

## Value

Character vector of learning state verbs.

## Examples

``` r
# All metacognitive verbs
get_learning_states("metacognitive")
#>  [1] "Plan"        "Monitor"     "Evaluate"    "Reflect"     "Regulate"   
#>  [6] "Adjust"      "Adapt"       "Check"       "Assess"      "Judge"      
#> [11] "Strategize"  "Prioritize"  "Set_goals"   "Track"       "Self_assess"
#> [16] "Calibrate"   "Diagnose"    "Forecast"    "Anticipate"  "Reconsider" 

# 10 random verbs from all categories
get_learning_states(n = 10, seed = 42)
#>  [1] "Reason"     "Edit"       "Satisfy"    "Rehearse"   "Worry"     
#>  [6] "Dedicate"   "Calendar"   "Initiate"   "Categorize" "Summarize" 

# 5 verbs from cognitive and behavioral
get_learning_states(c("cognitive", "behavioral"), n = 5)
#> [1] "Apply"      "Test"       "Generalize" "Write"      "Copy"      

# Self-regulated learning focus (metacognitive + motivational)
get_learning_states(c("metacognitive", "motivational"))
#>  [1] "Plan"        "Monitor"     "Evaluate"    "Reflect"     "Regulate"   
#>  [6] "Adjust"      "Adapt"       "Check"       "Assess"      "Judge"      
#> [11] "Strategize"  "Prioritize"  "Set_goals"   "Track"       "Self_assess"
#> [16] "Calibrate"   "Diagnose"    "Forecast"    "Anticipate"  "Reconsider" 
#> [21] "Focus"       "Persist"     "Explore"     "Create"      "Strive"     
#> [26] "Commit"      "Motivate"    "Endure"      "Overcome"    "Challenge"  
#> [31] "Aspire"      "Dedicate"    "Invest"      "Concentrate" "Attend"     
#> [36] "Sustain"     "Maintain"    "Initiate"    "Continue"    "Pursue"     
#> [41] "Drive"       "Hustle"      "Push"        "Achieve"     "Accomplish" 
#> [46] "Excel"       "Improve"     "Grow"        "Develop"     "Progress"   
```
