# Group Regulation Actions

The 9 standard group regulation (SSRL) actions used in collaborative
learning research.

## Usage

``` r
GROUP_REGULATION_ACTIONS
```

## Format

A character vector of 9 actions.

## Details

These actions represent Socially Shared Regulation of Learning (SSRL):

- adapt:

  Adjusting strategies based on feedback

- cohesion:

  Building and maintaining group unity

- consensus:

  Reaching group agreement

- coregulate:

  Mutually regulating each other's learning

- discuss:

  Engaging in group discussion

- emotion:

  Expressing or regulating emotions

- monitor:

  Tracking and evaluating progress

- plan:

  Planning group activities and strategies

- synthesis:

  Integrating and combining ideas

## Examples

``` r
GROUP_REGULATION_ACTIONS
#> [1] "adapt"      "cohesion"   "consensus"  "coregulate" "discuss"   
#> [6] "emotion"    "monitor"    "plan"       "synthesis" 

# Use in simulation
data <- simulate_long_data(
  n_groups = 10,
  states = GROUP_REGULATION_ACTIONS,
  seed = 42
)
```
