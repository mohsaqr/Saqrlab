# Saqrlab: A Modern Laboratory for Data Simulation

Saqrlab generates synthetic datasets with known ground-truth parameters
so you can test and validate statistical methods. It spans classic
designs (t-test, ANOVA, regression, factor analysis), modern
data-generating processes (IRT, survival, multilevel, growth curves,
hidden Markov, missing-data mechanisms), and a full Temporal Network
Analysis toolkit.

## Details

Two complementary tiers:

- **Explicit-parameter** simulators (e.g.
  [`simulate_ttest()`](https://pak.dynasite.org/Saqrlab/reference/simulate_ttest.md),
  [`simulate_irt()`](https://pak.dynasite.org/Saqrlab/reference/simulate_irt.md),
  [`simulate_mlm()`](https://pak.dynasite.org/Saqrlab/reference/simulate_mlm.md),
  and the
  [`simulate()`](https://pak.dynasite.org/Saqrlab/reference/simulate.md)
  dispatcher) return a
  [saqr_sim](https://pak.dynasite.org/Saqrlab/reference/saqr_sim.md)
  object with `$data` and the true `$params`, for parameter-recovery
  testing.

- **Random-parameter** generation via
  [`simulate_data()`](https://pak.dynasite.org/Saqrlab/reference/simulate_data.md)
  returns a bare data frame whose structure is invented from the seed,
  for stress/robustness testing.

Start with
[`list_simulators()`](https://pak.dynasite.org/Saqrlab/reference/list_simulators.md)
for the catalogue,
[`simulate()`](https://pak.dynasite.org/Saqrlab/reference/simulate.md)
to generate, and
[`validate_recovery()`](https://pak.dynasite.org/Saqrlab/reference/validate_recovery.md)
to score whether a fitted method recovered the truth. Browse the
per-category HTML reference in the package's `reference/` directory.

## See also

Useful links:

- <https://mohsaqr.github.io/Saqrlab/>

- <https://github.com/mohsaqr/Saqrlab>

- Report bugs at <https://github.com/mohsaqr/Saqrlab/issues>

## Author

**Maintainer**: Mohammed Saqr <saqr@saqr.me>
