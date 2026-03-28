# Session Handoff — 2026-03-28

## Completed This Session

### 1. saqr_sim S3 class (`R/saqr_sim.R`)
- Unified wrapper: every explicit-parameter simulation function returns `saqr_sim`
- S3 methods: print, summary, as.data.frame, [, dim, head, tail, str, names
- Backward-compatible: `$data` and `$params` still work, all existing tests pass unchanged
- 45 tests in `tests/testthat/test-saqr_sim.R`

### 2. Explicit-parameter statistical functions (`R/simulate_statistical.R`)
- `simulate_ttest(n_a, n_b, mean_a, mean_b, sd_a, sd_b)` — params include `cohens_d`
- `simulate_anova(n, means, sds)` — params include `eta_squared`; supports unequal groups
- `simulate_correlation(n, sigma, means)` — params include `is_correlation`
- `simulate_clusters(n, centers, sds, props)` — per-cluster sizes or proportions
- `simulate_prediction(n, coefs, cat_levels, cat_effects)` — params include `r_squared`
- 81 tests in `tests/testthat/test-simulate_statistical.R`

### 3. simulate_longitudinal() (`R/simulate_longitudinal.R`)
- VAR(1) panel data for mlVAR/ESM: explicit temporal (B), contemporaneous, between matrices
- ESM day/beep structure via `beeps_per_day` (day breaks reset carry-over)
- Complexity injection: na, outliers, heavy_tailed, heteroscedastic, tiny_n
- Auto-generated B guaranteed stationary via eigenvalue rescaling
- 42 tests in `tests/testthat/test-simulate_longitudinal.R`

### 4. Wrapped latent functions in saqr_sim (`R/simulate_latent.R`)
- simulate_lpa, simulate_lca, simulate_regression, simulate_fa, simulate_seq_clusters
- All return saqr_sim; 75 existing tests pass unchanged

### 5. Comprehensive manual (`docs/simulate_data-manual.md`)
- 944-line reference covering all 12 simulation functions
- Full signatures, return structure tables, working examples
- Parameter Recovery Cookbook with complete worked examples for every function type

### 6. Updated CLAUDE.md
- Rewrote to reflect post-split state (simulation-only package)
- Removed all references to network estimation code that moved to Nestimate

## Current State

### Tests
- **559 total assertions passing, 0 failures**
- New this session: 45 (saqr_sim) + 81 (statistical) + 42 (longitudinal) = 168 new tests
- All pre-existing tests unaffected

### Uncommitted Changes
All work from this session is **uncommitted**. `git status` shows:
- **Modified**: CHANGES.md, CLAUDE.md, FEATURES.md, HANDOFF.md, LEARNINGS.md, NAMESPACE, R/simulate_latent.R, man/Saqrlab-package.Rd, man/simulate_onehot_data.Rd
- **New (untracked)**: R/saqr_sim.R, R/simulate_longitudinal.R, R/simulate_statistical.R, man/saqr_sim.Rd, man/simulate_anova.Rd, man/simulate_clusters.Rd, man/simulate_correlation.Rd, man/simulate_longitudinal.Rd, man/simulate_prediction.Rd, man/simulate_ttest.Rd, tests/testthat/test-saqr_sim.R, tests/testthat/test-simulate_longitudinal.R, tests/testthat/test-simulate_statistical.R

### Architecture: two simulation tiers
1. **Explicit-parameter** (`simulate_ttest`, `simulate_anova`, etc.) → `saqr_sim` with `$params` containing ground truth. For parameter recovery testing.
2. **Random-parameter** (`simulate_data()`) → bare `data.frame` with attributes. Random structure from seed. For stress/robustness testing. NOT wrapped in saqr_sim (incompatible interface — tests use `d$group`, `names(d)` directly).

## Open Issues
1. Nestimate still needs `git init` + push to GitHub
2. Saqrlab DESCRIPTION needs `Imports: Nestimate` once Nestimate has a remote
3. `simulate_data()` stays as bare data.frame — wrapping breaks 169 existing tests that treat result as data.frame
4. Sequence functions (`simulate_sequences`, `simulate_tna_network`, etc.) not yet wrapped in saqr_sim
5. No complexity injection yet on latent functions (simulate_lpa, etc.) or statistical functions
6. docs/simulate_data-manual.md is local only — not yet pushed to GitHub

## Next Steps (prioritised)
1. Commit and push all session work to GitHub
2. Push `docs/simulate_data-manual.md` to GitHub for public access
3. Add complexity injection to explicit-parameter functions (simulate_ttest, simulate_lpa, etc.)
4. Add more longitudinal DGPs: growth curve, RI-CLPM, DSEM
5. Wire Nestimate as Saqrlab dependency
6. Consider wrapping sequence functions in saqr_sim

## Context
- Saqrlab: `/Users/mohammedsaqr/Documents/Github/Saqrlab/`
- Nestimate: `/Users/mohammedsaqr/Documents/Github/Nestimate/` (not on git)
- Manual: `docs/simulate_data-manual.md`
