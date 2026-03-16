# Session Handoff — 2026-03-16

## Completed
- Removed 20 computation files from Saqrlab R/ (moved to Nestimate in previous session)
- Fixed .onLoad in Saqrlab-package.R (removed .register_builtin_estimators() call)
- Rewrote simulate_data.R with complexity + batch generation features (169 tests pass)
- Added simulate_latent.R: simulate_lpa(), simulate_lca(), simulate_regression(), simulate_fa()
- Refactored all four to return list(data, params) — JSON-serializable, accessible to JS CLI apps
- 75 simulate_latent tests pass; 244 total simulation tests pass

## Current State
- simulate_data.R: 169 tests pass (complexity + batch)
- simulate_latent.R: 75 tests pass (LPA, LCA, regression, FA with ground-truth params)
- All functions return list(data, params) — no attr() anywhere
- docs/simulate_data-cookbook.md fully updated with all four functions + CLI export example
- NAMESPACE still exports functions now only in Nestimate (warnings on load_all, will error on check)
- Nestimate: standalone, not yet git repo, not yet wired as dependency

## Key Decisions
- list(data, params) return instead of data.frame + attr(): params visible to JSON/CLI/JS consumers
- simulate_fa uses raw eigenvalue clamp for Sigma (not .nearest_pd which normalises to correlation)
- .nearest_pd() used only for phi (correlation matrix); Sigma gets eigenvalue clamp only
- simulate_fa validates phi via isSymmetric() + tryCatch(chol(phi))
- psi auto-computed as 1 - communalities; error if any value <= 0 (over-factored model)
- sigma_implied stored in params for direct JS-side comparison with estimated covariance

## Open Issues
- NAMESPACE exports functions not in Saqrlab (moved to Nestimate) → warnings on load_all, will error on check
- Nestimate not yet a git repository
- Saqrlab DESCRIPTION not yet updated to depend on Nestimate

## Next Steps
1. Initialize git for Nestimate, commit everything
2. Add Nestimate to Saqrlab Imports / update DESCRIPTION
3. Update Saqrlab NAMESPACE to re-export Nestimate functions
4. Run devtools::check() to verify clean

## Context
- Nestimate: /Users/mohammedsaqr/Documents/Github/Nestimate/
- Saqrlab: /Users/mohammedsaqr/Documents/Github/Saqrlab/

---

# Previous Session Handoff — 2026-03-14

## Completed
- **Proximity timeline redesign**: Rewrote `.plot_proximity_timeline()` with smooth variable-width micro-segments, Okabe-Ito palette, direct endpoint labels, highlight parameter. Rewrote `.compute_proximity_mds()` with strength-based interleaved slot positioning. Removed `render_edges`, `smooth` params and `.build_edge_segments()`.
- **Package split**: Created Nestimate package at `/Users/mohammedsaqr/Documents/Github/Nestimate/` with 22 R computation files split from Saqrlab. 2137/2138 tests pass. Full project docs (CLAUDE.md, HANDOFF.md, LEARNINGS.md, CHANGES.md, FEATURES.md, ARCHITECTURE.md, STATUS.md).

## Current State
- **Saqrlab**: Still contains the original 22 computation files (not yet removed). 483 temporal_network tests pass. Uncommitted changes in R/temporal_network.R and tests/testthat/test-temporal_network.R.
- **Nestimate**: Standalone package, loads and tests pass, not yet a git repo. Contains the computation code that will eventually be removed from Saqrlab.

## Key Decisions
- Saqrlab retains: simulation, testing, verification, data resources (learning_states, global_names)
- Nestimate contains: estimation, bootstrap, permutation, HON, GIMME, MCML, mlVAR, temporal networks
- Saqrlab will depend on Nestimate (not yet wired up)

## Open Issues
- Nestimate not yet a git repository
- Saqrlab still has duplicate computation files (not yet removed)
- Saqrlab DESCRIPTION not yet updated to depend on Nestimate
- 1 marginal mlvar test failure in Nestimate (test helper data quality)
- Proximity timeline `.compute_proximity_mds()` uses strength-based interleaved slots — works but could be refined further

## Next Steps
1. Initialize git for Nestimate, commit everything
2. Remove 22 computation files from Saqrlab R/
3. Add Nestimate to Saqrlab Imports
4. Update Saqrlab NAMESPACE to re-export Nestimate functions
5. Commit Saqrlab changes

## Context
- Nestimate: `/Users/mohammedsaqr/Documents/Github/Nestimate/`
- Saqrlab: `/Users/mohammedsaqr/Documents/Github/Saqrlab/`
