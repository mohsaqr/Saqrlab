# Session Handoff — 2026-03-16

## Completed
- Removed 20 computation files from Saqrlab R/ (moved to Nestimate in previous session)
- Fixed .onLoad in Saqrlab-package.R (removed .register_builtin_estimators() call)
- Rewrote simulate_data.R with complexity + batch generation features
- 169 simulate_data tests pass

## Current State
- Saqrlab R/ now only contains simulation, testing, verification, data resource files
- simulate_data() fully enhanced: complexity parameter, n_batch parameter, type="batch" wildcard
- Package loads cleanly (after .onLoad fix)
- NAMESPACE still exports many functions that are now only in Nestimate (will cause warnings but not errors during devtools::load_all)
- Nestimate: standalone, not yet git repo, not yet wired as dependency

## Key Decisions
- `n_batch` placed after `...` in simulate_data() signature to prevent R partial matching of `n=` (common structural override)
- complexity = "clean" default preserves backward compatibility for all existing tests
- Auto-complexity uses Gumbel-max trick for weighted no-replacement sampling (vectorized, no for loop)

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
