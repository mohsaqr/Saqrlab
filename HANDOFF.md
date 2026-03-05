# Session Handoff — 2026-03-05

## Completed
- **Rewrite `build_mcml()`** — Uses Saqrlab's own `build_network()` for between/within estimation
  - Between: recodes states to cluster labels → collapses consecutive same-cluster repetitions → `build_network()` on collapsed data → zero diagonal by construction
  - Within: filters non-member states to NA → `build_network()` per cluster
  - Fallback for pre-computed inputs: off-diagonal matrix block aggregation + row-normalize (diagonal always zero)
  - `cograph::cluster_summary()` only used for `plot_mcml()` compatibility
  - 6 data input formats, 5 cluster formats, enriched `$edges` field
- **`bootstrap_mcml()`** — Same recode/collapse/filter approach via `bootstrap_network()`
  - Between + within bootstrap using identical estimation logic as `build_mcml()`
  - S3 class `mcml_bootstrap` with print, summary, plot
- **Documentation** — Two reference docs in `docs/`

## Current State
- `devtools::load_all()` succeeds
- **101 tests pass** in test-mcml.R
- **56 tests pass** in test-bootstrap_mcml.R
- `$between` and `$within` are plain matrices (not tna objects)
- Between-cluster matrix has zero diagonal (only actual cluster switches counted)

### Files Modified
- `R/mcml.R` — rewritten with native estimation pipeline + `.collapse_consecutive()` helper
- `tests/testthat/test-mcml.R` — 101 tests (updated for new return structure + zero diagonal)
- `NAMESPACE` — exports unchanged

### Files Created
- `tests/testthat/test-bootstrap_mcml.R` — 56 tests
- `docs/build_mcml-for-cograph.md` — short cograph developer reference
- `docs/mcml-technical-reference.md` — full technical reference

## Key Decisions
- **Native estimation over cluster_summary**: Between/within built by recode/collapse/filter + build_network, not post-hoc aggregation.
- **Collapse consecutive**: `.collapse_consecutive()` replaces consecutive same-cluster values with NA before counting, so between diagonal is always zero. Only actual cluster switches are counted.
- **cluster_summary kept for plotting only**: `$cluster_summary` field still populated for `cograph::plot_mcml()` compatibility.
- **Plain matrices**: `$between` and `$within` are matrices, not tna objects. Simpler, consistent interface.

## Next Steps
1. Git commit when ready

## Context
- Package: Saqrlab (R package), branch: main
- cograph in Imports (for plotting + edge list conversion)
