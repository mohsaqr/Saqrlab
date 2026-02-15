# Saqrlab Project Learnings

### 2026-02-15
- [tna input formats]: `tna::tna()` accepts both a data frame of sequences (rows=sequences, cols=time points) and a square numeric matrix (transition frequencies). Column prefix can be "T" (as in `group_regulation`) not just "V".
- [group_regulation_long]: Structure is `Actor` (int), `Achiever` (chr), `Group` (num), `Course` (chr), `Time` (POSIXct), `Action` (chr). Has 27,533 rows, 9 unique actions (adapt, cohesion, consensus, coregulate, discuss, emotion, monitor, plan, synthesis).
- [group_regulation]: Wide format with T1-T26 columns (not V1-V26). 2000 rows. NAs are trailing (variable-length sequences).
- [tna model structure]: `tna` object is a list with `weights` (transition matrix), `inits` (initial probs), `labels` (state names), `data` (encoded sequences). `attr(model, "type")` is "relative" for probability weights.
- [transition counting]: To build a frequency matrix from sequences, count consecutive pairs (action[t], action[t+1]) within each sequence. For long data, split by ID and order by time first. For wide data, iterate across columns per row. Use `table()` with factors for vectorized counting.
- [frequencies.R]: File contains two functions: `frequencies()` builds a transition frequency matrix (state x state), `convert_sequence_format()` converts to frequency counts, onehot, edgelist, or follows format. Both accept wide and long input.
- [dplyr programmatic columns]: Use `dplyr::across(dplyr::all_of(col_names))` inside `group_by()`, `count()`, `distinct()` for programmatic column selection. Works with character vectors of column names. Avoids `!!sym()` / rlang dependency.
- [pivot_longer default]: `tidyr::pivot_longer()` creates a "name" column by default when `names_to` is not specified. This column is harmless if not used in downstream operations.
- [R CMD check]: Package has pre-existing warning about `\usage` entries in some man pages (not from new code). Hidden files (.claude, .github) and non-portable file paths in docs/ cache trigger NOTEs.
- [determinant]: `determinant()` is from `base`, not `stats`. Cannot use `@importFrom stats determinant` — just call it directly.
- [glasso]: `glasso::glasso()` returns list with `w` (estimated covariance), `wi` (estimated precision/inverse covariance). Warm starts use `start = "warm"` with `w.init` and `wi.init`. Cold start is default.
- [EBICglasso validation]: Our native EBIC + glasso implementation matches `qgraph::EBICglasso` exactly on test data (same partial correlations, same edges, same sparsity).
- [compositional frequency data]: Per-sequence action frequencies are compositional (rows sum to ~constant sequence length). This causes all partial correlations to be strongly negative. Expected behavior, not a bug.
- [glasso vs glassoFast]: Both have zero deps and use Fortran. `glassoFast` is ~2x faster but NOT available as a WebR (WASM) binary. `glasso` IS on the WebR repo. `glassoFast` also lacks `penalize.diagonal` param (always penalizes diagonal; workaround: pass rho as a matrix with zero diagonal). Chose `glasso` for WebR compatibility.
- [engagement stslist]: `tna::engagement` is an stslist object. Converting to frequency via `convert_sequence_format(as.data.frame(engagement), format = "frequency")` produces a `%` column (void/padding marker). Non-syntactic column names like `%` need auto-cleaning.
- [cograph::tplot]: Does NOT accept `main` argument. Use `title()` after the plot call to add titles. Pass `directed = FALSE` for undirected networks.
- [pcor_network]: Function renamed from `estimate_pcor_network()` to `pcor_network()`. Auto-cleans non-syntactic column names, all-NA columns, NA rows, and zero-variance columns with informative messages. Uses `glasso` package (WebR-compatible).
