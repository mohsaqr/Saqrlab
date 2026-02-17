# 06 — Testing Process, Dos, Don'ts, and Pitfalls

This guide covers the testing strategy, critical pitfalls, performance tips, and rules for replicating the Saqrlab network analysis system.

---

## 1. Testing Strategy

### 1.1 Test Files

| Test file | Lines | Covers |
|---|---|---|
| `tests/testthat/test-build_network.R` | ~700 | `build_network()`, all 3 association methods, multilevel, predictability, S3 methods, `$data` field |
| `tests/testthat/test-estimate_network.R` | ~500 | Deprecated `estimate_network()`, all 6 methods, scaling, threshold, aliases, multilevel, edge extraction |
| `tests/testthat/test-bootstrap_network.R` | ~280 | `bootstrap_network()`, all 6 methods, inference types, CIs, summary, pruning, reproducibility |
| `tests/testthat/test-permutation_test.R` | ~275 | `permutation_test()`, all 3 transition + 2 association methods, paired mode, p.adjust, reproducibility |
| `tests/testthat/test-estimator_registry.R` | varies | Registry CRUD operations |
| `tests/testthat/test-frequencies.R` | varies | `frequencies()` and `convert_sequence_format()` |

### 1.2 Test Data Helpers

Every test file defines helper functions at the top to generate reproducible data:

**Wide sequence data** (for transition methods):
```r
.make_boot_wide <- function(n = 50, t = 10, states = c("A", "B", "C"),
                            seed = 42) {
  set.seed(seed)
  mat <- matrix(sample(states, n * t, replace = TRUE), nrow = n, ncol = t)
  colnames(mat) <- paste0("T", seq_len(t))
  as.data.frame(mat, stringsAsFactors = FALSE)
}
```

**Frequency data** (for association methods):
```r
.make_freq_data <- function(n = 100, p = 5, seed = 42) {
  set.seed(seed)
  mat <- matrix(rpois(n * p, lambda = 10), nrow = n, ncol = p)
  colnames(mat) <- paste0("state_", seq_len(p))
  df <- as.data.frame(mat)
  df$rid <- seq_len(n)  # some helpers include rid, some don't
  df
}
```

**Multilevel data**:
```r
.make_multilevel_data <- function(n_persons = 30, obs_per_person = 5,
                                  p = 5, seed = 42) {
  set.seed(seed)
  n_total <- n_persons * obs_per_person
  mat <- matrix(rpois(n_total * p, lambda = 10), nrow = n_total, ncol = p)
  colnames(mat) <- paste0("state_", seq_len(p))
  df <- as.data.frame(mat)
  df$person <- rep(seq_len(n_persons), each = obs_per_person)
  df$rid <- seq_len(n_total)
  df
}
```

**Key**: Use `rpois(lambda = 10)` for frequency data — Poisson counts are realistic for action frequencies and ensure all-positive, integer values. Use `sample(states, ..., replace = TRUE)` for sequence data.

### 1.3 Test Organization

Tests follow a consistent pattern:

1. **Input validation**: error messages for bad inputs
2. **Basic functionality**: each method works and returns correct structure
3. **Method aliases**: aliases resolve to canonical names
4. **Auto-cleaning**: messages for dropped columns/rows
5. **Method-specific fields**: extras present/absent as expected
6. **Cross-method consistency**: all methods produce structurally identical output
7. **S3 methods**: print output contains expected strings
8. **Composability**: params stored and replayable
9. **Reproducibility**: same seed → identical results

### 1.4 Small Iterations for Speed

All bootstrap and permutation tests use small iteration counts (20-50) to keep test runtime fast:

```r
boot <- bootstrap_network(wide, method = "relative", iter = 30L, seed = 1)
perm <- permutation_test(net1, net2, iter = 50L, seed = 42)
```

And small `nlambda` for glasso tests:

```r
net <- build_network(df, method = "glasso", params = list(nlambda = 20L))
```

### 1.5 `suppressWarnings()` for Deprecated Functions

testthat edition 3 treats warnings as failures. All calls to `estimate_network()` (which is deprecated) must be wrapped:

```r
net <- suppressWarnings(estimate_network(wide, method = "relative"))
```

The deprecation warning test itself explicitly expects the warning:

```r
test_that("estimate_network emits deprecation warning", {
  wide <- .make_wide_seq()
  expect_warning(estimate_network(wide, method = "relative"), "deprecated")
})
```

### 1.6 Integration Tests

Some tests use real data from the `tna` package:

```r
test_that("estimate_network works with tna::group_regulation (wide)", {
  skip_if_not_installed("tna")
  net <- suppressWarnings(
    estimate_network(tna::group_regulation, method = "relative")
  )
  expect_equal(net$n_nodes, 9)
  row_sums <- rowSums(net$matrix)
  expect_true(all(abs(row_sums - 1) < 1e-10))
})
```

**Key**: `skip_if_not_installed("tna")` makes these tests optional.

---

## 2. Critical Pitfalls

### 2.1 `byrow = TRUE` with tabulate-based encoding

**The bug**: Using `pair_idx = (from_int - 1) * n_states + to_int` with `tabulate()` produces a flat vector in row-major order. R's `matrix()` fills column-major by default.

**The fix**: Always use `byrow = TRUE`:
```r
matrix(counts, nrow = n_states, ncol = n_states, byrow = TRUE, ...)
```

**Without the fix**: From and to get silently swapped. The network matrix looks plausible (same shape, same non-zero pattern) but all directed edges point the wrong way. No error is thrown.

### 2.2 NA Handling: Pair First, Filter After

**The bug**: Stripping NAs from a sequence before creating consecutive pairs bridges gaps.

Example: Sequence `[A, NA, B, C]`
- **Wrong**: Strip NA → `[A, B, C]` → pairs: `A→B, B→C` (false `A→B` transition)
- **Correct**: Pairs first → `(A,NA), (NA,B), (B,C)` → filter NAs → `B→C` only

**The fix**:
```r
from_vec <- actions[-n]
to_vec   <- actions[-1L]
valid    <- !is.na(from_vec) & !is.na(to_vec)
from_vec <- from_vec[valid]
to_vec   <- to_vec[valid]
```

### 2.3 `sample.int(n)` vs `sample.int(n, k)` — RNG Divergence

**The bug**: `sample.int(n)` (full permutation) consumes `n` RNG values. `sample.int(n, k)` consumes only `k`. If one implementation uses `sample.int(n_total)` and another uses `sample.int(n_total, n_x)`, they consume different amounts of RNG state per iteration. After iteration 1, all subsequent permutations diverge.

**The fix**: Use `sample.int(n_total, n_x)` to match NCT's RNG consumption:
```r
idx_x <- sample.int(n_total, n_x)  # draws n_x values, consuming n_x RNG
```

**Impact**: Without this fix, p-value agreement with NCT drops from 100% to ~1.4%.

### 2.4 `determinant()` is from `base`, not `stats`

**The bug**: Writing `@importFrom stats determinant` causes a NAMESPACE error because `determinant()` is in the `base` package, not `stats`.

**The fix**: Just call `determinant()` directly. It's always available.

### 2.5 `glassoFast` Not Available on WebR

**The bug**: `glassoFast` is not compiled for WASM/WebR.

**The fix**: Use `glasso` package instead. It's on the WebR repo. Note that `glassoFast` also lacks the `penalize.diagonal` parameter.

### 2.6 Self-Pairs in Co-occurrence Are Double-Counted

**The bug**: The bidirectional approach counts `(A[i], A[j])` once as forward and once as reverse. For self-pairs (same state in both positions), this means diagonal entries are doubled.

**The fix**: `diag(cooc) <- diag(cooc) / 2` after accumulation.

### 2.7 `$data` Stores Cleaned Data, Not Raw Input

**The bug**: Downstream code (permutation test, bootstrap) assumes `$data` is ready to use. If it stored raw data, every iteration would need to re-clean (drop NAs, non-numeric columns, zero-variance columns), causing both performance issues and potential inconsistencies.

**The fix**: Estimators return `cleaned_data` which `build_network()` stores as `$data`. For data frame input to association methods, this is the cleaned numeric matrix. For matrix input, it's NULL. For transition methods, it's the original data frame (no cleaning needed).

### 2.8 Fixed Lambda Across Permutations is Wrong

**The bug**: Using the lambda selected for the original network across all permutation iterations drops significance agreement from 99% to 93%. Each permutation sample has different correlation structure → different optimal lambda.

**The fix**: Re-select lambda via EBIC on each permutation iteration. The `glassopath()` approach makes this fast — single Fortran call per iteration, then loop over the 3D array to find best EBIC.

### 2.9 `p.adjust()` Applied to Flat Vector, Then Reshaped

**The bug**: Applying `p.adjust()` to a matrix directly doesn't work correctly — it needs to see ALL p-values at once.

**The fix**: Flatten → adjust → reshape:
```r
p_values_flat <- p.adjust(p_values_flat, method = adjust)
p_mat <- matrix(p_values_flat, n_nodes, n_nodes, dimnames = list(nodes, nodes))
```

### 2.10 data.table NSE Symbols Need `globalVariables()`

**The bug**: data.table NSE symbols like `.SD`, `:=`, `.grp_key` trigger R CMD check NOTEs about "no visible binding for global variable".

**The fix**: Declare them in `utils::globalVariables()` in `R/Saqrlab-package.R`.

---

## 3. Performance Tips

### 3.1 `colSums()` on 2D Matrix vs `apply()` on 3D Array

The bootstrap fast path uses `colSums(trans_2d[idx, ])` — a single C-level call. The alternative (used by `tna::bootstrap()`) stores per-sequence counts in a 3D array and uses `apply(..., c(1,2), sum)` — R-level looping.

**Benchmark**: 2.8x faster for 2000 sequences, 500 iterations.

### 3.2 Running Accumulators vs Full Storage

Permutation test uses running counters:
```r
exceed_counts <- exceed_counts + (abs(perm_diff) >= abs(obs_flat))
sum_diffs     <- sum_diffs + perm_diff
sum_diffs_sq  <- sum_diffs_sq + perm_diff^2
```

This requires O(n²) memory instead of O(iter × n²). For 1000 iterations with 9 nodes: 81 numbers vs 81,000.

### 3.3 `glassopath()` vs Separate `glasso()` Calls

`glassopath()` solves the entire lambda path in one Fortran call with internal warm-starting. Separate `glasso()` calls (even with manual warm-starting in R) are 14.5x slower because of R↔Fortran overhead per call.

### 3.4 Pre-Clean Data Once

Association permutation test pools `$data` directly (already cleaned). Does NOT re-run `.prepare_association_input()` per iteration. The fast-path closures call `cor()` directly on the numeric matrix.

### 3.5 `vapply()` for Type-Safe Iteration

Bootstrap uses `vapply()` which returns a matrix directly (each iteration's vector becomes a column). No need for `lapply()` + `do.call(rbind, ...)`.

```r
boot_flat <- vapply(seq_len(iter), function(i) {
  ...
  as.vector(mat)
}, numeric(nbins))
# boot_flat is nbins × iter matrix
t(boot_flat)  # → iter × nbins
```

### 3.6 `data.table::data.table()` for Summary Construction

Summary data frames are built with `data.table::data.table()` for speed, then converted to `data.frame` at the end. The data.table filtering syntax (`dt[condition]`) is also used for row filtering.

---

## 4. Dos

1. **Always validate estimator output** has `matrix`, `nodes`, `directed`
2. **Always use `byrow = TRUE`** with tabulate-based flat indexing
3. **Always filter NAs after pairing**, not before
4. **Always store `cleaned_data`** in estimator return value
5. **Always use `set.seed()`** before the permutation/bootstrap loop, not inside it
6. **Always use `on.exit()`** to restore `par()` in plot methods
7. **Always test reproducibility** with the `seed` parameter
8. **Always use `sample.int(n_total, n_x)`** (partial draw) for permutation group assignment
9. **Always re-select lambda via EBIC** per permutation iteration for glasso
10. **Always use `as.vector(t(matrix))`** for row-major flattening in summaries
11. **Always use `stopifnot()`** for input validation at function entry
12. **Always use `tryCatch()`** around `solve()`, `glasso()`, and `glassopath()` — they can fail
13. **Always declare data.table NSE symbols** in `utils::globalVariables()`
14. **Always use `suppressWarnings()`** when testing deprecated functions in testthat edition 3
15. **Always use `skip_if_not_installed()`** for optional package tests

---

## 5. Don'ts

1. **Don't use `for` loops for counting** — use vectorized `tabulate()`
2. **Don't store full 3D arrays** when running accumulators suffice
3. **Don't re-clean data** that's already in `$data`
4. **Don't use `sample.int(n_total)`** when you only need `sample.int(n_total, k)`
5. **Don't use fixed lambda** across permutations for glasso
6. **Don't call `@importFrom stats determinant`** — it's from `base`
7. **Don't suppress NAs before transition pair extraction** — pair first, filter after
8. **Don't assume matrix input has `cleaned_data`** — it's NULL
9. **Don't add `dimnames` to bootstrap iteration matrices** — they're flattened immediately
10. **Don't use `glassoFast`** — not available on WebR
11. **Don't use `table()` for counting** — `tabulate()` is faster and doesn't create factors
12. **Don't forget `byrow = TRUE`** when constructing matrices from flat indices
13. **Don't forget `t()` before `as.vector()`** when building summary data frames
14. **Don't use `sapply()`** — use `vapply()` for type safety
15. **Don't use `matrix()` default** (column-major) with row-major flat indices

---

## 6. Verification Checklist

After implementing the system, verify:

### Correctness

- [ ] `build_network()` with all 6 methods produces correct output
- [ ] Method aliases resolve correctly (all 10 aliases)
- [ ] Transition matrices: rows sum to 1 for `"relative"`, integers for `"frequency"`
- [ ] Co-occurrence matrix is symmetric with halved diagonal
- [ ] Association matrices: zero diagonal, symmetric
- [ ] Multilevel: between aggregates to person means, within centers per-person
- [ ] `level = "both"` returns `netobject_ml` with both sub-networks
- [ ] Edge extraction: directed uses all non-diagonal non-zero; undirected uses upper triangle
- [ ] Scaling applied correctly (minmax, max, rank, normalize)
- [ ] Threshold zeroes entries with `abs(x) < threshold`

### Bootstrap

- [ ] Transition fast path matches full-estimator path (same results)
- [ ] Stability inference p-values in [0, 1]
- [ ] `ci_lower <= ci_upper` everywhere
- [ ] Summary has correct columns and filtering
- [ ] Pruned model has `<= edges` of original
- [ ] Same seed → identical results

### Permutation Test

- [ ] P-values in [0, 1]
- [ ] `diff = x$matrix - y$matrix`
- [ ] Same networks → diff = 0, effect_size = 0
- [ ] Paired mode works with equal-sized groups
- [ ] `p.adjust()` produces values `>=` raw p-values
- [ ] Same seed → identical results
- [ ] Transition path matches tna::permutation_test (r=1.0 for relative, co_occurrence)
- [ ] Association path matches NCT (r≥0.999 for cor, glasso)

### Registry

- [ ] 6 built-in estimators registered on package load
- [ ] Custom estimators can be registered, used, and removed
- [ ] `get_estimator()` on missing name gives helpful error
- [ ] `list_estimators()` returns sorted data frame

### Predictability

- [ ] glasso/pcor: R² = 1 - 1/Omega_jj, values in [0, 1]
- [ ] cor: multiple R² from correlation submatrix, values in [0, 1]
- [ ] Isolated nodes (no edges) → R² = 0
- [ ] Multilevel: returns list with `$between` and `$within`

---

## 7. Running Tests

```r
# Run all tests
devtools::test()

# Run specific test file
testthat::test_file("tests/testthat/test-build_network.R")

# Run with verbose output
devtools::test(reporter = "summary")
```

All tests should pass with zero failures and zero errors. Warnings may appear from deprecated function tests (expected).
