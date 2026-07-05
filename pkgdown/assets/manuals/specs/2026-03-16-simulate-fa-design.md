# Design Spec: `simulate_fa()` + List Return Pattern

**Date:** 2026-03-16
**File:** `R/simulate_latent.R`
**Tests:** `tests/testthat/test-simulate_latent.R`

---

## Problem

The existing `simulate_lpa()`, `simulate_lca()`, and `simulate_regression()` return a `data.frame` with ground-truth parameters stored in `attr(d, "true_params")`. This is opaque to non-R consumers — a JavaScript app that calls these scripts from the command line cannot access R attributes. The params are invisible.

Additionally, no function exists for factor analysis simulation with explicit loadings.

---

## Solution

### 1. Unified return structure: `list(data, params)`

All four simulation functions return a named list:

```r
list(
  data   = data.frame(...),   # the simulated observations
  params = list(...)          # ground-truth parameters, all plain R types
)
```

This is trivially serializable:

```bash
Rscript -e "
  source('R/simulate_latent.R')
  r <- simulate_fa(loadings, n = 500, seed = 1)
  jsonlite::write_json(r, 'output.json', auto_unbox = TRUE)
"
```

The JS app reads `output.data` for observations and `output.params` for ground truth — no R-specific knowledge required.

### 2. `simulate_fa(loadings, phi, psi, n, seed)`

New function in `R/simulate_latent.R`. Generates multivariate normal data from an explicit factor model for factor analysis parameter recovery testing.

---

## Function Specifications

### Updated: `simulate_lpa(means, sds, props, n, seed = NULL)`

**Returns:**
```r
list(
  data   = data.frame(y1, ..., yp, true_profile),
  params = list(means = means, sds = sds_mat, props = props)
)
```

### Updated: `simulate_lca(item_probs, class_probs, n, seed = NULL)`

**Returns:**
```r
list(
  data   = data.frame(item1, ..., itemm, true_class),
  params = list(item_probs = item_probs, class_probs = class_probs)
)
```

### Updated: `simulate_regression(coefs, predictor_sds, error_sd, n, seed = NULL)`

**Returns:**
```r
list(
  data   = data.frame(y, x1, ..., xk),
  params = list(coefs = coefs, predictor_sds = predictor_sds, error_sd = error_sd)
)
```

### New: `simulate_fa(loadings, phi = NULL, psi = NULL, n, seed = NULL)`

**Parameters:**

| Argument | Type | Default | Description |
|---|---|---|---|
| `loadings` | numeric matrix (p × m) | required | Factor loading matrix. Rows = observed variables, columns = factors. |
| `phi` | numeric matrix (m × m) | `diag(m)` | Factor correlation matrix. Must be symmetric positive-definite. |
| `psi` | numeric vector (length p) | auto-computed | Unique variances. Defaults to `1 - diag(loadings %*% phi %*% t(loadings))`. All must be > 0. |
| `n` | positive integer | required | Sample size. |
| `seed` | integer or NULL | `NULL` | Random seed. |

**Data generation model:**

```
Sigma = loadings %*% phi %*% t(loadings) + diag(psi)
X ~ MVN(0, Sigma)  via Cholesky decomposition
```

`.nearest_pd()` is defined in `simulate_data.R` (shared namespace). It normalises diagonals to 1 (correlation matrix output) so it is **only used for `phi`**, never for `Sigma`. For `Sigma`, use a raw eigenvalue clamp: `eig$vectors %*% diag(pmax(eig$values, 1e-6)) %*% t(eig$vectors)` — no diagonal normalisation.

**Returns:**
```r
list(
  data   = data.frame(y1, ..., yp),
  params = list(
    loadings      = loadings,   # p × m matrix
    phi           = phi,        # m × m matrix
    psi           = psi,        # length-p vector
    sigma_implied = Sigma       # p × p implied covariance matrix
  )
)
```

`sigma_implied` is included so the JS app can cross-check the covariance structure directly without recomputing it.

---

## Validation Rules

### `simulate_fa`

- `loadings` must be a numeric **matrix** (not a vector) with no NAs, ≥ 1 row, ≥ 1 column. Pass single-variable case as `matrix(0.8, nrow=1, ncol=1)`.
- `phi` symmetry checked with `isSymmetric(phi)`. Positive-definiteness verified via `tryCatch(chol(phi), error = function(e) stop("`phi` must be positive-definite: ", e$message))`. Off-diagonal range [−1, 1] is NOT checked separately (PD implies it for correlation matrices). `ncol(phi) == ncol(loadings)`.
- Computed `psi` (or user-supplied) must all be > 0; values ≤ 0 indicate communality ≥ 1 (over-factored model — error with informative message identifying which variables).
- `n` must be a positive integer. No lower bound enforced for `n` vs `p` — supported but documented in cookbook that `n < p` produces singular sample covariance (recovery tests should use `n > p`).
- `attr()` is not set on any return value. The `type` field is dropped entirely (redundant with the list structure — JS consumers identify type from the call site, not the output).

### Updated functions

No signature changes. Return type changes from `data.frame` to `list(data, params)`. All existing validation logic preserved.

---

## Test Scope

### Updated tests (45 existing)

All `attr(d, "true_params")` references change to `result$params`.
All `nrow(d)` / column access change to `nrow(result$data)` / `result$data$colname`.

### New tests for `simulate_fa` (~20)

| Test | What it checks |
|---|---|
| Returns list with `$data` and `$params` | Structure |
| `$data` has columns `y1`…`yp` | Columns |
| `nrow($data) == n` | Row count |
| `$params$loadings` identical to input | Ground truth preserved |
| `$params$phi` defaults to `diag(m)` | Default phi |
| `$params$psi` all > 0 for valid loadings | Positive uniquenesses |
| `$params$sigma_implied` matches `loadings %*% phi %*% t(loadings) + diag(psi)` | Sigma correctness |
| Implied covariance of generated data ≈ sigma_implied at large n | Recovery |
| Seed reproducibility (`identical()` on `$data`; `all.equal()` on `$params$sigma_implied`) | Determinism |
| Different seeds → different data | Non-determinism without seed |
| Error on non-matrix loadings | Validation |
| Error on phi dimension mismatch | Validation |
| Error on psi values ≤ 0 (user-supplied) | Validation |
| Error when auto-psi has values ≤ 0 (over-factored) | Validation |
| Single-factor model works | Edge case |
| Oblique model (phi ≠ I) generates correct structure | Oblique |

---

## Example Usage

```r
source("R/simulate_latent.R")

# Two orthogonal factors, 6 indicators
# nrow=6, ncol=2, byrow=FALSE (default): col1 = first 6 values, col2 = last 6
loadings <- matrix(c(0.8, 0.7, 0.6, 0,   0,   0,
                     0,   0,   0,   0.8, 0.7, 0.6),
                   nrow = 6, ncol = 2)
# Result: y1-y3 load on F1 (0.8, 0.7, 0.6), y4-y6 load on F2 (0.8, 0.7, 0.6)

result <- simulate_fa(loadings = loadings, n = 500, seed = 1)

result$data        # data.frame: y1..y6
result$params      # list: loadings, phi, psi, sigma_implied

# Recovery check: sample covariance ≈ sigma_implied at large n
result2 <- simulate_fa(loadings = loadings, n = 5000, seed = 1)
max(abs(cov(result2$data) - result2$params$sigma_implied))  # should be < 0.1

# Oblique model
phi <- matrix(c(1, 0.4, 0.4, 1), nrow = 2)
result3 <- simulate_fa(loadings = loadings, phi = phi, n = 500, seed = 1)
result3$params$phi   # factor correlation matrix

# JS-friendly export
jsonlite::write_json(result, "fa_data.json", auto_unbox = TRUE)
```
