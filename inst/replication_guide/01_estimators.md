# 01 — Building All 6 Estimators

All estimator implementations live in `R/estimators.R` (787 lines). This guide walks through every estimator and shared helper line by line.

---

## 1. Shared Helper: `.select_state_cols()` (lines 15-23)

Resolves which columns in a data frame contain state data.

```r
.select_state_cols <- function(data, id = NULL, cols = NULL) {
  if (!is.null(cols)) {
    cols                              # explicit columns win
  } else if (!is.null(id)) {
    setdiff(names(data), id)          # exclude ID columns
  } else {
    names(data)                       # all columns
  }
}
```

**Priority**: explicit `cols` > exclude `id` > all columns.

Used by: `.count_transitions_wide()`, `.count_cooccurrence_wide()`, `.precompute_per_sequence_wide()`.

---

## 2. Core Counting Engine: `.count_transitions()` (lines 42-59)

Dispatcher function. Determines format (auto-detect via presence of `action` column), then calls the appropriate implementation.

```r
.count_transitions <- function(data, format = "auto", action = "Action",
                               id = NULL, time = "Time", cols = NULL) {
  stopifnot(is.data.frame(data))

  if (format == "auto") {
    format <- if (action %in% names(data)) "long" else "wide"
  }

  if (format == "wide") {
    .count_transitions_wide(data, id = id, cols = cols)
  } else {
    .count_transitions_long(data, action = action, id = id, time = time)
  }
}
```

### 2.1 Wide Path: `.count_transitions_wide()` (lines 68-112)

**Algorithm**:

1. Resolve state columns via `.select_state_cols()`
2. Convert to matrix: `mat <- as.matrix(data[, state_cols, drop = FALSE])`
3. Create consecutive pairs via matrix slicing:
   - `from_vec <- as.vector(mat[, -nc, drop = FALSE])` — all columns except last
   - `to_vec <- as.vector(mat[, -1L, drop = FALSE])` — all columns except first
4. **Filter NAs AFTER pairing**: `valid <- !is.na(from_vec) & !is.na(to_vec)` — this is CRITICAL
5. Integer-encode: `from_int <- match(from_vec, all_states)`, same for `to_int`
6. Flat index: `pair_idx <- (from_int - 1L) * n_states + to_int`
7. Count: `counts <- tabulate(pair_idx, nbins = n_states * n_states)`
8. Reshape: `matrix(as.integer(counts), nrow = n_states, ncol = n_states, byrow = TRUE, ...)`

**CRITICAL — `byrow = TRUE`**: The flat index `(from - 1) * n + to` produces row-major ordering. R's `matrix()` fills column-major by default. Without `byrow = TRUE`, from and to get silently swapped. The network will look plausible but will be wrong.

**CRITICAL — NA handling**: NAs are filtered AFTER creating pairs. If you strip NAs before pairing (e.g., removing NA rows first), you create false transitions that bridge gaps. Example: sequence `[A, NA, B]` should produce zero transitions, but stripping the NA first creates the false transition `A → B`.

**Return**: Square integer matrix with sorted unique state names as row/column names.

### 2.2 Long Path: `.count_transitions_long()` (lines 121-218)

Uses `data.table` for fast grouped operations.

**Algorithm**:

1. Convert to data.table: `dt <- data.table::as.data.table(data)`
2. Order by ID + time: `data.table::setorderv(dt, order_cols)`
3. Build group key:
   - Single ID column: use directly
   - Multiple ID columns: composite key via `do.call(paste, c(.SD, sep = "\x1f"))` — uses unit separator to avoid collisions
4. Extract pairs per group using data.table `[.data.table]`:
   ```r
   pairs <- dt[, {
     a <- get(action_col)
     n <- length(a)
     if (n < 2L) list(from = character(0), to = character(0))
     else {
       f <- a[-n]; t <- a[-1L]
       ok <- !is.na(f) & !is.na(t)
       list(from = f[ok], to = t[ok])
     }
   }, by = grp_col]
   ```
5. Same tabulate + matrix pattern as wide path

**Same CRITICAL rules apply**: NAs filtered after pairing, `byrow = TRUE` in `matrix()`.

---

## 3. Transition Estimators

### 3.1 `.estimator_frequency()` (lines 226-244)

The simplest estimator. Calls `.count_transitions()` and returns the raw integer matrix.

```r
.estimator_frequency <- function(data, format = "auto", action = "Action",
                                  id = NULL, time = "Time", cols = NULL, ...) {
  freq_mat <- .count_transitions(data, format, action, id, time, cols)
  states <- rownames(freq_mat)
  list(
    matrix = freq_mat,
    nodes = states,
    directed = TRUE,
    cleaned_data = data,
    frequency_matrix = freq_mat
  )
}
```

**Key**: Returns `cleaned_data = data` (the original data frame, not a copy). For transition methods, no cleaning is needed — the data is used as-is.

### 3.2 `.estimator_relative()` (lines 249-275)

Calls `.count_transitions()`, then row-normalizes to get transition probabilities.

```r
freq_mat <- .count_transitions(...)
row_sums <- rowSums(freq_mat)
rel_mat <- freq_mat
storage.mode(rel_mat) <- "double"        # convert from integer to double
nonzero <- row_sums > 0
rel_mat[nonzero, ] <- rel_mat[nonzero, ] / row_sums[nonzero]
```

**Key detail**: `storage.mode(rel_mat) <- "double"` converts the integer matrix to double BEFORE division. Without this, integer division would truncate.

**Key detail**: Only normalizes rows with non-zero sums. Rows with all zeros (states that never appear as a source) remain all zeros.

Returns both `matrix` (relative) and `frequency_matrix` (raw counts).

### 3.3 `.estimator_co_occurrence()` (lines 285-312)

Positional co-occurrence: counts all pairs of column positions `(i, j)` where `i < j`. Unlike transitions, this is NOT about consecutive pairs — it's about all pairwise column combinations.

Dispatches to `.count_cooccurrence_wide()` or `.count_cooccurrence_long()`.

Returns `directed = FALSE` — co-occurrence networks are always undirected.

#### Wide co-occurrence: `.count_cooccurrence_wide()` (lines 325-374)

**Algorithm**:

1. Integer-encode entire matrix once: `int_mat <- matrix(match(mat, all_states), ...)`
2. Nested loop over all column pairs `(i, j)` where `i < j`:
   ```r
   for (i in seq_len(nc - 1L)) {
     col_i <- int_mat[, i]
     for (j in seq(i + 1L, nc)) {
       col_j <- int_mat[, j]
       valid <- !is.na(col_i) & !is.na(col_j)
       fi <- col_i[valid]
       tj <- col_j[valid]
       # Bidirectional: count both (fi→tj) and (tj→fi)
       idx_fwd <- (fi - 1L) * n_states + tj
       idx_rev <- (tj - 1L) * n_states + fi
       counts <- counts + tabulate(idx_fwd, nbins) + tabulate(idx_rev, nbins)
     }
   }
   ```
3. **Self-pair halving**: `diag(cooc) <- diag(cooc) / 2` — the bidirectional approach double-counts diagonal entries (where both columns have the same state). Must halve them.
4. Matrix uses `byrow = TRUE` and `as.numeric(counts)` (not integer).

#### Long co-occurrence: `.count_cooccurrence_long()` (lines 382-457)

Uses `utils::combn(n, 2)` per group to generate all position pairs, then counts both directions:

```r
pairs <- dt[!is.na(get(action_col)), {
  a <- get(action_col)
  n <- length(a)
  if (n < 2L) list(from = character(0), to = character(0))
  else {
    cp <- utils::combn(n, 2)
    f <- a[cp[1, ]]
    t <- a[cp[2, ]]
    list(from = c(f, t), to = c(t, f))  # both directions
  }
}, by = grp_col]
```

Same self-pair halving at the end.

---

## 4. Shared Association Helper: `.prepare_association_input()` (lines 469-560)

Universal input cleaner for all three association methods. Handles both data frame and matrix input.

### Data Frame Path (lines 472-527)

**Step-by-step cleaning**:

1. **Exclude columns**: remove `id_col` and `"rid"` (standard ID column name)
2. **Keep numeric only**: `vapply(data, is.numeric, logical(1))`
3. **Drop non-syntactic names**: `make.names(keep) == keep` — columns with names like `%` or `*` are dropped with a message
4. **Require ≥ 2 columns** after cleaning
5. **Convert to matrix**: `mat <- as.matrix(data[, keep, drop = FALSE])`
6. **Drop all-NA columns**: `apply(mat, 2, function(x) all(is.na(x)))`
7. **Drop NA rows**: `complete.cases(mat)` — rows with ANY NA are removed
8. **Require ≥ 3 complete rows**
9. **Drop zero-variance columns**: `apply(mat, 2, stats::var)` — columns with `var == 0`
10. **Require ≥ 2 columns** after all cleaning
11. **Compute correlation**: `S <- cor(mat, method = cor_method)`

**Returns**: `list(S = S, n = nrow(mat), mat = mat)` — `mat` is the cleaned numeric matrix, stored as `cleaned_data` in the estimator return.

### Matrix Path (lines 529-554)

For pre-computed correlation or covariance matrices:

1. **Validate**: square, symmetric (`isSymmetric(unname(data), tol = 1e-8)`)
2. **Require `n`**: sample size must be provided when input is a matrix
3. **Auto-detect type**: if all diagonal values ≈ 1, treat as correlation; otherwise treat as covariance
4. **Covariance → correlation**: `d <- sqrt(diag(data)); S <- data / outer(d, d)`
5. **Auto-name**: if no column names, generates `V1`, `V2`, ...
6. **mat = NULL**: no row-level data available from matrix input

**Returns**: `list(S = S, n = n, mat = NULL)`

---

## 5. Association Estimators

### 5.1 `.estimator_cor()` (lines 565-593)

The simplest association estimator.

```r
prepared <- .prepare_association_input(data, id_col, n, cor_method, input_type)
S <- prepared$S
net <- S
diag(net) <- 0                           # zero out diagonal (self-correlation = 1)
net[abs(net) < threshold] <- 0           # apply threshold

list(
  matrix = net, nodes = colnames(net), directed = FALSE,
  cleaned_data = prepared$mat,
  cor_matrix = S, n = prepared$n, p = ncol(S)
)
```

**Key**: The threshold is applied inside the estimator (for method-level thresholding). `build_network()` may apply an additional threshold afterward.

### 5.2 `.estimator_pcor()` (lines 598-638)

Unregularized partial correlations via matrix inversion.

```r
prepared <- .prepare_association_input(...)
S <- prepared$S

Wi <- tryCatch(
  solve(S),                              # invert correlation matrix
  error = function(e) {
    stop("Correlation matrix is singular (p >= n or collinear variables). ",
         "Use method = 'glasso' for regularised estimation.", call. = FALSE)
  }
)
colnames(Wi) <- rownames(Wi) <- colnames(S)

pcor <- .precision_to_pcor(Wi, threshold)
```

**Key**: `solve()` will fail when `p >= n` (more variables than observations) or with collinear variables. The error message suggests using `glasso` instead.

### 5.3 `.estimator_glasso()` (lines 739-786)

EBIC-selected graphical lasso regularization.

**Algorithm**:

1. **Prepare input**: `.prepare_association_input()`
2. **Validate params**: gamma ≥ 0, nlambda ≥ 2, lambda.min.ratio in (0, 1)
3. **Compute lambda path**: `.compute_lambda_path(S, nlambda, lambda.min.ratio)`
4. **Select best lambda**: `.select_ebic(S, lambda_path, n, gamma, penalize.diagonal)`
5. **Convert to partial correlations**: `.precision_to_pcor(selected$wi, threshold)`

Returns all intermediate results: `precision_matrix`, `cor_matrix`, `lambda_selected`, `ebic_selected`, `lambda_path`, `ebic_path`, `gamma`, `n`, `p`.

---

## 6. Shared Association Helpers

### 6.1 `.compute_lambda_path()` (lines 645-652)

Generates a log-spaced sequence of regularization parameters.

```r
.compute_lambda_path <- function(S, nlambda, lambda.min.ratio) {
  lambda_max <- max(abs(S[upper.tri(S)]))   # max off-diagonal correlation
  if (lambda_max <= 0) stop("All off-diagonal correlations are zero")
  lambda_min <- lambda_max * lambda.min.ratio
  exp(seq(log(lambda_max), log(lambda_min), length.out = nlambda))
}
```

**Key**: Uses `upper.tri(S)` to get off-diagonal elements only. The path goes from `lambda_max` (most regularized, sparsest) down to `lambda_min` (least regularized, densest).

### 6.2 `.select_ebic()` (lines 657-723)

Selects the best lambda via EBIC using warm-started glasso fits.

**Per-lambda iteration**:

```r
for (i in seq_along(lambda_path)) {
  fit <- tryCatch(
    glasso::glasso(
      s = S, rho = lam,
      penalize.diagonal = penalize_diagonal,
      start = if (is.null(w_prev)) "cold" else "warm",
      w.init = w_prev, wi.init = wi_prev,
      trace = FALSE
    ),
    error = function(e) NULL
  )
  # ... warm start update ...
  w_prev <- fit$w
  wi_prev <- fit$wi
```

**EBIC formula** (lines 699-702):

```r
log_det <- determinant(fit$wi, logarithm = TRUE)
loglik <- (n / 2) * (log_det_val - sum(diag(S %*% fit$wi)))
npar <- sum(abs(fit$wi[upper.tri(fit$wi)]) > 1e-10)
ebic <- -2 * loglik + npar * log(n) + 4 * npar * gamma * log(p)
```

**CRITICAL**: `determinant()` is from `base`, NOT `stats`. Do NOT write `@importFrom stats determinant` — it will cause a NAMESPACE error.

**Key**: `npar` counts the number of non-zero upper-triangular entries in the precision matrix (using threshold `1e-10`). Higher `gamma` penalizes complexity more, producing sparser networks.

### 6.3 `.precision_to_pcor()` (lines 728-734)

Converts a precision matrix to partial correlations.

```r
.precision_to_pcor <- function(Wi, threshold) {
  d <- sqrt(diag(Wi))
  pcor <- -Wi / outer(d, d)              # negate and standardize
  diag(pcor) <- 0                        # zero diagonal
  pcor[abs(pcor) < threshold] <- 0       # apply threshold
  pcor
}
```

**Key**: The negation (`-Wi / ...`) converts the precision matrix sign convention to partial correlation convention. Precision matrix entries are negative when variables are positively partially correlated.

---

## 7. Estimator Return Contract

Every estimator function MUST return a list with at least these three fields:

| Field | Type | Required |
|---|---|---|
| `matrix` | numeric matrix (square, named) | Yes |
| `nodes` | character vector (= row/colnames) | Yes |
| `directed` | logical scalar | Yes |

Additional fields are carried through to the `netobject` as extras:

| Field | Convention |
|---|---|
| `cleaned_data` | Cleaned data used for estimation (matrix or data.frame or NULL) |
| `frequency_matrix` | Raw counts (transition methods) |
| `cor_matrix` | Correlation matrix (association methods) |
| `precision_matrix` | Precision matrix (pcor, glasso) |
| `n` | Sample size (association methods) |
| `p` | Number of variables (association methods) |
| `lambda_selected` | Selected lambda (glasso) |
| `ebic_selected` | EBIC at selected lambda (glasso) |
| `lambda_path` | Full lambda path (glasso) |
| `ebic_path` | EBIC values per lambda (glasso) |
| `gamma` | EBIC gamma hyperparameter (glasso) |

`build_network()` validates that `matrix`, `nodes`, and `directed` are present. It copies `cleaned_data` to the `$data` field, and carries all other extras via:

```r
extras <- setdiff(names(est_result), c("matrix", "nodes", "directed"))
for (key in extras) {
  result[[key]] <- est_result[[key]]
}
```

---

## 8. Summary: What Each Estimator Does

| Estimator | Input → Processing → Output |
|---|---|
| `frequency` | sequences → `.count_transitions()` → integer count matrix |
| `relative` | sequences → `.count_transitions()` → row-normalize → probability matrix |
| `co_occurrence` | sequences → column-pair enumeration → bidirectional counting → symmetric matrix |
| `cor` | data/matrix → `.prepare_association_input()` → `cor()` → zero diagonal → threshold |
| `pcor` | data/matrix → `.prepare_association_input()` → `cor()` → `solve()` → `.precision_to_pcor()` |
| `glasso` | data/matrix → `.prepare_association_input()` → `cor()` → lambda path → EBIC selection → `.precision_to_pcor()` |
