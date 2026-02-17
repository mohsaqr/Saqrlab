# 02 — The Entry Point: `build_network()`

`build_network()` is the single public function users call to estimate any network. It lives in `R/build_network.R` (lines 111-224). This guide walks through every line.

---

## 1. Function Signature (lines 111-117)

```r
build_network <- function(data,
                          method,
                          params = list(),
                          scaling = NULL,
                          threshold = 0,
                          level = NULL,
                          id_col = NULL)
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| `data` | data.frame or matrix | required | Sequences (transition) or frequencies/correlations (association) |
| `method` | character(1) | required | Estimator name or alias |
| `params` | list | `list()` | Method-specific arguments (composability key) |
| `scaling` | character or NULL | NULL | Post-estimation scaling(s) |
| `threshold` | numeric(1) | 0 | Absolute-value threshold for zeroing edges |
| `level` | character or NULL | NULL | Multilevel decomposition: `"between"`, `"within"`, `"both"` |
| `id_col` | character or NULL | NULL | ID column for multilevel decomposition |

---

## 2. Input Validation (lines 118-146)

### Basic stopifnot (lines 118-120)

```r
stopifnot(is.character(method), length(method) == 1)
stopifnot(is.list(params))
stopifnot(is.numeric(threshold), length(threshold) == 1, threshold >= 0)
```

### Method alias resolution (line 123)

```r
method <- .resolve_method_alias(method)
```

Calls the alias table in `R/estimate_network.R` (lines 44-62). Maps e.g. `"tna"` → `"relative"`, `"ebicglasso"` → `"glasso"`. Unknown names pass through.

### Level parameter validation (lines 126-135)

```r
if (!is.null(level)) {
  level <- match.arg(level, c("between", "within", "both"))
  if (is.null(id_col)) {
    stop("'id_col' is required when 'level' is specified.", call. = FALSE)
  }
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame when 'level' is specified.", call. = FALSE)
  }
}
```

- `level` must be one of three values
- `id_col` is required when `level` is set
- Matrix input is not allowed with multilevel decomposition

### Scaling validation (lines 138-146)

```r
if (!is.null(scaling)) {
  valid_scaling <- c("minmax", "max", "rank", "normalize")
  bad <- setdiff(scaling, valid_scaling)
  if (length(bad) > 0) {
    stop("Unknown scaling method(s): ", paste(bad, collapse = ", "),
         ". Options: ", paste(valid_scaling, collapse = ", "), call. = FALSE)
  }
}
```

Scaling is validated before any estimation work begins.

---

## 3. Estimator Lookup (line 149)

```r
estimator <- get_estimator(method)
```

Retrieves from the global `.estimator_registry` environment. Returns `list(fn, description, directed)`. If the estimator doesn't exist, `get_estimator()` throws an informative error listing all available estimators.

---

## 4. Level = "both": Recursive Dispatch (lines 152-164)

```r
if (identical(level, "both")) {
  between <- build_network(
    data, method = method, params = params, scaling = scaling,
    threshold = threshold, level = "between", id_col = id_col
  )
  within_net <- build_network(
    data, method = method, params = params, scaling = scaling,
    threshold = threshold, level = "within", id_col = id_col
  )
  result <- list(between = between, within = within_net, method = method)
  class(result) <- "netobject_ml"
  return(result)
}
```

**Key**: `level = "both"` recursively calls `build_network()` twice — once with `"between"` and once with `"within"`. Returns early with `netobject_ml` class (not `netobject`).

---

## 5. Multilevel Decomposition (lines 167-169)

```r
if (!is.null(level) && !estimator$directed) {
  data <- .decompose_multilevel(data, id_col = id_col, level = level)
}
```

Only applied to association methods (`!estimator$directed`). Transition methods ignore the `level` parameter.

### `.decompose_multilevel()` (file: `R/estimate_network.R`, lines 155-210)

**Between-person** (`level = "between"`):
```r
mat$.id <- id_vals
agg <- aggregate(. ~ .id, data = mat, FUN = mean)
result <- agg[, names(agg) != ".id", drop = FALSE]
```
Aggregates to person means. Result has one row per person.

**Within-person** (`level = "within"`):
1. Drop persons with < 2 observations:
   ```r
   tab <- table(id_vals)
   multi <- names(tab[tab >= 2])
   keep_rows <- id_vals %in% multi
   ```
2. Person-mean center each variable:
   ```r
   for (j in seq_len(ncol(mat_m))) {
     mat_m[, j] <- mat_m[, j] - ave(mat_m[, j], id_vals, FUN = mean)
   }
   ```

**Key**: `ave(x, group, FUN = mean)` returns group means expanded to original length. Subtracting centers each variable within each person.

---

## 6. Estimator Call (line 172)

```r
est_result <- do.call(estimator$fn, c(list(data = data), params))
```

**Key**: `params` is splatted alongside `data = data`. This is the composability mechanism — any params the estimator accepts can be passed through without `build_network()` knowing about them.

For example, for glasso: `do.call(.estimator_glasso, list(data = data, gamma = 0.5, nlambda = 50))`.

---

## 7. Validate Estimator Output (lines 175-182)

```r
if (!is.list(est_result) ||
    is.null(est_result$matrix) ||
    is.null(est_result$nodes) ||
    is.null(est_result$directed)) {
  stop("Estimator '", method,
       "' must return a list with 'matrix', 'nodes', and 'directed'.",
       call. = FALSE)
}
```

This is the contract enforcement. Every estimator — built-in or custom — must return these three fields.

---

## 8. Apply Scaling (lines 189-191)

```r
if (!is.null(scaling)) {
  net_matrix <- .apply_scaling(net_matrix, scaling)
}
```

### `.apply_scaling()` (file: `R/estimate_network.R`, lines 69-108)

Applies each scaling in order:

| Scaling | Logic |
|---|---|
| `"minmax"` | Non-zero values rescaled to [0, 1]: `(x - min) / (max - min)` |
| `"max"` | Divide all values by `max(abs(mat))` |
| `"rank"` | Replace non-zero values with their ranks |
| `"normalize"` | Divide each row by its `rowSums(abs())` — rows sum to 1 in absolute value |

**Key**: All scaling methods preserve zeros. The `minmax` and `rank` scalings only operate on non-zero entries.

**Key**: Multiple scalings are applied in sequence. `c("rank", "minmax")` first ranks, then rescales ranks to [0, 1].

---

## 9. Apply Threshold (lines 194-196)

```r
if (threshold > 0) {
  net_matrix[abs(net_matrix) < threshold] <- 0
}
```

Sets entries with absolute value below the threshold to zero. This is in addition to any method-level thresholding (e.g., `threshold` param passed to the estimator).

---

## 10. Extract Edges (line 199)

```r
edges <- .extract_edges_from_matrix(net_matrix, directed = directed)
```

### `.extract_edges_from_matrix()` (file: `R/estimate_network.R`, lines 119-142)

```r
if (directed) {
  idx <- which(mat != 0 & row(mat) != col(mat), arr.ind = TRUE)
} else {
  idx <- which(upper.tri(mat) & mat != 0, arr.ind = TRUE)
}
```

- **Directed**: all non-zero, non-diagonal entries → separate A→B and B→A edges
- **Undirected**: upper triangle only → one edge per pair

Returns: `data.frame(from, to, weight)` with character node names.

---

## 11. Build `netobject` (lines 202-223)

```r
result <- list(
  data = est_result$cleaned_data,
  matrix = net_matrix,
  nodes = nodes,
  directed = directed,
  method = method,
  params = params,
  scaling = scaling,
  threshold = threshold,
  n_nodes = length(nodes),
  n_edges = nrow(edges),
  edges = edges,
  level = level
)

# Carry over method-specific extras
extras <- setdiff(names(est_result), c("matrix", "nodes", "directed"))
for (key in extras) {
  result[[key]] <- est_result[[key]]
}

structure(result, class = "netobject")
```

**Key**: `data` is set to `est_result$cleaned_data`, which is:
- For association methods (data frame input): the cleaned numeric matrix
- For transition methods: the original input data frame
- For matrix input: NULL

**Key**: All extra fields from the estimator (e.g., `precision_matrix`, `cor_matrix`, `lambda_selected`) are copied verbatim into the netobject. This means custom estimators can attach any metadata they want, and it will be preserved.

**Key**: `structure(result, class = "netobject")` sets the S3 class. This enables `print.netobject()`, `plot.netobject()`, and `predictability.netobject()` dispatch.

---

## 12. Complete Data Flow Summary

```
build_network(data, method = "glasso", params = list(gamma = 0.5, nlambda = 50))
  │
  ├─ method = .resolve_method_alias("glasso")  → "glasso"
  ├─ estimator = get_estimator("glasso")       → {fn, desc, directed=FALSE}
  │
  ├─ [level == "both"?]  → recursive dispatch → netobject_ml
  ├─ [level != NULL?]    → .decompose_multilevel(data, ...)
  │
  ├─ est_result = do.call(.estimator_glasso, list(data=data, gamma=0.5, nlambda=50))
  │   └─ .prepare_association_input()  → clean → cor() → S
  │   └─ .compute_lambda_path()        → lambda_path
  │   └─ .select_ebic()                → best precision matrix
  │   └─ .precision_to_pcor()          → partial correlations
  │   └─ returns: {matrix, nodes, directed, cleaned_data, precision_matrix, ...}
  │
  ├─ validate: matrix, nodes, directed present?
  ├─ .apply_scaling(matrix, scaling)
  ├─ threshold: abs(x) < threshold → 0
  ├─ .extract_edges_from_matrix()
  │
  └─ structure(list(...), class = "netobject")
```

---

## 13. The Deprecated Wrapper: `estimate_network()`

In `R/estimate_network.R` (lines 16-33), `estimate_network()` simply calls `.Deprecated("build_network")` then delegates:

```r
estimate_network <- function(data, method = "relative", ...) {
  .Deprecated("build_network")
  build_network(data = data, method = method, ...)
}
```

**Key difference**: `estimate_network()` defaults `method = "relative"` for backward compatibility. `build_network()` requires `method` to be specified (no default).

**Testing note**: In testthat edition 3, `.Deprecated()` emits a warning, which is treated as a test failure. All tests calling `estimate_network()` must use `suppressWarnings()`.
