# 05 — Estimator Registry & S3 Methods

This guide covers the global estimator registry (`R/estimator_registry.R`, 157 lines), all S3 methods for `netobject`, `netobject_ml`, `saqr_bootstrap`, and `saqr_permutation`, and the predictability system.

---

## 1. Estimator Registry

### 1.1 The Registry Object (line 4)

```r
.estimator_registry <- new.env(parent = emptyenv())
```

A clean environment with no parent scope. Using `parent = emptyenv()` prevents accidental inheritance from global scope.

### 1.2 CRUD Operations

#### `register_estimator()` (lines 41-51) — exported

```r
register_estimator <- function(name, fn, description, directed) {
  stopifnot(
    is.character(name), length(name) == 1, nzchar(name),
    is.function(fn),
    is.character(description), length(description) == 1,
    is.logical(directed), length(directed) == 1
  )
  entry <- list(fn = fn, description = description, directed = directed)
  assign(name, entry, envir = .estimator_registry)
  invisible(NULL)
}
```

Stores a list with three fields: `fn` (the estimator function), `description` (human-readable), `directed` (logical).

#### `get_estimator()` (lines 66-77) — exported

```r
get_estimator <- function(name) {
  stopifnot(is.character(name), length(name) == 1)
  if (!exists(name, envir = .estimator_registry, inherits = FALSE)) {
    available <- ls(envir = .estimator_registry)
    stop(sprintf("Estimator '%s' not found. Available: %s",
                 name, paste(sort(available), collapse = ", ")),
         call. = FALSE)
  }
  get(name, envir = .estimator_registry, inherits = FALSE)
}
```

**Key**: `inherits = FALSE` prevents looking up the parent environment chain. The error message lists all available estimators.

#### `list_estimators()` (lines 91-112) — exported

Returns a data frame with columns `name`, `description`, `directed`. Sorts by name.

#### `remove_estimator()` (lines 127-137) — exported

```r
remove_estimator <- function(name) {
  stopifnot(is.character(name), length(name) == 1)
  if (!exists(name, envir = .estimator_registry, inherits = FALSE)) {
    stop(sprintf("Estimator '%s' not found in registry.", name), call. = FALSE)
  }
  rm(list = name, envir = .estimator_registry)
  invisible(NULL)
}
```

### 1.3 Built-in Registration (lines 142-157)

```r
.register_builtin_estimators <- function() {
  register_estimator("relative", .estimator_relative,
                     "Row-normalized transition probabilities", directed = TRUE)
  register_estimator("frequency", .estimator_frequency,
                     "Raw transition frequency counts", directed = TRUE)
  register_estimator("co_occurrence", .estimator_co_occurrence,
                     "Co-occurrence within sequences", directed = FALSE)
  register_estimator("cor", .estimator_cor,
                     "Pairwise correlation network", directed = FALSE)
  register_estimator("pcor", .estimator_pcor,
                     "Unregularized partial correlations", directed = FALSE)
  register_estimator("glasso", .estimator_glasso,
                     "EBICglasso regularized partial correlations",
                     directed = FALSE)
}
```

Called from `.onLoad()` in `R/Saqrlab-package.R`:

```r
.onLoad <- function(libname, pkgname) {
  .register_builtin_estimators()
}
```

### 1.4 Custom Estimator Contract

A custom estimator must:

1. Accept `data` as its first argument and `...` for additional params
2. Return a list with at least: `matrix` (square numeric), `nodes` (character), `directed` (logical)
3. Optionally return `cleaned_data` and any extras

Example:

```r
custom_fn <- function(data, ...) {
  numeric_cols <- vapply(data, is.numeric, logical(1))
  mat <- as.matrix(data[, numeric_cols, drop = FALSE])
  S <- cor(mat)
  diag(S) <- 0
  list(matrix = S, nodes = colnames(S), directed = FALSE,
       cleaned_data = mat)
}
register_estimator("my_cor", custom_fn, "Custom correlation", directed = FALSE)
```

---

## 2. S3 Methods: `netobject`

### 2.1 `print.netobject()` (file: `R/build_network.R`, lines 235-283)

Displays a human-readable summary:

```
Partial Correlation Network (EBICglasso) [undirected]
  Data: 80 x 5
  Nodes: 5  |  Edges: 7
  Sample size: 80
  Gamma: 0.50  |  Lambda: 0.1234
```

**Key details**:
- Method labels are hardcoded in a named vector (lines 236-243)
- Directedness shown as `[directed]` or `[undirected]`
- Level shown as `[between-person]` or `[within-person]` if set
- Data dimensions shown if `$data` is not NULL
- Gamma/lambda shown only for `method == "glasso"`
- Scaling shown only if not NULL
- Threshold shown only if > 0
- Returns `invisible(x)`

### 2.2 `plot.netobject()` (lines 350-379)

Uses `cograph::splot()` for visualization.

```r
plot.netobject <- function(x, predictability = TRUE,
                           pie_color = "#377EB8", ...) {
  if (!requireNamespace("cograph", quietly = TRUE)) {
    stop("Package 'cograph' is required for plotting.")
  }

  node_cols <- .node_colors(x$n_nodes)

  dots <- list(
    x = x$matrix,
    directed = x$directed,
    node_fill = node_cols,
    edge_labels = TRUE,
    edge_label_size = 0.65,
    node_size = 8,
    theme = "colorblind",
    ...
  )

  if (predictability && x$method %in% c("glasso", "pcor", "cor")) {
    r2 <- predictability.netobject(x)
    dots$pie_values <- r2
    dots$pie_colors <- pie_color
  }

  do.call(cograph::splot, dots)
}
```

**Key**: Predictability pies only shown for association methods. Transition methods don't have precision matrices.

**Key**: `do.call()` used to pass args programmatically — this allows `...` to override defaults.

### 2.3 Node Colors: `.node_colors()` (lines 312-331)

Hand-picked palette of 15 distinct pastel colors, cycled via `rep_len(palette, p)`.

```r
palette <- c(
  "#A8D8EA",  # light blue
  "#FFCAB1",  # peach
  "#B5EAD7",  # mint
  "#E2B6CF",  # mauve
  "#FFDAC1",  # apricot
  "#C7CEEA",  # lavender
  "#F3E8C0",  # cream
  "#D4F0C0",  # pistachio
  "#F5C6D0",  # pink
  "#B8E0D2",  # seafoam
  "#EAC8A0",  # sand
  "#C8B8DB",  # lilac
  "#A0D2DB",  # teal
  "#F0D9A0",  # gold
  "#D8A8C8"   # orchid
)
```

---

## 3. S3 Methods: `netobject_ml`

### 3.1 `print.netobject_ml()` (lines 292-307)

```
Multilevel Network (method: glasso)
-- Between-person --
  Nodes: 5  |  Edges: 7
  Sample size: 30 (unique persons)
-- Within-person --
  Nodes: 5  |  Edges: 8
  Sample size: 150 (observations)
```

### 3.2 `plot.netobject_ml()` (lines 398-444)

Side-by-side plot using `par(mfrow = c(1, 2))` with `on.exit(graphics::par(old_par))` to restore.

**Key**: `on.exit()` ensures `par()` is restored even if an error occurs. This is essential for well-behaved plot methods.

---

## 4. S3 Methods: `saqr_bootstrap`

### 4.1 `print.saqr_bootstrap()` (file: `R/bootstrap_network.R`, lines 512-543)

```
Bootstrap: Transition Network (relative probabilities) [directed]
  Iterations: 500  |  CI level: 0.05  |  Inference: stability
  Nodes: 9  |  Original edges: 72  |  Significant edges: 45
  Consistency range: [0.75, 1.25]
```

### 4.2 `summary.saqr_bootstrap()` (lines 554-556)

Simply returns `object$summary` — the long-format data frame.

### 4.3 `plot.saqr_bootstrap()` (lines 570-598)

Side-by-side plot: Original network vs. Significant-only network. Same `par(mfrow = c(1, 2))` pattern with `on.exit()`.

---

## 5. S3 Methods: `saqr_permutation`

### 5.1 `print.saqr_permutation()` (file: `R/permutation_test.R`, lines 494-524)

```
Permutation Test: Transition Network (relative probabilities) [directed]
  Iterations: 1000  |  Alpha: 0.05  |  Paired  |  Adjust: BH
  Nodes: 9  |  Edges tested: 72  |  Significant: 12
```

### 5.2 `summary.saqr_permutation()` (lines 535-537)

Returns `object$summary`.

### 5.3 `plot.saqr_permutation()` (lines 551-579)

Side-by-side: Network X vs. Significant Differences. Same `par`/`on.exit` pattern.

---

## 6. Predictability System

### 6.1 Generic: `predictability()` (file: `R/build_network.R`, lines 489-491)

```r
predictability <- function(object, ...) {
  UseMethod("predictability")
}
```

### 6.2 `predictability.netobject()` (lines 496-521)

Two computation paths based on method:

**For glasso/pcor** (lines 497-500) — from precision matrix:

```r
if (object$method %in% c("glasso", "pcor")) {
  omega_diag <- diag(object$precision_matrix)
  r2 <- 1 - 1 / omega_diag
}
```

Formula: R²_j = 1 - 1/Omega_jj, where Omega is the precision matrix. This is the analytical formula from Haslbeck & Waldorp (2018).

**For cor** (lines 502-516) — multiple R² from correlation matrix:

```r
r2 <- vapply(seq_len(p), function(j) {
  neighbors <- which(net[j, ] != 0)
  if (length(neighbors) == 0L) return(0)
  if (length(neighbors) == 1L) return(S[neighbors, j]^2)
  r_vec <- S[neighbors, j]
  R_nn  <- S[neighbors, neighbors]
  tryCatch(
    as.numeric(crossprod(r_vec, solve(R_nn, r_vec))),
    error = function(e) 0
  )
}, numeric(1))
```

For each node j, computes the multiple R² of j regressed on its neighbors (nodes with non-zero edges):

- **0 neighbors**: R² = 0
- **1 neighbor**: R² = r² (squared bivariate correlation)
- **>1 neighbors**: R² = r'R_nn⁻¹r (multivariate R² formula using correlation submatrix)

**Post-processing** (lines 518-520):

```r
r2 <- pmin(pmax(r2, 0), 1)              # clamp to [0, 1]
names(r2) <- colnames(object$matrix)
```

**Validation**: Tested against mgm package — r = 0.999, mean |diff| = 0.008.

### 6.3 `predictability.netobject_ml()` (lines 526-531)

Returns a list:

```r
list(
  between = predictability(object$between),
  within  = predictability(object$within)
)
```

---

## 7. Global Variables Declaration

In `R/Saqrlab-package.R` (lines 5-11):

```r
utils::globalVariables(c(
  "density", "mean_val", "mean_value", "median_val", "sd_val", "sd_value",
  "value", "metric", "category", "model_type", "iteration", "mean", "sd",
  "sim_alpha", "sim_diag_c",
  ".grp_key", ".", ".SD", ":=",
  "cr_lower", "cr_upper", "weight_x", "weight_y", ".seq_grp"
))
```

These declarations suppress R CMD check NOTEs for data.table NSE symbols (`.SD`, `:=`, `.grp_key`, `.`) and column names used in data.table filtering expressions (`cr_lower`, `cr_upper`, `weight_x`, `weight_y`).

---

## 8. The `convert_sequence_format()` Function

In `R/frequencies.R` (lines 152-181). Converts sequence data between formats:

| Output format | Description |
|---|---|
| `"frequency"` | Count of each action per sequence |
| `"onehot"` | Binary presence/absence per sequence |
| `"edgelist"` | Consecutive transition pairs |
| `"follows"` | Each action paired with its predecessor |

Accepts both wide and long input. Uses `dplyr` and `tidyr` for pivoting.

**Key for association methods**: `convert_sequence_format(data, format = "frequency")` is how users prepare sequence data for association network estimation. The output has one row per sequence, one column per state, with integer counts.

---

## 9. The `frequencies()` Function

In `R/frequencies.R` (lines 69-85). Thin wrapper around `.count_transitions()`:

```r
frequencies <- function(data, action = "Action", id = NULL,
                        time = "Time", cols = NULL,
                        format = c("auto", "long", "wide")) {
  format <- match.arg(format)
  .count_transitions(data, format, action, id, time, cols)
}
```

Returns a square integer transition frequency matrix. Can be passed to `tna::tna()` directly.
