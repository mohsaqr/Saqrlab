# Saqrlab — HTML Function Reference

Self-contained HTML documentation for every exported function, organised by
category. Open **[`index.html`](index.html)** in a browser for the clickable
landing page (each function name jumps to its documented example).

Each `NN-category.html` is generated from the matching `NN-category.Rmd` by
running the package's live examples — the outputs and figures shown are real.

| Page | Category |
|------|----------|
| `01-statistical.html` | Statistical simulators |
| `02-latent.html` | Latent-variable models |
| `03-longitudinal-multilevel.html` | Longitudinal & multilevel |
| `04-irt.html` | Item Response Theory |
| `05-survival-hmm.html` | Survival & hidden Markov |
| `06-missing-data.html` | Missing-data mechanisms |
| `07-random-parameter.html` | Random-parameter generation |
| `08-tna-simulation.html` | TNA simulation |
| `09-sequences.html` | Sequences |
| `10-networks-graphs.html` | Networks & graphs |
| `11-reference-data.html` | Reference data & utilities |
| `12-comparison-fitting.html` | Comparison & model fitting |
| `13-visualization.html` | Visualization |
| `14-batch-grid-sampling.html` | Batch, grid & sampling |
| `15-lab-infrastructure.html` | Laboratory infrastructure |

## Regenerating

```r
# re-render one category after a code change:
rmarkdown::render("reference/01-statistical.Rmd")

# rebuild the index landing page:
Rscript reference/build_index.R
```
