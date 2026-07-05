# Generates reference/index.html with fully clickable cards and per-function
# deep links. Anchor for a function = tolower(name), except the two select_*
# helpers which share one section. Run: Rscript reference/build_index.R
setwd("/Users/mohammedsaqr/Documents/Github/Saqrlab/reference")

cats <- list(
  list(file="01-statistical.html", title="Statistical simulators",
       desc="Explicit-parameter datasets with known effect sizes for parameter-recovery testing.",
       fns=c("simulate_ttest","simulate_anova","simulate_correlation","simulate_clusters","simulate_prediction")),
  list(file="02-latent.html", title="Latent-variable models",
       desc="Latent profiles/classes, regression, factor analysis, and sequence clusters.",
       fns=c("simulate_lpa","simulate_lca","simulate_regression","simulate_fa","simulate_seq_clusters")),
  list(file="03-longitudinal-multilevel.html", title="Longitudinal & multilevel",
       desc="VAR(1) panel/ESM data, students-in-classes nesting, latent growth curves.",
       fns=c("simulate_longitudinal","simulate_mlm","simulate_growth")),
  list(file="04-irt.html", title="Item Response Theory",
       desc="1PL / 2PL / 3PL / GRM item responses with known item and ability parameters.",
       fns=c("simulate_irt")),
  list(file="05-survival-hmm.html", title="Survival & hidden Markov",
       desc="Time-to-event data with censoring, and hidden-Markov sequences with true latent paths.",
       fns=c("simulate_survival","simulate_hmm")),
  list(file="06-missing-data.html", title="Missing-data mechanisms",
       desc="First-class MCAR / MAR / MNAR injection with calibrated rates and metadata.",
       fns=c("inject_missingness")),
  list(file="07-random-parameter.html", title="Random-parameter generation",
       desc="Seed-driven datasets (15 types) with complexity injection and batch generation.",
       fns=c("simulate_data")),
  list(file="08-tna-simulation.html", title="TNA simulation",
       desc="Transition networks, group/hierarchical/multitype TNA, matrices and probabilities.",
       fns=c("simulate_tna_network","simulate_tna_networks","simulate_tna_matrix","simulate_tna_datasets",
             "simulate_group_tna_networks","simulate_htna","simulate_mlna","simulate_mtna","simulate_matrix",
             "generate_tna_matrix","generate_tna_networks","generate_tna_datasets","generate_group_tna_networks",
             "generate_probabilities","generate_sequence_data","sample_tna")),
  list(file="09-sequences.html", title="Sequences",
       desc="Markov-chain sequences and stability/instability sequence patterns.",
       fns=c("simulate_sequences","simulate_sequences_advanced")),
  list(file="10-networks-graphs.html", title="Networks & graphs",
       desc="igraph / statnet graphs, edge lists, one-hot encodings, hierarchical long data.",
       fns=c("simulate_igraph","simulate_network","simulate_edge_list","simulate_onehot_data","simulate_long_data")),
  list(file="11-reference-data.html", title="Reference data & utilities",
       desc="Learning-state vocabularies, global names, selection helpers, parameter validation.",
       fns=c("LEARNING_STATES","GLOBAL_NAMES","GROUP_REGULATION_ACTIONS","get_learning_states","get_global_names",
             "list_learning_categories","list_name_regions","select_states","smart_select_states","validate_sim_params")),
  list(file="12-comparison-fitting.html", title="Comparison & model fitting",
       desc="Fit TNA models and compare networks, centralities, edge recovery, reliability.",
       fns=c("fit_network_model","compare_networks","compare_centralities","compare_edge_recovery",
             "calculate_edge_recovery","compare_network_estimation","compare_estimation","compare_reliability",
             "compare_tna_models","cross_validate_tna")),
  list(file="13-visualization.html", title="Visualization",
       desc="Plot estimation comparisons, sampling distributions, and TNA comparisons.",
       fns=c("plot_network_estimation","plot_sampling_distribution","plot_tna_comparison")),
  list(file="14-batch-grid-sampling.html", title="Batch, grid & sampling",
       desc="Parameter grids, batch fitting, bootstrap/grid/network simulation, result summaries.",
       fns=c("generate_param_grid","create_param_grid","batch_apply","batch_fit_models","run_grid_simulation",
             "run_bootstrap_simulation","run_bootstrap_iteration","run_network_simulation","run_sampling_analysis",
             "summarize_grid_results","analyze_grid_results","evaluate_bootstrap","summarize_networks","summarize_simulation")),
  list(file="15-lab-infrastructure.html", title="Laboratory infrastructure",
       desc="Unified dispatcher, parameter-recovery scoring, scenario presets, tidy/export.",
       fns=c("simulate","list_simulators","validate_recovery","list_scenarios","get_scenario","run_scenario",
             "tidy_simulation_results","export_simulation","saqr_sim"))
)

anchor_for <- function(fn) {
  if (fn %in% c("select_states", "smart_select_states")) return("select_states-smart_select_states")
  tolower(fn)
}

chip <- function(file, fn) {
  sprintf('<a class="chip" href="%s#%s"><code>%s</code></a>', file, anchor_for(fn), fn)
}

card <- function(c) {
  chips <- paste(vapply(c$fns, function(f) chip(c$file, f), character(1)), collapse = "")
  sprintf(
'    <div class="card" data-href="%s">
      <h2><a href="%s">%s</a><span class="badge">%d</span></h2>
      <p class="desc">%s</p>
      <div class="fns">%s</div>
      <a class="open" href="%s">Open category &rarr;</a>
    </div>', c$file, c$file, c$title, length(c$fns), c$desc, chips, c$file)
}

cards <- paste(vapply(cats, card, character(1)), collapse = "\n\n")
total <- sum(vapply(cats, function(c) length(c$fns), integer(1)))

html <- sprintf('<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Saqrlab Function Reference (v0.4.0)</title>
<style>
  :root { --accent:#2c3e50; --accent2:#18bc9c; --bg:#f7f9fa; --card:#fff; --muted:#6c757d; }
  * { box-sizing:border-box; }
  body { margin:0; font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif; color:#2c3e50; background:var(--bg); line-height:1.5; }
  a { color:var(--accent2); }
  header { background:var(--accent); color:#fff; padding:2.4rem 1.5rem 2rem; }
  header .wrap { max-width:1100px; margin:0 auto; }
  header h1 { margin:0 0 .3rem; font-size:1.9rem; }
  header p { margin:.25rem 0; opacity:.92; }
  header code { background:rgba(255,255,255,.15); padding:.1rem .4rem; border-radius:4px; }
  main { max-width:1100px; margin:0 auto; padding:1.6rem 1.5rem 3rem; }
  .meta { color:var(--muted); margin:0 0 1.6rem; font-size:.95rem; }
  .grid { display:grid; grid-template-columns:repeat(auto-fill,minmax(340px,1fr)); gap:1.1rem; }
  .card { background:var(--card); border:1px solid #e3e8ec; border-radius:10px; padding:1.05rem 1.2rem 1.1rem; transition:box-shadow .15s,transform .15s,border-color .15s; display:flex; flex-direction:column; cursor:pointer; }
  .card:hover { box-shadow:0 8px 26px rgba(44,62,80,.14); transform:translateY(-2px); border-color:var(--accent2); }
  .card h2 { margin:0 0 .15rem; font-size:1.14rem; }
  .card h2 a { color:var(--accent); text-decoration:none; }
  .card h2 a:hover { color:var(--accent2); text-decoration:underline; }
  .card .desc { color:var(--muted); font-size:.9rem; margin:.1rem 0 .7rem; }
  .badge { display:inline-block; background:var(--accent2); color:#fff; font-size:.72rem; font-weight:600; border-radius:20px; padding:.05rem .55rem; margin-left:.45rem; vertical-align:middle; }
  .fns { margin:0 0 .7rem; display:flex; flex-wrap:wrap; gap:.3rem; }
  .chip { text-decoration:none; }
  .chip code { font-size:.8rem; background:#eef3f5; color:#2c3e50; padding:.12rem .4rem; border-radius:5px; display:inline-block; border:1px solid #e0e7ea; transition:background .12s,color .12s; }
  .chip:hover code { background:var(--accent2); color:#fff; border-color:var(--accent2); }
  .open { margin-top:auto; font-size:.85rem; font-weight:600; text-decoration:none; }
  .open:hover { text-decoration:underline; }
  footer { max-width:1100px; margin:0 auto; padding:0 1.5rem 2.5rem; color:var(--muted); font-size:.85rem; }
</style>
</head>
<body>
<header><div class="wrap">
  <h1>Saqrlab &mdash; Function Reference</h1>
  <p>A modern laboratory for data simulation. Version <code>0.4.0</code> &middot; %d functions across %d categories.</p>
  <p>Click any <strong>function name</strong> to jump straight to its documented example, or click a card to open the whole category.</p>
</div></header>
<main>
  <p class="meta">Every example was executed during rendering &mdash; the outputs and figures on each page are real.</p>
  <div class="grid">

%s

  </div>
</main>
<footer>Generated for Saqrlab 0.4.0. Each linked page is fully self-contained (figures and outputs embedded).</footer>
<script>
  // Make the whole card clickable, while keeping inner links (chips, title) working.
  document.querySelectorAll(".card").forEach(function(card){
    card.addEventListener("click", function(e){
      if (e.target.closest("a")) return;            // a real link was clicked
      window.location.href = card.getAttribute("data-href");
    });
  });
</script>
</body>
</html>
', total, length(cats), cards)

writeLines(html, "index.html")
cat("Wrote index.html with", total, "deep-linked function chips across", length(cats), "cards\n")
