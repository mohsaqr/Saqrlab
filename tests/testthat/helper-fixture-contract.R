# Shared helpers for the cross-package fixture contract guard.
# (testthat auto-sources helper-*.R before tests; the golden generator below
#  also sources this file so the fingerprint logic lives in exactly one place.)

# The six fixture types consumed by JStats / Carm / validation, and the column
# contract each consumer relies on. See test-fixture-contract.R for the why.
fixture_contract <- list(
  ttest           = c("group", "score"),
  anova           = c("group", "score"),
  correlation     = c("x1", "x2", "x3", "x4"),
  clusters        = c("x1", "x2", "x3", "x4", "x5", "true_cluster"),
  factor_analysis = c("x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8"),
  prediction      = c("y", "x1", "x2", "x3", "x4", "cat1", "cat2")
)

# Representative seeds within every type's real generation range (1..1000).
fixture_seeds <- c(1L, 2L, 3L, 100L)

# complexity="auto" batch cases mirroring saqrlab-batch-ref.R (SEED*mult + s).
fixture_batch_seeds <- c(
  ttest       = 202610L,
  anova       = 202620L,
  correlation = 202630L,
  prediction  = 202640L
)

# Deterministic, environment-stable fingerprint of a generated data.frame:
# captures dim, column names + order, column classes, and every value to
# 10 significant digits (mirrors the 12-digit JSON contract without float noise).
fingerprint_df <- function(d) {
  col_parts <- vapply(names(d), function(nm) {
    col <- d[[nm]]
    vals <- if (is.numeric(col)) {
      formatC(col, digits = 10, format = "g")
    } else {
      as.character(col)
    }
    paste0(nm, "<", class(col)[1], ">:", paste(vals, collapse = ","))
  }, character(1))
  paste(c(sprintf("dim=%dx%d", nrow(d), ncol(d)), col_parts), collapse = "|")
}

# Build the full golden map: key -> fingerprint, for both the default-complexity
# fixtures and the complexity="auto" batch fixtures.
build_fixture_fingerprints <- function() {
  out <- list()
  for (ty in names(fixture_contract)) {
    for (s in fixture_seeds) {
      out[[sprintf("%s_seed%03d", ty, s)]] <-
        fingerprint_df(simulate_data(ty, seed = s))
    }
  }
  for (ty in names(fixture_batch_seeds)) {
    out[[sprintf("%s_auto%d", ty, fixture_batch_seeds[[ty]])]] <-
      fingerprint_df(simulate_data(ty, seed = fixture_batch_seeds[[ty]], complexity = "auto"))
  }
  out
}

# Path to the committed golden reference.
fixture_golden_path <- function() {
  testthat::test_path("fixtures", "fixture-contract-golden.rds")
}
