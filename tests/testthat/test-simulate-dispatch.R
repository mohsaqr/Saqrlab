# Tests for the unified simulate() dispatcher and list_simulators().

test_that("simulate() dispatches to simulate_ttest with full fidelity", {
  via_dispatch <- simulate(
    "ttest", n_a = 50, n_b = 50, mean_a = 0, mean_b = 0.5, seed = 7
  )
  direct <- simulate_ttest(
    n_a = 50, n_b = 50, mean_a = 0, mean_b = 0.5, seed = 7
  )
  expect_s3_class(via_dispatch, "saqr_sim")
  expect_identical(via_dispatch$data, direct$data)
  expect_identical(via_dispatch$params, direct$params)
  expect_identical(via_dispatch$type, direct$type)
  expect_identical(via_dispatch$seed, direct$seed)
})

test_that("simulate() dispatches to simulate_anova with full fidelity", {
  via_dispatch <- simulate("anova", n = 30, means = c(10, 12, 15), seed = 3)
  direct <- simulate_anova(n = 30, means = c(10, 12, 15), seed = 3)
  expect_s3_class(via_dispatch, "saqr_sim")
  expect_identical(via_dispatch$data, direct$data)
  expect_identical(via_dispatch$params, direct$params)
})

test_that("simulate() dispatches to simulate_correlation with full fidelity", {
  sigma <- matrix(c(1, 0.6, 0.3,
                    0.6, 1, 0.5,
                    0.3, 0.5, 1), nrow = 3)
  via_dispatch <- simulate("correlation", n = 200, sigma = sigma, seed = 11)
  direct <- simulate_correlation(n = 200, sigma = sigma, seed = 11)
  expect_identical(via_dispatch$data, direct$data)
  expect_identical(via_dispatch$params, direct$params)
})

test_that("simulate() dispatches to simulate_regression with full fidelity", {
  coefs <- c(x1 = 2, x2 = -1)
  predictor_sds <- c(x1 = 1, x2 = 1)
  via_dispatch <- simulate(
    "regression", coefs = coefs, predictor_sds = predictor_sds,
    error_sd = 1, n = 100, seed = 5
  )
  direct <- simulate_regression(
    coefs = coefs, predictor_sds = predictor_sds, error_sd = 1, n = 100, seed = 5
  )
  expect_s3_class(via_dispatch, "saqr_sim")
  expect_identical(via_dispatch$data, direct$data)
  expect_identical(via_dispatch$params, direct$params)
})

test_that("seed forwarding yields reproducible output", {
  a <- simulate("ttest", n_a = 40, n_b = 40, mean_a = 1, mean_b = 2, seed = 99)
  b <- simulate("ttest", n_a = 40, n_b = 40, mean_a = 1, mean_b = 2, seed = 99)
  expect_identical(a$data, b$data)

  c <- simulate("ttest", n_a = 40, n_b = 40, mean_a = 1, mean_b = 2, seed = 100)
  expect_false(identical(a$data, c$data))
})

test_that("unknown type errors with a helpful message listing valid types", {
  expect_error(
    simulate("not_a_real_type", n = 10),
    "Unknown simulation type"
  )
  expect_error(
    simulate("not_a_real_type", n = 10),
    "ttest"
  )
  expect_error(
    simulate("not_a_real_type", n = 10),
    "hmm"
  )
})

test_that("list_simulators() returns a tidy data.frame, one row per simulator", {
  tbl <- list_simulators()
  expect_s3_class(tbl, "data.frame")
  expect_identical(names(tbl), c("type", "function", "description"))

  expected_types <- c(
    "ttest", "anova", "correlation", "clusters", "prediction", "regression",
    "lpa", "lca", "fa", "seq_clusters", "longitudinal", "mlm", "growth",
    "irt", "survival", "hmm"
  )
  expect_identical(tbl$type, expected_types)
  expect_identical(nrow(tbl), length(expected_types))
  expect_identical(tbl$`function`, paste0("simulate_", expected_types))
  expect_true(all(nzchar(tbl$description)))
})

test_that("every listed type dispatches without error to a saqr_sim object", {
  # Smoke-test the full map end to end for a handful of types with minimal args.
  r_lca <- simulate(
    "lca",
    item_probs = matrix(c(0.9, 0.1, 0.1, 0.9), nrow = 2, byrow = TRUE),
    class_probs = c(0.5, 0.5), n = 100, seed = 2
  )
  expect_s3_class(r_lca, "saqr_sim")
  expect_identical(r_lca$type, "lca")
})
