# Tests for scenario presets + tidy results / export

test_that("list_scenarios returns the expected tidy data.frame", {
  s <- list_scenarios()
  expect_s3_class(s, "data.frame")
  expect_identical(names(s), c("scenario", "description", "simulator", "n_cases"))
  expect_true(all(c("power_ttest", "mlm_icc_sweep", "irt_models",
                    "survival_censoring") %in% s$scenario))
  expect_type(s$n_cases, "integer")
  expect_true(all(s$n_cases >= 1L))
})


test_that("get_scenario returns the full list of cases", {
  cases <- get_scenario("power_ttest")
  expect_type(cases, "list")
  expect_length(cases, 6L)
  expect_true(all(vapply(cases, is.list, logical(1))))
})


test_that("get_scenario returns a single case as an arg-list", {
  one <- get_scenario("power_ttest", case = 1)
  expect_type(one, "list")
  expect_true(all(c("n_a", "n_b", "mean_a", "mean_b") %in% names(one)))
})


test_that("unknown scenario errors with valid names listed", {
  expect_error(get_scenario("nope"), "Unknown scenario")
  expect_error(get_scenario("nope"), "power_ttest")
  expect_error(run_scenario("nope"), "Unknown scenario")
})


test_that("get_scenario errors on out-of-range case", {
  expect_error(get_scenario("power_ttest", case = 99), "case must be between")
})


test_that("run_scenario returns a named list of saqr_sim objects", {
  s <- run_scenario("power_ttest", seed = 1)
  expect_type(s, "list")
  expect_length(s, 6L)
  expect_identical(names(s), paste0("case", 1:6))
  expect_true(all(vapply(s, inherits, logical(1), what = "saqr_sim")))
})


test_that("run_scenario reproduces under the same seed", {
  a <- run_scenario("power_ttest", seed = 42)
  b <- run_scenario("power_ttest", seed = 42)
  expect_equal(a$case1$data, b$case1$data)
  expect_equal(a$case3$data, b$case3$data)
})


test_that("run_scenario forwards distinct per-case seeds", {
  s <- run_scenario("power_ttest", seed = 10)
  expect_equal(s$case1$seed, 11L)
  expect_equal(s$case2$seed, 12L)
})


test_that("every scenario runs end to end with valid args", {
  all_names <- list_scenarios()$scenario
  runs <- lapply(all_names, function(nm) run_scenario(nm, seed = 7))
  ok <- vapply(runs, function(r) {
    all(vapply(r, inherits, logical(1), what = "saqr_sim"))
  }, logical(1))
  expect_true(all(ok))
})


test_that("tidy_simulation_results returns one tidy data.frame with .case", {
  s <- run_scenario("power_ttest", seed = 1)
  tidy <- tidy_simulation_results(s)
  expect_s3_class(tidy, "data.frame")
  expect_true(".case" %in% names(tidy))
  expect_setequal(unique(tidy$.case), paste0("case", 1:6))

  expected_rows <- sum(vapply(s, function(sim) nrow(sim$data), integer(1)))
  expect_equal(nrow(tidy), expected_rows)
})


test_that("tidy_simulation_results adds scalar params as columns", {
  s <- run_scenario("power_ttest", seed = 1)
  tidy <- tidy_simulation_results(s)
  # cohens_d is a scalar param of simulate_ttest
  expect_true("cohens_d" %in% names(tidy))
})


test_that("tidy_simulation_results accepts a single saqr_sim", {
  one <- simulate_ttest(n_a = 20, n_b = 20, mean_a = 0, mean_b = 1, seed = 1)
  tidy <- tidy_simulation_results(one)
  expect_s3_class(tidy, "data.frame")
  expect_true(all(tidy$.case == "case1"))
})


test_that("tidy_simulation_results unions columns across heterogeneous cases", {
  mixed <- list(
    a = simulate_ttest(n_a = 10, n_b = 10, mean_a = 0, mean_b = 1, seed = 1),
    b = simulate_anova(n = 10, means = c(1, 2, 3), seed = 2)
  )
  tidy <- tidy_simulation_results(mixed)
  expect_s3_class(tidy, "data.frame")
  expect_setequal(unique(tidy$.case), c("a", "b"))
})


test_that("export_simulation writes a readable CSV", {
  s <- run_scenario("power_ttest", seed = 1)
  path <- tempfile(fileext = ".csv")
  out <- export_simulation(s, path)
  expect_identical(out, path)
  expect_true(file.exists(path))

  back <- utils::read.csv(path, stringsAsFactors = FALSE)
  expect_s3_class(back, "data.frame")
  expect_true(".case" %in% names(back))
  expect_true("score" %in% names(back))
})


test_that("export_simulation accepts a plain data.frame", {
  df <- data.frame(a = 1:3, b = letters[1:3])
  path <- tempfile(fileext = ".csv")
  export_simulation(df, path)
  back <- utils::read.csv(path, stringsAsFactors = FALSE)
  expect_equal(nrow(back), 3L)
  expect_true(all(c("a", "b") %in% names(back)))
})
