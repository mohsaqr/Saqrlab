test_that("GLOBAL_NAMES is a named list", {
  expect_type(GLOBAL_NAMES, "list")
  expect_true(length(names(GLOBAL_NAMES)) > 0)
})

test_that("GLOBAL_NAMES has expected regions", {
  expected_regions <- c("arab", "east_asia", "south_asia", "nordic", "western_europe")
  expect_true(all(expected_regions %in% names(GLOBAL_NAMES)))
})

test_that("GLOBAL_NAMES has sufficient names", {
  total_names <- sum(sapply(GLOBAL_NAMES, length))
  expect_gte(total_names, 1000)
})

test_that("get_global_names returns correct number", {
  names <- get_global_names(10)
  expect_equal(length(names), 10)

  names <- get_global_names(50)
  expect_equal(length(names), 50)
})

test_that("get_global_names returns unique names", {
  names <- get_global_names(100)
  expect_equal(length(names), length(unique(names)))
})

test_that("get_global_names respects regions parameter", {
  arab_names <- get_global_names(10, regions = "arab")
  expect_equal(length(arab_names), 10)
  expect_true(all(arab_names %in% GLOBAL_NAMES$arab))
})

test_that("get_global_names handles region shortcuts", {
  # "africa" should expand to multiple African regions
  africa_names <- get_global_names(20, regions = "africa")
  expect_equal(length(africa_names), 20)
})

test_that("get_global_names handles multiple regions", {
  names <- get_global_names(20, regions = c("arab", "nordic"))

  # Names should come from either region
  valid_names <- c(GLOBAL_NAMES$arab, GLOBAL_NAMES$nordic)
  expect_true(all(names %in% valid_names))
})

test_that("get_global_names is reproducible with seed", {
  names1 <- get_global_names(20, seed = 42)
  names2 <- get_global_names(20, seed = 42)

  expect_equal(names1, names2)
})

test_that("list_name_regions returns data frame", {
  regions <- list_name_regions()

  expect_s3_class(regions, "data.frame")
  expect_true(nrow(regions) >= 20)
  expect_true("Region" %in% names(regions))
})

test_that("get_global_names handles 'all' regions", {
  names <- get_global_names(100, regions = "all")
  expect_equal(length(names), 100)
})
