test_that("LEARNING_STATES is a named list", {
  expect_type(LEARNING_STATES, "list")
  expect_true(length(names(LEARNING_STATES)) > 0)
})

test_that("LEARNING_STATES has expected categories", {
  expected <- c("metacognitive", "cognitive", "behavioral", "social",
                "motivational", "affective", "group_regulation", "lms")
  expect_true(all(expected %in% names(LEARNING_STATES)))
})

test_that("get_learning_states returns states from category", {
  states <- get_learning_states("metacognitive")

  expect_type(states, "character")
  expect_true(length(states) > 0)
  expect_true(all(states %in% LEARNING_STATES$metacognitive))
})

test_that("get_learning_states respects n parameter", {
  states <- get_learning_states("cognitive", n = 5)
  expect_equal(length(states), 5)
})

test_that("get_learning_states handles multiple categories", {
  states <- get_learning_states(c("metacognitive", "cognitive"), n = 10)

  expect_equal(length(states), 10)

  valid <- c(LEARNING_STATES$metacognitive, LEARNING_STATES$cognitive)
  expect_true(all(states %in% valid))
})

test_that("get_learning_states is reproducible with seed", {
  states1 <- get_learning_states("metacognitive", n = 5, seed = 42)
  states2 <- get_learning_states("metacognitive", n = 5, seed = 42)

  expect_equal(states1, states2)
})

test_that("list_learning_categories returns data frame", {
  cats <- list_learning_categories()

  expect_s3_class(cats, "data.frame")
  expect_true("Category" %in% names(cats))
  expect_true("metacognitive" %in% cats$Category)
  expect_true("cognitive" %in% cats$Category)
})

test_that("select_states returns correct number", {
  states <- select_states(n_states = 8)

  expect_equal(length(states), 8)
  expect_equal(length(unique(states)), 8)
})

test_that("select_states respects primary_categories", {
  states <- select_states(
    n_states = 5,
    primary_categories = "metacognitive"
  )

  expect_equal(length(states), 5)
  # Primary categories are preferred but not guaranteed exclusive
  all_states <- unlist(LEARNING_STATES)
  expect_true(all(states %in% all_states))
})

test_that("select_states is reproducible with seed", {
  states1 <- select_states(n_states = 6, seed = 42)
  states2 <- select_states(n_states = 6, seed = 42)

  expect_equal(states1, states2)
})
