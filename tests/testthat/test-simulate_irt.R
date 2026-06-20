# Tests for simulate_irt()
# Base-R parameter recovery; mirt recovery optional via skip_if_not_installed.

test_that("2PL has correct structure and value range", {
  fit <- simulate_irt(n_persons = 200, n_items = 15, model = "2PL", seed = 1)

  expect_s3_class(fit, "saqr_sim")
  expect_equal(fit$type, "irt")
  expect_equal(dim(fit$data), c(200L, 15L))
  expect_equal(names(fit$data), paste0("item", 1:15))

  vals <- unlist(as.data.frame(fit))
  expect_true(all(vals %in% c(0L, 1L)))

  expect_equal(fit$params$model, "2PL")
  expect_length(fit$params$a, 15L)
  expect_length(fit$params$b, 15L)
  expect_length(fit$params$theta, 200L)
  expect_equal(fit$params$n_categories, 2L)
  expect_true(all(fit$params$c == 0))
})

test_that("each model produces a valid response matrix (smoke + structure)", {
  f1 <- simulate_irt(n_persons = 150, n_items = 10, model = "1PL", seed = 11)
  f2 <- simulate_irt(n_persons = 150, n_items = 10, model = "2PL", seed = 12)
  f3 <- simulate_irt(n_persons = 150, n_items = 10, model = "3PL", seed = 13)
  fg <- simulate_irt(n_persons = 150, n_items = 10, model = "GRM",
                     n_categories = 4, seed = 14)

  expect_equal(dim(f1$data), c(150L, 10L))
  expect_true(all(unlist(as.data.frame(f1)) %in% c(0L, 1L)))
  expect_true(all(f1$params$a == 1))            # 1PL fixes discrimination

  expect_true(all(unlist(as.data.frame(f2)) %in% c(0L, 1L)))

  expect_true(all(unlist(as.data.frame(f3)) %in% c(0L, 1L)))
  expect_true(all(f3$params$c == 0.2))          # default guessing for 3PL

  expect_equal(dim(fg$data), c(150L, 10L))
  expect_true(all(unlist(as.data.frame(fg)) %in% 0:3))
  expect_equal(dim(fg$params$thresholds), c(10L, 3L))
})

test_that("reproducibility: same seed identical, different seed differs", {
  a <- simulate_irt(n_persons = 100, n_items = 8, model = "2PL", seed = 99)
  b <- simulate_irt(n_persons = 100, n_items = 8, model = "2PL", seed = 99)
  d <- simulate_irt(n_persons = 100, n_items = 8, model = "2PL", seed = 100)

  expect_identical(as.data.frame(a), as.data.frame(b))
  expect_identical(a$params$theta, b$params$theta)
  expect_false(identical(as.data.frame(a), as.data.frame(d)))
})

test_that("difficulty recovery: proportion-correct correlates negatively with b", {
  fit <- simulate_irt(n_persons = 4000, n_items = 30, model = "2PL", seed = 7)
  prop_correct <- colMeans(as.data.frame(fit))
  expect_lt(cor(prop_correct, fit$params$b), -0.8)
})

test_that("ability recovery: total score correlates positively with theta", {
  fit <- simulate_irt(n_persons = 4000, n_items = 30, model = "2PL", seed = 8)
  total <- rowSums(as.data.frame(fit))
  expect_gt(cor(total, fit$params$theta), 0.8)
})

test_that("1PL also recovers difficulty and ability", {
  fit <- simulate_irt(n_persons = 4000, n_items = 30, model = "1PL", seed = 21)
  prop_correct <- colMeans(as.data.frame(fit))
  total <- rowSums(as.data.frame(fit))
  expect_lt(cor(prop_correct, fit$params$b), -0.8)
  expect_gt(cor(total, fit$params$theta), 0.8)
})

test_that("GRM recovers ability and difficulty (ordered categories)", {
  fit <- simulate_irt(n_persons = 4000, n_items = 20, model = "GRM",
                      n_categories = 4, seed = 31)
  total <- rowSums(as.data.frame(fit))
  expect_gt(cor(total, fit$params$theta), 0.8)
  # Higher b -> lower mean category response.
  mean_cat <- colMeans(as.data.frame(fit))
  expect_lt(cor(mean_cat, fit$params$b), -0.8)
})

test_that("2PL: higher discrimination -> higher item-total point-biserial (soft)", {
  fit <- simulate_irt(n_persons = 5000, n_items = 30, model = "2PL", seed = 41)
  d <- as.data.frame(fit)
  total <- rowSums(d)
  # Point-biserial of each item against the total score.
  pb <- vapply(seq_along(d), function(j) {
    rest <- total - d[[j]]               # rest-score to avoid self-inclusion
    cor(d[[j]], rest)
  }, numeric(1))
  # Softer check: positive association between a and discrimination index.
  expect_gt(cor(pb, fit$params$a), 0.3)
})

test_that("supplied parameters are validated and used", {
  a <- runif(5, 0.5, 2)
  b <- seq(-2, 2, length.out = 5)
  theta <- rnorm(80)
  fit <- simulate_irt(n_persons = 80, n_items = 5, model = "2PL",
                      a = a, b = b, theta = theta, seed = 5)
  expect_equal(unname(fit$params$a), a)
  expect_equal(unname(fit$params$b), b)
  expect_equal(fit$params$theta, theta)

  expect_error(simulate_irt(n_persons = 80, n_items = 5, model = "2PL",
                            b = c(1, 2)))                 # wrong length
  expect_error(simulate_irt(n_persons = 80, n_items = 5, model = "GRM",
                            n_categories = 2))            # GRM needs > 2
  expect_error(simulate_irt(n_persons = 80, n_items = 5, model = "2PL",
                            n_categories = 3))            # dichotomous needs 2
})

test_that("mirt recovers item parameters when available", {
  skip_if_not_installed("mirt")
  fit <- simulate_irt(n_persons = 2000, n_items = 15, model = "2PL", seed = 61)
  mod <- mirt::mirt(as.data.frame(fit), 1, itemtype = "2PL", verbose = FALSE)
  est <- mirt::coef(mod, IRTpars = TRUE, simplify = TRUE)$items
  expect_gt(cor(est[, "b"], fit$params$b), 0.8)
  expect_gt(cor(est[, "a"], fit$params$a), 0.5)
})
