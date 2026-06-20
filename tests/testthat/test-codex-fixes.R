# Regression tests for correctness bugs found in the Codex line-by-line review
# (2026-06-20). Each guards a confirmed, reproduced bug.

test_that("inject_missingness MAR/MNAR calibration hits prop even when >0.5", {
  set.seed(1)
  d <- data.frame(y = rnorm(4000), x1 = rnorm(4000))
  for (p in c(0.5, 0.7, 0.9)) {
    mar  <- inject_missingness(d, "MAR",  prop = p, predictor = "x1", cols = "y", seed = 1)
    mnar <- inject_missingness(d, "MNAR", prop = p, cols = "y", seed = 1)
    expect_equal(attr(mar,  "missing_info")$realized_prop, p, tolerance = 0.05)
    expect_equal(attr(mnar, "missing_info")$realized_prop, p, tolerance = 0.05)
  }
})

test_that("inject_missingness MAR missingness still depends on the predictor", {
  set.seed(1)
  d <- data.frame(y = rnorm(4000), x1 = rnorm(4000))
  mar <- inject_missingness(d, "MAR", prop = 0.4, predictor = "x1", cols = "y", seed = 1)
  ind <- attr(mar, "missing_info")$indicator[, 1]
  expect_gt(cor(ind, d$x1), 0.2)   # clearly predictor-driven, not MCAR
})

test_that("compare_centralities does not crash with < 2 common states", {
  skip_if_not_installed("tna")
  mkseq <- function(st, s) { set.seed(s); as.data.frame(matrix(sample(st, 20 * 8, replace = TRUE), 20, 8)) }
  t1 <- tna::tna(mkseq(c("A", "B", "C"), 1))
  t2 <- tna::tna(mkseq(c("A", "X", "Y"), 2))   # only "A" shared
  expect_silent_result <- suppressWarnings(compare_centralities(t1, t2))
  expect_type(expect_silent_result, "list")
})

test_that("long_to_wide handles zero-row and NA-id input without crashing", {
  zero <- Saqrlab:::long_to_wide(
    data.frame(Actor = integer(0), Time = numeric(0), Action = character(0))
  )
  expect_s3_class(zero, "data.frame")
  expect_equal(nrow(zero), 0L)

  with_na <- Saqrlab:::long_to_wide(
    data.frame(Actor = c(1, NA, 1), Time = c(1, 2, 2), Action = c("a", "b", "c"),
               stringsAsFactors = FALSE)
  )
  expect_false(anyNA(with_na$Actor))   # NA ids dropped, not fabricated
})

test_that("safe_bind_rows does not corrupt factor vs character column conflicts", {
  a <- data.frame(g = factor(c("x", "y")), v = 1:2, stringsAsFactors = FALSE)
  b <- data.frame(g = c("z", "w"),         v = 3:4, stringsAsFactors = FALSE)
  out <- Saqrlab:::safe_bind_rows(list(a, b))
  expect_equal(nrow(out), 4L)
  expect_false(anyNA(out$g))            # factor+character merge kept all values
  expect_setequal(out$g, c("x", "y", "z", "w"))
})

test_that("simulate_hmm transition draw never yields NA states (FP-CDF guard)", {
  # near-degenerate rows whose cumsum can drift below 1 in floating point
  tr <- matrix(c(0.9999999, 0.0000001, 0.5, 0.5), 2, 2, byrow = TRUE)
  tr <- tr / rowSums(tr)
  sim <- simulate_hmm(n_sequences = 20, seq_length = 40, n_states = 2, n_symbols = 3,
                      trans = tr, seed = 1)
  expect_false(anyNA(sim$params$hidden_paths))
  expect_false(anyNA(sim$data$symbol))
})
