# ===========================================================================
# Tests for build_hon() — Higher-Order Network construction
# ===========================================================================

# --- Helper: simple trajectories with known structure ---
.make_hon_data <- function() {
  data.frame(
    T1 = c("A", "B", "C", "A", "D"),
    T2 = c("B", "A", "A", "B", "A"),
    T3 = c("C", "B", "B", "C", "B"),
    T4 = c("D", "C", "C", "D", "C"),
    T5 = c("A", "D", "D", "A", NA),
    T6 = c("B", "A", "A", "B", NA),
    T7 = c("C", "B", "B", "C", NA),
    stringsAsFactors = FALSE
  )
}

# ===========================================================================
# Section 1: Input validation
# ===========================================================================
test_that("build_hon rejects non-data.frame non-list input", {
  expect_error(build_hon(42), "data.frame or list")
})

test_that("build_hon rejects empty data.frame", {
  expect_error(build_hon(data.frame()), "at least one")
})

test_that("build_hon rejects max_order < 1", {
  expect_error(build_hon(.make_hon_data(), max_order = 0), "max_order")
})

test_that("build_hon rejects min_freq < 1", {
  expect_error(build_hon(.make_hon_data(), min_freq = 0), "min_freq")
})

test_that("build_hon accepts data.frame input", {
  result <- build_hon(.make_hon_data(), min_freq = 1L)
  expect_s3_class(result, "saqr_hon")
})

test_that("build_hon accepts list input", {
  trajs <- list(c("A", "B", "C"), c("B", "C", "A"))
  result <- build_hon(trajs, min_freq = 1L)
  expect_s3_class(result, "saqr_hon")
})

test_that("build_hon strips trailing NAs from data.frame rows", {
  df <- data.frame(T1 = c("A", "A"), T2 = c("B", "B"), T3 = c("C", NA),
                   stringsAsFactors = FALSE)
  result <- build_hon(df, min_freq = 1L)
  expect_s3_class(result, "saqr_hon")
})

test_that("build_hon collapse_repeats removes adjacent duplicates", {
  trajs <- list(c("A", "A", "B", "B", "C"))
  result <- build_hon(trajs, min_freq = 1L, collapse_repeats = TRUE)
  expect_true(result$n_edges > 0)
})

# ===========================================================================
# Section 2: Observation counting
# ===========================================================================
test_that("counts correct for single trajectory A->B->C", {
  trajs <- list(c("A", "B", "C"))
  count <- .hon_build_observations(trajs, max_order = 2L)
  expect_equal(count[[.hon_encode("A")]][["B"]], 1L)
  expect_equal(count[[.hon_encode("B")]][["C"]], 1L)
  expect_equal(count[[.hon_encode(c("A", "B"))]][["C"]], 1L)
})

test_that("counts accumulate across multiple trajectories", {
  trajs <- list(c("A", "B", "C"), c("A", "B", "D"), c("A", "B", "C"))
  count <- .hon_build_observations(trajs, max_order = 1L)
  expect_equal(count[[.hon_encode("A")]][["B"]], 3L)
  expect_equal(count[[.hon_encode("B")]][["C"]], 2L)
  expect_equal(count[[.hon_encode("B")]][["D"]], 1L)
})

test_that("max_order limits observation depth", {
  trajs <- list(c("A", "B", "C", "D"))
  count1 <- .hon_build_observations(trajs, max_order = 1L)
  expect_false(is.null(count1[[.hon_encode("A")]]))
  expect_true(is.null(count1[[.hon_encode(c("A", "B"))]]))
  count2 <- .hon_build_observations(trajs, max_order = 2L)
  expect_false(is.null(count2[[.hon_encode(c("A", "B"))]]))
  expect_true(is.null(count2[[.hon_encode(c("A", "B", "C"))]]))
})

test_that("short trajectories contribute only possible orders", {
  trajs <- list(c("X", "Y"))
  count <- .hon_build_observations(trajs, max_order = 5L)
  expect_equal(count[[.hon_encode("X")]][["Y"]], 1L)
  all_keys <- ls(count)
  expect_true(all(.hon_key_len(all_keys) == 1L))
})

# ===========================================================================
# Section 3: Distribution building
# ===========================================================================
test_that("distributions sum to 1 for each source", {
  trajs <- list(c("A", "B", "C"), c("A", "B", "D"), c("A", "B", "C"))
  count <- .hon_build_observations(trajs, max_order = 1L)
  distr <- .hon_build_distributions(count, min_freq = 1L)
  for (key in ls(distr)) {
    probs <- distr[[key]]
    if (length(probs) > 0L) {
      expect_equal(sum(probs), 1.0, tolerance = 1e-10)
    }
  }
})

test_that("min_freq filters low-count transitions", {
  trajs <- list(c("A", "B", "C"), c("A", "B", "D"), c("A", "B", "C"))
  count <- .hon_build_observations(trajs, max_order = 1L)
  distr <- .hon_build_distributions(count, min_freq = 2L)
  b_distr <- distr[[.hon_encode("B")]]
  expect_true(is.na(b_distr["D"]) || is.null(b_distr["D"]))
  expect_equal(unname(b_distr["C"]), 1.0)
})

test_that("sources with all counts below min_freq get empty distribution", {
  trajs <- list(c("X", "Y"))
  count <- .hon_build_observations(trajs, max_order = 1L)
  distr <- .hon_build_distributions(count, min_freq = 5L)
  expect_equal(length(distr[[.hon_encode("X")]]), 0L)
})

# ===========================================================================
# Section 4: KL-divergence and threshold
# ===========================================================================
test_that("KLD of identical distributions is 0", {
  d <- c(A = 0.5, B = 0.5)
  expect_equal(.hon_kld(d, d), 0.0)
})

test_that("KLD of peaked vs uniform is log2(2) = 1", {
  a <- c(X = 1.0)
  b <- c(X = 0.5, Y = 0.5)
  expect_equal(.hon_kld(a, b), 1.0)
})

test_that("KLD threshold decreases with more data", {
  count_env <- new.env(hash = TRUE, parent = emptyenv())
  count_env[["k"]] <- c(A = 5L, B = 5L)
  t1 <- .hon_kld_threshold(2L, "k", count_env)
  count_env[["k"]] <- c(A = 50L, B = 50L)
  t2 <- .hon_kld_threshold(2L, "k", count_env)
  expect_true(t2 < t1)
})

test_that("KLD threshold increases with order", {
  count_env <- new.env(hash = TRUE, parent = emptyenv())
  count_env[["k"]] <- c(A = 10L, B = 10L)
  t2 <- .hon_kld_threshold(2L, "k", count_env)
  t3 <- .hon_kld_threshold(3L, "k", count_env)
  expect_true(t3 > t2)
})

# ===========================================================================
# Section 5: End-to-end pipeline
# ===========================================================================
test_that("build_hon returns correct structure", {
  result <- build_hon(.make_hon_data(), max_order = 2L, min_freq = 1L)
  expect_true(is.matrix(result$matrix))
  expect_true(is.data.frame(result$edges))
  expect_true(is.character(result$nodes))
  expect_true(result$directed)
  expect_true(result$n_nodes > 0L)
  expect_true(result$n_edges > 0L)
})

test_that("build_hon nodes use pipe notation", {
  trajs <- list(c("A", "B", "C"))
  result <- build_hon(trajs, max_order = 1L, min_freq = 1L)
  # All nodes should contain "|"
  expect_true(all(grepl("|", result$nodes, fixed = TRUE)))
})

test_that("sequence_to_node produces correct notation", {
  expect_equal(.hon_sequence_to_node("A"), "A|")
  expect_equal(.hon_sequence_to_node(c("A", "B")), "B|A")
  expect_equal(.hon_sequence_to_node(c("X", "A", "B")), "B|A.X")
})

test_that("print and summary work without error", {
  result <- build_hon(.make_hon_data(), max_order = 2L, min_freq = 1L)
  expect_output(print(result), "Higher-Order Network")
  expect_output(summary(result), "Summary")
})

test_that("edge weights are probabilities (0, 1]", {
  result <- build_hon(.make_hon_data(), max_order = 2L, min_freq = 1L)
  expect_true(all(result$edges$weight > 0))
  expect_true(all(result$edges$weight <= 1))
})

test_that("encode/decode roundtrip preserves tuple", {
  tup <- c("alpha", "beta", "gamma")
  expect_equal(.hon_decode(.hon_encode(tup)), tup)
})

test_that("key_len returns correct lengths", {
  keys <- c(.hon_encode("A"), .hon_encode(c("A", "B")),
            .hon_encode(c("X", "Y", "Z")))
  expect_equal(.hon_key_len(keys), c(1L, 2L, 3L))
})

test_that("build_hon with max_order=1 produces only first-order nodes", {
  trajs <- list(c("A", "B", "C", "D"))
  result <- build_hon(trajs, max_order = 1L, min_freq = 1L)
  # All nodes should have order 1 (format "X|")
  node_orders <- vapply(result$nodes, function(nd) {
    parts <- strsplit(nd, "|", fixed = TRUE)[[1L]]
    if (length(parts) < 2L || parts[2L] == "") 1L
    else length(strsplit(parts[2L], ".", fixed = TRUE)[[1L]]) + 1L
  }, integer(1L))
  expect_true(all(node_orders == 1L))
})

test_that("adjacency matrix dimensions match n_nodes", {
  result <- build_hon(.make_hon_data(), max_order = 2L, min_freq = 1L)
  expect_equal(nrow(result$matrix), result$n_nodes)
  expect_equal(ncol(result$matrix), result$n_nodes)
})
