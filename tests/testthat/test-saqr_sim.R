# ===========================================================================
# saqr_sim S3 class
# ===========================================================================

test_that("saqr_sim: constructor builds valid object", {
  obj <- saqr_sim(
    data   = data.frame(x = 1:5, y = rnorm(5)),
    params = list(mu = 0, sd = 1),
    type   = "test",
    seed   = 42L
  )
  expect_s3_class(obj, "saqr_sim")
  expect_s3_class(obj, "list")
  expect_equal(obj$type, "test")
  expect_equal(obj$seed, 42L)
})

test_that("saqr_sim: $data and $params access works (backward compat)", {
  df <- data.frame(a = 1:3, b = 4:6)
  obj <- saqr_sim(data = df, params = list(k = 2), type = "test")
  expect_identical(obj$data, df)
  expect_identical(obj$params, list(k = 2))
})

test_that("saqr_sim: [ delegates to $data", {
  df <- data.frame(x = 1:10, y = 11:20)
  obj <- saqr_sim(data = df, params = list(), type = "test")
  expect_equal(obj[1:3, ], df[1:3, ])
  expect_equal(obj[, "x"], df[, "x"])
  expect_equal(obj[2, 1], df[2, 1])
})

test_that("saqr_sim: as.data.frame extracts $data", {
  df <- data.frame(v = rnorm(5))
  obj <- saqr_sim(data = df, params = list(), type = "test")
  expect_identical(as.data.frame(obj), df)
})

test_that("saqr_sim: dim/nrow/ncol delegate to $data", {
  df <- data.frame(a = 1:10, b = 11:20, c = 21:30)
  obj <- saqr_sim(data = df, params = list(), type = "test")
  expect_equal(dim(obj), c(10L, 3L))
  expect_equal(nrow(obj), 10L)
  expect_equal(ncol(obj), 3L)
})

test_that("saqr_sim: head/tail work", {
  df <- data.frame(x = 1:20)
  obj <- saqr_sim(data = df, params = list(), type = "test")
  expect_equal(nrow(head(obj)), 6L)
  expect_equal(nrow(head(obj, 3L)), 3L)
  expect_equal(nrow(tail(obj, 4L)), 4L)
})

test_that("saqr_sim: print runs without error", {
  obj <- saqr_sim(
    data   = data.frame(V1 = 1:3, V2 = 4:6),
    params = list(temporal = diag(2)),
    type   = "longitudinal",
    seed   = 1L
  )
  expect_output(print(obj), "saqr_sim")
  expect_output(print(obj), "longitudinal")
  expect_output(print(obj), "seed=1")
})

test_that("saqr_sim: names returns list element names for $ access", {
  obj <- saqr_sim(
    data   = data.frame(x = 1),
    params = list(a = 1),
    type   = "t",
    seed   = 1L
  )
  nms <- names(obj)
  expect_true(all(c("data", "params", "type", "seed") %in% nms))
})

test_that("saqr_sim: errors on invalid input", {
  expect_error(saqr_sim(data = "not a df", params = list(), type = "t"))
  expect_error(saqr_sim(data = data.frame(x = 1), params = "bad", type = "t"))
})


# ===========================================================================
# Backward compatibility: latent functions return saqr_sim
# ===========================================================================

test_that("simulate_lpa returns saqr_sim and backward compat works", {
  means <- matrix(c(0, 0, 10, 10), nrow = 2, ncol = 2)
  r <- simulate_lpa(means = means, sds = 0.5, props = c(0.5, 0.5), n = 50, seed = 1)

  # New: is saqr_sim

  expect_s3_class(r, "saqr_sim")
  expect_equal(r$type, "lpa")

  # Backward compat: $data, $params, names(r)
  expect_s3_class(r$data, "data.frame")
  expect_true("means" %in% names(r$params))
  expect_true(all(c("data", "params") %in% names(r)))

  # Backward compat: [ subscripting
  expect_equal(nrow(r[1:5, ]), 5L)

  # Backward compat: as.data.frame
  expect_identical(as.data.frame(r), r$data)
})

test_that("simulate_lca returns saqr_sim and backward compat works", {
  ip <- matrix(c(0.9, 0.1, 0.1, 0.9), nrow = 2, ncol = 2)
  r <- simulate_lca(item_probs = ip, class_probs = c(0.5, 0.5), n = 50, seed = 1)

  expect_s3_class(r, "saqr_sim")
  expect_equal(r$type, "lca")
  expect_s3_class(r$data, "data.frame")
  expect_true("item_probs" %in% names(r$params))
})

test_that("simulate_regression returns saqr_sim and backward compat works", {
  coefs <- c("(Intercept)" = 2, x1 = 3, x2 = -1)
  r <- simulate_regression(
    coefs = coefs, predictor_sds = c(x1 = 1, x2 = 1),
    error_sd = 0.5, n = 100, seed = 1
  )

  expect_s3_class(r, "saqr_sim")
  expect_equal(r$type, "regression")
  expect_s3_class(r$data, "data.frame")
  expect_identical(r$params$coefs, coefs)
})

test_that("simulate_fa returns saqr_sim and backward compat works", {
  loadings <- matrix(c(0.8, 0.7, 0.6, 0, 0, 0,
                        0, 0, 0, 0.8, 0.7, 0.6), nrow = 6, ncol = 2)
  r <- simulate_fa(loadings = loadings, n = 100, seed = 1)

  expect_s3_class(r, "saqr_sim")
  expect_equal(r$type, "fa")
  expect_s3_class(r$data, "data.frame")
  expect_true("sigma_implied" %in% names(r$params))
})

test_that("simulate_seq_clusters returns saqr_sim and backward compat works", {
  r <- simulate_seq_clusters(n = 50, n_clusters = 2, n_states = 3, seed = 1)

  expect_s3_class(r, "saqr_sim")
  expect_equal(r$type, "seq_clusters")
  expect_s3_class(r$data, "data.frame")
  expect_true("trans_list" %in% names(r$params))
})
