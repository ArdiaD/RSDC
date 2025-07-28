test_that("f_minvar returns valid output for synthetic data", {
  skip_if_not_installed("quadprog")

  set.seed(123)
  T <- 50
  K <- 3
  P <- choose(K, 2)

  y <- matrix(rnorm(T * K, 0, 1), ncol = K)
  colnames(y) <- paste0("Asset", 1:K)

  sigma_matrix <- matrix(1, T, K)
  colnames(sigma_matrix) <- paste0("Asset", 1:K)

  predicted_corr <- matrix(0.2, T, P)

  out <- f_minvar(
    sigma_matrix = sigma_matrix,
    value_cols = 1:K,
    predicted_corr = predicted_corr,
    y = y,
    rebalance = "daily",
    long_only = TRUE
  )

  expect_s3_class(out, "minvar_portfolio")
  expect_equal(dim(out$weights), c(T, K))
  expect_type(out$cov_matrices, "list")
  expect_length(out$cov_matrices, T)
  expect_true(is.numeric(out$volatility))
  expect_equal(dim(out$y), c(T, K))
  expect_equal(out$K, K)
})
