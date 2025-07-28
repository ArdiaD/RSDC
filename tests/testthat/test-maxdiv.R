test_that("f_maxdiv returns valid output on synthetic data", {
  skip_if_not_installed("Rsolnp")

  set.seed(42)
  T <- 40
  K <- 3
  P <- choose(K, 2)

  y <- matrix(rnorm(T * K), ncol = K)
  colnames(y) <- paste0("Asset", 1:K)

  sigma_matrix <- matrix(1, T, K)
  colnames(sigma_matrix) <- colnames(y)

  predicted_corr <- matrix(0.3, T, P)

  out <- f_maxdiv(
    sigma_matrix = sigma_matrix,
    value_cols = 1:K,
    predicted_corr = predicted_corr,
    y = y,
    long_only = TRUE,
    rebalance = "monthly"
  )

  expect_s3_class(out, "maxdiv_portfolio")
  expect_equal(dim(out$weights), c(T, K))
  expect_length(out$returns, T)
  expect_length(out$diversification_ratios, T)
  expect_true(is.numeric(out$mean_diversification))
  expect_equal(out$K, K)
  expect_equal(out$assets, 1:K)
  expect_true(is.numeric(out$volatility))
})
