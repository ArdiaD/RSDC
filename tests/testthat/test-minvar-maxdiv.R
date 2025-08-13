test_that("rsdc_minvar produces valid weights and volatility", {
  skip_if_not_installed("quadprog")

  T <- 20; K <- 3
  y <- toy_residuals(T, K)
  S <- toy_sigma_matrix(T, K)
  pairs <- toy_cor_pairs(K)

  # simple time-varying pairwise corr: small oscillations
  P <- length(pairs)
  pred_corr <- matrix(0, T, P)
  for (t in 1:T) pred_corr[t, ] <- c(0.2, 0.1, 0.05) + 0.01*sin(t/3)

  res <- rsdc_minvar(sigma_matrix = S,
                     value_cols = colnames(S),
                     predicted_corr = pred_corr,
                     y = y,
                     long_only = TRUE)

  expect_equal(dim(res$weights), c(T, K))
  # weights sum to ~1
  expect_true(all(abs(rowSums(res$weights) - 1) < 1e-6))
  # long-only honored
  expect_true(all(res$weights >= -1e-10))
  expect_true(is.finite(res$volatility))
})

test_that("rsdc_maxdiv produces valid outputs", {
  skip_if_not_installed("Rsolnp")

  T <- 20; K <- 3
  y <- toy_residuals(T, K)
  S <- toy_sigma_matrix(T, K)
  pairs <- toy_cor_pairs(K); P <- length(pairs)
  pred_corr <- matrix(rep(c(0.15, 0.10, 0.05), each = T), ncol = P)

  res <- rsdc_maxdiv(sigma_matrix = S,
                     value_cols = colnames(S),
                     predicted_corr = pred_corr,
                     y = y,
                     long_only = TRUE)

  expect_equal(dim(res$weights), c(T, K))
  expect_equal(length(res$returns), T)
  expect_equal(length(res$diversification_ratios), T)
  expect_true(all(abs(rowSums(res$weights) - 1) < 1e-6))
  expect_true(all(res$weights >= -1e-10))
  expect_true(is.finite(res$mean_diversification))
  expect_true(is.finite(res$volatility))
})

