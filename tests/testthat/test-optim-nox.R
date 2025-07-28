test_that("f_optim_noX returns valid structure", {
  set.seed(123)
  T <- 50
  K <- 3
  N <- 2
  residuals <- matrix(rnorm(T * K), ncol = K)

  res <- f_optim_noX(N = N, residuals = residuals, out_of_sample = FALSE, control = list(do_trace = FALSE))

  expect_type(res, "list")
  expect_named(res, c("transition_matrix", "means", "volatilities", "correlations", "covariances", "log_likelihood", "beta"))
  expect_equal(dim(res$transition_matrix), c(N, N))
  expect_equal(dim(res$means), c(N, K))
  expect_equal(dim(res$volatilities), c(N, K))
  expect_equal(dim(res$correlations), c(N, K * (K - 1) / 2))
  expect_equal(dim(res$covariances), c(K, K, N))
  expect_null(res$beta)
  expect_true(is.finite(res$log_likelihood))
})
