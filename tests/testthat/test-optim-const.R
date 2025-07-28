test_that("f_optim_const returns valid structure", {
  set.seed(123)
  T <- 50
  K <- 3
  residuals <- matrix(rnorm(T * K), ncol = K)

  res <- f_optim_const(residuals = residuals, control = list(do_trace = FALSE))

  expect_type(res, "list")
  expect_named(res, c("transition_matrix", "means", "volatilities", "correlations", "covariances", "log_likelihood", "beta"))
  expect_equal(dim(res$transition_matrix), c(1, 1))
  expect_equal(dim(res$means), c(1, K))
  expect_equal(dim(res$volatilities), c(1, K))
  expect_equal(dim(res$correlations), c(1, K * (K - 1) / 2))
  expect_equal(dim(res$covariances), c(K, K, 1))
  expect_null(res$beta)
  expect_true(is.finite(res$log_likelihood))
})
