test_that("f_optim returns valid structure with TVTP", {
  set.seed(123)
  T <- 50
  K <- 3
  N <- 2
  p <- 2
  residuals <- matrix(rnorm(T * K), ncol = K)
  X <- cbind(1, matrix(rnorm(T * (p - 1)), ncol = p - 1))  # include intercept

  res <- f_optim(N = N, residuals = residuals, X = X, out_of_sample = FALSE, control = list(do_trace = FALSE))

  expect_type(res, "list")
  expect_named(res, c("transition_matrix", "means", "volatilities", "correlations", "covariances", "log_likelihood", "beta"))
  expect_equal(dim(res$transition_matrix), c(N, N))
  expect_equal(dim(res$means), c(N, K))
  expect_equal(dim(res$volatilities), c(N, K))
  expect_equal(dim(res$correlations), c(N, K * (K - 1) / 2))
  expect_equal(dim(res$covariances), c(K, K, N))
  expect_equal(dim(res$beta), c(N, p))
  expect_true(is.finite(res$log_likelihood))
})
