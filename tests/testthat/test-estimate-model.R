test_that("estimate_model wrapper dispatches correctly", {
  set.seed(123)
  T <- 50
  K <- 2
  N <- 2
  residuals <- matrix(rnorm(T * K), ncol = K)
  X <- cbind(1, matrix(rnorm(T), ncol = 1))

  # TVTP
  res1 <- estimate_model("tvtp", residuals = residuals, N = N, X = X)
  expect_equal(dim(res1$transition_matrix), c(N, N))
  expect_equal(dim(res1$beta), c(N, ncol(X)))

  # noX
  res2 <- estimate_model("noX", residuals = residuals, N = N)
  expect_equal(dim(res2$transition_matrix), c(N, N))
  expect_null(res2$beta)

  # const
  res3 <- estimate_model("const", residuals = residuals)
  expect_equal(dim(res3$transition_matrix), c(1, 1))
  expect_null(res3$beta)
})
