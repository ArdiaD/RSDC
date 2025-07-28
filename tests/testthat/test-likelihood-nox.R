test_that("f_likelihood returns finite numeric value for fixed transition model", {
  set.seed(456)

  T <- 50
  K <- 2
  N <- 2

  y <- matrix(rnorm(T * K), ncol = K)
  n_pairs <- K * (K - 1) / 2
  n_trans <- N * (N - 1)
  n_rho <- N * n_pairs

  trans_params <- rep(0.8, n_trans)
  rho_params <- runif(n_rho, min = -0.3, max = 0.3)

  params <- c(trans_params, rho_params)

  val <- RSDC::f_likelihood(params = params, y = y, exog = NULL, K = K, N = N)

  expect_type(val, "double")
  expect_length(val, 1)
  expect_true(is.finite(val))
  expect_true(val > 0)
})
