test_that("f_likelihood returns finite numeric value for TVTP model", {
  set.seed(123)

  T <- 50
  K <- 3
  N <- 2
  p <- 2

  y <- matrix(rnorm(T * K), ncol = K)
  exog <- matrix(rnorm(T * p), ncol = p)

  n_pairs <- K * (K - 1) / 2
  n_beta <- N * p
  n_rho <- N * n_pairs

  beta_params <- rnorm(n_beta, mean = 0, sd = 0.1)
  rho_params <- runif(n_rho, min = -0.3, max = 0.3)

  params <- c(beta_params, rho_params)

  val <- RSDC::f_likelihood(params = params, y = y, exog = exog, K = K, N = N)

  expect_type(val, "double")
  expect_length(val, 1)
  expect_true(is.finite(val))
  expect_true(val > 0)
})

