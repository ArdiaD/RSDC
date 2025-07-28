test_that("f_simul_tvtp returns expected structure and dimensions", {
  set.seed(42)

  n <- 100
  K <- 3
  N <- 2
  p <- 2

  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  beta <- array(rnorm(N * N * p, mean = 0, sd = 0.2), dim = c(N, N, p))
  mu <- matrix(rnorm(N * K), nrow = N)
  sigma <- array(NA, dim = c(K, K, N))
  for (j in 1:N) {
    A <- matrix(rnorm(K * K), K)
    sigma[,,j] <- crossprod(A) + diag(0.5, K)  # ensure positive definite
  }

  result <- RSDC::f_simul_tvtp(n, X, beta, mu, sigma, N)

  expect_type(result, "list")
  expect_named(result, c("states", "observations", "transition_matrices"))

  expect_true(is.numeric(result$states))
  expect_equal(length(result$states), n)

  expect_true(is.matrix(result$observations))
  expect_equal(dim(result$observations), c(n, K))

  expect_true(is.array(result$transition_matrices))
  expect_equal(dim(result$transition_matrices), c(N, N, n))

  # Check each transition matrix is stochastic
  for (t in 2:n) {
    row_sums <- rowSums(result$transition_matrices[,,t])
    expect_true(all(abs(row_sums - 1) < 1e-6))
  }
})

test_that("f_simul_tvtp gives reproducible results with seed", {
  n <- 50
  K <- 2
  N <- 2
  p <- 1

  X <- matrix(rnorm(n * p), nrow = n)
  beta <- array(rnorm(N * N * p), dim = c(N, N, p))
  mu <- matrix(rnorm(N * K), nrow = N)
  sigma <- array(diag(K), dim = c(K, K, N))

  res1 <- RSDC::f_simul_tvtp(n, X, beta, mu, sigma, N, seed = 999)
  res2 <- RSDC::f_simul_tvtp(n, X, beta, mu, sigma, N, seed = 999)

  expect_equal(res1$states, res2$states)
  expect_equal(res1$observations, res2$observations)
  expect_equal(res1$transition_matrices, res2$transition_matrices)
})

