test_that("f_hamilton returns correct output structure and dimensions", {
  set.seed(42)

  T <- 50
  K <- 3
  N <- 2
  p <- 2

  # Simulated inputs
  y <- matrix(rnorm(T * K), ncol = K)
  X <- matrix(rnorm(T * p), ncol = p)
  beta <- matrix(rnorm(N * p), nrow = N)
  rho_matrix <- matrix(runif(N * choose(K, 2), -0.5, 0.5), nrow = N)

  # Ensure valid correlation matrices
  for (i in 1:N) {
    temp <- diag(K)
    temp[lower.tri(temp)] <- rho_matrix[i, ]
    temp[upper.tri(temp)] <- t(temp)[upper.tri(temp)]
    if (min(eigen(temp)$values) < 0.1) {
      rho_matrix[i, ] <- rho_matrix[i, ] * 0.5  # shrink values if matrix nearly singular
    }
  }

  result <- RSDC::f_hamilton(
    y = y,
    X = X,
    beta = beta,
    rho_matrix = rho_matrix,
    K = K,
    N = N
  )

  expect_type(result, "list")
  expect_named(result, c("filtered_probs", "smoothed_probs", "log_likelihood"))
  expect_equal(dim(result$filtered_probs), c(N, T))
  expect_equal(dim(result$smoothed_probs), c(N, T))
  expect_type(result$log_likelihood, "double")
  expect_true(is.finite(result$log_likelihood))
})

