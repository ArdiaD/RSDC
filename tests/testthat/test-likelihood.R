test_that("rsdc_likelihood matches Hamilton (no X)", {
  y <- toy_residuals(50, 3); K <- ncol(y); N <- 2
  rho <- toy_rho_matrix(K)
  P <- toy_P()

  # pack params: trans (2) + rho (N * K(K-1)/2 = 2*3=6)
  trans <- c(P[1,1], P[2,2])
  params <- c(trans, as.vector(t(rho)))
  nll <- rsdc_likelihood(params, y = y, exog = NULL, K = K, N = N)

  ref <- rsdc_hamilton(y = y, X = NULL, beta = NULL,
                       rho_matrix = rho, K = K, N = N, P = P)
  expect_equal(nll, -ref$log_likelihood, tolerance = 1e-1)
})

test_that("rsdc_likelihood matches Hamilton (TVTP)", {
  y <- toy_residuals(50, 3); K <- ncol(y); N <- 2
  rho <- toy_rho_matrix(K)
  X <- cbind(1, scale(seq_len(nrow(y))))
  p <- ncol(X)
  beta <- rbind(c(1.5, 0.1), c(0.8, -0.2))

  params <- c(as.vector(t(beta)), as.vector(t(rho)))
  nll <- rsdc_likelihood(params, y = y, exog = X, K = K, N = N)

  ref <- rsdc_hamilton(y = y, X = X, beta = beta,
                       rho_matrix = rho, K = K, N = N, P = NULL)
  expect_equal(nll, -ref$log_likelihood, tolerance = 1e-1)
})
