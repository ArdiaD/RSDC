test_that("rsdc_estimate(method='const') returns coherent structure", {
  skip_on_cran()
  y <- toy_residuals(T = 12, K = 2)   # small to keep it quick

  fit <- rsdc_estimate(method = "const", residuals = y)

  expect_true(is.list(fit))
  expect_equal(dim(fit$transition_matrix), c(1, 1))
  expect_equal(ncol(fit$correlations), choose(ncol(y), 2))  # for K=2 => 1
  expect_equal(dim(fit$covariances), c(2, 2, 1))
  expect_true(is.finite(fit$log_likelihood))
})

test_that("rsdc_estimate(method='noX') returns coherent structure", {
  skip_on_cran()
  y <- toy_residuals(T = 10, K = 2)

  fit <- rsdc_estimate(method = "noX", residuals = y, N = 2)

  expect_true(is.list(fit))
  expect_equal(dim(fit$transition_matrix), c(2, 2))
  expect_equal(dim(fit$correlations), c(2, choose(ncol(y), 2)))  # 2 x 1
  expect_equal(dim(fit$covariances), c(2, 2, 2))
  expect_true(is.finite(fit$log_likelihood))
})

test_that("rsdc_estimate(method='tvtp') returns coherent structure", {
  skip_on_cran()
  y <- toy_residuals(T = 10, K = 2)
  # X is not a helper; define inline
  X <- cbind(1, scale(seq_len(nrow(y))))

  fit <- rsdc_estimate(method = "tvtp", residuals = y, N = 2, X = X)

  expect_true(is.list(fit))
  expect_equal(dim(fit$transition_matrix), c(2, 2))
  expect_equal(dim(fit$correlations), c(2, choose(ncol(y), 2)))  # 2 x 1
  expect_equal(dim(fit$covariances), c(2, 2, 2))
  expect_equal(dim(fit$beta), c(2, ncol(X)))
  expect_true(is.finite(fit$log_likelihood))
})

test_that("rsdc_estimate errors when tvtp is called without X", {
  y <- toy_residuals(T = 8, K = 2)
  expect_error(
    rsdc_estimate(method = "tvtp", residuals = y, N = 2, X = NULL),
    regexp = "X must be provided for method = 'tvtp'"
  )
})

test_that("rsdc_estimate noX N=3: returned P and correlations are consistent with IS log-likelihood", {
  skip_on_cran()
  set.seed(42)
  y <- toy_residuals(T = 30, K = 2)
  K <- ncol(y); N <- 3

  fit <- rsdc_estimate("noX", residuals = y, N = N)

  P <- fit$transition_matrix
  # N=3: 2 free trans params per row (first two columns), then N*C rho params row-wise
  trans_params <- c(P[1, 1], P[1, 2],
                    P[2, 1], P[2, 2],
                    P[3, 1], P[3, 2])
  rho_vec <- as.vector(t(fit$correlations))
  params  <- c(trans_params, rho_vec)

  nll <- rsdc_likelihood(params, y = y, exog = NULL, K = K, N = N)
  expect_equal(-nll, fit$log_likelihood, tolerance = 1e-6)
})

