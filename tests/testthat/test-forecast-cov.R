test_that("f_forecast_cov returns correct output structure (out-of-sample = TRUE)", {
  set.seed(123)

  # Simulate simple data
  T <- 50
  K <- 2
  N <- 2
  p <- 2

  residuals <- matrix(rnorm(T * K), ncol = K)
  colnames(residuals) <- paste0("Asset", 1:K)

  X <- matrix(rnorm(T * p), ncol = p)
  sigma_matrix <- matrix(runif(T * K, 0.8, 1.2), ncol = K)
  colnames(sigma_matrix) <- colnames(residuals)
  value_cols <- 1:K

  # Estimate model
  fit <- RSDC::estimate_model(
    method = "tvtp",
    residuals = residuals,
    X = X,
    N = N,
    out_of_sample = TRUE,
    control = list(do_trace = FALSE)
  )

  # Forecast
  forecast <- RSDC::f_forecast_cov(
    method = "tvtp",
    N = N,
    residuals = residuals,
    X = X,
    final_params = fit,
    sigma_matrix = sigma_matrix,
    value_cols = value_cols,
    out_of_sample = TRUE
  )

  # Check structure of forecast output
  expect_type(forecast, "list")
  expect_true(all(c("cov_matrices", "predicted_correlations", "BIC", "y", "sigma_matrix") %in% names(forecast)))
  expect_equal(ncol(forecast$predicted_correlations), choose(K, 2))
  expect_equal(length(forecast$cov_matrices), nrow(forecast$sigma_matrix))
})

