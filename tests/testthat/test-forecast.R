test_that("rsdc_forecast returns expected shapes (const)", {
  K <- 3; T <- 25
  y <- toy_residuals(T, K)
  S <- toy_sigma_matrix(T, K)

  out <- rsdc_forecast(method = "const",
                       N = 1,
                       residuals = y,
                       X = NULL,
                       final_params = list(
                         correlations = matrix(c(0.2, 0.1, 0.05), nrow = 1),
                         log_likelihood = -100.0
                       ),
                       sigma_matrix = S,
                       value_cols = colnames(S),
                       out_of_sample = FALSE)

  expect_equal(dim(out$predicted_correlations), c(T, K*(K-1)/2))
  expect_equal(length(out$cov_matrices), T)
  expect_true(is.finite(out$BIC))
})

test_that("rsdc_forecast returns expected shapes (noX and tvtp)", {
  K <- 3; T <- 25; N <- 2
  y <- toy_residuals(T, K)
  S <- toy_sigma_matrix(T, K)
  rho <- toy_rho_matrix(K)

  # noX branch (uses provided P through final_params)
  P <- toy_P()
  out_noX <- rsdc_forecast(method = "noX",
                           N = N,
                           residuals = y,
                           X = NULL,
                           final_params = list(
                             correlations = rho,
                             transition_matrix = P,
                             log_likelihood = -200
                           ),
                           sigma_matrix = S,
                           value_cols = colnames(S),
                           out_of_sample = FALSE)

  expect_equal(dim(out_noX$predicted_correlations), c(T, K*(K-1)/2))
  expect_equal(dim(out_noX$smoothed_probs), c(N, T))
  expect_true(is.finite(out_noX$BIC))

  # tvtp branch (uses X and beta)
  X <- cbind(1, scale(seq_len(T)))
  beta <- rbind(c(1.2, 0.0), c(0.8, -0.1))

  out_tvtp <- rsdc_forecast(method = "tvtp",
                            N = N,
                            residuals = y,
                            X = X,
                            final_params = list(
                              correlations = rho,
                              beta = beta,
                              log_likelihood = -210
                            ),
                            sigma_matrix = S,
                            value_cols = colnames(S),
                            out_of_sample = FALSE)

  expect_equal(dim(out_tvtp$predicted_correlations), c(T, K*(K-1)/2))
  expect_equal(dim(out_tvtp$smoothed_probs), c(N, T))
  expect_true(is.finite(out_tvtp$BIC))
})


