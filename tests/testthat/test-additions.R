# Features added in 1.5-0: multi-start, Viterbi, multi-step forecast,
# correlation bands, broom tidiers, and the ggplot2 autoplot method.

test_that("multi-start keeps the best fit and records the start log-likelihoods", {
  skip_on_cran()
  set.seed(1); y <- scale(matrix(rnorm(300 * 2), 300, 2))
  f <- rsdc_estimate("noX", residuals = y, N = 2,
                     control = list(n_starts = 3L, itermax = 40))
  expect_s3_class(f, "rsdc_fit")
  expect_length(f$start_logliks, 3L)
  expect_equal(f$log_likelihood, max(f$start_logliks), tolerance = 1e-6)
})

test_that("rsdc_viterbi returns a valid most-likely regime path", {
  skip_on_cran()
  set.seed(2); y <- scale(matrix(rnorm(200 * 2), 200, 2))
  f <- rsdc_estimate("noX", residuals = y, N = 2, control = list(itermax = 40))
  v <- rsdc_viterbi(f)
  expect_length(v, 200L)
  expect_true(all(v %in% 1:2))
})

test_that("rsdc_forecast_ahead returns h-step regime and correlation forecasts", {
  skip_on_cran()
  set.seed(3); y <- scale(matrix(rnorm(200 * 2), 200, 2))
  f <- rsdc_estimate("noX", residuals = y, N = 2, control = list(itermax = 40))
  fa <- rsdc_forecast_ahead(f, horizon = 5L)
  expect_equal(dim(fa$regime_probs), c(5L, 2L))
  expect_equal(dim(fa$predicted_correlations), c(5L, 1L))
  expect_true(all(abs(rowSums(fa$regime_probs) - 1) < 1e-8))
})

test_that("rsdc_corr_bands returns ordered pointwise bands", {
  skip_on_cran()
  set.seed(4); y <- scale(matrix(rnorm(250 * 2), 250, 2))
  f <- rsdc_estimate("noX", residuals = y, N = 2,
                     control = list(itermax = 40, compute_se = TRUE))
  b <- rsdc_corr_bands(f, B = 50L, seed = 1)
  expect_length(b, 1L)                       # one asset pair for K = 2
  m <- b[[1]]
  expect_equal(colnames(m), c("fit", "lower", "upper"))
  expect_equal(nrow(m), nrow(y))
  expect_true(all(m[, "lower"] <= m[, "upper"]))
})

test_that("broom tidiers (tidy/glance/augment) work", {
  skip_on_cran()
  set.seed(5); y <- scale(matrix(rnorm(200 * 2), 200, 2))
  f <- rsdc_estimate("noX", residuals = y, N = 2, control = list(itermax = 40))
  td <- generics::tidy(f)
  expect_true(all(c("term", "estimate", "std.error", "statistic", "p.value") %in% names(td)))
  gl <- generics::glance(f)
  expect_equal(nrow(gl), 1L)
  expect_true(all(c("logLik", "AIC", "BIC", "nobs") %in% names(gl)))
  au <- generics::augment(f)
  expect_equal(nrow(au), 200L)
  expect_true(".state" %in% names(au))
})

test_that("autoplot returns a ggplot for switching models", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  set.seed(6); y <- scale(matrix(rnorm(200 * 2), 200, 2))
  f <- rsdc_estimate("noX", residuals = y, N = 2, control = list(itermax = 40))
  p <- ggplot2::autoplot(f)
  expect_s3_class(p, "ggplot")
})
