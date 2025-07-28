#' Estimate Regime-Switching Correlation Model with TVTP
#'
#' Fits a multivariate regime-switching model with time-varying transition probabilities (TVTP),
#' using differential evolution followed by L-BFGS-B optimization.
#'
#' @param N Integer. Number of regimes.
#' @param residuals A numeric matrix of dimension T × K. Typically standardized residuals or returns.
#' @param X A numeric matrix of dimension T × p. Exogenous variables used to drive the transition probabilities.
#' Must include a column of ones if an intercept is desired (not added automatically).
#' @param out_of_sample Logical. If \code{TRUE}, the first 70\% of the data is used for estimation,
#' and the remaining 30\% is reserved for out-of-sample forecasting.
#' @param control Optional list of control arguments (e.g., \code{do_trace = TRUE}) passed to internal DEoptim control.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{transition_matrix}{The average regime transition matrix evaluated at the mean of \code{X}.}
#'   \item{means}{A placeholder matrix of zeros (N × K) for regime-specific means (not estimated).}
#'   \item{volatilities}{A placeholder matrix of ones (N × K) for regime-specific volatilities (not estimated).}
#'   \item{correlations}{A matrix of regime-specific correlation coefficients (N × C),
#'   where C = K(K−1)/2 is the number of unique pairwise correlations.}
#'   \item{covariances}{An array of full correlation matrices for each regime (K × K × N).}
#'   \item{log_likelihood}{The log-likelihood of the fitted model.}
#'   \item{beta}{Matrix of estimated logistic regression coefficients for the TVTP specification (N × p).}
#' }
#'
#' @importFrom DEoptim DEoptim
#' @export
f_optim <- function(N, residuals, X, out_of_sample = FALSE, control = list()) {

  require(DEoptim)

  con <- list(do_trace = FALSE)
  con[names(control)] <- control

  # Parameter setup
  K <- ncol(residuals)
  p <- ncol(X)
  n_beta <- N * p
  n_rho <- N * (K * (K - 1) / 2)
  rho_indices <- (n_beta + 1):(n_beta + n_rho)

  # Bounds with numerical stability
  bounds <- list(
    lower = c(rep(-10, n_beta), rep(-1, n_rho)),
    upper = c(rep(10, n_beta), rep(1, n_rho))
  )

  # Out-of-sample handling (fixed indexing)
  if (out_of_sample) {
    is_cut <- round(0.7 * nrow(residuals))
    residuals_is <- residuals[1:is_cut, ]
    residuals_oos <- residuals[(is_cut + 1):nrow(residuals), ]
    X_is <- X[1:is_cut, ]
    X_oos <- X[(is_cut + 1):nrow(X), ]
    y_optim <- residuals_is
    X_optim <- X_is
  } else {
    y_optim <- residuals
    X_optim <- X
  }

  # DEoptim
  set.seed(123)
  results <- DEoptim(
    fn = f_likelihood,
    lower = bounds$lower,
    upper = bounds$upper,
    control = list(
      itermax = 500,
      NP = 10 * length(bounds$lower),
      strategy = 2,
      F = 0.8,
      CR = 0.9,
      trace = con$do_trace,
      parallelType = 1,
      parVar = c("f_hamilton", "f_likelihood", "plogis"),
      packages = c("mvtnorm"),
      reltol = 1e-8,
      steptol = 50
    ),
    y = y_optim,
    exog = X_optim,
    K = K,
    N = N
  )

  # Optim refinement (fixed residual reference)
  result_optim <- optim(
    par = results$optim$bestmem,
    fn = f_likelihood,
    method = "L-BFGS-B",
    lower = bounds$lower,
    upper = bounds$upper,
    y = y_optim,
    exog = X_optim,
    K = K,
    N = N,
    control = list(maxit = 1000)
  )

  # Parameter processing
  optim_params <- result_optim$par
  C_per_state <- K * (K - 1) / 2

  # Generalized state ordering
  corr_means <- sapply(1:N, function(i) {
    mean(optim_params[rho_indices][(1:C_per_state) + (i - 1)*C_per_state])
  })
  sorted_states <- order(corr_means)

  if (!all(sorted_states == 1:N)) {
    # Reorder parameters
    optim_params <- c(
      # Betas
      unlist(lapply(sorted_states, function(s) {
        optim_params[((s - 1)*p + 1):(s*p)]
      })),
      # Rhos
      unlist(lapply(sorted_states, function(s) {
        optim_params[rho_indices][((s - 1)*C_per_state + 1):(s*C_per_state)]
      }))
    )
  }

  # Parameter extraction
  beta <- matrix(optim_params[1:n_beta], nrow = N, ncol = p)
  rho_matrix <- matrix(optim_params[rho_indices], nrow = N, byrow = TRUE)

  # Covariance matrix reconstruction
  sigma_array <- array(dim = c(K, K, N))
  for (i in 1:N) {
    cor_mat <- diag(K)
    cor_mat[lower.tri(cor_mat)] <- rho_matrix[i, ]
    cor_mat[upper.tri(cor_mat)] <- t(cor_mat)[upper.tri(cor_mat)]
    sigma_array[,,i] <- cor_mat
  }

  # Generalized transition matrix
  avg_X <- colMeans(X)
  P <- matrix(nrow = N, ncol = N)
  for (i in 1:N) {
    for (j in 1:N) {
      if(i == j) {
        P[i,j] <- plogis(avg_X %*% beta[i, ])
      } else {
        P[i,j] <- (1 - plogis(avg_X %*% beta[i, ])) / (N - 1)
      }
    }
  }

  list(
    transition_matrix = P,
    means = matrix(0, N, K),
    volatilities = matrix(1, N, K),
    correlations = rho_matrix,
    covariances = sigma_array,
    log_likelihood = -result_optim$value,
    beta = beta
  )
}

#' Estimate Regime-Switching Correlation Model with Constant Transition Matrix
#'
#' Fits a multivariate regime-switching correlation model with a fixed (time-invariant) transition matrix.
#' This variant does not use any exogenous covariates; the transition probabilities are constant and estimated
#' as part of the optimization process.
#'
#' The estimation is done via differential evolution followed by local refinement using L-BFGS-B.
#'
#' @param N Integer. Number of regimes.
#' @param residuals A numeric matrix of dimension T × K. Typically standardized residuals or returns.
#' @param out_of_sample Logical. If \code{TRUE}, the first 70\% of the data is used for estimation,
#' and the remaining 30\% is reserved for out-of-sample forecasting.
#' @param control Optional list of control arguments (e.g., \code{do_trace = TRUE}) passed to internal DEoptim control.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{transition_matrix}{Estimated fixed transition matrix of dimension N × N.}
#'   \item{means}{A placeholder matrix of zeros (N × K) for regime-specific means (not estimated).}
#'   \item{volatilities}{A placeholder matrix of ones (N × K) for regime-specific volatilities (not estimated).}
#'   \item{correlations}{A matrix of regime-specific pairwise correlations (N × C),
#'   where C = K(K−1)/2 is the number of unique pairwise correlations.}
#'   \item{covariances}{An array of full correlation matrices for each regime (K × K × N).}
#'   \item{log_likelihood}{The log-likelihood of the fitted model.}
#'   \item{beta}{\code{NULL} (no exogenous variables used in this specification).}
#' }
#'
#' @importFrom DEoptim DEoptim
#' @export
f_optim_noX <- function(N, residuals, out_of_sample = FALSE, control = list()) {
  stopifnot(is.matrix(residuals))

  con <- list(do_trace = FALSE)
  con[names(control)] <- control

  K <- ncol(residuals)
  T <- nrow(residuals)
  n_trans <- N * (N - 1)
  C_per_state <- K * (K - 1) / 2
  n_rho <- N * C_per_state

  bounds <- list(
    lower = c(rep(0.01, n_trans), rep(-1, n_rho)),
    upper = c(rep(0.99, n_trans), rep(1, n_rho))
  )

  if (out_of_sample) {
    cut <- round(0.7 * T)
    y_optim <- residuals[1:cut, ]
  } else {
    y_optim <- residuals
  }

  set.seed(123)
  de_result <- DEoptim::DEoptim(
    fn = f_likelihood,
    lower = bounds$lower,
    upper = bounds$upper,
    control = list(itermax = 500, NP = 10 * length(bounds$lower),
                   strategy = 2, F = 0.8, CR = 0.9,
                   trace = con$do_trace, reltol = 1e-8, steptol = 50),
    y = y_optim, exog = NULL, K = K, N = N
  )

  optim_result <- optim(
    par = de_result$optim$bestmem,
    fn = f_likelihood,
    method = "L-BFGS-B",
    lower = bounds$lower,
    upper = bounds$upper,
    y = y_optim, exog = NULL, K = K, N = N,
    control = list(maxit = 1000)
  )

  final_par <- optim_result$par
  rho_matrix <- matrix(final_par[(n_trans + 1):(n_trans + n_rho)], nrow = N, byrow = TRUE)

  sigma_array <- array(dim = c(K, K, N))
  for (i in 1:N) {
    R <- diag(K)
    R[lower.tri(R)] <- rho_matrix[i, ]
    R[upper.tri(R)] <- t(R)[upper.tri(R)]
    sigma_array[,,i] <- R
  }

  # Build transition matrix
  trans_params <- final_par[1:n_trans]
  P <- matrix(0, N, N)
  if (N == 2) {
    P[1, 1] <- trans_params[1]
    P[1, 2] <- 1 - trans_params[1]
    P[2, 1] <- 1 - trans_params[2]
    P[2, 2] <- trans_params[2]
  }

  list(
    transition_matrix = P,
    means = matrix(0, N, K),
    volatilities = matrix(1, N, K),
    correlations = rho_matrix,
    covariances = sigma_array,
    log_likelihood = -optim_result$value,
    beta = NULL
  )
}

#' Estimate Constant Correlation Model
#'
#' Estimates a single correlation matrix under the assumption that asset correlations are constant over time.
#' This model does not include regime switching or exogenous covariates.
#'
#' Estimation is performed by maximizing the log-likelihood of a multivariate normal distribution
#' with constant correlation, using differential evolution followed by L-BFGS-B optimization.
#'
#' @param residuals A numeric matrix of dimension T × K. Typically standardized residuals or returns.
#' @param control Optional list of control parameters passed to the DEoptim optimizer (e.g., \code{do_trace = TRUE}).
#'
#' @return A list with the following components:
#' \describe{
#'   \item{transition_matrix}{A 1 × 1 matrix with value 1 (since there is no switching).}
#'   \item{means}{A placeholder matrix of zeros (1 × K) for the regime mean (not estimated).}
#'   \item{volatilities}{A placeholder matrix of ones (1 × K) for the regime volatility (not estimated).}
#'   \item{correlations}{A 1 × C matrix of estimated correlation coefficients, where C = K(K−1)/2.}
#'   \item{covariances}{A 3D array of dimension K × K × 1 containing the full correlation matrix.}
#'   \item{log_likelihood}{The log-likelihood value of the fitted constant-correlation model.}
#'   \item{beta}{\code{NULL} (no covariates used).}
#' }
#'
#' @importFrom DEoptim DEoptim
#' @export
f_optim_const <- function(residuals, control = list()) {
  stopifnot(is.matrix(residuals))

  con <- list(do_trace = FALSE)
  con[names(control)] <- control

  K <- ncol(residuals)
  T <- nrow(residuals)
  n_rho <- K * (K - 1) / 2

  bounds <- list(
    lower = rep(-1, n_rho),
    upper = rep(1, n_rho)
  )

  neg_loglik_const <- function(rho_vec, y, K) {
    R <- diag(K)
    R[lower.tri(R)] <- rho_vec
    R[upper.tri(R)] <- t(R)[upper.tri(R)]

    if (min(eigen(R)$values) < 1e-6) return(1e10)

    inv_R <- solve(R + diag(1e-8, K))
    log_det <- determinant(R, logarithm = TRUE)$modulus[1]
    quad_terms <- rowSums((y %*% inv_R) * y)
    ll <- -0.5 * sum(log_det + quad_terms + K * log(2 * pi))
    return(-ll)
  }

  set.seed(123)
  de_result <- DEoptim::DEoptim(
    fn = neg_loglik_const,
    lower = bounds$lower,
    upper = bounds$upper,
    control = list(itermax = 500, trace = con$do_trace),
    y = residuals, K = K
  )

  optim_result <- optim(
    par = de_result$optim$bestmem,
    fn = neg_loglik_const,
    method = "L-BFGS-B",
    lower = bounds$lower,
    upper = bounds$upper,
    y = residuals, K = K,
    control = list(maxit = 1000)
  )

  rho_vec <- optim_result$par
  R <- diag(K)
  R[lower.tri(R)] <- rho_vec
  R[upper.tri(R)] <- t(R)[upper.tri(R)]

  array_R <- array(R, dim = c(K, K, 1))

  list(
    transition_matrix = matrix(1, 1, 1),
    means = matrix(0, 1, K),
    volatilities = matrix(1, 1, K),
    correlations = matrix(rho_vec, nrow = 1),
    covariances = array_R,
    log_likelihood = -optim_result$value,
    beta = NULL
  )
}

#' Estimate Regime-Switching or Constant Correlation Model
#'
#' A unified wrapper function to estimate correlation models with or without regime switching,
#' and with or without exogenous transition variables. This function delegates to one of:
#' \itemize{
#'   \item \code{f_optim()} for regime switching with exogenous covariates (\code{method = "tvtp"})
#'   \item \code{f_optim_noX()} for regime switching with fixed transition probabilities (\code{method = "noX"})
#'   \item \code{f_optim_const()} for a constant correlation model with no regime switching (\code{method = "const"})
#' }
#'
#' @param method Character string. One of:
#'   \itemize{
#'     \item \code{"tvtp"}: Regime switching with time-varying transition probabilities driven by \code{X}.
#'     \item \code{"noX"}: Regime switching with constant transition probabilities (no covariates).
#'     \item \code{"const"}: Constant correlation model (no regime switching or covariates).
#'   }
#' @param residuals A numeric matrix of dimension T × K. Typically standardized residuals or returns.
#' @param N Integer. Number of regimes. Ignored if \code{method = "const"}.
#' @param X A numeric matrix of exogenous covariates (T × p). Required only for \code{method = "tvtp"}.
#' @param out_of_sample Logical. If \code{TRUE}, uses 70\% of the data for estimation and 30\% for out-of-sample analysis.
#' @param control Optional list of control options. Currently supports \code{do_trace = TRUE} for verbose output during optimization.
#'
#' @return A list containing the estimated model parameters. See the return values of
#' \code{\link{f_optim}}, \code{\link{f_optim_noX}}, and \code{\link{f_optim_const}} for details.
#'
#' @export
estimate_model <- function(method = c("tvtp", "noX", "const"), residuals, N = 2, X = NULL, out_of_sample = FALSE, control = list()) {
  method <- match.arg(method)

  con <- list(do_trace = FALSE)
  con[names(control)] <- control

  if (method == "tvtp") {
    if (is.null(X)) stop("X must be provided for method = 'tvtp'")
    return(f_optim(N = N, residuals = residuals, X = X, out_of_sample = out_of_sample, control = list(do_trace = con$do_trace)))
  }

  if (method == "noX") {
    return(f_optim_noX(N = N, residuals = residuals, out_of_sample = out_of_sample, control = list(do_trace = con$do_trace)))
  }

  if (method == "const") {
    return(f_optim_const(residuals = residuals, control = list(do_trace = con$do_trace)))
  }

  stop("Unknown method: ", method)
}




