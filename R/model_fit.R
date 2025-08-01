#' Hamilton Filter with Time-Varying or Fixed Transition Probabilities
#'
#' Implements the Hamilton filter for a multivariate time series under a regime-switching correlation model.
#' Supports both fixed transition probabilities and time-varying transition probabilities (TVTP)
#' modeled via logistic functions of exogenous covariates.
#'
#' The function computes the filtered and smoothed probabilities of the latent regime sequence, and
#' evaluates the log-likelihood under the specified model.
#'
#' @param y A numeric matrix of dimension T × K containing the observed time series (e.g., standardized residuals).
#' @param X An optional numeric matrix of dimension T × p of exogenous covariates used for TVTP.
#' Required if \code{beta} is provided.
#' @param beta An optional numeric matrix of dimension N × p. Logistic regression coefficients governing
#' the regime-specific persistence probabilities in the TVTP specification.
#' @param rho_matrix A numeric matrix of dimension N × C, where C = K(K−1)/2.
#' Each row contains the lower-triangular correlation parameters for one regime.
#' @param K Integer. Number of observed variables.
#' @param N Integer. Number of regimes.
#' @param P An optional N × N fixed transition matrix (used when \code{X} and \code{beta} are \code{NULL}).
#'
#' @return A list with components:
#' \describe{
#'   \item{filtered_probs}{A numeric matrix of dimension N × T.
#'   Contains the filtered regime probabilities \eqn{\Pr(S_t = j \mid \Omega_t)}.}
#'   \item{smoothed_probs}{A numeric matrix of dimension N × T.
#'   Contains the smoothed regime probabilities \eqn{\Pr(S_t = j \mid \Omega_T)}.}
#'   \item{log_likelihood}{The scalar log-likelihood value of the model given the data.}
#' }
#'
#' @details
#' When \code{X} and \code{beta} are both provided, the function builds a time-varying transition
#' matrix \eqn{\Pi_t} using a logistic specification for the diagonal elements, and equal off-diagonal
#' probabilities. When omitted, the transition matrix is either fixed (if \code{P} is supplied)
#' or assumed uniform (if \code{P} is \code{NULL}).
#'
#' The correlation matrices are reconstructed from the lower-triangular parameters in \code{rho_matrix}.
#' The likelihood is computed using multivariate normal densities with mean zero and unit variance.
#'
#' @export
rsdc_hamilton <- function(y, X = NULL, beta = NULL, rho_matrix, K, N, P = NULL) {

  if (!is.matrix(y)) stop("y must be a numeric matrix.")
  if (!is.null(X) && !is.matrix(X)) stop("X must be a numeric matrix or NULL.")
  if (!is.null(beta) && (!is.matrix(beta) || nrow(beta) != N)) stop("beta must be a matrix with N rows.")
  if (!is.matrix(rho_matrix) || nrow(rho_matrix) != N || ncol(rho_matrix) != K*(K - 1)/2) {
    stop("rho_matrix must be of dimension N × (K*(K-1)/2).")
  }
  if (!is.null(P) && (!is.matrix(P) || !all(dim(P) == c(N, N)))) {
    stop("P must be an N × N matrix if provided.")
  }

  T <- nrow(y)
  n_pairs <- K*(K - 1)/2

  # Build regime correlation matrices
  sigma <- array(dim = c(K, K, N))
  for (m in 1:N) {
    R <- diag(K)
    R[lower.tri(R)] <- rho_matrix[m, ]
    R[upper.tri(R)] <- t(R)[upper.tri(R)]

    # Check positive definiteness
    if (min(eigen(R)$values) < 1e-8) return(list(log_likelihood = -Inf))

    sigma[,,m] <- R
  }

  # Transition probability matrices
  P_mats <- array(dim = c(N, N, T))

  if (is.null(X) || is.null(beta)) {
    # No exogenous variables, use fixed transition matrix P
    if (is.null(P)) {
      # If no P matrix is provided, use equal probabilities
      for (t in 1:T) {
        P_mats[,,t] <- matrix(1/N, N, N)
      }
    } else {
      # Use provided P matrix
      for (t in 1:T) {
        P_mats[,,t] <- P
      }
    }
  } else {
    # Time-varying transition probabilities
    for (t in 1:T) {
      for (i in 1:N) {
        # Diagonal element using logistic function
        p_ii <- plogis(X[t,] %*% beta[i,])

        # Handle N=1 case
        if (N == 1) {
          P_mats[,,t] <- matrix(1)
          next
        }

        # Off-diagonal elements (equal probability)
        off_prob <- (1 - p_ii)/(N - 1)

        P_mats[i,,t] <- off_prob
        P_mats[i,i,t] <- p_ii
      }
    }
  }

  # Multivariate normal log-densities
  log_densities <- matrix(nrow = N, ncol = T)
  for (m in 1:N) {
    sigma_m <- sigma[,,m]
    inv_sigma <- try(solve(sigma_m + diag(1e-8, K)), silent = TRUE)
    if (inherits(inv_sigma, "try-error")) return(list(log_likelihood = -Inf))

    log_det <- determinant(sigma_m, logarithm = TRUE)$modulus[1]
    centered <- as.matrix(y)
    inv_sigma <- as.matrix(inv_sigma)
    centered <- apply(centered, 2, as.numeric)
    inv_sigma <- matrix(as.numeric(inv_sigma), nrow = nrow(inv_sigma))
    quad_form <- rowSums((centered %*% inv_sigma) * centered)
    log_densities[m,] <- -0.5*(log_det + quad_form + K*log(2*pi))
  }

  # Filtering and smoothing
  filtered <- matrix(0, N, T)
  predicted <- matrix(0, N, T)
  smoothed <- matrix(0, N, T)
  xi <- rep(1/N, N)
  log_lik <- 0

  for (t in 1:T) {
    P_t <- P_mats[,,t]
    predicted[,t] <- P_t %*% xi

    likelihood <- exp(log_densities[,t])
    filtered[,t] <- predicted[,t] * likelihood
    sum_filtered <- sum(filtered[,t])

    if (sum_filtered <= 0) return(list(log_likelihood = -Inf))

    filtered[,t] <- filtered[,t] / sum_filtered
    log_lik <- log_lik + log(sum_filtered)
    xi <- filtered[,t]
  }

  # Smoothing
  smoothed[,T] <- filtered[,T]
  for (t in (T - 1):1) {
    P_t1 <- P_mats[,, t + 1]
    temp <- filtered[,t] * (t(P_t1) %*% (smoothed[, t + 1]/predicted[, t + 1]))
    smoothed[,t] <- temp / sum(temp)
  }

  list(
    filtered_probs = filtered,
    smoothed_probs = smoothed,
    log_likelihood = log_lik
  )
}

#' Negative Log-Likelihood for Regime-Switching Correlation Model
#'
#' Computes the negative log-likelihood of a multivariate regime-switching correlation model,
#' with either fixed or time-varying transition probabilities. The model assumes
#' regime-specific correlation matrices and uses the Hamilton filter for likelihood evaluation.
#'
#' @param params A numeric vector of model parameters. This includes:
#' \itemize{
#'   \item If \code{exog} is \code{NULL}: \eqn{N(N-1)} transition probabilities (one per off-diagonal regime),
#'   followed by \eqn{N \cdot K(K-1)/2} correlation parameters.
#'   \item If \code{exog} is not \code{NULL}: \eqn{N \cdot p} logistic regression coefficients
#'   (where \code{p = ncol(exog)}), followed by \eqn{N \cdot K(K-1)/2} correlation parameters.
#' }
#' @param y A numeric matrix of dimension \eqn{T \times K}, where \eqn{T} is the number of time steps and
#' \eqn{K} the number of observed variables (e.g., asset returns).
#' @param exog Optional. A numeric matrix of dimension \eqn{T \times p} of exogenous variables.
#' If provided, a time-varying transition model is used.
#' @param K Integer. Number of observed variables (columns in \code{y}).
#' @param N Integer. Number of regimes.
#'
#' @return A numeric scalar: the negative log-likelihood of the model (to be minimized).
#'
#' @details
#' The function evaluates the log-likelihood of the data under a regime-switching
#' correlation model using the Hamilton filter. Two types of transition dynamics are supported:
#' \itemize{
#'   \item \strong{Fixed transition matrix} when \code{exog = NULL}. Transition probabilities are estimated directly.
#'   \item \strong{Time-varying transition probabilities (TVTP)} when \code{exog} is provided. In this case,
#'   regime persistence probabilities are modeled via a logistic function:
#'   \deqn{p_{ii,t} = \frac{1}{1 + \exp(-X_t^\top \beta_i)}}
#' }
#'
#' The Hamilton filter computes the filtered and smoothed regime probabilities internally.
#' The negative log-likelihood is returned for use in optimization routines.
#'
#' @seealso \code{\link{rsdc_hamilton}}, \code{\link{optim}}, \code{\link{DEoptim}}, \code{\link{plogis}}
#'
#' @export
rsdc_likelihood <- function(params, y, exog = NULL, K, N) {
  if (any(is.na(params)) || any(!is.finite(params))) return(1e10)

  # Parameter count calculation
  n_p <- ifelse(is.null(exog), N * (N - 1), 0)  # Include transition params only when no exog
  n_pairs <- K * (K - 1) / 2  # Number of correlation pairs
  n_rho <- N * n_pairs  # Correlation parameters
  n_beta <- if (!is.null(exog)) N * ncol(exog) else 0

  # Parameter extraction
  if (is.null(exog)) {
    # Without exogenous variables, extract transition parameters first
    trans_params <- params[1:n_p]
    rho_matrix <- matrix(params[(n_p + 1):(n_p + n_rho)], nrow = N)
    beta <- NULL

    # Build transition matrix
    P <- matrix(0, N, N)
    if (N == 2) {
      P <- matrix(c(trans_params[1], 1 - trans_params[1],
                    1 - trans_params[2], trans_params[2]),
                  nrow = N, byrow = TRUE)
    }
  } else {
    # With exogenous variables
    beta <- matrix(params[1:n_beta], nrow = N)
    rho_matrix <- matrix(params[(n_beta + 1):(n_beta + n_rho)], nrow = N)
    P <- NULL
  }

  # Check correlation validity
  if (any(abs(rho_matrix) >= 1)) return(1e10)

  # Run Hamilton filter
  result <- rsdc_hamilton(y, exog, beta, rho_matrix, K, N, P)

  if (!is.finite(result$log_likelihood)) return(1e10)

  return(-result$log_likelihood)
}
