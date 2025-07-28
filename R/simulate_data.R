#' Simulate Multivariate Regime-Switching Data with TVTP
#'
#' Simulates a multivariate time series from a regime-switching model
#' with time-varying transition probabilities (TVTP) driven by covariates.
#'
#' @param n Integer. Number of time steps to simulate.
#' @param X Numeric matrix (n × p) of covariates used to determine transition probabilities.
#'          Each row \code{X[t, ]} corresponds to covariates available at time \code{t}.
#' @param beta Numeric array (N × N × p). Coefficients for the multinomial logistic transition model.
#'             Each slice \code{beta[i, j, ]} corresponds to the transition from state \code{i} to \code{j}.
#' @param mu Numeric matrix (N × K). Regime-specific mean vectors.
#' @param sigma Numeric array (K × K × N). Regime-specific covariance matrices.
#' @param N Integer. Number of regimes.
#' @param seed Optional integer. If provided, sets the random seed for reproducibility.
#'
#' @return A list with components:
#' \describe{
#'   \item{states}{Integer vector (length \code{n}) indicating the active regime at each time.}
#'   \item{observations}{Numeric matrix (n × K) of simulated multivariate observations.}
#'   \item{transition_matrices}{Array (N × N × n) containing the transition matrix at each time \code{t}.}
#' }
#' @export
f_simul_tvtp <- function(n, X, beta, mu, sigma, N, seed = NULL) {
  if (!requireNamespace("mvtnorm", quietly = TRUE)) {
    stop("Please install 'mvtnorm'.")
  }

  if (!is.null(seed)) set.seed(seed)

  K <- ncol(mu)
  p <- ncol(X)

  if (!all(dim(beta) == c(N, N, p))) stop("beta must be an array of dim N × N × p")
  if (!all(dim(mu) == c(N, K))) stop("mu must be a matrix of dim N × K")
  if (!all(dim(sigma) == c(K, K, N))) stop("sigma must be an array of dim K × K × N")

  states <- numeric(n)
  observations <- matrix(NA, nrow = n, ncol = K)
  transition_matrices <- array(NA, dim = c(N, N, n))

  # Initial state
  states[1] <- sample(1:N, size = 1)
  observations[1, ] <- mvtnorm::rmvnorm(1, mean = mu[states[1], ], sigma = sigma[,,states[1]])

  # Iterate over time
  for (t in 2:n) {
    P_t <- matrix(NA, nrow = N, ncol = N)

    for (i in 1:N) {
      logits <- sapply(1:N, function(j) X[t - 1, , drop = FALSE] %*% beta[i, j, ])
      exp_logits <- exp(logits - max(logits))  # Numerical stability
      P_t[i, ] <- exp_logits / sum(exp_logits)
    }

    transition_matrices[,,t] <- P_t
    prev_state <- states[t - 1]
    states[t] <- sample(1:N, size = 1, prob = P_t[prev_state, ])
    observations[t, ] <- mvtnorm::rmvnorm(1, mean = mu[states[t], ], sigma = sigma[,,states[t]])
  }

  list(
    states = states,
    observations = observations,
    transition_matrices = transition_matrices
  )
}
