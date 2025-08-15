#' Simulate Multivariate Regime-Switching Data (TVTP)
#'
#' Simulates a multivariate time series from a regime-switching model with
#' \emph{time-varying transition probabilities} (TVTP) driven by covariates \code{X}.
#' Transition probabilities are generated via a multinomial logistic (softmax) link;
#' observations are drawn from regime-specific Gaussian distributions.
#'
#' @param n Integer. Number of time steps to simulate.
#' @param X Numeric matrix \eqn{n \times p} of covariates used to form the transition
#'   probabilities. Row \code{X[t, ]} corresponds to covariates available at time \code{t}.
#'   Only rows \code{1:(n-1)} are used to transition from \code{t-1} to \code{t}.
#' @param beta Numeric array \eqn{N \times N \times p}. Softmax coefficients for the
#'   multinomial transition model; \code{beta[i, j, ]} parameterizes the transition
#'   from state \eqn{i} to state \eqn{j}.
#' @param mu Numeric matrix \eqn{N \times K}. Regime-specific mean vectors.
#' @param sigma Numeric array \eqn{K \times K \times N}. Regime-specific covariance
#'   (here, correlation/variance) matrices; each \eqn{K \times K} slice must be
#'   symmetric positive definite.
#' @param N Integer. Number of regimes.
#' @param seed Optional integer. If supplied, sets the RNG seed for reproducibility.
#'
#' @returns A list with:
#' \describe{
#'   \item{states}{Integer vector of length \code{n}; the simulated regime index at each time.}
#'   \item{observations}{Numeric matrix \eqn{n \times K}; the simulated observations.}
#'   \item{transition_matrices}{Array \eqn{N \times N \times n}; the transition matrix \eqn{P_t}
#'         used at each time step (with \eqn{P_1} undefined by construction; see Details).}
#' }
#'
#' @details
#' \itemize{
#'   \item \strong{Initial state and first draw:} The initial regime \eqn{S_1} is sampled
#'         uniformly; the first observation \eqn{y_1} is drawn from
#'         \eqn{\mathcal{N}(\mu_{S_1}, \Sigma_{S_1})}.
#'   \item \strong{TVTP via softmax:} For \eqn{t \ge 2}, the row \eqn{i} of \eqn{P_t} is
#'         \deqn{P_t(i, j) = \frac{\exp\!\big(X_{t-1}^\top \beta_{i,j}\big)}
#'                           {\sum_{h=1}^N \exp\!\big(X_{t-1}^\top \beta_{i,h}\big)}\,,}
#'         computed with log-sum-exp stabilization.
#'   \item \strong{Sampling:} Given \eqn{S_{t-1}}, draw \eqn{S_t} from the categorical
#'         distribution with probabilities \eqn{P_t(S_{t-1}, \cdot)} and
#'         \eqn{y_t \sim \mathcal{N}(\mu_{S_t}, \Sigma_{S_t})}.
#' }
#'
#' @examples
#' set.seed(123)
#' n <- 200; K <- 3; N <- 2; p <- 2
#' # Covariates: intercept + trend
#' X <- cbind(1, scale(seq_len(n)))
#'
#' # TVTP coefficients: encourage persistence on the diagonal
#' beta <- array(0, dim = c(N, N, p))
#' beta[1, 1, ] <- c(1.2,  0.0)
#' beta[2, 2, ] <- c(1.0, -0.1)
#'
#' # Regime means and covariances
#' mu <- rbind(c(0, 0, 0),
#'             c(0, 0, 0))
#' rho <- rbind(c(0.10, 0.05, 0.00),
#'              c(0.60, 0.40, 0.30))  # lower-tri for K=3: (2,1), (3,1), (3,2)
#' Sig <- array(0, dim = c(K, K, N))
#' for (m in 1:N) {
#'   R <- diag(K); R[lower.tri(R)] <- rho[m, ]; R[upper.tri(R)] <- t(R)[upper.tri(R)]
#'   Sig[, , m] <- R
#'
#' sim <- rsdc_simulate(n = n, X = X, beta = beta, mu = mu, sigma = Sig, N = N, seed = 99)
#' table(sim$states)
#' dim(sim$observations)
#' }
#'
#' @seealso \code{\link{rsdc_hamilton}} (filter/evaluation),
#'   \code{\link{rsdc_estimate}} (estimators),
#'   \code{\link{rsdc_forecast}} (forecasting)
#'
#' @references
#'   Hamilton, J. D. (1989). A new approach to the economic analysis of nonstationary time series
#'   and the business cycle. \emph{Econometrica}, 57(2), 357-384.
#'
#' @note Requires \pkg{mvtnorm} for multivariate normal sampling (called as \code{mvtnorm::rmvnorm}).
#'
#' @export
rsdc_simulate <- function(n, X, beta, mu, sigma, N, seed = NULL) {
  if (!requireNamespace("mvtnorm", quietly = TRUE)) {
    stop("Please install 'mvtnorm'.")
  }

  if (!is.null(seed)) set.seed(seed)

  K <- ncol(mu)
  p <- ncol(X)

  if (!all(dim(beta) == c(N, N, p))) stop("beta must be an array of dim N x N x p")
  if (!all(dim(mu) == c(N, K))) stop("mu must be a matrix of dim N x K")
  if (!all(dim(sigma) == c(K, K, N))) stop("sigma must be an array of dim K x K x N")

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
