#' Hamilton Filter (Fixed P or TVTP)
#'
#' Runs the Hamilton (1989) filter for a multivariate regime-switching \emph{correlation} model.
#' Supports either a fixed (time-invariant) transition matrix \eqn{P} or time-varying transition
#' probabilities (TVTP) built from exogenous covariates \code{X} via a logistic link.
#' Returns filtered/smoothed regime probabilities and the log-likelihood.
#'
#' @param y Numeric matrix \eqn{T \times K} of observations (e.g., standardized residuals/returns).
#'   Columns are treated as mean-zero with unit variance; only the correlation structure is modeled.
#' @param X Optional numeric matrix \eqn{T \times p} of covariates for TVTP. Required if \code{beta} is supplied.
#' @param beta Optional numeric matrix \eqn{N \times p}. TVTP coefficients; row \eqn{i} governs
#'   persistence of regime \eqn{i} via \code{plogis(X[t, ] \%*\% beta[i, ])}.
#' @param rho_matrix Numeric matrix \eqn{N \times C} of regime correlation parameters, where
#'   \eqn{C = K(K-1)/2}. Each row is the lower-triangular part (by \code{lower.tri}) of a regime's
#'   correlation matrix.
#' @param K Integer. Number of observed series (columns of \code{y}).
#' @param N Integer. Number of regimes.
#' @param P Optional \eqn{N \times N} fixed transition matrix. Used only when \code{X} or \code{beta} is \code{NULL}.
#'
#' @returns A list with:
#' \describe{
#'   \item{filtered_probs}{\eqn{N \times T} matrix of filtered probabilities
#'     \eqn{\Pr(S_t = j \mid \Omega_t)}.}
#'   \item{smoothed_probs}{\eqn{N \times T} matrix of smoothed probabilities
#'     \eqn{\Pr(S_t = j \mid \Omega_T)}.}
#'   \item{log_likelihood}{Scalar log-likelihood of the model given \code{y}.}
#' }
#'
#' @details
#' \itemize{
#'   \item \strong{Correlation rebuild:} For regime \eqn{m}, a correlation matrix \eqn{R_m} is reconstructed
#'         from \code{rho_matrix[m, ]} (lower-triangular fill + symmetrization). Non-PD proposals are penalized.
#'   \item \strong{Transition dynamics:}
#'         \itemize{
#'           \item \emph{Fixed P:} If \code{X} or \code{beta} is missing, a constant \eqn{P} is used
#'                 (user-provided via \code{P}; otherwise uniform \eqn{1/N} rows).
#'           \item \emph{TVTP:} With \code{X} and \code{beta}, diagonal entries use
#'                 \code{plogis(X[t, ] \%*\% beta[i, ])}. Off-diagonals are equal and sum to \eqn{1 - p_{ii,t}}.
#'                 For \eqn{N=1}, \eqn{P_t = [1]}.
#'         }
#'   \item \strong{Numerical safeguards:} A small ridge is added before inversion; if filtering
#'         degenerates at a time step, \code{log_likelihood = -Inf} is returned.
#' }
#'
#' @examples
#' set.seed(1)
#' T <- 50; K <- 3; N <- 2
#' y <- scale(matrix(rnorm(T * K), T, K), center = TRUE, scale = TRUE)
#'
#' # Example rho: two regimes with different average correlations
#' rho <- rbind(c(0.10, 0.05, 0.00),
#'              c(0.60, 0.40, 0.30))  # lower-tri order for K=3
#'
#' # Fixed-P filtering
#' Pfix <- matrix(c(0.9, 0.1,
#'                  0.2, 0.8), nrow = 2, byrow = TRUE)
#' out_fix <- rsdc_hamilton(y = y, X = NULL, beta = NULL,
#'                          rho_matrix = rho, K = K, N = N, P = Pfix)
#' str(out_fix$filtered_probs)
#'
#' # TVTP filtering (include an intercept yourself)
#' X <- cbind(1, scale(seq_len(T)))
#' beta <- rbind(c(1.2, 0.0),
#'               c(0.8, -0.1))
#' out_tvtp <- rsdc_hamilton(y = y, X = X, beta = beta,
#'                           rho_matrix = rho, K = K, N = N)
#' out_tvtp$log_likelihood
#'
#' @seealso \code{\link{rsdc_likelihood}} and \code{\link{rsdc_estimate}}.
#'
#' @references
#' \insertRef{hamilton1989}{RSDC}
#'
#' @note TVTP uses a logistic link on the diagonal; off-diagonals are equal by construction.
#'
#' @importFrom stats plogis
#' @export
rsdc_hamilton <- function(y, X = NULL, beta = NULL, rho_matrix, K, N, P = NULL) {

  if (!is.matrix(y)) stop("y must be a numeric matrix.")
  if (!is.null(X) && !is.matrix(X)) stop("X must be a numeric matrix or NULL.")
  if (!is.null(beta) && (!is.matrix(beta) || nrow(beta) != N)) stop("beta must be a matrix with N rows.")
  if (!is.matrix(rho_matrix) || nrow(rho_matrix) != N || ncol(rho_matrix) != K*(K - 1)/2) {
    stop("rho_matrix must be of dimension N x (K*(K-1)/2).")
  }
  if (!is.null(P) && (!is.matrix(P) || !all(dim(P) == c(N, N)))) {
    stop("P must be an N x N matrix if provided.")
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

#' Negative Log-Likelihood for Regime-Switching Correlation Models
#'
#' Computes the negative log-likelihood for a multivariate \emph{correlation-only}
#' regime-switching model, with either a fixed (time-invariant) transition matrix or
#' time-varying transition probabilities (TVTP) driven by exogenous covariates.
#' Likelihood evaluation uses the Hamilton (1989) filter.
#'
#' @param params Numeric vector of model parameters packed as:
#' \itemize{
#'   \item \strong{No exogenous covariates} (\code{exog = NULL}): first
#'         \eqn{N(N-1)} transition parameters (for the fixed transition matrix),
#'         followed by \eqn{N \times K(K-1)/2} correlation parameters, stacked
#'         \emph{row-wise by regime} in \code{lower.tri} order.
#'   \item \strong{With exogenous covariates} (\code{exog} provided): first
#'         \eqn{N \times p} TVTP coefficients (\code{beta}, row \eqn{i} corresponds
#'         to regime \eqn{i}), followed by \eqn{N \times K(K-1)/2} correlation parameters,
#'         stacked \emph{row-wise by regime} in \code{lower.tri} order.
#' }
#' @param y Numeric matrix \eqn{T \times K} of observations (e.g., standardized residuals).
#'          Columns are treated as mean-zero, unit-variance; only correlation is modeled.
#' @param exog Optional numeric matrix \eqn{T \times p} of exogenous covariates.
#'        If supplied, a TVTP specification is used.
#' @param K Integer. Number of observed series (columns of \code{y}).
#' @param N Integer. Number of regimes.
#'
#' @returns Numeric scalar: the \emph{negative} log-likelihood to be minimized by an optimizer.
#'
#' @details
#' \itemize{
#'   \item \strong{Transition dynamics:}
#'         \itemize{
#'           \item \emph{Fixed P (no \code{exog}):} \code{params} begins with transition
#'                 parameters. For \eqn{N=2}, the implementation maps them to
#'                 \eqn{P=\begin{pmatrix} p_{11} & 1-p_{11}\\ 1-p_{22} & p_{22}\end{pmatrix}}.
#'           \item \emph{TVTP:} with \code{exog}, diagonal persistence is
#'                 \eqn{p_{ii,t} = \mathrm{logit}^{-1}(X_t^\top \beta_i)}; off-diagonals are equal
#'                 and sum to \eqn{1-p_{ii,t}}.
#'         }
#'   \item \strong{Correlation build:} per regime, the lower-triangular vector is filled into
#'         a symmetric correlation matrix. Non-positive-definite proposals or \eqn{|\rho|\ge 1}
#'         are penalized via a large objective value.
#'   \item \strong{Evaluation:} delegates to \code{\link{rsdc_hamilton}}; if the filter returns
#'         \code{log_likelihood = -Inf}, a large penalty is returned.
#' }
#'
#' @examples
#' # Small toy example (N = 2, K = 3), fixed P (no exog)
#' set.seed(1)
#' T <- 40; K <- 3; N <- 2
#' y <- scale(matrix(rnorm(T * K), T, K), center = TRUE, scale = TRUE)
#'
#' # Pack parameters: trans (p11, p22), then rho by regime (lower-tri order)
#' p11 <- 0.9; p22 <- 0.8
#' rho1 <- c(0.10, 0.05, 0.00)  # (2,1), (3,1), (3,2)
#' rho2 <- c(0.60, 0.40, 0.30)
#' params <- c(p11, p22, rho1, rho2)
#'
#' nll <- rsdc_likelihood(params, y = y, exog = NULL, K = K, N = N)
#' nll
#'
#' # TVTP example: add X and beta (pack beta row-wise, then rho)
#' X <- cbind(1, scale(seq_len(T)))
#' beta <- rbind(c(1.2, 0.0),
#'               c(0.8, -0.1))
#' params_tvtp <- c(as.vector(t(beta)), rho1, rho2)
#' nll_tvtp <- rsdc_likelihood(params_tvtp, y = y, exog = X, K = K, N = N)
#' nll_tvtp
#'
#' @seealso \code{\link{rsdc_hamilton}} (filter),
#' \code{\link[stats]{optim}}, and \code{\link[DEoptim]{DEoptim}}
#'
#' @note The function is written for use inside optimizers; it performs inexpensive validation
#'       and returns large penalties for invalid parameterizations instead of stopping with errors.
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
    rho_matrix <- matrix(params[(n_p + 1):(n_p + n_rho)], nrow = N, byrow = TRUE)
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
