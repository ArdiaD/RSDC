#' Regime-Switching Correlation (TVTP) — DEoptim + L-BFGS-B
#'
#' Estimates a multivariate correlation model with \emph{time-varying transition probabilities} (TVTP).
#' Regime persistence \eqn{p_{ii,t}} is modeled via a logistic link on exogenous covariates \code{X};
#' off-diagonal transition probabilities are set equal to \eqn{(1 - p_{ii,t})/(N-1)}. Correlations are
#' regime-specific and reconstructed from lower-triangular parameters. Optimization uses global search
#' (\pkg{DEoptim}) followed by local refinement (L-BFGS-B).
#'
#' @param N Integer. Number of regimes.
#' @param residuals Numeric matrix \eqn{T \times K}. Typically standardized residuals or returns
#'   (unit variance per column is assumed).
#' @param X Numeric matrix \eqn{T \times p} of exogenous covariates driving TVTP.
#'   Include an intercept column yourself if desired (no automatic intercept).
#' @param out_of_sample Logical. If \code{TRUE}, estimation uses the first 70% of rows; the remainder is
#'   held out (indices are fixed by a 70/30 split).
#' @param control Optional list for basic verbosity control, currently supporting
#'   \code{do_trace = TRUE} to print optimizer progress and \code{seed} for
#'   reproducibility (default = 123). Algorithmic hyperparameters are set internally.
#'
#' @returns A list with:
#' \describe{
#'   \item{transition_matrix}{Average \eqn{N \times N} transition matrix evaluated at \code{colMeans(X)}.}
#'   \item{means}{\eqn{N \times K} zeros (placeholders; means are not estimated).}
#'   \item{volatilities}{\eqn{N \times K} ones (placeholders; vols are not estimated).}
#'   \item{correlations}{\eqn{N \times C} matrix of lower-triangular correlations, \eqn{C = K(K-1)/2}.}
#'   \item{covariances}{Array \eqn{K \times K \times N} of regime correlation matrices.}
#'   \item{log_likelihood}{Maximized log-likelihood (numeric scalar).}
#'   \item{beta}{\eqn{N \times p} matrix of TVTP coefficients.}
#' }
#'
#' @details
#' \itemize{
#'   \item \strong{TVTP parameterization:} diagonal entries use \code{plogis(X \%*\% beta[i, ])};
#'         off-diagonals are equal and sum to one.
#'   \item \strong{Bounds:} \eqn{\beta \in [-10, 10]} elementwise; correlations \eqn{\in (-1, 1)}.
#'   \item \strong{State ordering:} regimes are reordered by ascending mean correlation for identifiability.
#'   \item \strong{Numerical safeguards:} tiny ridge is added in inverses; non-PD proposals are penalized.
#' }
#'
#' @section Assumptions:
#' Inputs in \code{residuals} are treated as mean-zero with unit variance; only the correlation structure is estimated.
#'
#' @seealso \code{\link{f_optim_noX}}, \code{\link{f_optim_const}}, \code{\link{rsdc_estimate}}, \code{\link{rsdc_hamilton}}
#' @references
#' \insertRef{R-DEoptim}{RSDC}
#' \insertRef{DEoptim-JSS}{RSDC}
#' \insertRef{hamilton1989}{RSDC}
#'
#' @note Uses \code{DEoptim::DEoptim()} then \code{stats::optim(method = "L-BFGS-B")}.
#'
#' @importFrom stats optim plogis
f_optim <- function(N, residuals, X, out_of_sample = FALSE, control = list()) {

  con <- list(seed = 123, do_trace = FALSE)
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
  set.seed(con$seed)
  results <- DEoptim::DEoptim(
    fn = rsdc_likelihood,
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
      parVar = c("rsdc_hamilton", "rsdc_likelihood", "plogis"),
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
    fn = rsdc_likelihood,
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

#' Regime-Switching Correlation (Fixed Transition Matrix, no X)
#'
#' Estimates a multivariate correlation model with a \emph{time-invariant} (fixed) regime
#' transition matrix. No exogenous covariates are used; both the fixed transition probabilities
#' and the regime-specific correlation parameters are estimated by global search
#' (\pkg{DEoptim}) followed by local refinement (L-BFGS-B).
#'
#' @param N Integer. Number of regimes.
#' @param residuals Numeric matrix \eqn{T \times K}. Typically standardized residuals or returns
#'   (columns are treated as mean-zero with unit variance).
#' @param out_of_sample Logical. If \code{TRUE}, estimation uses the first 70\% of rows; the
#'   remainder is held out. (Split index is a fixed 70/30 cut.)
#' @param control Optional list for basic verbosity control, currently supporting
#'   \code{do_trace = TRUE} to print optimizer progress and \code{seed} for
#'   reproducibility (default = 123). Algorithmic hyperparameters are set internally.
#'
#' @returns A list with:
#' \describe{
#'   \item{transition_matrix}{Estimated \eqn{N \times N} fixed transition matrix.}
#'   \item{means}{\eqn{N \times K} zeros (placeholders; means are not estimated).}
#'   \item{volatilities}{\eqn{N \times K} ones (placeholders; vols are not estimated).}
#'   \item{correlations}{\eqn{N \times C} matrix of lower-triangular correlations,
#'     with \eqn{C = K(K-1)/2}.}
#'   \item{covariances}{Array \eqn{K \times K \times N} of regime correlation matrices.}
#'   \item{log_likelihood}{Maximized log-likelihood (numeric scalar).}
#'   \item{beta}{\code{NULL} (no covariates in this specification).}
#' }
#'
#' @details
#' \itemize{
#'   \item \strong{Parameterization:} For \eqn{N=2}, the transition parameters are
#'         \eqn{\{p_{11}, p_{22}\}}; off-diagonals are \eqn{1 - p_{ii}}.
#'   \item \strong{Bounds:} \eqn{p_{ii} \in [0.01, 0.99]} and \eqn{\rho \in (-1, 1)}
#'     for numerical stability and identifiability.
#'   \item \strong{Correlation build:} Each regime’s correlation matrix is reconstructed from
#'     its lower-triangular vector and symmetrized; non-PD proposals are penalized in the
#'     likelihood routine.
#' }
#'
#' @section Assumptions:
#' Inputs in \code{residuals} are treated as mean-zero with unit variance; only the correlation
#' structure and fixed transition probabilities are estimated.
#'
#' @seealso \code{\link{f_optim}} (TVTP), \code{\link{f_optim_const}} (constant correlation),
#'   \code{\link{rsdc_estimate}} (wrapper), \code{\link{rsdc_hamilton}}, \code{\link{rsdc_likelihood}}
#'
#' @references
#' \insertRef{R-DEoptim}{RSDC}
#' \insertRef{DEoptim-JSS}{RSDC}
#' \insertRef{hamilton1989}{RSDC}
#'
#' @note Uses \code{DEoptim::DEoptim()} then \code{stats::optim(method = "L-BFGS-B")}.
#'
#' @importFrom stats optim
f_optim_noX <- function(N, residuals, out_of_sample = FALSE, control = list()) {
  stopifnot(is.matrix(residuals))

  con <- list(seed = 123, do_trace = FALSE)
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

  set.seed(con$seed)
  de_result <- DEoptim::DEoptim(
    fn = rsdc_likelihood,
    lower = bounds$lower,
    upper = bounds$upper,
    control = list(itermax = 500, NP = 10 * length(bounds$lower),
                   strategy = 2, F = 0.8, CR = 0.9,
                   trace = con$do_trace, reltol = 1e-8, steptol = 50),
    y = y_optim, exog = NULL, K = K, N = N
  )

  optim_result <- optim(
    par = de_result$optim$bestmem,
    fn = rsdc_likelihood,
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

#' Constant-Correlation Model (Single Regime)
#'
#' Estimates a single (time-invariant) correlation matrix assuming asset correlations
#' are constant over time. No regime switching or exogenous covariates are used.
#' Optimization uses a global stage (\pkg{DEoptim}) followed by local refinement
#' (L-BFGS-B via \code{stats::optim}).
#'
#' @param residuals Numeric matrix \eqn{T \times K}. Typically standardized residuals or
#'   returns (columns treated as mean-zero with unit variance).
#' @param out_of_sample Logical. If \code{TRUE}, estimation uses the first 70\% of rows; the
#'   remainder is held out. (Split index is a fixed 70/30 cut.)
#' @param control Optional list for basic verbosity control, currently supporting
#'   \code{do_trace = TRUE} to print optimizer progress and \code{seed} for
#'   reproducibility (default = 123). Algorithmic hyperparameters are set internally.
#'
#' @returns A list with:
#' \describe{
#'   \item{transition_matrix}{\eqn{1 \times 1} matrix equal to 1 (no switching).}
#'   \item{means}{\eqn{1 \times K} zeros (placeholders; means are not estimated).}
#'   \item{volatilities}{\eqn{1 \times K} ones (placeholders; vols are not estimated).}
#'   \item{correlations}{\eqn{1 \times C} vector of lower-triangular correlations,
#'     where \eqn{C = K(K-1)/2}.}
#'   \item{covariances}{Array \eqn{K \times K \times 1} with the full correlation matrix.}
#'   \item{log_likelihood}{Maximized log-likelihood (numeric scalar).}
#'   \item{beta}{\code{NULL} (no covariates).}
#' }
#'
#' @details
#' \itemize{
#'   \item \strong{Parameterization:} the free parameters are the \eqn{C = K(K-1)/2}
#'         pairwise correlations stacked in \code{lower.tri} order.
#'   \item \strong{Bounds / PD checks:} elementwise bounds \eqn{\rho \in (-1, 1)};
#'         non–positive-definite proposals are penalized in the objective.
#'   \item \strong{Numerical safeguards:} a tiny ridge is added before inverting
#'         the correlation matrix to improve stability.
#'   \item \strong{Scaling:} since inputs are treated as unit-variance, only the
#'         correlation structure is estimated.
#' }
#'
#' @seealso \code{\link{rsdc_estimate}} (wrapper),
#'   \code{\link{f_optim}} (TVTP), \code{\link{f_optim_noX}} (fixed transition matrix)
#'
#' @references
#' \insertRef{R-DEoptim}{RSDC}
#' \insertRef{DEoptim-JSS}{RSDC}
#' \insertRef{hamilton1989}{RSDC}
#' \insertRef{pelletier2006regime}{RSDC}
#'
#' @note Uses \code{DEoptim::DEoptim()} then \code{stats::optim(method = "L-BFGS-B")}.
#'
#' @importFrom stats optim
f_optim_const <- function(residuals, out_of_sample = FALSE, control = list()) {
  stopifnot(is.matrix(residuals))

  con <- list(seed = 123, do_trace = FALSE)
  con[names(control)] <- control

  K <- ncol(residuals)
  n_rho <- K * (K - 1) / 2

  neg_loglik_const <- function(rho_vec, y, K) {
    R <- diag(K)
    R[lower.tri(R)] <- rho_vec
    R[upper.tri(R)] <- t(R)[upper.tri(R)]
    # PD check
    eig <- tryCatch(eigen(R, only.values = TRUE)$values, error = function(e) NA_real_)
    if (any(!is.finite(eig)) || min(eig) < 1e-8) return(1e10)

    inv_R  <- solve(R + diag(1e-8, K))
    logdet <- determinant(R, logarithm = TRUE)$modulus[1]
    quad   <- rowSums((y %*% inv_R) * y)
    # negative log-likelihood
    return(-(-0.5 * sum(logdet + quad + K * log(2 * pi))))
  }

  # Split (if requested)
  if (out_of_sample) {
    cut <- round(0.7 * nrow(residuals))
    y_is  <- residuals[1:cut, , drop = FALSE]
    y_oos <- residuals[(cut + 1):nrow(residuals), , drop = FALSE]
  } else {
    y_is  <- residuals
    y_oos <- NULL
  }

  # Fit on IS only
  set.seed(con$seed)
  de_result <- DEoptim::DEoptim(
    fn = neg_loglik_const,
    lower = rep(-1, n_rho),
    upper = rep(1, n_rho),
    control = list(itermax = 500, trace = con$do_trace),
    y = y_is, K = K
  )

  optim_result <- optim(
    par = de_result$optim$bestmem,
    fn = neg_loglik_const,
    method = "L-BFGS-B",
    lower = rep(-1, n_rho),
    upper = rep(1, n_rho),
    y = y_is, K = K,
    control = list(maxit = 1000)
  )

  rho_vec <- optim_result$par
  R <- diag(K)
  R[lower.tri(R)] <- rho_vec
  R[upper.tri(R)] <- t(R)[upper.tri(R)]
  array_R <- array(R, dim = c(K, K, 1))

  # Single returned likelihood depending on the flag
  ll_return <- if (out_of_sample) {
    -neg_loglik_const(rho_vec, y_is, K)  # OOS ll
  } else {
    -optim_result$value                   # IS ll
  }

  list(
    transition_matrix = matrix(1, 1, 1),
    means        = matrix(0, 1, K),
    volatilities = matrix(1, 1, K),
    correlations = matrix(rho_vec, nrow = 1),
    covariances  = array_R,
    log_likelihood = ll_return,
    beta = NULL
  )
}

#' Estimate Regime-Switching or Constant Correlation Model (Wrapper)
#'
#' Unified front-end that dispatches to one of three estimators:
#' \itemize{
#'   \item \code{f_optim()} — TVTP specification (\code{method = "tvtp"}).
#'   \item \code{f_optim_noX()} — fixed transition matrix (\code{method = "noX"}).
#'   \item \code{f_optim_const()} — constant correlation, single regime (\code{method = "const"}).
#' }
#'
#' @param method Character. One of \code{"tvtp"}, \code{"noX"}, \code{"const"}.
#' @param residuals Numeric matrix \eqn{T \times K}. Typically standardized residuals/returns.
#' @param N Integer. Number of regimes. Ignored when \code{method = "const"}.
#' @param X Numeric matrix \eqn{T \times p} of exogenous covariates (required for \code{"tvtp"}).
#' @param out_of_sample Logical. If \code{TRUE}, a fixed 70/30 split is applied prior to estimation.
#' @param control Optional list. Currently forwards \code{do_trace = FALSE} and \code{seed = 123} to the backends.
#'
#' @return
#' \describe{
#'   \item{\code{transition_matrix}}{Estimated transition matrix (\eqn{1 \times 1} for \code{"const"}).}
#'   \item{\code{correlations}}{Regime lower-triangular correlations.}
#'   \item{\code{covariances}}{Array of full correlation matrices.}
#'   \item{\code{log_likelihood}}{Maximized log-likelihood.}
#'   \item{\code{beta}}{TVTP coefficients (only for \code{"tvtp"}).}
#' }
#'
#' @details
#' \itemize{
#'   \item \strong{Method selection:} \code{match.arg()} validates \code{method}.
#'   \item \strong{Inputs:} \code{"tvtp"} requires non-NULL \code{X}; \code{N} is ignored for \code{"const"}.
#'   \item \strong{Split:} If \code{out_of_sample = TRUE}, the first 70\% is used for fitting.
#' }
#'
#' @references
#' \insertRef{R-DEoptim}{RSDC} \cr
#' \insertRef{DEoptim-JSS}{RSDC} \cr
#' \insertRef{hamilton1989}{RSDC} \cr
#' \insertRef{pelletier2006regime}{RSDC}
#'
#' @examples
#' \donttest{
#' y <- scale(matrix(rnorm(100 * 3), 100, 3))
#' rsdc_estimate("const", residuals = y)
#' rsdc_estimate("noX", residuals = y, N = 2)
#' X <- cbind(1, scale(seq_len(nrow(y))))
#' rsdc_estimate("tvtp", residuals = y, N = 2, X = X)
#' }
#'
#' @seealso \code{\link{f_optim}}, \code{\link{f_optim_noX}}, \code{\link{f_optim_const}},
#'   \code{\link{rsdc_hamilton}}, \code{\link{rsdc_likelihood}}
#'
#' @export
rsdc_estimate <- function(method = c("tvtp", "noX", "const"),
                          residuals, N = 2, X = NULL,
                          out_of_sample = FALSE, control = list()) {
  method <- match.arg(method)

  con <- list(seed = 123, do_trace = FALSE)
  con[names(control)] <- control

  if (method == "tvtp") {
    if (is.null(X)) stop("X must be provided for method = 'tvtp'")
    return(f_optim(N = N, residuals = residuals, X = X, out_of_sample = out_of_sample, control = list(seed = con$seed, do_trace = con$do_trace)))
  }

  if (method == "noX") {
    return(f_optim_noX(N = N, residuals = residuals, out_of_sample = out_of_sample, control = list(seed = con$seed, do_trace = con$do_trace)))
  }

  if (method == "const") {
    return(f_optim_const(residuals = residuals, out_of_sample = out_of_sample,
                         control = list(seed = con$seed, do_trace = con$do_trace)))
  }

  stop("Unknown method: ", method)
}
