#' Regime-Switching Correlation (TVTP) — DEoptim + L-BFGS-B
#'
#' Estimates a multivariate correlation model with \emph{time-varying transition probabilities} (TVTP).
#' Regime persistence is modeled via a logistic link for \eqn{N=2} or a softmax for \eqn{N \ge 3}
#' on exogenous covariates \code{X}. For \eqn{N=2}, off-diagonals equal \eqn{(1-p_{ii,t})/(N-1)};
#' for \eqn{N \ge 3}, all entries are independently determined by the softmax. Correlations are
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
#'   \item{transition_matrix}{Representative \eqn{N \times N} transition matrix evaluated at
#'     the in-sample covariate means (\code{colMeans} of the estimation window; the first 70\%
#'     when \code{out_of_sample = TRUE}). This is a point summary at the mean covariate, \emph{not} the
#'     time-average \eqn{\frac{1}{T}\sum_t P_t}; under the nonlinear logit/softmax link the
#'     two differ (Jensen's inequality). The full per-period dynamics are governed by \code{beta}.}
#'   \item{means}{\eqn{N \times K} zeros (placeholders; means are not estimated).}
#'   \item{volatilities}{\eqn{N \times K} ones (placeholders; vols are not estimated).}
#'   \item{correlations}{\eqn{N \times C} matrix of lower-triangular correlations, \eqn{C = K(K-1)/2}.}
#'   \item{covariances}{Array \eqn{K \times K \times N} of regime correlation matrices.}
#'   \item{log_likelihood}{In-sample log-likelihood (numeric scalar); full-sample when \code{out_of_sample = FALSE}.}
#'   \item{log_likelihood_oos}{OOS log-likelihood evaluated on the held-out 30\% (numeric scalar),
#'     or \code{NULL} when \code{out_of_sample = FALSE}.}
#'   \item{beta}{\eqn{N \times p} matrix of TVTP coefficients for \eqn{N=2};
#'     \eqn{N \times (N-1)p} for \eqn{N \ge 3}, packed row-wise.}
#' }
#'
#' @details
#' \itemize{
#'   \item \strong{TVTP parameterization:} for \eqn{N=2}, diagonal entries use
#'         \code{plogis(X \%*\% beta[i, ])} and off-diagonals are equal;
#'         for \eqn{N \ge 3}, a softmax over \eqn{(N-1)} free beta vectors per row,
#'         with all entries independently determined.
#'   \item \strong{Bounds:} \eqn{\beta \in [-10, 10]} elementwise; correlations \eqn{\in (-1, 1)}.
#'   \item \strong{State ordering:} regimes are reordered by ascending mean correlation for identifiability.
#'   \item \strong{Numerical safeguards:} tiny ridge is added in inverses; non-PD proposals are penalized.
#' }
#'
#' @section Assumptions:
#' Inputs in \code{residuals} are treated as mean-zero with unit variance; only the correlation structure is estimated.
#'
#' @seealso \code{\link{rsdc_estimate}}, \code{\link{rsdc_hamilton}}
#' @references
#' \insertRef{DEoptim-JSS}{RSDC}
#'
#' \insertRef{hamilton1989}{RSDC}
#'
#' @note Uses \code{DEoptim::DEoptim()} then \code{stats::optim(method = "L-BFGS-B")}.
#'
#' @importFrom stats optim plogis
#' @noRd
f_optim <- function(N, residuals, X, out_of_sample = FALSE, control = list()) {
  if (N > 3) stop("Only N = 2 or N = 3 is supported.")
  stopifnot(is.matrix(residuals))

  con <- list(seed = 123, do_trace = FALSE)
  con[names(control)] <- control

  # Parameter setup
  K <- ncol(residuals)
  p <- ncol(X)
  # N=2: one logistic beta vector per regime; N=3: (N-1) softmax beta vectors per regime row
  n_beta <- if (N == 2) N * p else N * (N - 1L) * p
  n_rho <- N * (K * (K - 1) / 2)
  rho_indices <- (n_beta + 1):(n_beta + n_rho)

  # Bounds (beta: [-10, 10] for all N; rho: (-1, 1))
  bounds <- list(
    lower = c(rep(-10, n_beta), rep(-1, n_rho)),
    upper = c(rep(10, n_beta), rep(1, n_rho))
  )

  # Out-of-sample handling (fixed indexing)
  if (out_of_sample) {
    is_cut        <- round(0.7 * nrow(residuals))
    residuals_is  <- residuals[1:is_cut, , drop = FALSE]
    residuals_oos <- residuals[(is_cut + 1):nrow(residuals), , drop = FALSE]
    X_is  <- X[1:is_cut, , drop = FALSE]
    X_oos <- X[(is_cut + 1):nrow(X), , drop = FALSE]
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
      parallelType = 0,
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

  # beta_dim: parameters per regime row in the flat vector
  beta_dim <- if (N == 2) p else (N - 1L) * p

  if (!all(sorted_states == 1:N)) {
    if (N == 2L) {
      beta_reordered <- unlist(lapply(sorted_states, function(s)
        optim_params[((s - 1L) * beta_dim + 1L):(s * beta_dim)]))
    } else {
      # N=3: beta columns must also be permuted and shifted for the new reference.
      # new_beta[new_row_i, new_col_j] = old_beta[s_i, s_j] - old_beta[s_i, s_N]
      # where s_i = sorted_states[i], s_j = sorted_states[j], s_N = sorted_states[N].
      ref_col <- sorted_states[N]
      beta_reordered <- unlist(lapply(sorted_states, function(s) {
        old_betas <- lapply(seq_len(N), function(col) {
          if (col < N)
            optim_params[((s - 1L) * beta_dim + (col - 1L) * p + 1L):((s - 1L) * beta_dim + col * p)]
          else rep(0, p)
        })
        ref_vec <- old_betas[[ref_col]]
        unlist(lapply(seq_len(N - 1L), function(j)
          old_betas[[sorted_states[j]]] - ref_vec))
      }))
    }
    optim_params <- c(
      beta_reordered,
      unlist(lapply(sorted_states, function(s)
        optim_params[rho_indices][((s - 1L) * C_per_state + 1L):(s * C_per_state)]))
    )
  }

  # Parameter extraction
  beta <- matrix(optim_params[1:n_beta], nrow = N, byrow = TRUE)
  rho_matrix <- matrix(optim_params[rho_indices], nrow = N, byrow = TRUE)

  # Covariance matrix reconstruction
  sigma_array <- array(dim = c(K, K, N))
  for (i in 1:N) {
    cor_mat <- diag(K)
    cor_mat[lower.tri(cor_mat)] <- rho_matrix[i, ]
    cor_mat[upper.tri(cor_mat)] <- t(cor_mat)[upper.tri(cor_mat)]
    sigma_array[,,i] <- cor_mat
  }

  # Transition matrix evaluated at the in-sample covariate means. X_optim is the
  # estimation window (X_is when out_of_sample = TRUE, full X otherwise), so this
  # summary does not leak held-out OOS covariates that beta was never fit on.
  avg_X <- colMeans(X_optim)
  P <- matrix(nrow = N, ncol = N)
  if (N == 2) {
    for (i in 1:N) {
      p_ii <- plogis(avg_X %*% beta[i, ])
      P[i, i] <- p_ii
      P[i, -i] <- (1 - p_ii) / (N - 1)
    }
  } else {
    # N=3: softmax
    for (i in 1:N) {
      logits <- c(sapply(seq_len(N - 1L), function(j)
        sum(avg_X * beta[i, ((j - 1L) * p + 1L):(j * p)])), 0)
      exp_l <- exp(logits - max(logits))
      P[i, ] <- exp_l / sum(exp_l)
    }
  }

  ll_oos <- if (out_of_sample)
    -rsdc_likelihood(optim_params, y = residuals_oos, exog = X_oos, K = K, N = N)
  else NULL

  list(
    transition_matrix  = P,
    means              = matrix(0, N, K),
    volatilities       = matrix(1, N, K),
    correlations       = rho_matrix,
    covariances        = sigma_array,
    log_likelihood     = -result_optim$value,
    log_likelihood_oos = ll_oos,
    beta               = beta
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
#’   \item{log_likelihood}{In-sample log-likelihood (numeric scalar); full-sample when \code{out_of_sample = FALSE}.}
#’   \item{log_likelihood_oos}{OOS log-likelihood evaluated on the held-out 30\% (numeric scalar),
#’     or \code{NULL} when \code{out_of_sample = FALSE}.}
#’   \item{beta}{\code{NULL} (no covariates in this specification).}
#’ }
#’
#’ @details
#’ \itemize{
#’   \item \strong{Parameterization:} For \eqn{N=2}, the transition parameters are
#’         \eqn{\{p_{11}, p_{22}\}}; off-diagonals are \eqn{1 - p_{ii}}.
#’   \item \strong{Bounds:} \eqn{p_{ii} \in [0.01, 0.99]} and \eqn{\rho \in (-1, 1)}
#’     for numerical stability and identifiability.
#’   \item \strong{Correlation build:} Each regime’s correlation matrix is reconstructed from
#’     its lower-triangular vector and symmetrized; non-PD proposals are penalized in the
#’     likelihood routine.
#’   \item \strong{State ordering:} regimes are reordered by ascending mean correlation for identifiability.
#’ }
#'
#' @section Assumptions:
#' Inputs in \code{residuals} are treated as mean-zero with unit variance; only the correlation
#' structure and fixed transition probabilities are estimated.
#'
#' @seealso \code{\link{rsdc_estimate}} (wrapper), \code{\link{rsdc_hamilton}}, \code{\link{rsdc_likelihood}}
#'
#' @references
#' \insertRef{DEoptim-JSS}{RSDC}
#'
#' \insertRef{hamilton1989}{RSDC}
#'
#' @note Uses \code{DEoptim::DEoptim()} then \code{stats::optim(method = "L-BFGS-B")}.
#'
#' @importFrom stats optim
#' @noRd
f_optim_noX <- function(N, residuals, out_of_sample = FALSE, control = list()) {
  if (N > 3) stop("Only N = 2 or N = 3 is supported.")
  stopifnot(is.matrix(residuals))

  con <- list(seed = 123, do_trace = FALSE)
  con[names(control)] <- control

  K <- ncol(residuals)
  n_obs <- nrow(residuals)
  n_trans <- N * (N - 1)
  C_per_state <- K * (K - 1) / 2
  n_rho <- N * C_per_state

  # Each free transition prob in [0.01, 0.99] for all N. For N=3 the implied third
  # prob (1 - p_i1 - p_i2) can be negative for some proposals; rsdc_likelihood rejects
  # those rows (returns 1e10). This keeps persistent diagonals (p_ii up to 0.99)
  # reachable instead of capping them at 0.485.
  trans_upper <- 0.99
  bounds <- list(
    lower = c(rep(0.01, n_trans), rep(-1, n_rho)),
    upper = c(rep(trans_upper, n_trans), rep(1, n_rho))
  )

  if (out_of_sample) {
    cut     <- round(0.7 * n_obs)
    y_optim <- residuals[1:cut, , drop = FALSE]
    y_oos   <- residuals[(cut + 1):n_obs, , drop = FALSE]
  } else {
    y_optim <- residuals
    y_oos   <- NULL
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

  # Regime reordering for identifiability (ascending mean correlation)
  trans_per_regime <- N - 1L
  rho_params <- final_par[(n_trans + 1):(n_trans + n_rho)]
  corr_means <- sapply(1:N, function(i) mean(rho_params[(1:C_per_state) + (i - 1L) * C_per_state]))
  sorted_states <- order(corr_means)
  if (!all(sorted_states == 1:N)) {
    orig_trans <- final_par[1:n_trans]
    if (N == 2L) {
      trans_reordered <- unlist(lapply(sorted_states, function(s)
        orig_trans[((s - 1L) * trans_per_regime + 1L):(s * trans_per_regime)]))
    } else {
      # N=3: reconstruct full P row, then select columns in new state order
      trans_reordered <- unlist(lapply(sorted_states, function(s) {
        free     <- orig_trans[((s - 1L) * trans_per_regime + 1L):(s * trans_per_regime)]
        full_row <- c(free, 1 - sum(free))
        full_row[sorted_states][seq_len(N - 1L)]
      }))
    }
    final_par <- c(
      trans_reordered,
      unlist(lapply(sorted_states, function(s)
        rho_params[((s - 1L) * C_per_state + 1L):(s * C_per_state)]))
    )
  }

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
  } else if (N == 3) {
    for (i in 1:N) {
      p_i1 <- trans_params[2L * (i - 1L) + 1L]
      p_i2 <- trans_params[2L * (i - 1L) + 2L]
      P[i, ] <- c(p_i1, p_i2, 1 - p_i1 - p_i2)
    }
  }

  ll_oos <- if (!is.null(y_oos))
    -rsdc_likelihood(final_par, y = y_oos, exog = NULL, K = K, N = N)
  else NULL

  list(
    transition_matrix  = P,
    means              = matrix(0, N, K),
    volatilities       = matrix(1, N, K),
    correlations       = rho_matrix,
    covariances        = sigma_array,
    log_likelihood     = -optim_result$value,
    log_likelihood_oos = ll_oos,
    beta               = NULL
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
#'   \item{log_likelihood}{In-sample log-likelihood (numeric scalar); full-sample when \code{out_of_sample = FALSE}.}
#'   \item{log_likelihood_oos}{OOS log-likelihood evaluated on the held-out 30\% (numeric scalar),
#'     or \code{NULL} when \code{out_of_sample = FALSE}.}
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
#' @seealso \code{\link{rsdc_estimate}} (wrapper).
#'
#' @references
#' \insertRef{DEoptim-JSS}{RSDC}
#'
#' \insertRef{hamilton1989}{RSDC}
#'
#' \insertRef{pelletier2006regime}{RSDC}
#'
#' @note Uses \code{DEoptim::DEoptim()} then \code{stats::optim(method = "L-BFGS-B")}.
#'
#' @importFrom stats optim
#' @noRd
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

    R_rdg  <- R + diag(1e-8, K)            # ridge used consistently for inverse and log-det
    inv_R  <- solve(R_rdg)
    logdet <- determinant(R_rdg, logarithm = TRUE)$modulus[1]
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

  # IS log-likelihood (full-sample when out_of_sample = FALSE)
  ll_return <- -optim_result$value

  # OOS log-likelihood: evaluate IS-fitted params on held-out data
  ll_oos <- if (!is.null(y_oos))
    -neg_loglik_const(rho_vec, y = y_oos, K = K)
  else NULL

  list(
    transition_matrix  = matrix(1, 1, 1),
    means              = matrix(0, 1, K),
    volatilities       = matrix(1, 1, K),
    correlations       = matrix(rho_vec, nrow = 1),
    covariances        = array_R,
    log_likelihood     = ll_return,
    log_likelihood_oos = ll_oos,
    beta               = NULL
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
#'   \item{\code{log_likelihood}}{In-sample log-likelihood; full-sample when \code{out_of_sample = FALSE}.}
#'   \item{\code{log_likelihood_oos}}{OOS log-likelihood on held-out 30%, or \code{NULL} when \code{out_of_sample = FALSE}.}
#'   \item{\code{beta}}{TVTP coefficients (only for \code{"tvtp"}).}
#' }
#'
#' @details
#' \itemize{
#'   \item \strong{Method selection:} \code{match.arg()} validates \code{method}.
#'   \item \strong{Inputs:} \code{"tvtp"} requires non-NULL \code{X}; \code{N} is ignored for \code{"const"}.
#'   \item \strong{Split:} If \code{out_of_sample = TRUE}, the first 70% is used for fitting.
#' }
#'
#' @references
#' \insertRef{DEoptim-JSS}{RSDC} \cr
#'
#' \insertRef{hamilton1989}{RSDC} \cr
#'
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
#' @seealso \code{\link{rsdc_hamilton}} and \code{\link{rsdc_likelihood}}.
#'
#' @export
rsdc_estimate <- function(method = c("tvtp", "noX", "const"),
                          residuals, N = 2, X = NULL,
                          out_of_sample = FALSE, control = list()) {
  method <- match.arg(method)
  if (!is.null(N) && N < 2) stop("N must be at least 2. Use method='const' for a single-regime model.")

  # Correlation-only model: columns are assumed mean-zero, unit-variance. Warn (do not
  # auto-transform) when inputs deviate materially so raw returns are not silently misspecified.
  if (is.matrix(residuals)) {
    col_sd <- apply(residuals, 2, stats::sd)
    if (any(is.finite(col_sd) & (col_sd < 0.5 | col_sd > 2)))
      warning("residuals columns are not ~unit-variance (sd outside [0.5, 2]); the model treats ",
              "inputs as standardized. Consider scale() so only the correlation structure is fit.")
  }

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
