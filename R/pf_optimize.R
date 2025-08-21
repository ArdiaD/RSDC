#' Minimum-Variance Portfolio (Rolling Weights)
#'
#' Computes rolling minimum-variance (MV) portfolio weights from a sequence of
#' per-period covariance matrices implied by forecasted volatilities and
#' pairwise correlations. Supports long-only or unconstrained MV. If the QP
#' solver fails at a time step, the routine falls back to equal weights.
#'
#' @param sigma_matrix Numeric matrix \eqn{T \times K} of forecasted \emph{volatilities}
#'   (standard deviations), one column per asset.
#' @param value_cols Character or integer vector giving the columns in \code{sigma_matrix}
#'   to use as assets (order defines the asset order).
#' @param predicted_corr Numeric matrix \eqn{T \times P} of pairwise correlations, where
#'   \eqn{P = \binom{K}{2}} and the columns correspond to \code{combn(K, 2)} order.
#' @param y Numeric matrix \eqn{T \times K} of asset returns aligned with \code{sigma_matrix}.
#'   Used only to compute the realized portfolio volatility.
#' @param long_only Logical. If \code{TRUE} (default), imposes long-only MV with the full-investment
#'   constraint \eqn{\sum_i w_i = 1} and \eqn{w_i \ge 0}. If \code{FALSE}, solves unconstrained MV
#'   with only \eqn{\sum_i w_i = 1}.
#'
#' @returns An object of class \code{"minvar_portfolio"}:
#' \describe{
#'   \item{weights}{\eqn{T \times K} matrix of MV weights (one row per time).}
#'   \item{cov_matrices}{List of length \eqn{T} with the per-period \eqn{K \times K} covariance matrices.}
#'   \item{volatility}{Realized standard deviation of portfolio returns (see Note on units).}
#'   \item{y}{The input \code{y} matrix (coerced to \eqn{T \times K}).}
#'   \item{K}{Number of assets.}
#' }
#'
#' @details
#' \itemize{
#'   \item \strong{Covariance build:} For each \eqn{t}, a correlation matrix \eqn{R_t}
#'         is reconstructed ... Let \eqn{D_t = \mathrm{diag}(\sigma_{t,1},\dots,\sigma_{t,K})}
#'         and \eqn{\Sigma_t = D_t R_t D_t}.
#'   \item \strong{Optimization:} Minimize \eqn{\tfrac{1}{2} w^\top \Sigma_t w} subject to
#'         \eqn{\mathbf{1}^\top w = 1} and, if \code{long_only}, \eqn{w \ge 0}
#'         (solved with \code{quadprog::solve.QP}).
#'   \item \strong{Failure handling:} If the QP fails at time \(t\), weights default to equal
#'         allocation \(w_i = 1/K\).
#' }
#'
#' @section Note on units:
#' The realized portfolio return at time \(t\) is computed as \code{sum(y[t, ] * weights[t, ]) / 100}.
#' This assumes \code{y} is expressed in \% (percentage) units. If your \code{y} is already in decimals,
#' remove the \code{/ 100} in the implementation or convert inputs accordingly.
#'
#' @examples
#' \dontrun{
#' # Toy example with K = 3
#' T <- 50; K <- 3
#' set.seed(42)
#' vols <- matrix(0.2 + 0.05*abs(sin(seq_len(T)/7)), T, K)
#' colnames(vols) <- paste0("A", 1:K)
#' # simple, stationary correlations
#' pred_corr <- cbind(rep(0.20, T), rep(0.10, T), rep(0.05, T))  # order: (2,1), (3,1), (3,2)
#' rets <- matrix(rnorm(T*K, sd = 0.01), T, K); colnames(rets) <- colnames(vols)
#'
#' mv <- rsdc_minvar(sigma_matrix  = vols,
#'                   value_cols    = colnames(vols),
#'                   predicted_corr= pred_corr,
#'                   y             = rets,
#'                   long_only     = TRUE)
#' head(mv$weights)
#' mv$volatility
#' }
#'
#' @seealso \code{\link{rsdc_maxdiv}} (maximum diversification),
#'   \code{\link[quadprog]{solve.QP}}
#'
#' @importFrom utils combn
#' @importFrom stats sd
#' @export
rsdc_minvar <- function(sigma_matrix, value_cols, predicted_corr, y,
                     long_only = TRUE) {

  if (!requireNamespace("quadprog", quietly = TRUE)) {
    stop("Package 'quadprog' is required. Install via install.packages('quadprog').")
  }

  K <- length(value_cols)
  y <- matrix(as.numeric(y), ncol = K)
  colnames(y) <- value_cols

  stopifnot(
    ncol(predicted_corr) == choose(K, 2),
    ncol(y) == K,
    ncol(sigma_matrix) >= K,
    nrow(sigma_matrix) == nrow(predicted_corr),
    nrow(sigma_matrix) == nrow(y)
  )

  T <- nrow(y)
  cor_pairs <- combn(K, 2, simplify = FALSE)
  portfolio_weights <- matrix(NA, nrow = T, ncol = K)
  portfolio_returns <- numeric(T)
  cov_matrices <- vector("list", T)

  for (t in 1:T) {
    R <- diag(K)
    for (i in seq_along(cor_pairs)) {
      p <- cor_pairs[[i]]
      R[p[1], p[2]] <- R[p[2], p[1]] <- predicted_corr[t, i]
    }

    D <- diag(as.numeric(sigma_matrix[t, value_cols]))
    cov_t <- D %*% R %*% D
    cov_matrices[[t]] <- cov_t

    tryCatch({
      Dmat <- 2 * cov_t
      dvec <- rep(0, K)

      if (long_only) {
        Amat <- cbind(rep(1, K), diag(K))
        bvec <- c(1, rep(0, K))
      } else {
        Amat <- matrix(rep(1, K), ncol = 1)
        bvec <- 1
      }

      qp <- quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq = 1)
      weights <- qp$solution
    }, error = function(e) {
      weights <<- rep(1 / K, K)
    })

    portfolio_weights[t, ] <- weights
    portfolio_returns[t] <- sum(y[t, ] * weights) / 100
  }

  structure(list(
    cov_matrices = cov_matrices,
    weights = portfolio_weights,
    volatility = sd(portfolio_returns, na.rm = TRUE),
    y = y,
    K = K
  ), class = "minvar_portfolio")
}

#' Maximum-Diversification Portfolio (Rolling Weights)
#'
#' Computes rolling maximum-diversification (MaxDiv) portfolio weights from a sequence
#' of per-period covariance matrices implied by forecasted volatilities and correlations.
#' Falls back to equal weights if the nonlinear solver fails.
#'
#' @param sigma_matrix Numeric matrix \eqn{T \times K} of forecasted standard deviations.
#' @param value_cols Character/integer vector naming columns in \code{sigma_matrix} (asset order).
#' @param predicted_corr Numeric matrix \eqn{T \times \binom{K}{2}} of pairwise correlations
#'   in \code{combn(K, 2)} column order.
#' @param y Numeric matrix \eqn{T \times K} of asset returns (for realized stats).
#' @param long_only Logical. If \code{TRUE}, impose \eqn{w \ge 0} and \eqn{\sum_i w_i = 1};
#'   otherwise bounds are \eqn{-1 \le w_i \le 1} with \eqn{\sum_i w_i = 1}.
#'
#' @return
#' \describe{
#'   \item{\code{weights}}{\eqn{T \times K} matrix of weights.}
#'   \item{\code{returns}}{Vector of realized portfolio returns \code{sum(y[t,] * weights[t,])}.}
#'   \item{\code{diversification_ratios}}{Vector of realized diversification ratios.}
#'   \item{\code{mean_diversification}}{Average diversification ratio.}
#'   \item{\code{K}}{Number of assets.}
#'   \item{\code{assets}}{Asset names.}
#'   \item{\code{volatility}}{Standard deviation of realized portfolio returns.}
#' }
#'
#' @details
#' \itemize{
#'   \item \strong{Covariance build:} For each \eqn{t}, reconstruct \eqn{R_t} from the
#'         pairwise vector; set \eqn{D_t=\mathrm{diag}(\sigma_{t,1},\dots,\sigma_{t,K})} and
#'         \eqn{\Sigma_t = D_t R_t D_t}.
#'   \item \strong{Objective (MaxDiv):} maximize
#'         \eqn{\mathrm{DR}(w) = \frac{\sum_i w_i \sigma_{t,i}}{\sqrt{w^\top \Sigma_t w}}}
#'         subject to \eqn{\sum_i w_i = 1} and bounds on \eqn{w}. Implemented by minimizing
#'         the negative ratio.
#'   \item \strong{Solver:} \code{Rsolnp::solnp} with equality \eqn{\sum_i w_i = 1} and
#'         bounds by \code{long_only}; on error, weights default to \eqn{1/K}.
#' }
#'
#' @seealso \code{\link{rsdc_minvar}}, \code{\link[Rsolnp]{solnp}}
#' @importFrom utils combn
#' @importFrom stats sd
#'
#' @examples
#' \dontrun{
#' # Toy example with K = 3
#' if (requireNamespace("Rsolnp", quietly = TRUE)) {
#'   T <- 50; K <- 3
#'   set.seed(42)
#'   vols <- matrix(0.2 + 0.05*abs(sin(seq_len(T)/7)), T, K)
#'   colnames(vols) <- paste0("A", 1:K)
#'   # simple, stationary correlations (order: (2,1), (3,1), (3,2))
#'   pred_corr <- cbind(rep(0.20, T), rep(0.10, T), rep(0.05, T))
#'   rets <- matrix(rnorm(T*K, sd = 0.01), T, K); colnames(rets) <- colnames(vols)
#'
#'   mx <- rsdc_maxdiv(sigma_matrix   = vols,
#'                     value_cols     = colnames(vols),
#'                     predicted_corr = pred_corr,
#'                     y              = rets,
#'                     long_only      = TRUE)
#'   head(mx$weights)
#'   mx$mean_diversification
#' }
#' }
#'
#' @export
rsdc_maxdiv <- function(sigma_matrix, value_cols, predicted_corr, y,
                     long_only = TRUE) {
  if (!requireNamespace("Rsolnp", quietly = TRUE)) {
    stop("Package 'Rsolnp' is required. Install via install.packages('Rsolnp').")
  }

  K <- length(value_cols)
  y <- matrix(as.numeric(y), ncol = K)
  colnames(y) <- value_cols

  stopifnot(
    ncol(predicted_corr) == choose(K, 2),
    ncol(y) == K,
    ncol(sigma_matrix) >= K,
    nrow(sigma_matrix) == nrow(predicted_corr),
    nrow(sigma_matrix) == nrow(y)
  )

  cor_pairs <- combn(K, 2, simplify = FALSE)
  cov_matrices <- vector("list", nrow(sigma_matrix))
  portfolio_weights <- matrix(NA, nrow = nrow(sigma_matrix), ncol = K)
  portfolio_returns <- numeric(nrow(sigma_matrix))
  diversification_ratios <- numeric(nrow(sigma_matrix))

  LB <- if (long_only) rep(0, K) else rep(-1, K)
  UB <- rep(1, K)

  for (t in 1:nrow(sigma_matrix)) {
    R <- diag(K)
    for (i in seq_along(cor_pairs)) {
      p <- cor_pairs[[i]]
      R[p[1], p[2]] <- R[p[2], p[1]] <- predicted_corr[t, i]
    }

    sigma_vals <- as.numeric(sigma_matrix[t, value_cols])
    cov_t <- diag(sigma_vals) %*% R %*% diag(sigma_vals)
    cov_matrices[[t]] <- cov_t

    vol_vec <- sigma_vals

    div_obj <- function(w) {
      port_vol <- sqrt(t(w) %*% cov_t %*% w)
      weighted_avg_vol <- sum(w * vol_vec)
      - weighted_avg_vol / port_vol
    }

    tryCatch({
      sol <- Rsolnp::solnp(
        pars = rep(1 / K, K),
        fun = div_obj,
        eqfun = function(w) sum(w),
        eqB = 1,
        LB = LB,
        UB = UB,
        control = list(trace = 0, tol = 1e-8)
      )
      weights <- sol$pars
    }, error = function(e) {
      weights <<- rep(1 / K, K)
    })

    portfolio_weights[t, ] <- weights
    portfolio_returns[t] <- sum(y[t, ] * weights)

    dr <- sum(weights * vol_vec) / sqrt(t(weights) %*% cov_t %*% weights)
    diversification_ratios[t] <- ifelse(is.nan(dr), 1, dr)
  }

  structure(list(
    weights = portfolio_weights,
    returns = portfolio_returns,
    diversification_ratios = diversification_ratios,
    mean_diversification = mean(diversification_ratios, na.rm = TRUE),
    K = K,
    assets = value_cols,
    volatility = sd(portfolio_returns)
  ), class = "maxdiv_portfolio")
}
