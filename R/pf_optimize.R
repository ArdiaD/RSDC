#' Minimum Variance Portfolio Weights
#'
#' Computes dynamic minimum-variance portfolio weights given a sequence of covariance matrices,
#' with optional rebalancing frequency for volatility scaling. Rebalancing is assumed at every row unless
#' a custom list of rebalance dates is provided.
#'
#' @param sigma_matrix Matrix (T × K) of forecasted volatilities.
#' @param value_cols Character or integer vector identifying volatility columns in \code{sigma_matrix}.
#' @param predicted_corr Matrix (T × P) of predicted pairwise correlations.
#' @param y Matrix (T × K) of asset returns.
#' @param rebalance Character. Rebalancing frequency used for volatility scaling. One of \code{"daily"}, \code{"monthly"}, \code{"quarterly"}, \code{"annually"}, or \code{"custom"}. Default is \code{"daily"}.
#' @param dates Optional vector of custom rebalancing dates (only used if \code{rebalance = "custom"}).
#' @param long_only Logical. If TRUE, applies long-only constraints (default: TRUE).
#'
#' @return A list of class \code{"minvar_portfolio"} containing:
#' \describe{
#'   \item{weights}{T × K matrix of portfolio weights.}
#'   \item{cov_matrices}{List of covariance matrices at each time $t$.}
#'   \item{volatility}{Annualized standard deviation of portfolio returns.}
#'   \item{y}{Input return matrix.}
#'   \item{K}{Number of assets.}
#' }
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

#' Maximum Diversification Portfolio Weights
#'
#' Computes dynamic max-diversification portfolio weights based on forecasted covariances.
#'
#' @param sigma_matrix Matrix (T × K) of forecasted volatilities.
#' @param value_cols Character or integer vector identifying volatility columns in \code{sigma_matrix}.
#' @param predicted_corr Matrix (T × P) of predicted pairwise correlations.
#' @param y Matrix (T × K) of asset returns.
#' @param long_only Logical. If TRUE, applies long-only constraints (default: TRUE).
#' @param rebalance Character. Rebalancing frequency used for volatility scaling.
#'   One of \code{"daily"}, \code{"monthly"}, \code{"quarterly"}, \code{"yearly"}, or \code{"custom"}.
#' @param dates Optional. Vector of dates (class \code{Date}) required only if \code{rebalance = "custom"}.
#'
#' @return A list of class \code{"maxdiv_portfolio"} containing:
#' \describe{
#'   \item{weights}{T × K matrix of portfolio weights.}
#'   \item{returns}{Vector of portfolio returns.}
#'   \item{diversification_ratios}{Vector of realized diversification ratios.}
#'   \item{mean_diversification}{Mean diversification ratio.}
#'   \item{K}{Number of assets.}
#'   \item{assets}{Asset names.}
#' }
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
