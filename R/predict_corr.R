#' Forecast Covariance Matrices from Regime-Switching Model
#'
#' Computes predicted correlation and covariance matrices from a fitted
#' regime-switching model. Supports constant, fixed-transition, and
#' time-varying transition probability (TVTP) specifications.
#'
#' @param method Character. One of \code{"tvtp"} (time-varying transition probabilities),
#'   \code{"noX"} (fixed transition matrix), or \code{"const"} (constant correlation).
#' @param N Integer. Number of regimes (ignored if \code{method = "const"}).
#' @param residuals Numeric matrix of dimension \eqn{T \times K}, where \eqn{T} is the number
#'   of time steps and \eqn{K} is the number of assets.
#' @param X Optional matrix of exogenous variables, of dimension \eqn{T \times p}. Required if \code{method = "tvtp"}.
#' @param final_params A list of model parameters, typically returned by \code{\link{estimate_model}}.
#'   It must include either:
#'   \itemize{
#'     \item \code{beta} and \code{correlations} (for \code{"tvtp"}),
#'     \item or \code{transition_matrix} and \code{correlations} (for \code{"noX"}),
#'     \item or just correlation info (for \code{"const"}).
#'   }
#' @param sigma_matrix Matrix of predicted standard deviations (T × K).
#' @param value_cols Character or integer vector: names or indices of the volatility columns in \code{sigma_matrix}.
#' @param out_of_sample Logical. If \code{TRUE}, use only the out-of-sample portion (last 30%) for forecasting.
#'
#' @return A named list with components:
#' \describe{
#'   \item{smoothed_probs}{Matrix (N × T) of smoothed regime probabilities (not available if \code{method = "const"}).}
#'   \item{sigma_matrix}{Filtered matrix of standard deviations used for forecasting.}
#'   \item{cov_matrices}{List of predicted covariance matrices for each time point (length T).}
#'   \item{predicted_correlations}{Matrix (T × P) of predicted pairwise correlations, with \eqn{P = K(K-1)/2}.}
#'   \item{BIC}{Bayesian Information Criterion value of the fitted model.}
#'   \item{y}{Matrix of residuals used for forecasting (in-sample or out-of-sample).}
#' }
#'
#' @seealso \code{\link{rsdc_hamilton}}, \code{\link{estimate_model}}, \code{\link{f_minvar}}, \code{\link{f_maxdiv}}
#'
#' @export
rsdc_forecast <- function(method = c("tvtp", "noX", "const"),
                         N,
                         residuals,
                         X = NULL,
                         final_params,
                         sigma_matrix,
                         value_cols,
                         out_of_sample = FALSE) {


  ### !!! METTRE UNE POSSIBILITÉ DE CHOISIR LES NOMS QU'ON VEUT POUR FINAL PARAMS !!! ###

  method <- match.arg(method)
  K <- ncol(residuals)
  n_obs <- nrow(residuals)

  if (out_of_sample) {
    is_cut <- round(0.7 * n_obs)
    residuals_is <- residuals[1:is_cut, ]
    residuals_oos <- residuals[(is_cut + 1):n_obs, ]
    sigma_matrix_oos <- sigma_matrix[(is_cut + 1):n_obs, ]
    X_is <- if (!is.null(X)) X[1:is_cut, ] else NULL
    X_oos <- if (!is.null(X)) X[(is_cut + 1):n_obs, ] else NULL
  } else {
    residuals_oos <- residuals
    sigma_matrix_oos <- sigma_matrix
    X_oos <- X
  }

  # Default outputs
  smoothed_probs <- NULL
  predicted_correl <- NULL

  if (method == "const") {
    # Use constant empirical correlation
    corr_matrix <- cor(residuals)
    predicted_correl <- matrix(rep(corr_matrix[lower.tri(corr_matrix)], nrow(residuals_oos)),
                               ncol = K * (K - 1) / 2, byrow = TRUE)
  } else {
    # Use Hamilton filter to get smoothed_probs
    result_filter <- RSDC::rsdc_hamilton(
      y = if (out_of_sample) residuals_is else residuals_oos,
      X = if (method == "tvtp") if (out_of_sample) X_is else X_oos else NULL,
      beta = if (method == "tvtp") final_params$beta else NULL,
      rho_matrix = final_params$correlations,
      K = K,
      N = N,
      P = if (method == "noX") final_params$transition_matrix else NULL
    )

    smoothed_probs <- result_filter$smoothed_probs

    # Weighted correlation forecast
    cor_pairs <- combn(K, 2, simplify = FALSE)
    predicted_correl <- matrix(0, nrow = nrow(residuals_oos), ncol = length(cor_pairs))

    for (t in 1:nrow(predicted_correl)) {
      for (s in 1:N) {
        predicted_correl[t, ] <- predicted_correl[t, ] +
          smoothed_probs[s, t] * final_params$correlations[s, ]
      }
    }
  }

  # Covariance matrix construction
  cor_pairs <- combn(K, 2, simplify = FALSE)
  cov_matrices <- lapply(1:nrow(sigma_matrix_oos), function(t) {
    sigma_vals <- as.numeric(sigma_matrix_oos[t, value_cols])
    R <- diag(K)

    for (i in seq_along(cor_pairs)) {
      pair <- cor_pairs[[i]]
      R[pair[1], pair[2]] <- R[pair[2], pair[1]] <- predicted_correl[t, i]
    }

    D <- diag(sigma_vals)
    D %*% R %*% D
  })

  # BIC calculation
  k <- switch(
    method,
    "tvtp" = N * ncol(X) + N * K * (K - 1) / 2,
    "noX"  = N * (N - 1) + N * K * (K - 1) / 2,
    "const" = K * (K - 1) / 2
  )
  n <- if (out_of_sample) nrow(residuals_oos) else n_obs
  BIC <- log(n) * k - 2 * final_params$log_likelihood

  list(
    smoothed_probs = smoothed_probs,
    sigma_matrix = sigma_matrix_oos,
    cov_matrices = cov_matrices,
    predicted_correlations = predicted_correl,
    BIC = BIC,
    y = residuals_oos
  )
}
