#' Forecast Covariance/Correlation Paths from an RSDC Model
#'
#' Generates per-period correlation and covariance matrices from a fitted model:
#' \code{"const"} (constant correlation), \code{"noX"} (fixed transition matrix), or
#' \code{"tvtp"} (time-varying transition probabilities).
#'
#' @param method Character. One of \code{"tvtp"}, \code{"noX"}, \code{"const"}.
#' @param N Integer. Number of regimes (ignored for \code{"const"}).
#' @param residuals Numeric matrix \eqn{T \times K} used to compute correlations or run the filter.
#' @param X Optional numeric matrix \eqn{T \times p} (required for \code{"tvtp"}).
#' @param final_params List with fitted parameters (e.g., from \code{\link{rsdc_estimate}}):
#'   must include \code{correlations}, and either \code{transition_matrix} (\code{"noX"})
#'   or \code{beta} (\code{"tvtp"}); include \code{log_likelihood} for BIC computation.
#' @param sigma_matrix Numeric matrix \eqn{T \times K} of forecasted standard deviations.
#' @param value_cols Character/integer vector of columns in \code{sigma_matrix} that define asset order.
#' @param out_of_sample Logical. If \code{TRUE}, use a fixed 70/30 split; otherwise use the whole sample.
#' @param control Optional list; supports \code{threshold} (in (0,1), default \code{0.7}).
#'
#' @return
#' \describe{
#'   \item{\code{smoothed_probs}}{\eqn{N \times T^\ast} smoothed probabilities (\code{"noX"}/\code{"tvtp"} only).}
#'   \item{\code{sigma_matrix}}{\eqn{T^\ast \times K} slice aligned to the forecast horizon.}
#'   \item{\code{cov_matrices}}{List of \eqn{K \times K} covariance matrices \eqn{\Sigma_t = D_t R_t D_t}.}
#'   \item{\code{predicted_correlations}}{\eqn{T^\ast \times \binom{K}{2}} pairwise correlations in \code{combn(K, 2)} order.}
#'   \item{\code{BIC}}{Bayesian Information Criterion \eqn{\mathrm{BIC} = \log(n)\,k - 2\,\ell}.}
#'   \item{\code{y}}{\eqn{T^\ast \times K} residual matrix aligned to the forecast horizon.}
#' }
#'
#' @details
#' \itemize{
#'   \item \strong{Forecast horizon:} If \code{out_of_sample = TRUE}, filter on the first \code{threshold}
#'         fraction and forecast on the remainder.
#'   \item \strong{Correlation paths:}
#'         \itemize{
#'           \item \code{"const"} — empirical correlation of \code{residuals}, repeated across time.
#'           \item \code{"noX"}/\code{"tvtp"} — smoothed-probability weighted average of regime correlations.
#'         }
#'   \item \strong{Covariance build:} Reconstruct \eqn{R_t} from the pairwise vector (columns ordered by
#'         \code{combn(K, 2)}), set \eqn{D_t = \mathrm{diag}(\sigma_{t,1},\dots,\sigma_{t,K})}, and
#'         \eqn{\Sigma_t = D_t R_t D_t}.
#'   \item \strong{BIC:} Parameter count \eqn{k} is
#'         \code{N * ncol(X) + N * K * (K - 1) / 2} for \code{"tvtp"},
#'         \code{N * (N - 1) + N * K * (K - 1) / 2} for \code{"noX"},
#'         and \code{K * (K - 1) / 2} for \code{"const"}.
#' }
#'
#' @examples
#' set.seed(123)
#' T <- 60; K <- 3; N <- 2
#' y <- scale(matrix(rnorm(T*K), T, K))
#' vols <- matrix(0.2 + 0.05*abs(sin(seq_len(T)/7)), T, K)
#' rho <- rbind(c(0.10, 0.05, 0.00), c(0.60, 0.40, 0.30))
#' Pfix <- matrix(c(0.9, 0.1, 0.2, 0.8), 2, 2, byrow = TRUE)
#' rsdc_forecast("noX", N, y, NULL,
#'               list(correlations = rho, transition_matrix = Pfix, log_likelihood = -200),
#'               vols, 1:K)
#'
#' @seealso \code{\link{rsdc_hamilton}}, \code{\link{rsdc_estimate}},
#'   \code{\link{rsdc_minvar}}, \code{\link{rsdc_maxdiv}}
#'
#' @importFrom stats cor
#' @importFrom utils combn
#' @export
rsdc_forecast <- function(method = c("tvtp", "noX", "const"),
                         N,
                         residuals,
                         X = NULL,
                         final_params,
                         sigma_matrix,
                         value_cols,
                         out_of_sample = FALSE,
                         control = list()) {

  con <- list(threshold = 0.7)
  con[names(control)] <- control
  thresh <- con$threshold

  method <- match.arg(method)
  K <- ncol(residuals)
  n_obs <- nrow(residuals)

  if (out_of_sample) {
    is_cut <- round(thresh * n_obs)
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
