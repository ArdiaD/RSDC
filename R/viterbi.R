#' Most likely regime path (Viterbi decoding)
#'
#' Returns the single most likely sequence of regimes
#' \eqn{(S_1, \dots, S_T)} under a fitted \code{"rsdc_fit"} model, via the Viterbi
#' algorithm. This is the maximum a posteriori \emph{joint} state path and
#' complements the marginal filtered/smoothed regime probabilities from
#' \code{\link{rsdc_forecast}}: the smoothed probabilities are decoded pointwise,
#' whereas Viterbi enforces a globally coherent transition path.
#'
#' @param object An \code{"rsdc_fit"} object from \code{\link{rsdc_estimate}}.
#' @param residuals Optional \eqn{T \times K} residual matrix; defaults to the
#'   matrix stored on the fit.
#' @param X Optional \eqn{T \times p} covariate matrix (\code{"tvtp"} only);
#'   defaults to the matrix stored on the fit.
#'
#' @returns An integer vector of length \eqn{T} giving the most likely regime at
#'   each time, in the model's ascending-correlation regime labeling
#'   (1 = lowest-correlation state).
#'
#' @seealso \code{\link{rsdc_forecast}} for marginal regime probabilities.
#' @examples
#' \donttest{
#' y <- scale(matrix(rnorm(200 * 2), 200, 2))
#' fit <- rsdc_estimate("noX", residuals = y, N = 2)
#' table(rsdc_viterbi(fit))
#' }
#' @export
rsdc_viterbi <- function(object, residuals = NULL, X = NULL) {
  if (!inherits(object, "rsdc_fit")) stop("object must be an 'rsdc_fit'.")
  method <- object$method; N <- object$N; K <- object$K
  if (is.null(residuals)) residuals <- object$residuals
  if (is.null(residuals)) stop("No residuals available; pass 'residuals'.")
  if (is.null(X) && method == "tvtp") X <- object$X
  Tn <- nrow(residuals)
  if (method == "const" || N == 1L) return(rep(1L, Tn))

  rho_matrix <- object$correlations
  beta <- object$beta

  # Per-regime Gaussian log-densities (N x T).
  logdens <- matrix(-Inf, N, Tn)
  for (s in seq_len(N)) {
    R <- diag(K)
    R[lower.tri(R)] <- rho_matrix[s, ]
    R[upper.tri(R)] <- t(R)[upper.tri(R)]
    R <- R + diag(1e-8, K)
    inv <- tryCatch(solve(R), error = function(e) NULL)
    if (is.null(inv)) next
    ld  <- determinant(R, logarithm = TRUE)$modulus[1]
    quad <- rowSums((residuals %*% inv) * residuals)
    logdens[s, ] <- -0.5 * (ld + quad + K * log(2 * pi))
  }

  # Transition matrix into time t (contemporaneous covariate for "tvtp").
  P_t <- function(t) {
    if (method == "tvtp") .rsdc_P_at(beta, X[t, ], N) else object$transition_matrix
  }
  logP <- function(t) {
    lp <- log(P_t(t)); lp[!is.finite(lp)] <- -1e10; lp
  }

  delta <- matrix(-Inf, N, Tn)
  psi   <- matrix(0L,   N, Tn)
  delta[, 1] <- log(1 / N) + logdens[, 1]
  for (t in 2:Tn) {
    lp <- logP(t)                      # lp[i, s] = log Pr(S_t = s | S_{t-1} = i)
    for (s in seq_len(N)) {
      val <- delta[, t - 1] + lp[, s]
      psi[s, t]   <- which.max(val)
      delta[s, t] <- max(val) + logdens[s, t]
    }
  }

  path <- integer(Tn)
  path[Tn] <- which.max(delta[, Tn])
  for (t in seq.int(Tn - 1L, 1L)) path[t] <- psi[path[t + 1L], t + 1L]
  path
}
