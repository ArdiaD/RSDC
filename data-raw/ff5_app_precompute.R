# Precompute the heavy results of the ff5-application vignette.
# The vignette displays this pipeline verbatim (eval = FALSE) and renders
# tables/figures from the slim ff5_app_precomputed.rds written here, so the
# package check stays fast while every number in the vignette is reproducible
# by sourcing this script (~1 h on 10 cores).
#
# Run from the package root: source("data-raw/ff5_app_precompute.R")

suppressMessages({ library(RSDC); library(rugarch) })
stopifnot(packageVersion("RSDC") >= "1.7.0")
CORES <- 14L

## ---- 1. Data and GARCH standardization ---------------------------------------
data("ff5ind", package = "RSDC")
industries <- c("Manuf", "Enrgy", "HiTec", "Hlth", "Utils")
ret <- as.matrix(ff5ind[, industries])

spec <- ugarchspec(
  variance.model = list(model = "gjrGARCH", garchOrder = c(1, 1)),
  mean.model     = list(armaOrder = c(0, 0), include.mean = FALSE))
z <- s <- matrix(NA_real_, nrow(ret), ncol(ret),
                 dimnames = list(NULL, industries))
for (i in seq_along(industries)) {
  g <- ugarchfit(spec, ret[, i], solver = "hybrid")
  z[, i] <- as.numeric(residuals(g, standardize = TRUE))
  s[, i] <- as.numeric(sigma(g))
}
for (i in seq_along(industries)) z[, i] <- z[, i] / sd(z[, i])

X_mccc <- cbind(intercept = 1, MCCC   = as.numeric(scale(ff5ind$mccc)))
X_vix  <- cbind(intercept = 1, logVIX = as.numeric(scale(log(ff5ind$vix))))
X_both <- cbind(X_mccc, logVIX = X_vix[, "logVIX"])

## ---- 2. Estimation: all specifications, multi-start protocol ------------------
ctrl <- list(compute_se = FALSE, n_starts = 4, cores = CORES)
fit1 <- function(label, method, N, X = NULL) {
  t0 <- Sys.time()
  f <- rsdc_estimate(method, residuals = z, N = N, X = X, control = ctrl)
  cat(sprintf("%-11s logLik %10.1f | npar %3d | spread %5.1f | %4.1f min\n",
              label, as.numeric(logLik(f)), f$npar,
              diff(range(f$start_logliks)),
              round(as.numeric(Sys.time() - t0, units = "mins"), 1)))
  f
}
fits <- list(
  "const"      = fit1("const",      "const", 1),
  "noX(N=2)"   = fit1("noX(N=2)",   "noX",   2),
  "noX(N=3)"   = fit1("noX(N=3)",   "noX",   3),
  "tvtp2-MCCC" = fit1("tvtp2-MCCC", "tvtp",  2, X_mccc),
  "tvtp2-VIX"  = fit1("tvtp2-VIX",  "tvtp",  2, X_vix),
  "tvtp2-both" = fit1("tvtp2-both", "tvtp",  2, X_both),
  "tvtp3-MCCC" = fit1("tvtp3-MCCC", "tvtp",  3, X_mccc),
  "tvtp3-VIX"  = fit1("tvtp3-VIX",  "tvtp",  3, X_vix),
  "tvtp3-both" = fit1("tvtp3-both", "tvtp",  3, X_both))

comparison <- data.frame(
  model  = names(fits),
  logLik = sapply(fits, function(f) as.numeric(logLik(f))),
  npar   = sapply(fits, function(f) f$npar),
  AIC    = sapply(fits, AIC),
  BIC    = sapply(fits, BIC),
  row.names = NULL)
print(comparison, digits = 8)
retained <- fits[["tvtp3-both"]]

## ---- 3. Smoothed correlation path (retained model) ----------------------------
fc <- rsdc_forecast("tvtp", N = 3, residuals = z, X = X_both,
                    final_params = retained, sigma_matrix = s,
                    value_cols = industries)
mean_corr <- rowMeans(fc$predicted_correlations)
vit <- rsdc_viterbi(retained)

## ---- 4. Portfolios: regime-based MV / MD vs equal weight ----------------------
mv <- rsdc_minvar(sigma_matrix = s, value_cols = industries,
                  predicted_corr = fc$predicted_correlations, y = ret,
                  long_only = TRUE, lag = TRUE)
md <- rsdc_maxdiv(sigma_matrix = s, value_cols = industries,
                  predicted_corr = fc$predicted_correlations, y = ret,
                  long_only = TRUE, lag = TRUE)
r_mv <- rowSums(rbind(NA, mv$weights[-nrow(mv$weights), ]) * ret)
r_md <- rowSums(rbind(NA, md$weights[-nrow(md$weights), ]) * ret)
r_ew <- rowMeans(ret)
pf <- data.frame(DATE = ff5ind$DATE, minvar = r_mv, maxdiv = r_md, equal = r_ew)
pf_stats <- sapply(pf[-1, -1], function(r)
  c(ann_vol = sd(r) * sqrt(252), ann_ret = mean(r) * 252))

## ---- 5. Bootstrap inference on the retained model -----------------------------
t0 <- Sys.time()
bs <- rsdc_bootstrap(retained, B = 199, X = X_both, seed = 123, cores = CORES)
cat(sprintf("bootstrap: %d effective replicates | %.1f min\n",
            bs$B, as.numeric(Sys.time() - t0, units = "mins")))

## ---- 6. Slim export for the vignette -------------------------------------------
out <- list(
  built        = Sys.time(),
  comparison   = comparison,
  logliks      = setNames(comparison$logLik, comparison$model),
  npars        = setNames(comparison$npar,  comparison$model),
  dates        = ff5ind$DATE,
  mean_corr    = mean_corr,
  smoothed     = retained$smoothed_probs,      # N x T
  viterbi      = vit,
  regime_means = rowMeans(retained$correlations),
  correlations = retained$correlations,        # N x C regime correlations
  par          = retained$par,                 # warm-start anchor for patches
  beta         = retained$beta,
  P_bar        = retained$transition_matrix,
  start_logliks = lapply(fits, function(f) f$start_logliks),
  pf           = pf,
  pf_stats     = pf_stats,
  weights_mv   = mv$weights,
  weights_md   = md$weights,
  mv_fallback_share = mean(apply(abs(mv$weights - 1 / length(industries)) < 1e-12,
                                 1, all)),
  boot_ci      = bs$ci,
  boot_se      = bs$se,
  boot_B       = bs$B)
saveRDS(out, "vignettes/ff5_app_precomputed.rds", compress = "xz")
cat("vignettes/ff5_app_precomputed.rds:",
    round(file.size("vignettes/ff5_app_precomputed.rds") / 1024), "KB\n")
