# =============================================================================
# Monte Carlo Study — RSDC Package  (parallelized)
# =============================================================================
#
# Five independent cases, K=2:
#   Case 1 — const, N=1
#   Case 2 — noX,   N=2
#   Case 3 — tvtp,  N=2
#   Case 4 — noX,   N=3   
#   Case 5 — tvtp,  N=3 
#
# Design: M=100 replications x T in {500, 1000, 2000}
# Replications are parallelized via parallel::mclapply with dynamic
#   scheduling (mc.preschedule = FALSE) for load balancing across uneven DEoptim
#   runtimes (forking; macOS/Linux. On Windows mclapply runs serially).
# Outputs:
#   - PNG plots          -> inst/simulation/plots/
#   - CSV summaries       -> inst/simulation/results/
#   - RDS results bundle  -> inst/simulation/results/monte_carlo_results.rds
#
# Source from the package root in RStudio:
#   source("inst/simulation/monte_carlo_study.R")
# =============================================================================
# Resolve the package root (the folder containing DESCRIPTION) without hardcoding a
# machine-specific path. rprojroot ships with devtools, so no extra dependency is
# needed; if the root cannot be found we fall back to the current directory.
root <- tryCatch(rprojroot::find_root(rprojroot::has_file("DESCRIPTION")),
                 error = function(e) getwd())
setwd(root)
if (!file.exists("DESCRIPTION"))
  stop("Could not locate the RSDC package root. Set the working directory to the ",
       "folder containing DESCRIPTION, then re-source this script.")

devtools::load_all()

# ── Global configuration ──────────────────────────────────────────────────────
M        <- 100L
T_grid   <- c(500L, 1000L, 2000L)
K_grid   <- c(2L, 3L, 4L)            # cross-section sizes swept in the N=3 cases (4 & 5)
# ── Parallelization ───────────────────────────────────────────────────────────
# Parallelism lives at the REPLICATION level: each (case[, K], T) cell maps its M
# replications across workers via parallel::mclapply (forking; serial on Windows).
# DEoptim runs serially inside every worker (parallelType = 0 in the package), so
# there is no nested parallelism to oversubscribe the CPU. Replications self-seed
# (set.seed(base_seed + m)), so results are bit-identical regardless of the worker
# count or scheduling order — only the wall-clock time changes.
#
# Workers: fixed at 14 (On the 18-core target machine).
MC_CORES <- 14L
# Pin BLAS to a single thread so matrix ops inside the forked workers don't spawn
# their own thread pools on top of the M-way fork (no-op if RhpcBLASctl is absent).
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) RhpcBLASctl::blas_set_num_threads(1L)

plots_dir   <- file.path("inst", "simulation", "plots")
results_dir <- file.path("inst", "simulation", "results")
dir.create(plots_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

saved_files <- character(0)   # collect output paths for the final confirmation

# ── True parameters ───────────────────────────────────────────────────────────
C1_RHO <- 0.5
C1_K   <- 2L

C2_P11  <- 0.9;  C2_P22  <- 0.8
C2_RHO1 <- 0.2;  C2_RHO2 <- 0.7
C2_N    <- 2L;   C2_K    <- 2L
C2_BETA <- matrix(c(qlogis(C2_P11), qlogis(C2_P22)), nrow = 2L, ncol = 1L)

# Lower persistence (P[i,i] ~ 0.70 at x=0) + larger slopes (+/-0.8) so X has a
# detectable, identifiable effect on the transition probabilities.
C3_BETA <- matrix(c(qlogis(0.70), -0.8,
                    qlogis(0.70),  0.8), nrow = 2L, byrow = TRUE)
C3_RHO1 <- 0.2;  C3_RHO2 <- 0.7
C3_N    <- 2L;   C3_K    <- 2L;  C3_PHI <- 0.7
# Representative switching probabilities at the covariate mean (scaled X has mean 0,
# so p_ii = plogis(intercept) = plogis(qlogis(0.70)) = 0.70 for both regimes).
C3_P11  <- plogis(C3_BETA[1L, 1L]);  C3_P22 <- plogis(C3_BETA[2L, 1L])

# ── Case 4 (noX, N=3) true parameters ───────────────────────────────────────
# The simulator consumes transition probabilities through the intercept-only
# softmax beta (reference = last regime). For an intercept-only X (p=1) and a
# target row P[i, ], beta[i, ] = (log(P[i,1]/P[i,3]), log(P[i,2]/P[i,3])).
# Regimes are listed in ASCENDING-rho order so the estimator's relabeling (by
# ascending mean correlation) maps estimate i <-> true i.
# The transition structure is K-invariant; correlations are EQUICORRELATION
# matrices (every pair = C4_RHO[r] in regime r), so the same per-regime rho is
# reused as K is swept (the K-sweep happens in the Case-4 execution block).
C4_N <- 3L
C4_P <- matrix(c(0.90, 0.07, 0.03,
                 0.08, 0.85, 0.07,
                 0.05, 0.15, 0.80), nrow = 3L, byrow = TRUE)
C4_RHO  <- c(0.2, 0.5, 0.8)          # per-regime equicorrelation value
C4_BETA <- t(apply(C4_P, 1L, function(pr) c(log(pr[1] / pr[3]), log(pr[2] / pr[3]))))  # 3 x 2

# ── Case 5 (tvtp, N=3) true parameters ──────────────────────────────────────
# Baseline transition matrix at the covariate mean (x=0), ascending-rho order.
# Intercepts come from the softmax inverse of C5_P0; AR(1)-covariate slopes are
# the even columns of the N x (N-1)*p = 3 x 4 beta. Packing per row:
#   beta[i, ] = (int_logit1, slope_logit1, int_logit2, slope_logit2).
# Correlations are EQUICORRELATION matrices reused across the K-sweep.
C5_N <- 3L;  C5_PHI <- 0.7
C5_P0 <- matrix(c(0.80, 0.15, 0.05,
                  0.10, 0.80, 0.10,
                  0.05, 0.15, 0.80), nrow = 3L, byrow = TRUE)
C5_RHO   <- c(0.2, 0.5, 0.8)          # per-regime equicorrelation value
C5_SLOPE <- matrix(c( 0.5, 0.0,
                      0.0, 0.5,
                     -0.5, 0.0), nrow = 3L, byrow = TRUE)   # slope on scaled X per free logit
C5_BETA  <- t(vapply(1:3, function(i) {
  ints <- c(log(C5_P0[i, 1] / C5_P0[i, 3]), log(C5_P0[i, 2] / C5_P0[i, 3]))
  c(ints[1], C5_SLOPE[i, 1], ints[2], C5_SLOPE[i, 2])
}, numeric(4L)))                                            # 3 x 4

# ── Shared helpers ────────────────────────────────────────────────────────────
make_corr <- function(K, rho_vals) {
  R               <- diag(K)
  R[lower.tri(R)] <- rho_vals
  R[upper.tri(R)] <- t(R)[upper.tri(R)]
  R
}

make_summary <- function(est_mat, true_vals) {
  means <- colMeans(est_mat)
  bias  <- means - true_vals
  rmse  <- sqrt(colMeans(sweep(est_mat, 2L, true_vals)^2))
  mcse  <- apply(est_mat, 2L, sd)
  data.frame(param = names(true_vals), true = unname(true_vals),
             mean_hat = means, bias = bias, rmse = rmse, mcse = mcse,
             row.names = NULL)
}

print_summary <- function(smry) { cat("\n"); print(smry, digits = 4L, row.names = FALSE); invisible(smry) }

# write a summary table to CSV and remember the path
save_csv <- function(df, fname) {
  path <- file.path(results_dir, fname)
  utils::write.csv(df, path, row.names = FALSE)
  saved_files <<- c(saved_files, path)
  invisible(path)
}

# draw on the active (RStudio) device, then save an identical PNG copy to disk
draw_and_save <- function(plot_fun, fname, width = 1400, height = 900, res = 130) {
  plot_fun()                                       # 1) display in RStudio
  path <- file.path(plots_dir, fname)
  grDevices::png(path, width = width, height = height, res = res)
  plot_fun()                                       # 2) save identical copy to disk
  grDevices::dev.off()
  saved_files <<- c(saved_files, path)
  invisible(path)
}

# AR(1) covariate (no internal seeding; caller controls the RNG via set.seed)
ar1_X <- function(T, phi) {
  sd_stat <- 1 / sqrt(1 - phi^2)
  x <- numeric(T); x[1] <- rnorm(1L, sd = sd_stat)
  for (t in 2:T) x[t] <- phi * x[t - 1L] + rnorm(1L)
  cbind(1, as.numeric(scale(x)))
}

# collect mclapply output into a data.frame, replacing any crashed worker with `fail`
collect <- function(lst, fail) {
  bad <- !vapply(lst, is.data.frame, logical(1L))
  if (any(bad)) lst[bad] <- rep(list(fail), sum(bad))
  do.call(rbind, lst)
}

# ── Case 1: const ──────────────────────────────────────────────────────────────
cat("================================================================\n")
cat("Case 1 — const, K=2 | True: rho=0.50\n")
cat("================================================================\n")

Sigma1 <- array(make_corr(C1_K, C1_RHO), dim = c(C1_K, C1_K, 1L))
mu1    <- matrix(0, nrow = 1L, ncol = C1_K)
FAIL1  <- data.frame(converged = FALSE, rho_hat = NA_real_)

run_case1 <- function(m, TT, base_seed) {
  tryCatch({
    set.seed(base_seed + m)
    sim <- rsdc_simulate(n = TT, X = matrix(1, TT, 1L), beta = matrix(nrow = 1L, ncol = 0L),
                         mu = mu1, sigma = Sigma1, N = 1L, seed = NULL)
    fit <- rsdc_estimate("const", residuals = sim$observations)
    data.frame(converged = TRUE, rho_hat = fit$correlations[1L, 1L])
  }, error = function(e) FAIL1)
}

results_case1 <- vector("list", length(T_grid)); names(results_case1) <- paste0("T", T_grid)
for (TT in T_grid) {
  cat(sprintf("\n--- T = %d  (M=%d, cores=%d) ---\n", TT, M, MC_CORES))
  t0  <- proc.time()["elapsed"]
  res <- collect(parallel::mclapply(seq_len(M), run_case1, TT = TT,
                                    base_seed = 20000L + TT * 1000L, mc.cores = MC_CORES, mc.preschedule = FALSE), FAIL1)
  elapsed <- proc.time()["elapsed"] - t0
  ok <- res$converged
  cat(sprintf("  Converged: %d / %d  (%.1f s)\n", sum(ok), M, elapsed))

  if (sum(ok) > 0L) {
    est_mat <- matrix(res$rho_hat[ok], ncol = 1L, dimnames = list(NULL, "rho"))
    smry <- make_summary(est_mat, c(rho = C1_RHO)); print_summary(smry)
    smry$T <- TT; save_csv(smry, sprintf("case1_const_T%d_summary.csv", TT))
    draw_and_save(function() {
      boxplot(est_mat[, "rho"], main = sprintf("Case 1 - const | T=%d", TT),
              col = "lightblue", ylab = "rho"); abline(h = C1_RHO, col = "red", lwd = 2)
    }, sprintf("case1_const_T%d.png", TT))
  }
  results_case1[[paste0("T", TT)]] <- list(T = TT, M = M, n_converged = sum(ok),
                                           elapsed_s = elapsed, reps = res)
}

# ── Case 2: noX ────────────────────────────────────────────────────────────────
cat("\n================================================================\n")
cat("Case 2 — noX, N=2, K=2 | True: p11=0.90, p22=0.80, rho1=0.20, rho2=0.70\n")
cat("================================================================\n")

Sigma2 <- array(0, dim = c(C2_K, C2_K, C2_N))
Sigma2[,, 1L] <- make_corr(C2_K, C2_RHO1); Sigma2[,, 2L] <- make_corr(C2_K, C2_RHO2)
mu2   <- matrix(0, nrow = C2_N, ncol = C2_K)
FAIL2 <- data.frame(converged = FALSE, p11_hat = NA_real_, p22_hat = NA_real_,
                    rho1_hat = NA_real_, rho2_hat = NA_real_)

run_case2 <- function(m, TT, base_seed) {
  tryCatch({
    set.seed(base_seed + m)
    sim <- rsdc_simulate(n = TT, X = matrix(1, TT, 1L), beta = C2_BETA,
                         mu = mu2, sigma = Sigma2, N = C2_N, seed = NULL)
    fit <- rsdc_estimate("noX", residuals = sim$observations, N = C2_N)
    data.frame(converged = TRUE,
               p11_hat = fit$transition_matrix[1L, 1L], p22_hat = fit$transition_matrix[2L, 2L],
               rho1_hat = fit$correlations[1L, 1L],     rho2_hat = fit$correlations[2L, 1L])
  }, error = function(e) FAIL2)
}

results_case2 <- vector("list", length(T_grid)); names(results_case2) <- paste0("T", T_grid)
for (TT in T_grid) {
  cat(sprintf("\n--- T = %d  (M=%d, cores=%d) ---\n", TT, M, MC_CORES))
  t0  <- proc.time()["elapsed"]
  res <- collect(parallel::mclapply(seq_len(M), run_case2, TT = TT,
                                    base_seed = 10000L + TT * 1000L, mc.cores = MC_CORES, mc.preschedule = FALSE), FAIL2)
  elapsed <- proc.time()["elapsed"] - t0
  ok <- res$converged
  cat(sprintf("  Converged: %d / %d  (%.1f s)\n", sum(ok), M, elapsed))

  if (sum(ok) > 0L) {
    est_mat <- as.matrix(res[ok, c("p11_hat", "p22_hat", "rho1_hat", "rho2_hat")])
    colnames(est_mat) <- c("p11", "p22", "rho1", "rho2")
    smry <- make_summary(est_mat, c(p11 = C2_P11, p22 = C2_P22, rho1 = C2_RHO1, rho2 = C2_RHO2))
    print_summary(smry); smry$T <- TT; save_csv(smry, sprintf("case2_noX_T%d_summary.csv", TT))
    draw_and_save(function() {
      par(mfrow = c(1, 4))
      boxplot(est_mat[, "p11"],  main = sprintf("p11 | T=%d", TT),  col = "lightgreen", ylab = "p11");  abline(h = C2_P11,  col = "red", lwd = 2)
      boxplot(est_mat[, "p22"],  main = sprintf("p22 | T=%d", TT),  col = "lightgreen", ylab = "p22");  abline(h = C2_P22,  col = "red", lwd = 2)
      boxplot(est_mat[, "rho1"], main = sprintf("rho1 | T=%d", TT), col = "lightgreen", ylab = "rho1"); abline(h = C2_RHO1, col = "red", lwd = 2)
      boxplot(est_mat[, "rho2"], main = sprintf("rho2 | T=%d", TT), col = "lightgreen", ylab = "rho2"); abline(h = C2_RHO2, col = "red", lwd = 2)
      par(mfrow = c(1, 1))
    }, sprintf("case2_noX_T%d.png", TT), width = 1800, height = 600)
  }
  results_case2[[paste0("T", TT)]] <- list(T = TT, M = M, n_converged = sum(ok),
                                           elapsed_s = elapsed, reps = res)
}

# ── Case 3: tvtp ───────────────────────────────────────────────────────────────
cat("\n================================================================\n")
cat("Case 3 — tvtp, N=2, K=2\n")
cat(sprintf("True: beta[1,]=(%.3f,-0.800), beta[2,]=(%.3f,+0.800), rho1=0.20, rho2=0.70\n",
            qlogis(0.70), qlogis(0.70)))
cat("================================================================\n")

Sigma3 <- array(0, dim = c(C3_K, C3_K, C3_N))
Sigma3[,, 1L] <- make_corr(C3_K, C3_RHO1); Sigma3[,, 2L] <- make_corr(C3_K, C3_RHO2)
mu3   <- matrix(0, nrow = C3_N, ncol = C3_K)
FAIL3 <- data.frame(converged = FALSE, label_ok = FALSE, bound_hit = NA,
                    p11_hat = NA_real_, p22_hat = NA_real_,
                    b11_hat = NA_real_, b12_hat = NA_real_, b21_hat = NA_real_,
                    b22_hat = NA_real_, rho1_hat = NA_real_, rho2_hat = NA_real_)

run_case3 <- function(m, TT, base_seed) {
  tryCatch({
    set.seed(base_seed + m)
    X   <- ar1_X(TT, C3_PHI)
    sim <- rsdc_simulate(n = TT, X = X, beta = C3_BETA, mu = mu3, sigma = Sigma3, N = C3_N, seed = NULL)
    fit <- rsdc_estimate("tvtp", residuals = sim$observations, N = C3_N, X = X)
    r1 <- fit$correlations[1L, 1L]; r2 <- fit$correlations[2L, 1L]
    data.frame(converged = TRUE, label_ok = r2 > r1, bound_hit = any(abs(fit$beta) >= 9.5),
               p11_hat = fit$transition_matrix[1L, 1L], p22_hat = fit$transition_matrix[2L, 2L],
               b11_hat = fit$beta[1L, 1L], b12_hat = fit$beta[1L, 2L],
               b21_hat = fit$beta[2L, 1L], b22_hat = fit$beta[2L, 2L],
               rho1_hat = r1, rho2_hat = r2)
  }, error = function(e) FAIL3)
}

results_case3 <- vector("list", length(T_grid)); names(results_case3) <- paste0("T", T_grid)
for (TT in T_grid) {
  cat(sprintf("\n--- T = %d  (M=%d, cores=%d) ---\n", TT, M, MC_CORES))
  t0  <- proc.time()["elapsed"]
  res <- collect(parallel::mclapply(seq_len(M), run_case3, TT = TT,
                                    base_seed = 31000L + TT * 1000L, mc.cores = MC_CORES, mc.preschedule = FALSE), FAIL3)
  elapsed <- proc.time()["elapsed"] - t0
  ok      <- res$converged
  lbl_ok  <- ok & res$label_ok
  bnd_hit <- ok & res$bound_hit            # ok=FALSE short-circuits the NA in failed rows
  usable  <- lbl_ok & !bnd_hit
  cat(sprintf("  Converged: %d / %d  |  Label OK: %d  |  Bound-hit excluded: %d  |  Usable: %d  (%.1f s)\n",
              sum(ok), M, sum(lbl_ok), sum(bnd_hit), sum(usable), elapsed))

  if (sum(usable) > 0L) {
    est_mat <- as.matrix(res[usable, c("p11_hat", "p22_hat", "b11_hat", "b12_hat", "b21_hat", "b22_hat", "rho1_hat", "rho2_hat")])
    colnames(est_mat) <- c("p11", "p22", "b11", "b12", "b21", "b22", "rho1", "rho2")
    true_vals <- c(p11 = C3_P11, p22 = C3_P22,
                   b11 = C3_BETA[1, 1], b12 = C3_BETA[1, 2], b21 = C3_BETA[2, 1],
                   b22 = C3_BETA[2, 2], rho1 = C3_RHO1, rho2 = C3_RHO2)
    smry <- make_summary(est_mat, true_vals); print_summary(smry)
    smry$T <- TT; save_csv(smry, sprintf("case3_tvtp_T%d_summary.csv", TT))
    draw_and_save(function() {
      par(mfrow = c(2, 4))
      boxplot(est_mat[, "p11"],  main = sprintf("p11 | T=%d", TT),  col = "lightyellow", ylab = "p11");  abline(h = C3_P11,        col = "red", lwd = 2)
      boxplot(est_mat[, "p22"],  main = sprintf("p22 | T=%d", TT),  col = "lightyellow", ylab = "p22");  abline(h = C3_P22,        col = "red", lwd = 2)
      boxplot(est_mat[, "b11"],  main = sprintf("b11 | T=%d", TT),  col = "lightyellow", ylab = "b11");  abline(h = C3_BETA[1, 1], col = "red", lwd = 2)
      boxplot(est_mat[, "b12"],  main = sprintf("b12 | T=%d", TT),  col = "lightyellow", ylab = "b12");  abline(h = C3_BETA[1, 2], col = "red", lwd = 2)
      boxplot(est_mat[, "b21"],  main = sprintf("b21 | T=%d", TT),  col = "lightyellow", ylab = "b21");  abline(h = C3_BETA[2, 1], col = "red", lwd = 2)
      boxplot(est_mat[, "b22"],  main = sprintf("b22 | T=%d", TT),  col = "lightyellow", ylab = "b22");  abline(h = C3_BETA[2, 2], col = "red", lwd = 2)
      boxplot(est_mat[, "rho1"], main = sprintf("rho1 | T=%d", TT), col = "lightyellow", ylab = "rho1"); abline(h = C3_RHO1,      col = "red", lwd = 2)
      boxplot(est_mat[, "rho2"], main = sprintf("rho2 | T=%d", TT), col = "lightyellow", ylab = "rho2"); abline(h = C3_RHO2,      col = "red", lwd = 2)
      par(mfrow = c(1, 1))
    }, sprintf("case3_tvtp_T%d.png", TT), width = 2000, height = 1000)
  }
  results_case3[[paste0("T", TT)]] <- list(T = TT, M = M, n_converged = sum(ok), n_label_ok = sum(lbl_ok),
                                           n_bound_hit = sum(bnd_hit), n_usable = sum(usable),
                                           elapsed_s = elapsed, reps = res)
}

# ── Case 4: noX, N=3 ─────────────────────────────────────────────────────────
cat("\n================================================================\n")
cat(sprintf("Case 4 — noX, N=3 | K in {%s}\n", paste(K_grid, collapse = ", ")))
cat(sprintf("True: diag(P)=(%.2f, %.2f, %.2f), equicorr rho=(%.2f, %.2f, %.2f)\n",
            C4_P[1, 1], C4_P[2, 2], C4_P[3, 3], C4_RHO[1], C4_RHO[2], C4_RHO[3]))
cat("================================================================\n")

FAIL4 <- data.frame(converged = FALSE,
                    p11_hat = NA_real_, p22_hat = NA_real_, p33_hat = NA_real_,
                    rho1_hat = NA_real_, rho2_hat = NA_real_, rho3_hat = NA_real_)

# KK enters the worker: regime correlations are K x K equicorrelation matrices;
# transition beta (C4_BETA) is K-invariant. Per-regime rho is summarised by the
# mean of the estimated lower-triangular correlations (all equal in truth).
run_case4 <- function(m, TT, KK, base_seed) {
  tryCatch({
    set.seed(base_seed + m)
    Sigma <- array(0, dim = c(KK, KK, C4_N))
    for (i in seq_len(C4_N)) Sigma[,, i] <- make_corr(KK, rep(C4_RHO[i], KK * (KK - 1) / 2))
    mu  <- matrix(0, nrow = C4_N, ncol = KK)
    sim <- rsdc_simulate(n = TT, X = matrix(1, TT, 1L), beta = C4_BETA,
                         mu = mu, sigma = Sigma, N = C4_N, seed = NULL)
    fit <- rsdc_estimate("noX", residuals = sim$observations, N = C4_N)
    data.frame(converged = TRUE,
               p11_hat = fit$transition_matrix[1L, 1L],
               p22_hat = fit$transition_matrix[2L, 2L],
               p33_hat = fit$transition_matrix[3L, 3L],
               rho1_hat = mean(fit$correlations[1L, ]),
               rho2_hat = mean(fit$correlations[2L, ]),
               rho3_hat = mean(fit$correlations[3L, ]))
  }, error = function(e) FAIL4)
}

results_case4 <- vector("list", length(K_grid)); names(results_case4) <- paste0("K", K_grid)
for (KK in K_grid) {
  cat(sprintf("\n········· Case 4 — K = %d ·········\n", KK))
  res_K <- vector("list", length(T_grid)); names(res_K) <- paste0("T", T_grid)
  for (TT in T_grid) {
    cat(sprintf("\n--- K = %d, T = %d  (M=%d, cores=%d) ---\n", KK, TT, M, MC_CORES))
    t0  <- proc.time()["elapsed"]
    res <- collect(parallel::mclapply(seq_len(M), run_case4, TT = TT, KK = KK,
                                      base_seed = 40000000L + KK * 3000000L + TT * 1000L, mc.cores = MC_CORES, mc.preschedule = FALSE), FAIL4)
    elapsed <- proc.time()["elapsed"] - t0
    ok <- res$converged
    cat(sprintf("  Converged: %d / %d  (%.1f s)\n", sum(ok), M, elapsed))

    if (sum(ok) > 0L) {
      est_mat <- as.matrix(res[ok, c("p11_hat", "p22_hat", "p33_hat", "rho1_hat", "rho2_hat", "rho3_hat")])
      colnames(est_mat) <- c("p11", "p22", "p33", "rho1", "rho2", "rho3")
      smry <- make_summary(est_mat, c(p11 = C4_P[1, 1], p22 = C4_P[2, 2], p33 = C4_P[3, 3],
                                      rho1 = C4_RHO[1], rho2 = C4_RHO[2], rho3 = C4_RHO[3]))
      print_summary(smry); smry$K <- KK; smry$T <- TT
      save_csv(smry, sprintf("case4_noX_N3_K%d_T%d_summary.csv", KK, TT))
      draw_and_save(function() {
        par(mfrow = c(2, 3))
        boxplot(est_mat[, "p11"],  main = sprintf("p11 | K=%d,T=%d", KK, TT),  col = "lightgreen", ylab = "p11");  abline(h = C4_P[1, 1], col = "red", lwd = 2)
        boxplot(est_mat[, "p22"],  main = sprintf("p22 | K=%d,T=%d", KK, TT),  col = "lightgreen", ylab = "p22");  abline(h = C4_P[2, 2], col = "red", lwd = 2)
        boxplot(est_mat[, "p33"],  main = sprintf("p33 | K=%d,T=%d", KK, TT),  col = "lightgreen", ylab = "p33");  abline(h = C4_P[3, 3], col = "red", lwd = 2)
        boxplot(est_mat[, "rho1"], main = sprintf("rho1 | K=%d,T=%d", KK, TT), col = "lightgreen", ylab = "rho1"); abline(h = C4_RHO[1], col = "red", lwd = 2)
        boxplot(est_mat[, "rho2"], main = sprintf("rho2 | K=%d,T=%d", KK, TT), col = "lightgreen", ylab = "rho2"); abline(h = C4_RHO[2], col = "red", lwd = 2)
        boxplot(est_mat[, "rho3"], main = sprintf("rho3 | K=%d,T=%d", KK, TT), col = "lightgreen", ylab = "rho3"); abline(h = C4_RHO[3], col = "red", lwd = 2)
        par(mfrow = c(1, 1))
      }, sprintf("case4_noX_N3_K%d_T%d.png", KK, TT), width = 1800, height = 1000)
    }
    res_K[[paste0("T", TT)]] <- list(K = KK, T = TT, M = M, n_converged = sum(ok),
                                     elapsed_s = elapsed, reps = res)
  }
  results_case4[[paste0("K", KK)]] <- res_K
}


# ── Case 5 parameters redefined for K=2 — WELL-IDENTIFIED design ───────────────
C5_N   <- 3L
C5_PHI <- 0.7
M5     <- 100                        # TEST (set to 100L for the full study)
KK     <- 2L
T_grid <- c(1000L, 2000L)                  # large T => N=3 recovered cleanly

# WELL-SEPARATED correlations (key to convergence) + PERSISTENT regimes
C5_RHO   <- c(-0.6, 0.1, 0.8)
C5_P0    <- matrix(c(0.90, 0.07, 0.03,
                     0.05, 0.90, 0.05,
                     0.03, 0.07, 0.90), nrow = 3L, byrow = TRUE)
C5_SLOPE <- matrix(c( 0.5, 0.0,
                      0.0, 0.5,
                      -0.5, 0.0), nrow = 3L, byrow = TRUE)
# beta = N x (N-1)*p = 3 x 4 (p = ncol(X) = 2), packed as (int1, slope1, int2, slope2);
# derived from C5_P0 (reference = regime 3) => consistent with the DGP.
C5_BETA <- t(vapply(1:3, function(i) {
  ints <- c(log(C5_P0[i, 1] / C5_P0[i, 3]), log(C5_P0[i, 2] / C5_P0[i, 3]))
  c(ints[1], C5_SLOPE[i, 1], ints[2], C5_SLOPE[i, 2])
}, numeric(4L)))

# ── Case 5: tvtp, N=3 ─────────────────────────────────────────────────────────
cat("\n================================================================\n")
cat(sprintf("Case 5 — tvtp, N=3 | K = %d | M=%d (well-identified design)\n", KK, M5))
cat(sprintf("True (at x=0): diag(P0)=(%.2f, %.2f, %.2f), rho=(%.2f, %.2f, %.2f)\n",
            C5_P0[1, 1], C5_P0[2, 2], C5_P0[3, 3], C5_RHO[1], C5_RHO[2], C5_RHO[3]))
cat("================================================================\n")

FAIL5 <- data.frame(converged = FALSE, label_ok = FALSE, bound_hit = NA,
                    p11_hat = NA_real_, p22_hat = NA_real_, p33_hat = NA_real_,
                    b11_hat = NA_real_, b12_hat = NA_real_, b13_hat = NA_real_, b14_hat = NA_real_,
                    b21_hat = NA_real_, b22_hat = NA_real_, b23_hat = NA_real_, b24_hat = NA_real_,
                    b31_hat = NA_real_, b32_hat = NA_real_, b33_hat = NA_real_, b34_hat = NA_real_,
                    rho1_hat = NA_real_, rho2_hat = NA_real_, rho3_hat = NA_real_)

run_case5 <- function(m, TT, KK, base_seed) {
  tryCatch({
    set.seed(base_seed + m)
    X     <- ar1_X(TT, C5_PHI)
    Sigma <- array(0, dim = c(KK, KK, C5_N))
    for (i in seq_len(C5_N)) Sigma[,, i] <- make_corr(KK, rep(C5_RHO[i], KK * (KK - 1) / 2))
    mu    <- matrix(0, nrow = C5_N, ncol = KK)
    sim <- rsdc_simulate(n = TT, X = X, beta = C5_BETA, mu = mu, sigma = Sigma, N = C5_N, seed = NULL)
    fit <- rsdc_estimate("tvtp", residuals = sim$observations, N = C5_N, X = X)
    r1 <- mean(fit$correlations[1L, ]); r2 <- mean(fit$correlations[2L, ]); r3 <-
      mean(fit$correlations[3L, ])
    data.frame(converged = TRUE, label_ok = (r1 < r2) & (r2 < r3),
               bound_hit = any(abs(fit$beta) >= 9.5),
               p11_hat = fit$transition_matrix[1L, 1L],
               p22_hat = fit$transition_matrix[2L, 2L],
               p33_hat = fit$transition_matrix[3L, 3L],
               b11_hat = fit$beta[1L, 1L], b12_hat = fit$beta[1L, 2L], b13_hat = fit$beta[1L, 3L],
               b14_hat = fit$beta[1L, 4L],
               b21_hat = fit$beta[2L, 1L], b22_hat = fit$beta[2L, 2L], b23_hat = fit$beta[2L, 3L],
               b24_hat = fit$beta[2L, 4L],
               b31_hat = fit$beta[3L, 1L], b32_hat = fit$beta[3L, 2L], b33_hat = fit$beta[3L, 3L],
               b34_hat = fit$beta[3L, 4L],
               rho1_hat = r1, rho2_hat = r2, rho3_hat = r3)
  }, error = function(e) FAIL5)
}

results_case5 <- vector("list", 1L); names(results_case5) <- paste0("K", KK)
res_K <- vector("list", length(T_grid)); names(res_K) <- paste0("T", T_grid)

for (TT in T_grid) {
  cat(sprintf("\n--- K = %d, T = %d  (M=%d, cores=%d) ---\n", KK, TT, M5, MC_CORES))
  t0  <- proc.time()["elapsed"]
  res <- collect(parallel::mclapply(seq_len(M5), run_case5, TT = TT, KK = KK,
                                    base_seed = 70000000L + KK * 3000000L + TT * 1000L,
                                    mc.cores = MC_CORES, mc.preschedule = FALSE), FAIL5)
  elapsed <- proc.time()["elapsed"] - t0
  ok      <- res$converged
  lbl_ok  <- ok & res$label_ok
  bnd_hit <- ok & res$bound_hit
  usable  <- lbl_ok & !bnd_hit
  cat(sprintf("  Converged: %d / %d  |  Label OK: %d  |  Bound-hit excluded: %d  |  Usable: %d  (%.1f
  s)\n",
              sum(ok), M5, sum(lbl_ok), sum(bnd_hit), sum(usable), elapsed))
  
  if (sum(usable) > 0L) {
    cols <- c("p11_hat", "p22_hat", "p33_hat",
              "b11_hat", "b12_hat", "b13_hat", "b14_hat",
              "b21_hat", "b22_hat", "b23_hat", "b24_hat",
              "b31_hat", "b32_hat", "b33_hat", "b34_hat",
              "rho1_hat", "rho2_hat", "rho3_hat")
    est_mat <- as.matrix(res[usable, cols]); colnames(est_mat) <- sub("_hat$", "", cols)
    true_vals <- c(p11 = C5_P0[1, 1], p22 = C5_P0[2, 2], p33 = C5_P0[3, 3],
                   b11 = C5_BETA[1, 1], b12 = C5_BETA[1, 2], b13 = C5_BETA[1, 3], b14 = C5_BETA[1, 4],
                   b21 = C5_BETA[2, 1], b22 = C5_BETA[2, 2], b23 = C5_BETA[2, 3], b24 = C5_BETA[2, 4],
                   b31 = C5_BETA[3, 1], b32 = C5_BETA[3, 2], b33 = C5_BETA[3, 3], b34 = C5_BETA[3, 4],
                   rho1 = C5_RHO[1], rho2 = C5_RHO[2], rho3 = C5_RHO[3])
    smry <- make_summary(est_mat, true_vals); print_summary(smry)
    smry$K <- KK; smry$T <- TT; save_csv(smry, sprintf("case5_tvtp_N3_K%d_T%d_summary.csv", KK, TT))
    draw_and_save(function() {
      par(mfrow = c(3, 6))
      for (nm in colnames(est_mat)) {
        boxplot(est_mat[, nm], main = sprintf("%s | K=%d,T=%d", nm, KK, TT), col = "lightyellow", ylab
                = nm)
        abline(h = true_vals[[nm]], col = "red", lwd = 2)
      }
      par(mfrow = c(1, 1))
    }, sprintf("case5_tvtp_N3_K%d_T%d.png", KK, TT), width = 2400, height = 1200)
  }
  res_K[[paste0("T", TT)]] <- list(K = KK, T = TT, M = M5, n_converged = sum(ok), n_label_ok =
                                     sum(lbl_ok),
                                   n_bound_hit = sum(bnd_hit), n_usable = sum(usable),
                                   elapsed_s = elapsed, reps = res)
}

results_case5[[paste0("K", KK)]] <- res_K

# ── Save RDS bundle + final confirmation ───────────────────────────────────────
rds_path <- file.path(results_dir, "monte_carlo_results.rds")
saveRDS(list(created = Sys.time(), M = M, T_grid = T_grid, mc_cores = MC_CORES,
             case1 = results_case1, case2 = results_case2, case3 = results_case3,
             case4 = results_case4, case5 = results_case5), rds_path)
saved_files <- c(saved_files, rds_path)

cat("\n================================================================\n")
cat(sprintf("DONE. %d files written:\n", length(saved_files)))
for (f in saved_files) cat("  ", normalizePath(f, mustWork = FALSE), "\n", sep = "")
cat("================================================================\n")
