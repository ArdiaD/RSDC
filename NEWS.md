# Changes in Version 1.6-0 (DA,BS,RN)
- Bug fix: the estimation backends no longer call `set.seed()` when a warm
  start skips the global search. Previously every warm-started refit reset the
  global RNG stream, so `rsdc_bootstrap()` — which interleaves simulation and
  warm-started re-estimation — produced nearly identical replicates from the
  second one on, drastically understating bootstrap standard errors and
  intervals (also via `vcov(type = "bootstrap")` and
  `confint(type = "bootstrap")`).

# Changes in Version 1.5-0 (DA,BS,RN)
- Bug fix (audit): `rsdc_forecast()` now returns a scalar `NA` BIC (instead of a
  silent zero-length `numeric(0)`) when `final_params$log_likelihood` is missing.
- Bug fix (audit): a warm start (`control$start`) lying outside the default
  L-BFGS-B search box — as produced by the identifiability relabeling for
  `N >= 3` (re-referenced softmax betas, former row-complement transition
  entries) — is no longer silently projected onto the box; the box is widened
  element-wise to keep the supplied start feasible (affects multi-start refits
  and `rsdc_bootstrap()` warm starts).
- CRAN hygiene: stray local session files under `inst/simulation/results/` are
  now build-ignored, and the declared `Rdpack` import is used
  (`importFrom(Rdpack, reprompt)`), clearing both `R CMD check` NOTEs.
- `rsdc_minvar()`/`rsdc_maxdiv()`: when `value_cols` is a character vector and
  `y` carries column names, misaligned columns are now detected — reordered to
  `value_cols` (with a warning) when the names match as a set, otherwise a
  warning states that positional order is assumed. Previously the names were
  silently overwritten.
- `rsdc_forecast_ahead()` now validates that `X_future` has the same number of
  columns as the estimation covariate matrix `X`.
- `plot.rsdc_fit()` error message corrected: stored regime probabilities do not
  depend on `compute_se`.
- Documentation: `rsdc_corr_bands()` documents the truncation of invalid
  Gaussian draws near parameter bounds.
- Removed dead `requireNamespace("mvtnorm")` guards (`mvtnorm` is a hard
  dependency in `Imports`).
- Multi-start estimation: `control$n_starts` repeats the global+local search from
  several seeds and keeps the highest-likelihood fit, storing `start_logliks`.
  Since `DEoptim` is a stochastic global search, this is primarily a stability
  diagnostic across seeds (is the optimum reproducible?), not a replacement for it.
- `rsdc_forecast_ahead()`: genuine multi-step-ahead forecasts of the regime
  distribution and implied correlations, propagating the terminal filtered state
  through the Markov chain (`X_future` for the time-varying case).
- `rsdc_corr_bands()`: pointwise uncertainty bands for the predicted correlation
  path, by drawing parameters from the asymptotic sampling distribution and
  re-running the filter.
- `rsdc_viterbi()`: most likely (MAP) regime path via the Viterbi algorithm.
- broom tidiers `tidy()`, `glance()`, `augment()` and a `ggplot2::autoplot()`
  method for `rsdc_fit` (ggplot2 is a soft dependency in `Suggests`).
- Documentation: the "Getting started" vignette gains a section demonstrating the
  new features above (multi-start, bands, Viterbi, multi-step forecasts, tidy
  output), and the broom/autoplot methods gain runnable examples.

# Changes in Version 1.4-0 (DA,BS,RN)
- Performance: the Hamilton-filter log-likelihood used during estimation is now
  evaluated in C++ (Rcpp/RcppArmadillo), matching the R reference to ~1e-8 and
  cutting estimation time by roughly an order of magnitude. The pure-R
  `rsdc_hamilton()` is retained (and used as the equivalence reference).
- Optimizer control: `control` now also accepts `cores` (parallel `DEoptim` via
  `parallelType = 1`) and `start` (a warm-start parameter vector that skips the
  global search and goes straight to local refinement).
- Parametric bootstrap: `rsdc_bootstrap()` simulates from the fitted model,
  re-estimates each replicate warm-started at the MLE, and returns bootstrap
  standard errors / percentile intervals that respect the parameter bounds.
- Robust standard errors: `vcov()`/`confint()` gain a `type` argument —
  `"hessian"` (default), `"opg"` (outer product of gradients), or `"sandwich"`
  (QML). Per-observation scores are stored on the fitted object.
- `summary()` now reports delta-method standard errors on the natural scale of
  the transition probabilities, plus regime diagnostics (stay probability,
  expected duration, ergodic distribution).
- Added a `plot()` method for fitted models (smoothed/filtered regime
  probabilities).
- Fitted objects now store the estimation residuals/covariates and the
  filtered/smoothed regime probabilities.
- Input-validation hardening (audit): `rsdc_estimate()`/`rsdc_hamilton()` now
  reject single-series input (`K < 2`) and a covariate matrix `X` whose row
  count does not match the data; `rsdc_hamilton()` validates that a user-supplied
  fixed transition matrix `P` is row-stochastic (finite, non-negative, rows
  summing to 1).
- `confint()` now uses the same safe square-root as `summary()`, so non-positive
  variance estimates map to `NA` instead of producing inconsistent `NaN`
  intervals.
- For the fixed-transition `noX` model with `N >= 3`, a transition row that the
  optimizer leaves infeasible (free probabilities summing above 1) is still
  projected onto the simplex, but the fit is now flagged as non-converged and its
  standard errors are suppressed (the projected point is not a valid basis for
  inference).

# Changes in Version 1.3-0 (DA,BS,RN)
- `rsdc_estimate()` now returns an object of class `"rsdc_fit"` with standard S3
  methods: `print`, `summary`, `coef`, `logLik`, `nobs`, `vcov`, `confint`,
  `predict`, and `simulate`. `AIC()`/`BIC()` work out of the box via `logLik`.
- Standard errors: the estimator now returns the observed-information
  variance-covariance (`vcov`) and standard errors (`summary`/`confint`), computed
  from the numerical Hessian of the negative log-likelihood at the MLE.
- Arbitrary number of regimes: `N >= 4` is now supported for `"noX"` and `"tvtp"`
  (previously capped at `N = 3`).
- `control` now forwards optimizer settings (`itermax`, `NP`, `parallelType`,
  `steptol`, `maxit`, and `compute_se`) to `DEoptim`/`optim`, enabling faster runs.
- Robustness: warns when a transition probability is pinned to its bound and when
  a `"tvtp"` covariate matrix `X` has no intercept/constant column.
- Added a "Getting started" vignette.
- Moved the 596 KB source workbook out of `inst/extdata` into build-ignored
  `data-raw/` to shrink the package tarball.

# Changes in Version 1.2-0 (BS,DA,RN)
- Added `mccc` dataset: daily Media Climate Change Concerns (Aggregate) index,
  forward-filled from monthly and aligned row-for-row to `greenbrown`. Supplies
  the exogenous TVTP covariate used in the vignette, so the empirical workflow
  is fully reproducible from packaged data alone (no external Excel files).
- Added lightweight raw source `inst/extdata/mccc-monthly.csv` and builder
  `data-raw/mccc.R`.
- Fixed `greenbrown` documentation: date range is 2014-01-02 to 2022-12-30.

# Changes in Version 1.1-2 (BS,DA)
- DESCRIPTION improved following the CRAN guidelines
- Seed is now a control parameter
- Removed examples for unexported functions
- Fixed do-not-run examples
- Added green-brown-ptf data in extdata
- Data added under the right .rda format
- Data Documentation added

# Changes in Version 1.1-1 (DA)
- URLs fixed

# Changes in Version 1.1-0 (DA)
- First release public version
