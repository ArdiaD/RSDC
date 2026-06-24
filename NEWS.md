# Changes in Version 1.4-0 (DA,BS,RN)
- Performance: the Hamilton-filter log-likelihood used during estimation is now
  evaluated in C++ (Rcpp/RcppArmadillo), matching the R reference to ~1e-8 and
  cutting estimation time by roughly an order of magnitude. The pure-R
  `rsdc_hamilton()` is retained (and used as the equivalence reference).
- Optimiser control: `control` now also accepts `cores` (parallel `DEoptim` via
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

# Changes in Version 1.3-0 (DA,BS,RN)
- `rsdc_estimate()` now returns an object of class `"rsdc_fit"` with standard S3
  methods: `print`, `summary`, `coef`, `logLik`, `nobs`, `vcov`, `confint`,
  `predict`, and `simulate`. `AIC()`/`BIC()` work out of the box via `logLik`.
- Standard errors: the estimator now returns the observed-information
  variance-covariance (`vcov`) and standard errors (`summary`/`confint`), computed
  from the numerical Hessian of the negative log-likelihood at the MLE.
- Arbitrary number of regimes: `N >= 4` is now supported for `"noX"` and `"tvtp"`
  (previously capped at `N = 3`).
- `control` now forwards optimiser settings (`itermax`, `NP`, `parallelType`,
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
- DESCRIPTION improved following CRAN's guidelines
- Seed is now a control parameter
- Removed examples for unexported functions
- Fixed donot run examples
- Added green-brown-ptf data in extdata
- Data added under the right .rda format
- Data Documentation added

# Changes in Version 1.1-1 (DA)
- URls fixed

# Changes in Version 1.1-0 (DA)
- First release public version
