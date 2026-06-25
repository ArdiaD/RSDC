## Submission

This is an update of the package already on CRAN (current version 1.1-2,
published 2025-09-03) to version 1.5-0. The update adds the C++ Hamilton filter,
the `mccc`/`greenbrown` datasets, the `"rsdc_fit"` S3 interface (print/summary/
coef/logLik/nobs/vcov/confint/predict/simulate/plot), numerical and bootstrap
standard errors, multi-step forecasts, parameter-uncertainty bands, Viterbi
decoding, and broom/ggplot2 methods. See NEWS.md for the full history.

## Test environments
* Local: macOS Sequoia 15.6, R 4.5.1
* win-builder: R release and R devel
* R-hub: Windows, Ubuntu (R devel)

## R CMD check results
0 errors | 0 warnings | 1 note

## Notes
* "Namespace in Imports field not imported from: 'Rdpack'". Rdpack is imported
  only to provide the Rd macros (e.g. \insertRef) used for references in the
  documentation; it is declared under RdMacros in DESCRIPTION, so the import is
  required even though no R object is imported from it.

## Reverse dependencies
There are no reverse dependencies on CRAN.

## Notes for the maintainer
* The package contains compiled code (Rcpp/RcppArmadillo). A local clang build
  may emit a `-Wfixed-enum-extension` warning that originates in R's own headers
  (Boolean.h), not in package code; it does not appear on CRAN's build machines.
* Examples that estimate a model are wrapped in \donttest{} to keep check times
  short; the full suite (including those paths) runs locally and on the CI
  environments above.
