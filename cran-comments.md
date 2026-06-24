## Submission

This is the first CRAN submission of RSDC (regime-switching dynamic correlation
models). The version number is 1.5-0 because the package was developed and
versioned iteratively before this initial public release (see NEWS.md); the
functionality is stable and the API is settled.

## Test environments
* Local: macOS Sequoia 15.6, R 4.5.1
* win-builder: R release and R devel
* R-hub: Windows, Ubuntu (R devel)

## R CMD check results
0 errors | 0 warnings | 1 note

## Notes
* New submission.
* "Namespace in Imports field not imported from: 'Rdpack'". Rdpack is imported
  only to provide the Rd macros (e.g. \insertRef) used for references in the
  documentation; it is declared under RdMacros in DESCRIPTION, so the import is
  required even though no R object is imported from it.

## Notes for the maintainer
* The package contains compiled code (Rcpp/RcppArmadillo). A local clang build
  may emit a `-Wfixed-enum-extension` warning that originates in R's own headers
  (Boolean.h), not in package code; it does not appear on CRAN's build machines.
* Examples that estimate a model are wrapped in \donttest{} to keep check times
  short; the full suite (including those paths) runs locally and on the CI
  environments above.
