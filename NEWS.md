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
