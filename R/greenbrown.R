#' Green vs Brown portfolio dataset
#'
#' Daily returns for a green and a brown portfolios constructed following
#' the equal-weighted 10-90 percentile approach.
#'
#' @format A data frame with 2266 rows and three columns:
#' \describe{
#'   \item{DATE}{Dates ranging from 2014-01-02 to 2024-12-30.}
#'   \item{return_green}{Numeric returns for the green portfolio.}
#'   \item{return_brown}{Numeric returns for the brown portfolio.}
#' }
#'
#' @source Originally data in \code{inst/extdata/green-brown-ptf.xlsx}.
#'
#' @examples
#' data("greenbrown")
#' str(greenbrown)
#' head(greenbrown)
#'
#' @keywords datasets
#' @docType data
#' @usage data(greenbrown)
"greenbrown"
