#' Example portfolio data (green vs brown)
#'
#' This is an Excel file shipped with the package in \code{inst/extdata/}.
#' It illustrates how to read external data (e.g., green and brown portfolio returns)
#' for use in empirical illustrations of the package functions.
#'
#' @format An Excel file with daily returns of green and brown portfolio
#' Each file corresponds to green and brown portfolio constructed in a different way.
#'
#' @details
#' To access the file path:
#' \preformatted{
#' file <- system.file("extdata", "green-brown-ptf.xlsx", package = "RSDC")
#' }
#' Then read it with your preferred reader, e.g.:
#' \preformatted{
#' data <- readxl::read_excel(file)
#' }
#'
#' @seealso \link[readxl:read_excel]{readxl::read_excel}
#'
#' @keywords datasets
#' @name green_brown_ptf
#' @docType data
NULL
