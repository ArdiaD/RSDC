#' @name RSDC
#' @aliases RSDC-package
#' @title RSDC: Regime-Switching Correlation Models for Portfolio Analysis
#'
#' @description
#' The \pkg{RSDC} package provides a comprehensive framework for modeling,
#' estimating, and forecasting correlation structures in multivariate time
#' series under regime-switching dynamics. It supports both fixed transition
#' probabilities and \emph{time-varying transition probabilities} (TVTP)
#' driven by exogenous variables.
#'
#' The methodology is particularly suited to empirical asset pricing and
#' portfolio management applications, enabling users to incorporate macroeconomic,
#' financial, or climate-related predictors into the regime dynamics. The package
#' integrates the full workflow — from model estimation to covariance matrix
#' reconstruction and portfolio optimization — in a single, reproducible pipeline.
#'
#' @section Main Features:
#' \itemize{
#'   \item \strong{Model estimation and filtering:}
#'     \code{\link{rsdc_hamilton}} (Hamilton filter),
#'     \code{\link{rsdc_likelihood}} (likelihood computation),
#'     \code{\link{rsdc_estimate}} (parameter estimation).
#'   \item \strong{Correlation and covariance forecasting:}
#'     \code{\link{rsdc_forecast}}.
#'   \item \strong{Portfolio construction:}
#'     \code{\link{rsdc_minvar}} (minimum-variance portfolios),
#'     \code{\link{rsdc_maxdiv}} (maximum-diversification portfolios).
#'   \item \strong{Simulation:}
#'     \code{\link{rsdc_simulate}} (simulate TVTP regime-switching series).
#' }
#'
#' @section Authors:
#' David Ardia and Benjamin Seguin
#'
#' @references
#' \insertRef{engle2002dynamic}{RSDC} \cr
#'
#' \insertRef{hamilton1989}{RSDC} \cr
#'
#' \insertRef{pelletier2006regime}{RSDC}
#'
#' @importFrom DEoptim DEoptim
#' @importFrom stats optim plogis cor sd
#' @importFrom utils combn
NULL

