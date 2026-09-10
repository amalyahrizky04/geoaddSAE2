#' Simulated area-level small area estimation data
#'
#' A synthetic area-level dataset used in the package examples and tests. It
#' mimics a typical Geoadditive SAE demonstration setup: a linear covariate,
#' a nonlinear covariate, spatial coordinates, and a direct estimator with
#' known sampling variance -- generated from a true model with a linear
#' effect, a nonlinear covariate effect, and a smooth spatial effect.
#'
#' @format A data frame with 100 rows and 9 columns:
#' \describe{
#'   \item{area}{Area identifier.}
#'   \item{y}{Direct estimator (response).}
#'   \item{x1}{Covariate with a linear effect on \code{y}.}
#'   \item{x2}{Another covariate with a linear effect on \code{y}.}
#'   \item{x3}{Covariate with a nonlinear effect on \code{y}.}
#'   \item{lat}{Latitude-like spatial coordinate.}
#'   \item{lon}{Longitude-like spatial coordinate.}
#'   \item{vardir}{Known sampling variance of the direct estimator.}
#'   \item{theta}{True small area parameter used to generate the simulated data.}
#' }
#' @examples
#' data(simulated_sae)
#' head(simulated_sae)
"simulated_sae"
