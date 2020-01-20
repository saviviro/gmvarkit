#' Euro area and U.S. long-term government bond yields and Euro-U.S. dollar exchange rate.
#'
#' A dataset containing time series of the difference between the monthly Euro area and U.S.
#' long-term government bond yields and monthly average Euro - U.S. dollar exchange rate. The data
#' covers the time period January 1989 - December 2009 with monthly frequency. This is the same data
#' (in non-scaled form) that is used by Kalliovirta et. al. (2016).
#'
#' @format A numeric matrix of class \code{'ts'} with 252 rows and 2 columns with one time series in each column:
#' \describe{
#'   \item{First column:}{The difference between the monthly Euro area and U.S. long-term government bond yields
#'   (10 year maturity, i_euro - i_us), from January 1989 to December 2009. calculated by the ECB and the
#'   Federal Reserve Board; prior to 2001, the Euro area data refer to the "EU11" countries, and afterwards
#'   with changing composition eventually to the "EU17" by the end of the data period.}
#'   \item{Second column:}{Monthly average Euro - U.S. dollar exchange rate, from January 1989 to December 2009.
#'   Based on the ECU - USD exchange rate prior to 1999.}
#' }
#'
#' @inherit GMVAR references
#' @source OECD Statistics
"eurusd"
