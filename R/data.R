#' Euro area and U.S. long-term government bond yields and Euro-U.S. dollar exchange rate
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


#' U.S. data containing log-differences of industrial production index, consumer price index, and M1, and an interest rate
#' variable
#'
#' A dataset containing the monthly U.S. data covering the period from February 1959 to December 2019
#' (731 observations) and consisting of four variables: the log-difference of industrial production index (IP),
#' the log-difference of consumer price index (CPI), the log-difference of M1 monetary aggregate (M1), and an interest rate
#' variable (RATE). The log-differences are multiplied by hundred. The interest rate variable is the effective
#' federal funds rate from February 1959 to August 2008 after which we replaced it with the Wu and Xia
#' (2016) shadow rate, which is not constrained by the zero lower bound and also
#' quantifies unconventional monetary policy measures.
#'
#'  The Wu and Xia (2016) shadow rate data was retrieved from the Federal Reserve Bank of Atlanta's website
#'  (\url{https://www.atlantafed.org/cqer/research/wu-xia-shadow-federal-funds-rate}) and the rest of the data was
#'  retrieved from the Federal Reserve Bank of St. Louis database.
#'
#' @format A numeric matrix of class \code{'ts'} with 731 rows and 4 columns with one time series in each column:
#' \describe{
#'   \item{IP:}{The log-difference multiplied by hundred, \url{https://fred.stlouisfed.org/series/INDPRO}}
#'   \item{CPI:}{The log-difference multiplied by hundred, \url{https://fred.stlouisfed.org/series/CPIAUCSL}}
#'   \item{M1:}{The log-difference multiplied by hundred, \url{https://fred.stlouisfed.org/series/M1SL}}
#'   \item{RATE:}{From 1959 February to 2008 August \url{https://fred.stlouisfed.org/series/FEDFUNDS} and
#'                from 2008 September onwards \url{https://www.atlantafed.org/cqer/research/wu-xia-shadow-federal-funds-rate}}
#' }
#'
#' @references
#'  \itemize{
#'    \item J.C. Wu and F.D. Xia. 2016. Measuring the Macroeconomic Impact of Monetary Policy at the Zero Lower Bound.
#'          \emph{Journal of Money, Credit and Banking}, 48(2-3): 253-291.
#'  }
#' @source The Federal Reserve Bank of St. Louis database and the Federal Reserve Bank of Atlanta's website
"usamone_prec"


#' U.S. data containing log-differences of industrial production index, consumer price index, and M1, and an interest rate
#' variable
#'
#' A dataset containing the monthly U.S. data covering the period from February 1959 to December 2020
#' (743 observations) and consisting of four variables: the log-difference of industrial production index (IP),
#' the log-difference of consumer price index (CPI), the log-difference of M1 monetary aggregate (M1), and an interest rate
#' variable (RATE). The log-differences are multiplied by hundred. The interest rate variable is the effective
#' federal funds rate from February 1959 to August 2008 after which we replaced it with the Wu and Xia
#' (2016) shadow rate, which is not constrained by the zero lower bound and also
#' quantifies unconventional monetary policy measures.
#'
#'  The Wu and Xia (2016) shadow rate data was retrieved from the Federal Reserve Bank of Atlanta's website
#'  (\url{https://www.atlantafed.org/cqer/research/wu-xia-shadow-federal-funds-rate}) and the rest of the data was
#'  retrieved from the Federal Reserve Bank of St. Louis database.
#'
#' @format A numeric matrix of class \code{'ts'} with 743 rows and 4 columns with one time series in each column:
#' \describe{
#'   \item{IP:}{The log-difference multiplied by hundred, \url{https://fred.stlouisfed.org/series/INDPRO}}
#'   \item{CPI:}{The log-difference multiplied by hundred, \url{https://fred.stlouisfed.org/series/CPIAUCSL}}
#'   \item{M1:}{The log-difference multiplied by hundred, \url{https://fred.stlouisfed.org/series/M1SL}}
#'   \item{RATE:}{From 1959 February to 2008 August \url{https://fred.stlouisfed.org/series/FEDFUNDS} and
#'                from 2008 September onwards \url{https://www.atlantafed.org/cqer/research/wu-xia-shadow-federal-funds-rate}}
#' }
#'
#' @inherit usamone_prec references source
"usamone"
