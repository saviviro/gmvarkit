#' U.S. real GDP percent change and GDP implicit price deflator percent change.
#'
#' A dataset containing a quarterly U.S. time series with two components:
#' the percentage change of real GDP and the percentage change of GDP implicit price deflator,
#' covering the period from 1959Q1 - 2019Q4.
#'
#' @format A numeric matrix of class \code{'ts'} with 244 rows and 2 columns with one time series in each column:
#' \describe{
#'   \item{First column (GDP):}{The quarterly percent change of real U.S. GDP, from 1959Q1 to 2019Q4, \url{https://fred.stlouisfed.org/series/GDPC1}.}
#'   \item{Second column (GDPDEF):}{The quarterly percent change of U.S. GDP implicit price deflator, from 1959Q1 to 2019Q4, \url{https://fred.stlouisfed.org/series/GDPDEF}.}
#' }
#'
#' @source The Federal Reserve Bank of St. Louis database
"gdpdef"


#' A quarterly U.S. data covering the period from 1954Q3 to 2021Q4 (270 observations) and consisting four variables:
#' cyclical component of the log of real GDP, the log-difference of GDP implicit price deflator, the log-diffence of producer
#' price index (all commodities), and an interest rate variable. The interest rate variable is the effective federal funds
#' rate from 1954Q3 to 2008Q2 and after that the Wu and Xia (2016) shadow rate, which is not constrained by the zero lower
#' bound and also quantifies unconventional monetary policy measures. The log-differences of the GDP deflator and producer price
#' index are multiplied by hundred. This data is the one that was used in Virolainen (2022).
#'
#' The cyclical component of the log of real GDP was obtained by applying a one-sided Hodrick-Prescott (HP) filter with the
#' standard smoothing parameter lambda=1600. The one-sided filter was obtained from the two-sided HP filter by applying the
#' filter up to horizon t, taking the last observation, and repeating this procedure for the full sample t=1,...,T.
#' In order to allow the series to start from any phase of the cycle, we applied the one-sided filter to the full available
#' sample from 1947Q1 to 2021Q1 before extracting our sample period from it. We computed the two-sided HP filters with the R
#' package lpirfs (Adämmer, 2021)
#'
#' @format A numeric matrix of class \code{'ts'} with 270 rows and 4 columns with one time series in each column:
#' \describe{
#'   \item{First column (GDP):}{The cyclical component of the log of real GDP, \url{https://fred.stlouisfed.org/series/GDPC1}.}
#'   \item{Second column (GDPDEF):}{The log-difference of GDP implicit price deflator, \url{https://fred.stlouisfed.org/series/GDPDEF}.}
#'   \item{Third column (PPI):}{The log-difference of producer price index (all commodities), \url{https://fred.stlouisfed.org/series/PPIACO}.}
#'   \item{Third column (RATE):}{The Federal funds rate from 1954Q3 to 2008Q2 and after that the Wu and Xia (2016) shadow rate,
#'    \url{https://fred.stlouisfed.org/series/FEDFUNDS}, \url{https://www.atlantafed.org/cqer/research/wu-xia-shadow-federal-funds-rate}.}
#' }
#'
#' @source The Federal Reserve Bank of St. Louis database and the Federal Reserve Bank of Atlanta's website
#' @references
#'  \itemize{
#'    \item Adämmer P. 2021. lprfs: Local Projections Impulse Response Functions. R package version: 0.2.0,
#'      \url{https://CRAN.R-project.org/package=lpirfs}.
#'    \item Virolainen S. 2022. Structural Gaussian mixture vector autoregressive model with application to the asymmetric
#'      effects of monetary policy shocks. Unpublished working paper, available as arXiv:2007.04713.
#'    \item Wu J. and Xia F. 2016. Measuring the macroeconomic impact of monetary policy at the zero lower bound.
#'      \emph{Journal of Money, Credit and Banking}, 48(2-3): 253-291.
#'  }
"usamone"


#' A monthly Euro area data covering the period from January 1999 to December 2021 (276 observations) and consisting four variables:
#' cyclical component of log industrial production index, the log-difference of harmonized consumer price index, the log-difference
#' of Brent crude oil prices (Europe), and an interest rate variable. The interest rate variable is the Euro overnight index average
#' rate (EONIA) from January 1999 to October 2008, and after that the Wu and Xia (2016) shadow rate, which is not constrained by the zero lower
#' bound and also quantifies unconventional monetary policy measures. The log-difference of the harmonized consumer price index is
#' multiplied by hundred and the log-difference of oil price by ten. This data is the one that was used in Virolainen (2022).
#'
#' The cyclical component of the log of industrial production index was obtained by applying the linear projection filter proposed
#' by Hamilton (2018) using the parameter values h=24 and p=12. In order to obtain as accurate estimates as possible, we applied the
#' filter to the full available sample from January 1991 to December 2021 before extracting our sample period from it.
#' package lpirfs (Adämmer, 2021).
#'
#' @format A numeric matrix of class \code{'ts'} with 276 rows and 4 columns with one time series in each column:
#' \describe{
#'   \item{First column (IPI):}{The cyclical component of the log of industrial production index,
#'     url is https://sdw.ecb.europa.eu/quickview.do?SERIES_KEY=132.STS.M.I8.Y.PROD.NS0010.4.000.}
#'   \item{Second column (HCPI):}{The log-difference of harmonized consumer price index,
#'     url is https://sdw.ecb.europa.eu/quickview.do?SERIES_KEY=122.ICP.M.U2.Y.000000.3.INX.}
#'   \item{Third column (OIL):}{The log-difference of Brent crude oil price (Europe),
#'     \url{https://fred.stlouisfed.org/series/MCOILBRENTEU}.}
#'   \item{Third column (RATE):}{The EONIA from January 1999 to October 2008 and after that the Wu and Xia (2016) shadow rate,
#'     urls are https://sdw.ecb.europa.eu/quickview.do?SERIES_KEY=143.FM.M.U2.EUR.4F.MM.EONIA.HSTA and
#'     \url{https://sites.google.com/view/jingcynthiawu/shadow-rates}.}
#' }
#'
#' @source The Federal Reserve Bank of St. Louis database and the Federal Reserve Bank of Atlanta's website
#' @references
#'  \itemize{
#'    \item Virolainen S. 2022. Gaussian and Student's t mixture vector autoregressive model with application to the
#'      asymmetric effects of monetary policy shocks in the Euro area. Unpublished working
#'      paper, available as arXiv:2109.13648.
#'    \item Wu J. and Xia F. 2016. Measuring the macroeconomic impact of monetary policy at the zero lower bound.
#'      \emph{Journal of Money, Credit and Banking}, 48(2-3): 253-291.
#'  }
"euromone"
