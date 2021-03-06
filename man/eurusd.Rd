% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{eurusd}
\alias{eurusd}
\title{Euro area and U.S. long-term government bond yields and Euro-U.S. dollar exchange rate}
\format{
A numeric matrix of class \code{'ts'} with 252 rows and 2 columns with one time series in each column:
\describe{
  \item{First column:}{The difference between the monthly Euro area and U.S. long-term government bond yields
  (10 year maturity, i_euro - i_us), from January 1989 to December 2009. calculated by the ECB and the
  Federal Reserve Board; prior to 2001, the Euro area data refer to the "EU11" countries, and afterwards
  with changing composition eventually to the "EU17" by the end of the data period.}
  \item{Second column:}{Monthly average Euro - U.S. dollar exchange rate, from January 1989 to December 2009.
  Based on the ECU - USD exchange rate prior to 1999.}
}
}
\source{
OECD Statistics
}
\usage{
eurusd
}
\description{
A dataset containing time series of the difference between the monthly Euro area and U.S.
long-term government bond yields and monthly average Euro - U.S. dollar exchange rate. The data
covers the time period January 1989 - December 2009 with monthly frequency. This is the same data
(in non-scaled form) that is used by Kalliovirta et. al. (2016).
}
\references{
\itemize{
   \item Kalliovirta L., Meitz M. and Saikkonen P. 2016. Gaussian mixture vector autoregression.
         \emph{Journal of Econometrics}, \strong{192}, 485-498.
   \item Kalliovirta L. and Saikkonen P. 2010. Reliable Residuals for Multivariate Nonlinear
         Time Series Models. \emph{Unpublished Revision of HECER Discussion Paper No. 247}.
   \item Virolainen S. 2020. Structural Gaussian mixture vector autoregressive model. Unpublished working
     paper, available as arXiv:2007.04713.
 }
}
\keyword{datasets}
