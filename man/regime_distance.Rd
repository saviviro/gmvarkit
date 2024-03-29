% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parameterReforms.R
\name{regime_distance}
\alias{regime_distance}
\title{Calculate "distance" between two (scaled) regimes
 \strong{\eqn{\upsilon_{m}}}\eqn{ = (\phi_{m,0},}\strong{\eqn{\phi_{m}}}\eqn{,\sigma_{m})}}
\usage{
regime_distance(regime_pars1, regime_pars2)
}
\arguments{
\item{regime_pars1}{a length \eqn{pd^2+d+d(d+1)/2} vector
\strong{\eqn{\upsilon_{m}}}\eqn{ = (\phi_{m,0},}\strong{\eqn{\phi_{m}}}\eqn{,\sigma_{m})}.}

\item{regime_pars2}{a length \eqn{pd^2+d+d(d+1)/2} vector
\strong{\eqn{\upsilon_{m}}}\eqn{ = (\phi_{m,0},}\strong{\eqn{\phi_{m}}}\eqn{,\sigma_{m})}.}
}
\value{
Returns "distance" between \code{regime_pars1} and \code{regime_pars2}. Values are scaled
  before calculating the "distance". Read the source code for more details.
}
\description{
\code{regime_distance} calculates "distance" between two scaled regimes. This is used in
  the genetic algorithm.
}
\section{Warning}{

 No argument checks!
}

\references{
\itemize{
   \item Kalliovirta L., Meitz M. and Saikkonen P. 2016. Gaussian mixture vector autoregression.
         \emph{Journal of Econometrics}, \strong{192}, 485-498.
   \item Virolainen S. (forthcoming). A statistically identified structural vector autoregression with endogenously
          switching volatility regime. \emph{Journal of Business & Economic Statistics}.
   \item Virolainen S. 2022. Gaussian and Student's t mixture vector autoregressive model with application to the
     asymmetric effects of monetary policy shocks in the Euro area. Unpublished working
     paper, available as arXiv:2109.13648.
 }
 @keywords internal
}
\keyword{internal}
