% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GIRFandGFEVD.R, R/plotMethods.R,
%   R/printMethods.R
\name{GIRF}
\alias{GIRF}
\alias{plot.girf}
\alias{print.girf}
\title{Estimate generalized impulse response function for
 structural (and reduced form) GMVAR, StMVAR, and G-StMVAR models.}
\usage{
GIRF(
  gsmvar,
  which_shocks,
  shock_size = 1,
  N = 30,
  R1 = 250,
  R2 = 250,
  init_regimes = 1:sum(gsmvar$model$M),
  init_values = NULL,
  which_cumulative = numeric(0),
  scale = NULL,
  scale_type = c("instant", "peak"),
  scale_horizon = N,
  ci = c(0.95, 0.8),
  include_mixweights = TRUE,
  ncores = 2,
  plot_res = TRUE,
  seeds = NULL,
  ...
)

\method{plot}{girf}(x, add_grid = FALSE, margs, ...)

\method{print}{girf}(x, ..., digits = 2, N_to_print)
}
\arguments{
\item{gsmvar}{an object of class \code{'gsmvar'}, typically created with \code{fitGSMVAR} or \code{GSMVAR}.}

\item{which_shocks}{a numeric vector of length at most \eqn{d}
(\code{=ncol(data)}) and elements in \eqn{1,...,d} specifying the
structural shocks for which the GIRF should be estimated.}

\item{shock_size}{a non-zero scalar value specifying the common size for all scalar
components of the structural shock. Note that the conditional covariance
matrix of the structural shock is an identity matrix and that the
(generalized) impulse responses may not be symmetric to the sign and size
of the shock.}

\item{N}{a positive integer specifying the horizon how far ahead should the
generalized impulse responses be calculated.}

\item{R1}{the number of repetitions used to estimate GIRF for each initial
value.}

\item{R2}{the number of initial values to be drawn from a stationary
distribution of the process or of a specific regime? The confidence bounds
will be sample quantiles of the GIRFs based on different initial values.
Ignored if the argument \code{init_value} is specified.}

\item{init_regimes}{a numeric vector of length at most \eqn{M} and elements
in \eqn{1,...,M} specifying the regimes from which the initial values
should be generated from. The initial values will be generated from a
mixture distribution with the mixture components being the stationary
distributions of the specific regimes and the (proportional) mixing weights
given by the mixing weight parameters of those regimes. Note that if
\code{init_regimes=1:M}, the initial values are generated from the
stationary distribution of the process and if \code{init_regimes=m}, the
initial value are generated from the stationary distribution of the
\eqn{m}th regime. Ignored if the argument \code{init_values} is specified.}

\item{init_values}{a size \eqn{(p\times d)} matrix specifying the initial values, where d is the number
of time series in the system. The \strong{last} row will be used as initial values for the first lag,
the second last row for second lag etc. If not specified, initial values will be drawn according to
mixture distribution specifed by the argument \code{init_regimes}.}

\item{which_cumulative}{a numeric vector with values in \eqn{1,...,d}
(\code{d=ncol(data)}) specifying which the variables for which the impulse
responses should be cumulative. Default is none.}

\item{scale}{should the GIRFs to some of the shocks be scaled so that they
correspond to a specific magnitude of instantaneous or peak response
of some specific variable (see the argument \code{scale_type})?
Provide a length three vector where the shock of interest
is given in the first element (an integer in \eqn{1,...,d}), the variable of
interest is given in the second element (an integer in \eqn{1,...,d}), and
the magnitude of its instantaneous or peak response in the third element
(a non-zero real number). If the GIRFs of multiple shocks should be scaled, provide
a matrix which has one column for each of the shocks with the columns being
the length three vectors described above.}

\item{scale_type}{If argument \code{scale} is specified, should the GIRFs be
scaled to match an instantaneous response (\code{"instant"}) or peak response
(\code{"peak"}). If \code{"peak"}, the scale is based on the largest magnitude
of peak response in absolute value. Ignored if \code{scale} is not specified.}

\item{scale_horizon}{If \code{scale_type == "peak"} what the maximum horizon up
to which peak response is expected? Scaling won't based on values after this.}

\item{ci}{a numeric vector with elements in \eqn{(0, 1)} specifying the
confidence levels of the confidence intervals.}

\item{include_mixweights}{should the generalized impulse response be
calculated for the mixing weights as well? \code{TRUE} or \code{FALSE}.}

\item{ncores}{the number CPU cores to be used in parallel computing. Only
single core computing is supported if an initial value is specified (and
the GIRF won't thus be estimated multiple times).}

\item{plot_res}{\code{TRUE} if the results should be plotted, \code{FALSE} if
not.}

\item{seeds}{a length \code{R2} vector containing the random number generator
seed for estimation of each GIRF. A single number of an initial value is
specified. or \code{NULL} for not initializing the seed. Exists for
creating reproducible results.}

\item{...}{arguments passed to \code{grid} which plots grid to the figure.}

\item{x}{object of class \code{'girf'} generated by the function \code{GIRF}.}

\item{add_grid}{should grid be added to the plots?}

\item{margs}{numeric vector of length four that adjusts the
\code{[bottom_marginal, left_marginal, top_marginal, right_marginal]}
as the relative sizes of the marginals to the figures of the responses.}

\item{digits}{the number of decimals to print}

\item{N_to_print}{an integer specifying the horizon how far to print the estimates and
confidence intervals. The default is that all the values are printed.}
}
\value{
Returns a class \code{'girf'} list with the GIRFs in the first
  element (\code{$girf_res}) and the used arguments the rest. The first
  element containing the GIRFs is a list with the \eqn{m}th element
  containing the point estimates for the GIRF in \code{$point_est} (the first
  element) and confidence intervals in \code{$conf_ints} (the second
  element). The first row is for the GIRF at impact \eqn{(n=0)}, the second
  for \eqn{n=1}, the third for \eqn{n=2}, and so on.

  The element \code{$all_girfs} is a list containing results from all the individual GIRFs
  obtained from the MC repetitions. Each element is for one shock and results are in
  array of the form \code{[horizon, variables, MC-repetitions]}.
}
\description{
\code{GIRF} estimates generalized impulse response function for
 structural (and reduced form) GMVAR, StMVAR, and G-StMVAR models.
}
\details{
The model DOES NOT need to be structural in order for this function to be
  applicable. When an identified structural GMVAR, StMVAR, or G-StMVAR model is
  provided in the argument \code{gsmvar}, the identification imposed by the model
  is used. When a reduced form model is provided in the argument \code{gsmvar},
  lower triangular Cholesky identification is used to identify the shocks.

  The confidence bounds reflect uncertainty about the initial state (but
  currently not about the parameter estimates) if initial values are not
  specified. If initial values are specified, there won't currently be
  confidence intervals. See the cited paper by Virolainen (2022) for details
  about the algorithm.

  Note that if the argument \code{scale} is used, the scaled responses of
  the mixing weights might be more than one in absolute value.
}
\section{Functions}{
\itemize{
\item \code{plot(girf)}: plot method

\item \code{print(girf)}: print method

}}
\examples{
 \donttest{
 # These are long-running examples that use parallel computing.
 # It takes approximately 30 seconds to run all the below examples.

 ## StMVAR(1, 2), d=2 model identified recursively by lower-triangular
 ## Cholesky decomposition (i.e., reduced form model is specified):
 params12t <- c(0.55, 0.11, 0.34, 0.05, -0.01, 0.72, 0.58, 0.01, 0.06, 0.17,
   0.25, 0.34, 0.05, -0.01, 0.72, 0.50, -0.01, 0.20, 0.60, 3.00, 12.00)
 mod12t <- GSMVAR(gdpdef, p=1, M=2, params=params12t, model="StMVAR")

 # Estimating the GIRFs of both structural shocks with initial values
 # drawn from the stationary distribution of the process,
 # 12 periods ahead, confidence levels 0.95 and 0.8:
 girf0 <- GIRF(mod12t, N=12, R1=100, R2=100)
 girf0
 plot(girf0)
 ## NOTE: Small R1 and R2 is used here to shorten the estimation time.
 ## Larger R1 and R2 should be considered in empirical applications!

 ## Structural GMVAR(2, 2), d=2 model identified with sign-constraints:
 params22s <- c(0.36, 0.121, 0.484, 0.072, 0.223, 0.059, -0.151, 0.395,
  0.406, -0.005, 0.083, 0.299, 0.218, 0.02, -0.119, 0.722, 0.093, 0.032,
   0.044, 0.191, 0.057, 0.172, -0.46, 0.016, 3.518, 5.154, 0.58)
 W_22 <- matrix(c(1, 1, -1, 1), nrow=2, byrow=FALSE)
 mod22s <- GSMVAR(gdpdef, p=2, M=2, params=params22s,
  structural_pars=list(W=W_22))
 mod22s
 # Alternatively, use:
 #fit22s <- fitGSMVAR(gdpdef, p=2, M=2, structural_pars=list(W=W_22),
 #                   ncalls=20, seeds=1:20)
 # To obtain an estimated version of the same model.

 # Estimating the GIRFs of both structural shocks with initial values
 # drawn from the stationary distribution of the process,
 # 12 periods ahead, confidence levels 0.95 and 0.8:
 girf1 <- GIRF(mod22s, N=12, R1=100, R2=100)
 girf1
 plot(girf1)

 # Estimating the GIRF of the second shock only, 12 periods ahead
 # and shock size 1, initial values drawn from the stationary distribution
 # of the first regime, confidence level 0.9:
 girf2 <- GIRF(mod22s, which_shocks=2, shock_size=1, N=12, init_regimes=1,
               ci=0.9, R1=100, R2=100)

 # Estimating the GIRFs of both structural shocks, negative one standard
 # error shock, N=20 periods ahead, estimation based on 200 Monte Carlo
 # simulations, and fixed initial values given by the last p observations
 # of the data:
 girf3 <- GIRF(mod22s, shock_size=-1, N=20, R1=200,
               init_values=mod22s$data)
 }
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
\seealso{
\code{\link{GFEVD}}, \code{\link{linear_IRF}}, \code{\link{fitGSMVAR}}, \code{\link{GSMVAR}},
  \code{\link{gsmvar_to_sgsmvar}}, \code{\link{reorder_W_columns}},
  \code{\link{swap_W_signs}}, \code{\link{simulate.gsmvar}},
  \code{\link{predict.gsmvar}}, \code{\link{profile_logliks}},
  \code{\link{quantile_residual_tests}}, \code{\link{LR_test}},
  \code{\link{Wald_test}}
}
