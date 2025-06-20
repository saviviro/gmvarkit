% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loglikelihood.R
\name{loglikelihood_int}
\alias{loglikelihood_int}
\title{Compute log-likelihood of a GMVAR, StMVAR, and G-StMVAR models}
\usage{
loglikelihood_int(
  data,
  p,
  M,
  params,
  model = c("GMVAR", "StMVAR", "G-StMVAR"),
  conditional = TRUE,
  parametrization = c("intercept", "mean"),
  constraints = NULL,
  same_means = NULL,
  weight_constraints = NULL,
  structural_pars = NULL,
  to_return = c("loglik", "mw", "mw_tplus1", "loglik_and_mw", "terms", "regime_cmeans",
    "regime_ccovs", "total_cmeans", "total_ccovs", "arch_scalars",
    "loglik_mw_archscalars"),
  check_params = TRUE,
  minval = NULL,
  stat_tol = 0.001,
  posdef_tol = 1e-08,
  df_tol = 1e-08
)
}
\arguments{
\item{data}{a matrix or class \code{'ts'} object with \code{d>1} columns. Each column is taken to represent
a univariate time series. \code{NA} values are not supported.}

\item{p}{a positive integer specifying the autoregressive order of the model.}

\item{M}{\describe{
  \item{For \strong{GMVAR} and \strong{StMVAR} models:}{a positive integer specifying the number of mixture components.}
  \item{For \strong{G-StMVAR} models:}{a size \eqn{(2\times 1)} integer vector specifying the number of \emph{GMVAR type} components \code{M1}
   in the first element and \emph{StMVAR type} components \code{M2} in the second element. The total number of mixture components
   is \code{M=M1+M2}.}
}}

\item{params}{a real valued vector specifying the parameter values.
  \describe{
    \item{\strong{For unconstrained models:}}{
      Should be size \eqn{((M(pd^2+d+d(d+1)/2+2)-M1-1)x1)} and have the form
      \strong{\eqn{\theta}}\eqn{ = }(\strong{\eqn{\upsilon}}\eqn{_{1}},
      ...,\strong{\eqn{\upsilon}}\eqn{_{M}}, \eqn{\alpha_{1},...,\alpha_{M-1},}\strong{\eqn{\nu}}\eqn{)}, where
      \itemize{
        \item \strong{\eqn{\upsilon}}\eqn{_{m}} \eqn{ = (\phi_{m,0},}\strong{\eqn{\phi}}\eqn{_{m}}\eqn{,\sigma_{m})}
        \item \strong{\eqn{\phi}}\eqn{_{m}}\eqn{ = (vec(A_{m,1}),...,vec(A_{m,p})}
        \item and \eqn{\sigma_{m} = vech(\Omega_{m})}, m=1,...,M,
        \item \strong{\eqn{\nu}}\eqn{=(\nu_{M1+1},...,\nu_{M})}
        \item \eqn{M1} is the number of GMVAR type regimes.
      }
    }
    \item{\strong{For constrained models:}}{
      Should be size \eqn{((M(d+d(d+1)/2+2)+q-M1-1)x1)} and have the form
      \strong{\eqn{\theta}}\eqn{ = (\phi_{1,0},...,\phi_{M,0},}\strong{\eqn{\psi},}
      \eqn{\sigma_{1},...,\sigma_{M},\alpha_{1},...,\alpha_{M-1},}\strong{\eqn{\nu}}), where
      \itemize{
        \item \strong{\eqn{\psi}} \eqn{(qx1)} satisfies (\strong{\eqn{\phi}}\eqn{_{1}}\eqn{,...,}
        \strong{\eqn{\phi}}\eqn{_{M}) =} \strong{\eqn{C \psi}} where \strong{\eqn{C}} is a \eqn{(Mpd^2xq)}
        constraint matrix.
      }
    }
    \item{\strong{For same_means models:}}{
      Should have the form
      \strong{\eqn{\theta}}\eqn{ = (}\strong{\eqn{\mu},}\strong{\eqn{\psi},}
      \eqn{\sigma_{1},...,\sigma_{M},\alpha_{1},...,\alpha_{M-1},}\strong{\eqn{\nu}}\eqn{)}, where
      \itemize{
        \item \strong{\eqn{\mu}}\eqn{= (\mu_{1},...,\mu_{g})} where
          \eqn{\mu_{i}} is the mean parameter for group \eqn{i} and
          \eqn{g} is the number of groups.
        \item If AR constraints are employed, \strong{\eqn{\psi}} is as for constrained
          models, and if AR constraints are not employed, \strong{\eqn{\psi}}\eqn{ = }
          (\strong{\eqn{\phi}}\eqn{_{1}}\eqn{,...,}\strong{\eqn{\phi}}\eqn{_{M})}.
      }
    }
     \item{\strong{For models with weight_constraints:}}{Drop \eqn{\alpha_1,...,\alpha_{M-1}} from
      the parameter vector.}
    \item{\strong{For structural models:}}{
      Reduced form models can be directly used as recursively identified structural models. If the structural model is
      identified by conditional heteroskedasticity, the parameter vector should have the form
      \strong{\eqn{\theta}}\eqn{ = (\phi_{1,0},...,\phi_{M,0},}\strong{\eqn{\phi}}\eqn{_{1},...,}\strong{\eqn{\phi}}\eqn{_{M},
      vec(W),}\strong{\eqn{\lambda}}\eqn{_{2},...,}\strong{\eqn{\lambda}}\eqn{_{M},\alpha_{1},...,\alpha_{M-1},}\strong{\eqn{\nu}}\eqn{)},
       where
      \itemize{
        \item\strong{\eqn{\lambda}}\eqn{_{m}=(\lambda_{m1},...,\lambda_{md})} contains the eigenvalues of the \eqn{m}th mixture component.
      }
      \describe{
        \item{\strong{If AR parameters are constrained: }}{Replace \strong{\eqn{\phi}}\eqn{_{1}}\eqn{,...,}
        \strong{\eqn{\phi}}\eqn{_{M}} with \strong{\eqn{\psi}} \eqn{(qx1)} that satisfies (\strong{\eqn{\phi}}\eqn{_{1}}\eqn{,...,}
        \strong{\eqn{\phi}}\eqn{_{M}) =} \strong{\eqn{C \psi}}, as above.}
        \item{\strong{If same_means: }}{Replace \eqn{(\phi_{1,0},...,\phi_{M,0})} with \eqn{(\mu_{1},...,\mu_{g})},
          as above.}
        \item{\strong{If \eqn{W} is constrained:}}{Remove the zeros from \eqn{vec(W)} and make sure the other entries satisfy
         the sign constraints.}
        \item{\strong{If \eqn{\lambda_{mi}} are constrained via \code{C_lambda}:}}{Replace \strong{\eqn{\lambda}}\eqn{_{2},...,}
        \strong{\eqn{\lambda}}\eqn{_{M}} with \strong{\eqn{\gamma}} \eqn{(rx1)} that satisfies (\strong{\eqn{\lambda}}\eqn{_{2}}
        \eqn{,...,} \strong{\eqn{\lambda}}\eqn{_{M}) =} \strong{\eqn{C_{\lambda} \gamma}} where \eqn{C_{\lambda}} is
         a \eqn{(d(M-1) x r)} constraint matrix.}
        \item{\strong{If \eqn{\lambda_{mi}} are constrained via \code{fixed_lambdas}:}}{Drop \strong{\eqn{\lambda}}\eqn{_{2},...,}
        \strong{\eqn{\lambda}}\eqn{_{M}} from the parameter vector.}
      }
    }
  }
  Above, \eqn{\phi_{m,0}} is the intercept parameter, \eqn{A_{m,i}} denotes the \eqn{i}th coefficient matrix of the \eqn{m}th
  mixture component, \eqn{\Omega_{m}} denotes the error term covariance matrix of the \eqn{m}:th mixture component, and
  \eqn{\alpha_{m}} is the mixing weight parameter. The \eqn{W} and \eqn{\lambda_{mi}} are structural parameters replacing the
  error term covariance matrices (see Virolainen, 2022). If \eqn{M=1}, \eqn{\alpha_{m}} and \eqn{\lambda_{mi}} are dropped.
  If \code{parametrization=="mean"}, just replace each \eqn{\phi_{m,0}} with regimewise mean \eqn{\mu_{m}}.
  \eqn{vec()} is vectorization operator that stacks columns of a given matrix into a vector. \eqn{vech()} stacks columns
  of a given matrix from the principal diagonal downwards (including elements on the diagonal) into a vector.

  In the \strong{GMVAR model}, \eqn{M1=M} and \strong{\eqn{\nu}} is dropped from the parameter vector. In the \strong{StMVAR} model,
  \eqn{M1=0}. In the \strong{G-StMVAR} model, the first \code{M1} regimes are \emph{GMVAR type} and the rest \code{M2} regimes are
  \emph{StMVAR type}. In \strong{StMVAR} and \strong{G-StMVAR} models, the degrees of freedom parameters in \strong{\eqn{\nu}} should
  be strictly larger than two.

  The notation is similar to the cited literature.}

\item{model}{is "GMVAR", "StMVAR", or "G-StMVAR" model considered? In the G-StMVAR model, the first \code{M1} components
are GMVAR type and the rest \code{M2} components are StMVAR type.}

\item{conditional}{a logical argument specifying whether the conditional or exact log-likelihood function
should be used.}

\item{parametrization}{\code{"intercept"} or \code{"mean"} determining whether the model is parametrized with intercept
parameters \eqn{\phi_{m,0}} or regime means \eqn{\mu_{m}}, m=1,...,M.}

\item{constraints}{a size \eqn{(Mpd^2 \times  q)} constraint matrix \strong{\eqn{C}} specifying general linear constraints
to the autoregressive parameters. We consider constraints of form
(\strong{\eqn{\phi}}\eqn{_{1}}\eqn{,...,}\strong{\eqn{\phi}}\eqn{_{M}) = }\strong{\eqn{C \psi}},
where \strong{\eqn{\phi}}\eqn{_{m}}\eqn{ = (vec(A_{m,1}),...,vec(A_{m,p}) (pd^2 x 1), m=1,...,M},
contains the coefficient matrices and \strong{\eqn{\psi}} \eqn{(q x 1)} contains the related parameters.
For example, to restrict the AR-parameters to be the same for all regimes, set \strong{\eqn{C}}=
[\code{I:...:I}]\strong{'} \eqn{(Mpd^2 \times  pd^2)} where \code{I = diag(p*d^2)}.
Ignore (or set to \code{NULL}) if linear constraints should \strong{not} be employed.}

\item{same_means}{Restrict the mean parameters of some regimes to be the same? Provide a list of numeric vectors
such that each numeric vector contains the regimes that should share the common mean parameters. For instance, if
\code{M=3}, the argument \code{list(1, 2:3)} restricts the mean parameters of the second and third regime to be
the same but the first regime has freely estimated (unconditional) mean. Ignore or set to \code{NULL} if mean parameters
should not be restricted to be the same among any regimes. \strong{This constraint is available only for mean parametrized models;
that is, when \code{parametrization="mean"}.}}

\item{weight_constraints}{a numeric vector of length \eqn{M-1} specifying fixed parameter values for the mixing weight parameters
\eqn{\alpha_m, \ m=1,...,M-1}. Each element should be strictly between zero and one, and the sum of all the elements should
be strictly less than one.}

\item{structural_pars}{If \code{NULL} a reduced form model is considered. Reduced models can be used directly as recursively
identified structural models. For a structural model identified by conditional heteroskedasticity, should be a list containing
at least the first one of the following elements:
\itemize{
  \item \code{W} - a \eqn{(dxd)} matrix with its entries imposing constraints on \eqn{W}: \code{NA} indicating that the element is
    unconstrained, a positive value indicating strict positive sign constraint, a negative value indicating strict
    negative sign constraint, and zero indicating that the element is constrained to zero.
  \item \code{C_lambda} - a \eqn{(d(M-1) \times  r)} constraint matrix that satisfies (\strong{\eqn{\lambda}}\eqn{_{2}}\eqn{,...,}
    \strong{\eqn{\lambda}}\eqn{_{M}) =} \strong{\eqn{C_{\lambda} \gamma}} where \strong{\eqn{\gamma}} is the new \eqn{(r \times 1)}
    parameter subject to which the model is estimated (similarly to AR parameter constraints). The entries of \code{C_lambda}
    must be either \strong{positive} or \strong{zero}. Ignore (or set to \code{NULL}) if the eigenvalues \eqn{\lambda_{mi}}
    should not be constrained.
  \item \code{fixed_lambdas} - a length \eqn{d(M-1)} numeric vector (\strong{\eqn{\lambda}}\eqn{_{2}}\eqn{,...,}
    \strong{\eqn{\lambda}}\eqn{_{M})} with elements strictly larger than zero specifying the fixed parameter values for the
    parameters \eqn{\lambda_{mi}} should be constrained to. This constraint is alternative \code{C_lambda}.
    Ignore (or set to \code{NULL}) if the eigenvalues \eqn{\lambda_{mi}} should not be constrained.
}
See Virolainen (forthcoming) for the conditions required to identify the shocks and for the B-matrix as well (it is \eqn{W} times
a time-varying diagonal matrix with positive diagonal entries).}

\item{to_return}{should the returned object be the log-likelihood value, which is default, or something else?
See the section "Return" for all the options.}

\item{check_params}{should it be checked that the parameter vector satisfies the model assumptions? Can be skipped to save
computation time if it does for sure.}

\item{minval}{the value that will be returned if the parameter vector does not lie in the parameter space
(excluding the identification condition).}

\item{stat_tol}{numerical tolerance for stationarity of the AR parameters: if the "bold A" matrix of any regime
has eigenvalues larger that \code{1 - stat_tol} the model is classified as non-stationary. Note that if the
tolerance is too small, numerical evaluation of the log-likelihood might fail and cause error.}

\item{posdef_tol}{numerical tolerance for positive definiteness of the error term covariance matrices: if
the error term covariance matrix of any regime has eigenvalues smaller than this, the model is classified
as not satisfying positive definiteness assumption. Note that if the tolerance is too small, numerical
evaluation of the log-likelihood might fail and cause error.}

\item{df_tol}{the parameter vector is considered to be outside the parameter space if all degrees of
freedom parameters are not larger than \code{2 + df_tol}.}
}
\value{
\describe{
  \item{By default:}{log-likelihood value of the specified GMVAR, StMVAR, or G-StMVAR model,}
  \item{If \code{to_return=="mw"}:}{a size \eqn{((n_obs-p)\times  M) }matrix containing the mixing weights: for m:th component in m:th column.}
  \item{If \code{to_return=="mw_tplus1"}:}{a size \eqn{((n_obs-p+1) \times  M)} matrix containing the mixing weights: for m:th component in m:th column.
    The last row is for \eqn{\alpha_{m,T+1}}.}
  \item{If \code{to_return=="terms"}:}{a size \eqn{((n_obs-p) \times 1)} numeric vector containing the terms \eqn{l_{t}}.}
  \item{if \code{to_return=="loglik_and_mw"}:}{a list of two elements. The first element contains the log-likelihood value and the
    second element contains the mixing weights.}
  \item{If \code{to_return=="regime_cmeans"}:}{an \code{[T-p, d, M]} array containing the regimewise conditional means
   (the first p values are used as the initial values).}
   \item{If \code{to_return=="regime_ccovs"}:}{an \code{[d, d, T-p, M]} array containing the regimewise conditional
   covariance matrices (the first p values are used as the initial values). The index \code{[ , , t, m]} gives the time
   \code{t} conditional covariance matrix for the regime \code{m}.}
  \item{If \code{to_return=="total_cmeans"}:}{a \code{[T-p, d]} matrix containing the conditional means of the process
   (the first p values are used as the initial values).}
  \item{If \code{to_return=="total_ccov"}:}{an \code{[d, d, T-p]} array containing the conditional covariance matrices of the process
   (the first p values are used as the initial values).}
  \item{If \code{to_return=="arch_scalars"}:}{a \code{[T-p, M]} matrix containing the regimewise arch scalars
   multiplying error term covariance matrix in the conditional covariance matrix of the regime. For GMVAR type regimes, these
   are all ones (the first p values are used as the initial values).}
  \item{if \code{to_return=="loglik_mw_archscalars"}:}{a list of three elements. The first element contains the log-likelihood value, the
    second element contains the mixing weights, the third element contains the arch scalars
    (this is used in \code{quantile_residuals_int}).}
 }
}
\description{
\code{loglikelihood_int} computes log-likelihoodof a GMVAR, StMVAR, and G-StMVAR models.
}
\details{
\code{loglikelihood_int} takes use of the function \code{dmvn} from the package \code{mvnfast}.
}
\references{
\itemize{
   \item Kalliovirta L., Meitz M. and Saikkonen P. 2016. Gaussian mixture vector autoregression.
         \emph{Journal of Econometrics}, \strong{192}, 485-498.
   \item Lütkepohl H. 2005. New Introduction to Multiple Time Series Analysis,
         \emph{Springer}.
   \item McElroy T. 2017. Computation of vector ARMA autocovariances.
         \emph{Statistics and Probability Letters}, \strong{124}, 92-96.
   \item Virolainen S. 2025. A statistically identified structural vector autoregression with endogenously
         switching volatility regime. \emph{Journal of Business & Economic Statistics}, \strong{43}, 1, 44-54.
   \item Virolainen S. 2022. Gaussian and Student's t mixture vector autoregressive model with application to the
     asymmetric effects of monetary policy shocks in the Euro area. Unpublished working
     paper, available as arXiv:2109.13648.
 }
}
\keyword{internal}
