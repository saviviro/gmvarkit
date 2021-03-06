% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pickParams.R
\name{pick_Am}
\alias{pick_Am}
\title{Pick coefficient matrices}
\usage{
pick_Am(p, M, d, params, m, structural_pars = NULL)
}
\arguments{
\item{p}{a positive integer specifying the autoregressive order of the model.}

\item{M}{a positive integer specifying the number of mixture components.}

\item{d}{number of time series in the system, i.e. the dimension.}

\item{params}{a real valued vector specifying the parameter values.
\describe{
  \item{\strong{For unconstrained models:}}{
    Should be size \eqn{((M(pd^2+d+d(d+1)/2+1)-1)x1)} and have the form
    \strong{\eqn{\theta}}\eqn{ = }(\strong{\eqn{\upsilon}}\eqn{_{1}},
    ...,\strong{\eqn{\upsilon}}\eqn{_{M}}, \eqn{\alpha_{1},...,\alpha_{M-1}}), where
    \itemize{
      \item \strong{\eqn{\upsilon}}\eqn{_{m}} \eqn{ = (\phi_{m,0},}\strong{\eqn{\phi}}\eqn{_{m}}\eqn{,\sigma_{m})}
      \item \strong{\eqn{\phi}}\eqn{_{m}}\eqn{ = (vec(A_{m,1}),...,vec(A_{m,p})}
      \item and \eqn{\sigma_{m} = vech(\Omega_{m})}, m=1,...,M.
    }
  }
  \item{\strong{For constrained models:}}{
    Should be size \eqn{((M(d+d(d+1)/2+1)+q-1)x1)} and have the form
    \strong{\eqn{\theta}}\eqn{ = (\phi_{1,0},...,\phi_{M,0},}\strong{\eqn{\psi},}
    \eqn{\sigma_{1},...,\sigma_{M},\alpha_{1},...,\alpha_{M-1})}, where
    \itemize{
      \item \strong{\eqn{\psi}} \eqn{(qx1)} satisfies (\strong{\eqn{\phi}}\eqn{_{1}}\eqn{,...,}
      \strong{\eqn{\phi}}\eqn{_{M}) =} \strong{\eqn{C \psi}} where \strong{\eqn{C}} is \eqn{(Mpd^2xq)}
      constraint matrix.
    }
  }
  \item{\strong{For same_means models:}}{
    Should have the form
    \strong{\eqn{\theta}}\eqn{ = (}\strong{\eqn{\mu},}\strong{\eqn{\psi},}
    \eqn{\sigma_{1},...,\sigma_{M},\alpha_{1},...,\alpha_{M-1})}, where
    \itemize{
      \item \strong{\eqn{\mu}}\eqn{= (\mu_{1},...,\mu_{g})} where
        \eqn{\mu_{i}} is the mean parameter for group \eqn{i} and
        \eqn{g} is the number of groups.
      \item If AR constraints are employed, \strong{\eqn{\psi}} is as for constrained
        models, and if AR constraints are not employed, \strong{\eqn{\psi}}\eqn{ = }
        (\strong{\eqn{\phi}}\eqn{_{1}}\eqn{,...,}\strong{\eqn{\phi}}\eqn{_{M})}.
    }
  }
  \item{\strong{For structural GMVAR model:}}{
    Should have the form
    \strong{\eqn{\theta}}\eqn{ = (\phi_{1,0},...,\phi_{M,0},}\strong{\eqn{\phi}}\eqn{_{1},...,}\strong{\eqn{\phi}}\eqn{_{M},
    vec(W),}\strong{\eqn{\lambda}}\eqn{_{2},...,}\strong{\eqn{\lambda}}\eqn{_{M},\alpha_{1},...,\alpha_{M-1})}, where
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
      \item{\strong{If \eqn{\lambda_{mi}} are constrained:}}{Replace \strong{\eqn{\lambda}}\eqn{_{2},...,}\strong{\eqn{\lambda}}\eqn{_{M}}
       with \strong{\eqn{\gamma}} \eqn{(rx1)} that satisfies (\strong{\eqn{\lambda}}\eqn{_{2}}\eqn{,...,}
      \strong{\eqn{\lambda}}\eqn{_{M}) =} \strong{\eqn{C_{\lambda} \gamma}} where \eqn{C_{\lambda}} is a \eqn{(d(M-1) x r)}
       constraint matrix.}
    }
  }
}
Above, \eqn{\phi_{m,0}} is the intercept parameter, \eqn{A_{m,i}} denotes the \eqn{i}th coefficient matrix of the \eqn{m}th
mixture component, \eqn{\Omega_{m}} denotes the error term covariance matrix of the \eqn{m}:th mixture component, and
\eqn{\alpha_{m}} is the mixing weight parameter. The \eqn{W} and \eqn{\lambda_{mi}} are structural parameters replacing the
error term covariance matrices (see Virolainen, 2020). If \eqn{M=1}, \eqn{\alpha_{m}} and \eqn{\lambda_{mi}} are dropped.
If \code{parametrization=="mean"}, just replace each \eqn{\phi_{m,0}} with regimewise mean \eqn{\mu_{m}}.
\eqn{vec()} is vectorization operator that stacks columns of a given matrix into a vector. \eqn{vech()} stacks columns
of a given matrix from the principal diagonal downwards (including elements on the diagonal) into a vector.
The notation is in line with the cited article by \emph{Kalliovirta, Meitz and Saikkonen (2016)} introducing the GMVAR model.}

\item{m}{which component?}

\item{structural_pars}{If \code{NULL} a reduced form model is considered. For structural model, should be a list containing
the following elements:
\itemize{
  \item \code{W} - a \eqn{(dxd)} matrix with its entries imposing constraints on \eqn{W}: \code{NA} indicating that the element is
    unconstrained, a positive value indicating strict positive sign constraint, a negative value indicating strict
    negative sign constraint, and zero indicating that the element is constrained to zero.
  \item \code{C_lambda} - a \eqn{(d(M-1) x r)} constraint matrix that satisfies (\strong{\eqn{\lambda}}\eqn{_{2}}\eqn{,...,}
    \strong{\eqn{\lambda}}\eqn{_{M}) =} \strong{\eqn{C_{\lambda} \gamma}} where \strong{\eqn{\gamma}} is the new \eqn{(r x 1)}
    parameter subject to which the model is estimated (similarly to AR parameter constraints). The entries of \code{C_lambda}
    must be either \strong{positive} or \strong{zero}. Ignore (or set to \code{NULL}) if the eigenvalues \eqn{\lambda_{mi}}
    should not be constrained.
}
See Virolainen (2020) for the conditions required to identify the shocks and for the B-matrix as well (it is \eqn{W} times
a time-varying diagonal matrix with positive diagonal entries).}
}
\value{
Returns a 3D array containing the coefficient matrices of the given component.
 A coefficient matrix \eqn{A_{m,i}} can be obtained by choosing \code{[, , i]}.
}
\description{
\code{pick_Am} picks the coefficient matrices \eqn{A_{m,i} (i=1,..,p)}
  from the given parameter vector so that they are arranged in a 3D array with the
  third dimension indicating each lag.
}
\details{
Does not support constrained parameter vectors.
}
\section{Warning}{

 No argument checks!
}

\references{
\itemize{
   \item Kalliovirta L., Meitz M. and Saikkonen P. 2016. Gaussian mixture vector autoregression.
         \emph{Journal of Econometrics}, \strong{192}, 485-498.
   \item Virolainen S. 2020. Structural Gaussian mixture vector autoregressive model. Unpublished working
     paper, available as arXiv:2007.04713.
 }
}
