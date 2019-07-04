#' @title Reform data
#'
#' @description \code{reform_data} reforms the data into a form that is
#'   easier to use when calculating log-likelihood values etc.
#'
#' @inheritParams loglikelihood_int
#' @return Returns the data reformed into a \eqn{((n_obs-p+1)x(dp))} matrix. The i:th row
#'   of the matrix contains the vector \eqn{(y_{i-1}',...,y_{i-p}')} \eqn{((dp)x1)}, where
#'   \eqn{y_{i}=(y_{1i},...,y_{di})} \eqn{(dx1)}.
#' @section Warning:
#'  No argument checks!

reform_data <- function(data, p) {
  d <- ncol(data)
  n_obs <- nrow(data)
  T_obs <- n_obs - p
  matrix(vapply(1:p, function(i1) as.vector(data[(p-i1+1):(T_obs + p - i1 + 1),]), numeric((n_obs-p+1)*d)), nrow=n_obs-p+1, byrow=FALSE)
}


#' @title Reform constrained parameter vector into the "standard" form
#'
#' @description \code{reform_constrained_pars} reforms constrained parameter vector
#'   into the form that corresponds to unconstrained "regular" parameter vectors.
#'
#' @inheritParams loglikelihood_int
#' @inheritParams is_stationary
#' @param change_na change NA parameter values of constrained models to -9.999?
#' @return Returns "regular model" parameter vector corresponding to the constraints.
#' @section Warning:
#'  No argument checks!
#' @inherit in_paramspace_int references


reform_constrained_pars <- function(p, M, d, params, constraints=NULL, change_na=FALSE) {
  if(is.null(constraints)) {
    return(params)
  }
  q <- ncol(constraints)
  psi <- params[(M*d+1):(M*d+q)]
  if(change_na==TRUE) {
    if(length(psi[is.na(psi)]) > 0) warning("Replaced some NA values with -9.999")
    psi[is.na(psi)] <- -9.999
  }
  psi_expanded <- constraints%*%psi
  pars <- as.vector(vapply(1:M, function(m) c(params[((m - 1)*d + 1):(m*d)], psi_expanded[((m - 1)*p*d^2 + 1):(m*p*d^2)],
                                              params[(M*d + q + (m-1)*d*(d + 1)/2 + 1):(M*d + q + m*d*(d + 1)/2)]),
                    numeric(p*d^2 + d + d*(d + 1)/2)))
  if(M==1) {
    return(pars)
  } else {
    return(c(pars, params[(M*d + q + M*d*(d + 1)/2 + 1):(M*d + q + M*d*(d + 1)/2 + M - 1)]))
  }
}


#' @title Form the \eqn{((dp)x(dp))}  "bold A" matrices related to the VAR processes
#'
#' @description \code{form_boldA} creates the "bold A" coefficient matrices related to
#'   VAR processes.
#'
#' @inheritParams pick_allA
#' @param all_A 4D array containing all coefficient matrices \eqn{A_{m,i}}, obtained from \code{pick_allA}.
#' @return Returns 3D array containing the \eqn{((dp)x(dp))} "bold A" matrices related to each component VAR-process.
#'  The matrix \strong{\eqn{A_{m}}} can be obtained by choosing \code{[, , m]}
#' @section Warning:
#'  No argument checks!
#' @inherit is_stationary references

form_boldA <- function(p, M, d, all_A) {
  I_all <- diag(nrow=d*(p-1))
  ZER_all <- matrix(0, nrow=d*(p-1), ncol=d)
  array(vapply(1:M, function(m) rbind(matrix(all_A[, , 1:p, m], nrow=d, byrow=FALSE), cbind(I_all, ZER_all)), numeric((d*p)^2)), dim=c(d*p, d*p, M))
}



#' @title Sort components in parameter vector by mixing weights into a decreasing order
#'
#' @description \code{sort_components} sorts mixture components in the parameter vector by
#'   by mixing weights into a decreasing order.
#'
#' @inheritParams is_stationary
#' @details Constrained parameter vectors are not supported!
#' @return Returns sorted parameter vector with \eqn{\alpha_{1}>...>\alpha_{m}}, that has form
#'   \strong{\eqn{\theta}}\eqn{ = }(\strong{\eqn{\upsilon_{1}}},...,\strong{\eqn{\upsilon_{M}}},
#'   \eqn{\alpha_{1},...,\alpha_{M-1}}), where:
#'  \itemize{
#'    \item \strong{\eqn{\upsilon_{m}}} \eqn{ = (\phi_{m,0},}\strong{\eqn{\phi_{m}}}\eqn{,\sigma_{m})}
#'    \item \strong{\eqn{\phi_{m}}}\eqn{ = (vec(A_{m,1}),...,vec(A_{m,1})}
#'    \item and \eqn{\sigma_{m} = vech(\Omega_{m})}, m=1,...,M.
#'  }
#'  Above \eqn{\phi_{m,0}} is the intercept parameter, \eqn{A_{m,i}} denotes the \eqn{i}:th coefficient matrix of the \eqn{m}:th
#'  component, \eqn{\Omega_{m}} denotes the error term covariance matrix of the \eqn{m}:th component and \eqn{\alpha_{m}} is the
#'  mixing weight parameter.
#'  \eqn{vec()} is vectorization operator that stack columns of the given matrix into a vector. \eqn{vech()} stacks columns
#'  of the given matrix from the principal diagonal downwards (including elements on the diagonal) to form a vector.
#'  The notations are in line with the cited article by KMS (2016)
#' @section Warning:
#'  No argument checks!
#' @inherit in_paramspace_int references

sort_components <- function(p, M, d, params) {
  alphas <- pick_alphas(p=p, M=M, d=d, params=params)
  ord <- order(alphas, decreasing=TRUE, method="radix")
  if(all(ord == 1:M)) {
    return(params)
  } else {
    q <- d + p*d^2 + d*(d + 1)/2
    qm1 <- (1:M-1)*q
    qm <- qm1[ord]
    pars <- vapply(1:M, function(m) params[(qm[m]+1):(qm[m]+q)], numeric(q))
    c(pars, alphas[ord][-M])
  }
}


#' @title Change parametrization of the parameter vector
#'
#' @description \code{change_parametrization()} changes the parametrization of the given parameter
#'   vector to \code{change_to}.
#'
#' @inheritParams loglikelihood_int
#' @inheritParams form_boldA
#' @param change_to either "intercept" or "mean" specifying to which parametrization it should be switched to.
#'   If set to \code{"intercept"}, it's assumed that \code{params} is mean-parametrized, and if set to \code{"mean"}
#'   it's assumed that \code{params} is intercept-parametrized.
#' @return Returns parameter vector described in \code{params}, but with parametrization changed from intercept to mean
#'   (when \code{change_to==mean}) or from mean to intercept (when \code{change_to==intercept}).
#' @section Warning:
#'  No argument checks!
#' @inherit is_stationary references

change_parametrization <- function(p, M, d, params, constraints=NULL, change_to=c("intercept", "mean")) {
  re_params <- params
  params <- reform_constrained_pars(p=p, M=M, d=d, params=params, constraints=constraints) # Parameters in regular form
  change_to <- match.arg(change_to)
  Id <- diag(nrow=d)
  all_A <- pick_allA(p=p, M=M, d=d, params=params)
  all_phi0_or_mu <- pick_phi0(p=p, M=M, d=d, params=params)

  if(is.null(constraints)) {
    qm1 <- (1:M-1)*(d + p*d^2 + d*(d + 1)/2)
  } else {
    qm1 <- (1:M-1)*d
  }

  if(change_to == "mean") { # params has original parametrization with intercept
    for(m in 1:M) {
      re_params[(qm1[m] + 1):(qm1[m] + d)] <- solve(Id - rowSums(all_A[, , , m, drop=FALSE], dims=2), all_phi0_or_mu[,m]) # Insert mu_m
    }
  } else { # mean parameters instead of phi0
    for(m in 1:M) {
      re_params[(qm1[m] + 1):(qm1[m] + d)] <- (Id - rowSums(all_A[, , , m, drop=FALSE], dims=2))%*%all_phi0_or_mu[,m] # Insert phi_{m,0}
    }
  }
  re_params
}


#' @title Change regime parameters \strong{\eqn{\upsilon_{m}}}\eqn{ = (\phi_{m,0},}\strong{\eqn{\phi_{m}}}\eqn{,\sigma_{m})}
#'   of the given parameter vector.
#'
#' @description \code{change_regime} changes the given regime parameters (excluding mixing weights parameter)
#'   to the given new parameters.
#'
#' @inheritParams pick_regime
#' @param regime_pars a size \eqn{((pd^2+d+d(d+1)/2)x1)} vector
#'   \strong{\eqn{\upsilon_{m}}}\eqn{ = (\phi_{m,0},}\strong{\eqn{\phi_{m}}}\eqn{,\sigma_{m})}.
#' @return Returns parameter vector with \code{m}:th regime changed to \code{regime_pars}.
#' @section Warning:
#'  No argument checks!
#' @inherit in_paramspace_int references


change_regime <- function(p, M, d, params, m, regime_pars) {
  qm1 <- (m-1)*(d + p*d^2 + d*(d + 1)/2)
  params[(qm1 + 1):(qm1 + d + p*d^2 + d*(d + 1)/2)] <- regime_pars
  params
}


#' @title Calculate "distance" between two (scaled) regimes \strong{\eqn{\upsilon_{m}}}\eqn{ = (\phi_{m,0},}\strong{\eqn{\phi_{m}}}\eqn{,\sigma_{m})}
#'
#' @description \code{regime_distance} calculates "distance" between two scaled regimes. This is used in
#'   the genetic algorithm.
#'
#' @param regime_pars1 a length \eqn{pd^2+d+d(d+1)/2} vector
#'   \strong{\eqn{\upsilon_{m}}}\eqn{ = (\phi_{m,0},}\strong{\eqn{\phi_{m}}}\eqn{,\sigma_{m})}.
#' @param regime_pars2 a length \eqn{pd^2+d+d(d+1)/2} vector
#'   \strong{\eqn{\upsilon_{m}}}\eqn{ = (\phi_{m,0},}\strong{\eqn{\phi_{m}}}\eqn{,\sigma_{m})}.
#' @return Returns "distance" between \code{regime_pars1} and \code{regime_pars2}. Values are scaled
#'   before calculating the "distance". Read the source code for more details.
#' @section Warning:
#'  No argument checks!
#' @inherit in_paramspace_int references

regime_distance <- function(regime_pars1, regime_pars2) {
  dist_fun <- function(x) {
    x <- abs(x)
    if(x < 1) {
      return(1)
    } else {
      return(10^ceiling(abs(log10(x))))
    }
  }
  scales1 <- vapply(regime_pars1, dist_fun, numeric(1))
  scales2 <- vapply(regime_pars2, dist_fun, numeric(1))
  sqrt(crossprod(regime_pars1/scales1 - regime_pars2/scales2))
}

