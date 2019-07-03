
#' @title Pick coefficient matrix
#'
#' @description \code{pick_Ami} picks the coefficient matrix \eqn{A_{m,i}} from the given parameter vector
#'
#' @param d number of time series in the system, i.e. the dimension.
#' @param m which component?
#' @param i which lag in 1,...,p?
#' @param unvec if \code{FALSE} then vectorized version of \eqn{A_{m,i}} will be returned instead of matrix.
#'   Default if \code{TRUE}.
#' @inheritParams loglikelihood_int
#' @return Returns the i:th lag coefficient matrix of m:th component, \eqn{A_{m,i}}.
#' @section Warning:
#'  No argument checks!
#' @inherit is_stationary references

pick_Ami <- function(p, M, d, params, m, i, unvec=TRUE) {
  qm1 <- (m-1)*(d + p*d^2 + d*(d + 1)/2)
  Ami <- params[(qm1 + d + (i - 1)*d^2 + 1):(qm1 + d + i*d^2)]
  if(unvec==TRUE) {
    return(unvec(d=d, a=Ami))
  } else {
    return(Ami)
  }
}


#' @title Pick coefficient matrices
#'
#' @description \code{pick_Am} picks the coefficient matrices \eqn{A_{m,i} (i=1,..,p)}
#'   from the given parameter vector so that they are arranged in a 3D array with the
#'   third dimension indicating each lag.
#'
#' @inheritParams pick_Ami
#' @return Returns a 3D array containing the coefficient matrices of the given component. Coefficient matrix
#'  \eqn{A_{m,i}} can be obtained by choosing \code{[, , i]}.
#' @section Warning:
#'  No argument checks!
#' @inherit is_stationary references

pick_Am <- function(p, M, d, params, m) {
  qm1 <- (m-1)*(d + p*d^2 + d*(d + 1)/2)
  array(vapply(1:p, function(i1) params[(qm1 + d + (i1 - 1)*d^2 + 1):(qm1 + d + i1*d^2)], numeric(d^2)), dim=c(d, d, p))
}


#' @title Pick coefficient all matrices
#'
#' @description \code{pick_allA} picks all coefficient matrices \eqn{A_{m,i} (i=1,..,p, m=1,..,M)}
#'   from the given parameter vector so that they are arranged in a 4D array with the fourth dimension
#'   indicating each component and third dimension indicating each lag.
#'
#' @inheritParams is_stationary
#' @return Returns a 4D array containing the coefficient matrices of the all components. Coefficient matrix
#'  \eqn{A_{m,i}} can be obtained by choosing \code{[, , i, m]}.
#' @section Warning:
#'  No argument checks!
#' @inherit is_stationary references

pick_allA <- function(p, M, d, params) {
  qm1 <- (1:M-1)*(d + p*d^2 + d*(d + 1)/2)
  tmp <- vapply(1:M, function(m) {
    vapply(1:p, function(i1) params[(qm1[m] + d + (i1 - 1)*d^2 + 1):(qm1[m] + d + i1*d^2)], numeric(d^2))
  }, numeric(p*d^2))
  array(tmp, dim=c(d, d, p, M))
}


#' @title Pick \eqn{\phi_{m,0}} or \eqn{\mu_{m}}, m=1,..,M vectors from the given parameter vector
#'
#' @description \code{pick_phi0} picks the intercept or mean parameters from the given parameter vector.
#'
#' @inheritParams is_stationary
#' @return Returns a \eqn{(dxM)} matrix containing \eqn{\phi_{m,0}} in the m:th column or
#'   \eqn{\mu_{m}} if the parameter vector is mean-parametrized, \eqn{, m=1,..,M}.
#' @section Warning:
#'  No argument checks!
#' @inherit is_stationary references

pick_phi0 <- function(p, M, d, params) {
  qm1 <- (1:M-1)*(d + p*d^2 + d*(d + 1)/2)
  vapply(1:M, function(m) params[(qm1[m] + 1):(qm1[m] + d)], numeric(d))
}


#' @title Pick all \eqn{\phi_{m,0}} or \eqn{\mu_{m}} and \eqn{A_{m}} parameter values from the given parameter vector.
#'
#' @description \code{pick_all_phi0_A} picks the intercept or mean parameters and vectorized coefficient
#'   matrices from the given parameter vector.
#'
#' @inheritParams is_stationary
#' @return Returns a \eqn{((pd^2+d)xM)} matrix containing \eqn{(\phi_{m,0}, vec(A_{m}))} in the m:th column,
#'  or \eqn{(\mu_{m}, vec(A_{m}))} if the parameter vector is mean-parametrized, m=1,..,M.
#' @section Warning:
#'  No argument checks!
#' @inherit is_stationary references

pick_all_phi0_A <- function(p, M, d, params) {
  q0 <- d + p*d^2
  qm1 <- (1:M-1)*(d + p*d^2 + d*(d + 1)/2)
  vapply(1:M, function(m) params[(qm1[m] + 1):(qm1[m] + q0)], numeric(q0))
}



#' @title Pick covariance matrices
#'
#' @description \code{pick_Omegas} picks the covariance matrices \eqn{\Omega_{m} (m=1,..,M)}
#'  from the given parameter vector so that they are arranged in a 3D array with the third
#'  dimension indicating each component.
#'
#' @inheritParams is_stationary
#' @return Returns a 3D array containing the covariance matrices of the given model. Coefficient matrix
#'  \eqn{\Omega_{m}} can be obtained by choosing \code{[, , m]}.
#' @section Warning:
#'  No argument checks!
#' @inherit in_paramspace_int references

pick_Omegas <- function(p, M, d, params) {
  qm1 <- (1:M-1)*(d + p*d^2 + d*(d + 1)/2)
  tmp <- vapply(1:M, function(m) unvech(d=d, a=params[(qm1[m] + d + p*d^2 + 1):(qm1[m] + d + p*d^2 + d*(d + 1)/2)]), numeric(d^2))
  array(tmp, dim=c(d, d, M))
}


#' @title Pick mixing weight parameters \eqn{\alpha_{m}, m=1,...,M} from the given parameter vector.
#'
#' @description \code{pick_alphas} picks the mixing weight parameters from the given parameter vector.
#'
#' @inheritParams is_stationary
#' @return Returns length M vector containing the mixing weight parameters \eqn{alpha_{m}, m=1,...,M},
#'   including non-parametrized \eqn{alpha_{M}}.
#' @section Warning:
#'  No argument checks!
#' @inherit in_paramspace_int references


pick_alphas <- function(p, M, d, params) {
  if(M==1) {
    return(1)
  } else {
    qM <- M*(d + p*d^2 + d*(d + 1)/2)
    alphas <- params[(qM + 1):(qM + M - 1)]
    return(c(alphas, 1-sum(alphas)))
  }
}



#' @title Pick regime parameters \strong{\eqn{\upsilon_{m}}}\eqn{ = (\phi_{m,0},}\strong{\eqn{\phi_{m}}}\eqn{,\sigma_{m})}
#'   from the given parameter vector.
#'
#' @description \code{pick_regime} picks the regime-parameters from the given parameter vector.
#'
#' @inheritParams pick_Am
#' @return Returns length \eqn{pd^2+d+d(d+1)/2} vector containing
#'  \strong{\eqn{\upsilon_{m}}}\eqn{ = (\phi_{m,0},}\strong{\eqn{\phi_{m}}}\eqn{,\sigma_{m})}, where
#'  \strong{\eqn{\phi_{m}}}\eqn{ = (vec(A_{m,1}),...,vec(A_{m,1})} and \eqn{\sigma_{m} = vech(\Omega_{m})}.
#' @section Warning:
#'  No argument checks!
#' @inherit is_stationary references

pick_regime <- function(p, M, d, params, m) {
  qm1 <- (m-1)*(d + p*d^2 + d*(d + 1)/2)
  params[(qm1 + 1):(qm1 + d + p*d^2 + d*(d + 1)/2)]
}


#' @title Calculate absolute values of the eigenvalues of the "bold A" matrices containing the AR coefficients
#'
#' @description \code{get_boldA_eigens} calculates absolute values of the eigenvalues of
#'   the "bold A" matrices containing the AR coefficients for each mixture component
#'
#' @inheritParams simulateGMVAR
#' @return Returns a list with \eqn{M} elements - one for each regime. Each element contains
#'  the absolute values (or modulus) of the eigenvalues of the "bold A" matrix containing
#'  the AR coefficients.
#' @inherit is_stationary references
#' @examples
#' params222 <- c(-11.904, 154.684, 1.314, 0.145, 0.094, 1.292, -0.389,
#'  -0.070, -0.109, -0.281, 0.920, -0.025, 4.839, 11.633, 124.983, 1.248,
#'   0.077, -0.040, 1.266, -0.272, -0.074, 0.034, -0.313, 5.855, 3.570,
#'   9.838, 0.740)
#' mod222 <- GMVAR(d=2, p=2, M=2, params=params222, parametrization="mean")
#' get_boldA_eigens(mod222)
#' @export

get_boldA_eigens <- function(gmvar) {
  check_gmvar(gmvar)
  p <- gmvar$model$p
  M <- gmvar$model$M
  d <- gmvar$model$d
  params <- reform_constrained_pars(p=p, M=M, d=d, params=gmvar$params, constraints=gmvar$model$constraints)
  all_A <- pick_allA(p=p, M=M, d=d, params=params)
  all_boldA <- form_boldA(p=p, M=M, d=d, all_A=all_A)
  lapply(1:M, function(m) abs(eigen(all_boldA[, , m], symmetric=FALSE, only.values=TRUE)$'values'))
}


#' @title Calculate the eigenvalues of the "Omega" error term covariance matrices
#'
#' @description \code{get_omega_eigens} the eigenvalues of the "Omega" error term covariance matrices
#'   for each mixture component
#'
#' @inheritParams simulateGMVAR
#' @return Returns a list with \eqn{M} elements - one for each regime. Each element contains
#'  the eigenvalues of the "Omega" error term covariance matrix.
#' @inherit is_stationary references
#' @examples
#' params222 <- c(-11.904, 154.684, 1.314, 0.145, 0.094, 1.292, -0.389,
#'  -0.070, -0.109, -0.281, 0.920, -0.025, 4.839, 11.633, 124.983, 1.248,
#'   0.077, -0.040, 1.266, -0.272, -0.074, 0.034, -0.313, 5.855, 3.570,
#'   9.838, 0.740)
#' mod222 <- GMVAR(d=2, p=2, M=2, params=params222, parametrization="mean")
#' get_omega_eigens(mod222)
#' @export

get_omega_eigens <- function(gmvar) {
  check_gmvar(gmvar)
  p <- gmvar$model$p
  M <- gmvar$model$M
  d <- gmvar$model$d
  params <- reform_constrained_pars(p=p, M=M, d=d, params=gmvar$params, constraints=gmvar$model$constraints)
  all_Omega <- pick_Omegas(p=p, M=M, d=d, params=params)
  lapply(1:M, function(m) eigen(all_Omega[, , m], symmetric=TRUE, only.values=TRUE)$'values')
}

