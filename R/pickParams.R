
#' @title Pick coefficient matrix
#'
#' @description \code{pick_Ami} picks the coefficient matrix \eqn{A_{m,i}} from the given parameter vector.
#'
#' @inheritParams loglikelihood_int
#' @param d number of time series in the system, i.e. the dimension.
#' @param m which component?
#' @param i which lag in 1,...,p?
#' @param unvec if \code{FALSE} then vectorized version of \eqn{A_{m,i}} will be returned instead of matrix.
#'   Default if \code{TRUE}.
#' @details Does not support constrained parameter vectors.
#' @return Returns the i:th lag coefficient matrix of m:th component, \eqn{A_{m,i}}.
#' @section Warning:
#'  No argument checks!
#' @inherit is_stationary references
#' @keywords internal

pick_Ami <- function(p, M, d, params, m, i, structural_pars=NULL, unvec=TRUE) {
  if(is.null(structural_pars)) {
    qm1 <- (m - 1)*(d + p*d^2 + d*(d + 1)/2)
    Ami <- params[(qm1 + d + (i - 1)*d^2 + 1):(qm1 + d + i*d^2)]
  } else {
    qm1 <- d*M + d^2*p*(m - 1)
    Ami <- params[(qm1 + d^2*(i - 1) + 1):(qm1 + d^2*i)]
  }
  if(unvec == TRUE) {
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
#' @return Returns a 3D array containing the coefficient matrices of the given component.
#'  A coefficient matrix \eqn{A_{m,i}} can be obtained by choosing \code{[, , i]}.
#' @inherit pick_Ami details
#' @inheritSection pick_Ami Warning
#' @inherit is_stationary references
#' @keywords internal

pick_Am <- function(p, M, d, params, m, structural_pars=NULL) {
  if(is.null(structural_pars)) {
    qm1 <- (m - 1)*(d + p*d^2 + d*(d + 1)/2)
    lowers <- qm1 + d + (1:p - 1)*d^2 + 1
    wd <- d^2 # How many params
    tmp <- matrix(0:(wd - 1), nrow=wd, ncol=length(lowers), byrow=FALSE)
    coefs <- params[tmp + matrix(rep(lowers, times=wd), nrow=wd, byrow=TRUE)]
  } else {
    coefs <- params[(d*M + d^2*p*(m - 1) + 1):(d*M + d^2*p*m)]
  }
  array(coefs, dim=c(d, d, p))
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
#' @inherit pick_Ami details
#' @inheritSection pick_Ami Warning
#' @inherit is_stationary references
#' @keywords internal

pick_allA <- function(p, M, d, params, structural_pars=NULL) {
  if(is.null(structural_pars)) {
    qm1 <- (1:M - 1)*(d + p*d^2 + d*(d + 1)/2)
    lowers <- qm1 + d + 1
    wd <- d^2*p
    tmp <- matrix(0:(wd - 1), nrow=wd, ncol=length(lowers), byrow=FALSE)
    coefs <- params[tmp + matrix(rep(lowers, times=wd), nrow=wd, byrow=TRUE)]
  } else {
    coefs <- params[(d*M + 1):(d*M + d^2*p*M)]
  }
  array(coefs, dim=c(d, d, p, M))
}


#' @title Pick \eqn{\phi_{m,0}} or \eqn{\mu_{m}}, m=1,..,M vectors
#'
#' @description \code{pick_phi0} picks the intercept or mean parameters from the given parameter vector.
#'
#' @inheritParams is_stationary
#' @return Returns a \eqn{(dxM)} matrix containing \eqn{\phi_{m,0}} in the m:th column or
#'   \eqn{\mu_{m}} if the parameter vector is mean-parametrized, \eqn{, m=1,..,M}.
#' @inherit pick_Ami details
#' @inheritSection pick_Ami Warning
#' @inherit is_stationary references
#' @keywords internal

pick_phi0 <- function(p, M, d, params, structural_pars=NULL) {
  if(is.null(structural_pars)) {
    qm1 <- (1:M - 1)*(d + p*d^2 + d*(d + 1)/2)
    tmp <- matrix(1:d, nrow=d, ncol=length(qm1), byrow=FALSE)
    coefs <- params[tmp + matrix(rep(qm1, times=d), nrow=d, byrow=TRUE)]
  } else {
    coefs <- params[1:(d*M)]
  }
  matrix(coefs, nrow=d, byrow=FALSE)
}


#' @title Pick all \eqn{\phi_{m,0}} or \eqn{\mu_{m}} and \eqn{A_{m,1},...,A_{m,p}} parameter values
#'
#' @description \code{pick_all_phi0_A} picks the intercept or mean parameters and vectorized coefficient
#'   matrices from the given parameter vector.
#'
#' @inheritParams is_stationary
#' @return Returns a \eqn{((pd^2+d)xM)} matrix containing \eqn{(\phi_{m,0}, vec(A_{m,1}),...,vec(A_{m,p}))} in the m:th column,
#'  or \eqn{(\mu_{m}, vec(A_{m,1}),...,vec(A_{m,p}))} if the parameter vector is mean-parametrized, m=1,..,M.
#' @inherit pick_Ami details
#' @inheritSection pick_Ami Warning
#' @inherit is_stationary references
#' @keywords internal

pick_all_phi0_A <- function(p, M, d, params, structural_pars=NULL) {
  if(is.null(structural_pars)) {
    q0 <- d + p*d^2
    qm1 <- (1:M - 1)*(d + p*d^2 + d*(d + 1)/2)
    ret <- vapply(1:M, function(m) params[(qm1[m] + 1):(qm1[m] + q0)], numeric(q0))
  } else {
    all_phi0 <- matrix(params[1:(d*M)], nrow=d, byrow=FALSE)
    all_A <- matrix(params[(d*M + 1):(d*M + d^2*p*M)], nrow=d^2*p, byrow=FALSE)
    ret <- rbind(all_phi0, all_A)
  }
  ret
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
#' @inherit pick_Ami details
#' @inheritSection pick_Ami Warning
#' @inherit in_paramspace_int references
#' @keywords internal

pick_Omegas <- function(p, M, d, params, structural_pars=NULL) {
  Omegas <- array(dim=c(d, d, M))
  if(is.null(structural_pars)) {
    qm1 <- (1:M - 1)*(d + p*d^2 + d*(d + 1)/2)
    for(m in 1:M) {
      Omegas[, , m] <- unvech(d=d, a=params[(qm1[m] + d + p*d^2 + 1):(qm1[m] + d + p*d^2 + d*(d + 1)/2)])
    }
  } else {
    W <- unvec(d=d, a=params[(d*M*(1 + d*p) + 1):(d*M*(1 + d*p) + d^2)])
    Omegas[, , 1] <- tcrossprod(W)
    if(M > 1) {
      for(m in 2:M) {
        lambdas <- params[(d*M*(1 + d*p) + d^2 + d*(m - 2) + 1):(d*M*(1 + d*p) + d^2 + d*(m - 1))]
        Omegas[, , m] <- W%*%tcrossprod(diag(lambdas), W)
      }
    }
  }
  Omegas
}


#' @title Pick mixing weight parameters \eqn{\alpha_{m}, m=1,...,M}
#'
#' @description \code{pick_alphas} picks the mixing weight parameters from the given parameter vector.
#'
#' @inheritParams is_stationary
#' @return Returns a length M vector containing the mixing weight parameters \eqn{alpha_{m}, m=1,...,M},
#'   including non-parametrized \eqn{alpha_{M}}.
#' @inheritSection pick_Ami Warning
#' @inherit in_paramspace_int references
#' @keywords internal

pick_alphas <- function(p, M, d, params) {
  if(M == 1) {
    return(1)
  } else {
    alphas <- params[(length(params) - M + 2):length(params)]
    return(c(alphas, 1 - sum(alphas)))
  }
}


#' @title Pick the structural parameter matrix W
#'
#' @description \code{pick_W} picks the structural parameter matrix W from a parameter vector
#'
#' @inheritParams is_stationary
#' @details Constrained parameter vectors are not supported. Not even constraints in \eqn{W}!
#' @return Returns a \eqn{(d x d)} matrix \eqn{W} from a parameter vector of a SGMVAR model.
#'   Returns \code{NULL} for reduced form models.
#' @inheritSection pick_Ami Warning
#' @inherit in_paramspace_int references
#' @keywords internal

pick_W <- function(p, M, d, params, structural_pars=NULL) {
  if(is.null(structural_pars)) return(NULL)
  unvec(d=d, a=params[(M*d + d^2*p*M + 1):(M*d + d^2*p*M + d^2)])
}


#' @title Pick the structural parameters eigenvalue 'lambdas'
#'
#' @description \code{pick_lambdas} picks the structural parameters eigenvalue 'lambdas from a parameter vector
#'
#' @inheritParams is_stationary
#' @return Returns a length \eqn{(d*(M - 1))} vector \eqn{(\lambda_{2},...,\lambda_{M})}
#'  (see the argument \code{params}) from a parameter vector of a SGMVAR model.
#'   Returns \code{numeric(0)} for reduced form models or when \eqn{M=1}.
#' @inherit pick_W details
#' @inheritSection pick_Ami Warning
#' @inherit in_paramspace_int references
#' @keywords internal

pick_lambdas <- function(p, M, d, params, structural_pars=NULL) {
  if(is.null(structural_pars) || M == 1) return(numeric(0))
  params[(M*d + d^2*p*M + d^2 + 1):((M*d + d^2*p*M + d^2 + d*(M - 1)))]
}


#' @title Pick regime parameters \strong{\eqn{\upsilon_{m}}}\eqn{ = (\phi_{m,0},}\strong{\eqn{\phi_{m}}}\eqn{,\sigma_{m})}
#'
#' @description \code{pick_regime} picks the regime-parameters from the given parameter vector.
#'
#' @inheritParams pick_Am
#' @details Models with AR, mean, or lambda parameter constraints are currently not supported.
#' @return
#'   \describe{
#'     \item{For reduced form models:}{returns length \eqn{pd^2+d+d(d+1)/2} vector containing
#'       \strong{\eqn{\upsilon_{m}}}\eqn{ = (\phi_{m,0},}\strong{\eqn{\phi_{m}}}\eqn{,\sigma_{m})}, where
#'       \strong{\eqn{\phi_{m}}}\eqn{ = (vec(A_{m,1}),...,vec(A_{m,1})} and \eqn{\sigma_{m} = vech(\Omega_{m})}.}
#'     \item{For structural models:}{returns the length \eqn{pd^2 + d} vector \eqn{(\phi_{m,0},}\strong{\eqn{\phi_{m}}}\eqn{)}.}
#'   }
#' @inheritSection pick_Ami Warning
#' @inherit is_stationary references
#' @keywords internal

pick_regime <- function(p, M, d, params, m, structural_pars=NULL) {
  if(is.null(structural_pars)) {
    qm1 <- (m - 1)*(d + p*d^2 + d*(d + 1)/2)
    return(params[(qm1 + 1):(qm1 + d + p*d^2 + d*(d + 1)/2)])
  } else {
    n_zeros <- sum(structural_pars$W == 0, na.rm=TRUE)
    phi0 <- pick_phi0(p=p, M=M, d=d, params=params, structural_pars=structural_pars)[,m]
    all_Am <- pick_Am(p=p, M=M, d=d, params=params, m=m, structural_pars=structural_pars)
    return(c(phi0, as.vector(all_Am)))
  }
}


#' @title Calculate absolute values of the eigenvalues of the "bold A" matrices containing the AR coefficients
#'
#' @description \code{get_boldA_eigens} calculates absolute values of the eigenvalues of
#'   the "bold A" matrices containing the AR coefficients for each mixture component.
#'
#' @inheritParams simulateGMVAR
#' @return Returns a matrix with \eqn{d*p} rows and \eqn{M} columns - one column for each regime.
#'  The \eqn{m}th column contains the absolute values (or modulus) of the eigenvalues of the "bold A" matrix containing
#'  the AR coefficients correspinding to regime \eqn{m}.
#' @inherit is_stationary references
#' @examples
#' # GMVAR(2, 2), d=2 model
#' params22 <- c(0.36, 0.121, 0.223, 0.059, -0.151, 0.395, 0.406, -0.005,
#'  0.083, 0.299, 0.215, 0.002, 0.03, 0.484, 0.072, 0.218, 0.02, -0.119,
#'   0.722, 0.093, 0.032, 0.044, 0.191, 1.101, -0.004, 0.105, 0.58)
#' mod22 <- GMVAR(p=2, M=2, d=2, params=params22)
#' get_boldA_eigens(mod22)
#' @export

get_boldA_eigens <- function(gmvar) {
  check_gmvar(gmvar)
  p <- gmvar$model$p
  M <- gmvar$model$M
  d <- gmvar$model$d
  params <- reform_constrained_pars(p=p, M=M, d=d, params=gmvar$params,
                                    constraints=gmvar$model$constraints,
                                    same_means=gmvar$model$same_means,
                                    structural_pars=gmvar$model$structural_pars)
  structural_pars <- get_unconstrained_structural_pars(structural_pars=gmvar$model$structural_pars)
  all_A <- pick_allA(p=p, M=M, d=d, params=params, structural_pars=structural_pars)
  all_boldA <- form_boldA(p=p, M=M, d=d, all_A=all_A)
  matrix(vapply(1:M, function(m) abs(eigen(all_boldA[, , m], symmetric=FALSE, only.values=TRUE)$'values'), numeric(d*p)),
         nrow=d*p, ncol=M, byrow=FALSE)
}


#' @title Calculate the eigenvalues of the "Omega" error term covariance matrices
#'
#' @description \code{get_omega_eigens} calculates the eigenvalues of the "Omega" error
#'  term covariance matrices for each mixture component.
#'
#' @inheritParams simulateGMVAR
#' @return Returns a matrix with \eqn{d} rows and \eqn{M} columns - one column for each regime.
#'  The \eqn{m}th column contains the eigenvalues of the "Omega" error term covariance matrix
#'  of the \eqn{m}th regime.
#' @inherit is_stationary references
#' @examples
#' # GMVAR(2, 2), d=2 model
#' params22 <- c(0.36, 0.121, 0.223, 0.059, -0.151, 0.395, 0.406, -0.005,
#'  0.083, 0.299, 0.215, 0.002, 0.03, 0.484, 0.072, 0.218, 0.02, -0.119,
#'   0.722, 0.093, 0.032, 0.044, 0.191, 1.101, -0.004, 0.105, 0.58)
#' mod22 <- GMVAR(p=2, M=2, d=2, params=params22)
#' get_omega_eigens(mod22)
#' @export

get_omega_eigens <- function(gmvar) {
  check_gmvar(gmvar)
  p <- gmvar$model$p
  M <- gmvar$model$M
  d <- gmvar$model$d
  params <- reform_constrained_pars(p=p, M=M, d=d, params=gmvar$params,
                                    constraints=gmvar$model$constraints,
                                    same_means=gmvar$model$same_means,
                                    structural_pars=gmvar$model$structural_pars)
  structural_pars <- get_unconstrained_structural_pars(structural_pars=gmvar$model$structural_pars)
  all_Omega <- pick_Omegas(p=p, M=M, d=d, params=params, structural_pars=structural_pars)
  matrix(vapply(1:M, function(m) eigen(all_Omega[, , m], symmetric=TRUE, only.values=TRUE)$'values', numeric(d)),
         nrow=d, ncol=M, byrow=FALSE)
}

