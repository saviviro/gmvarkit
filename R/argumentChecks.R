#' @title Check the stationary condition of a given GMVAR, StMVAR, or G-StMVAR model
#'
#' @description \code{is_stationary} checks the stationarity condition of a GMVAR, StMVAR, or G-StMVAR model.
#'
#' @inheritParams loglikelihood_int
#' @param d the number of time series in the system.
#' @param params a real valued vector specifying the parameter values.
#'   \describe{
#'     \item{\strong{For reduced form models:}}{
#'       Should be size \eqn{((M(pd^2+d+d(d+1)/2+2)-M1-1)x1)} and have the form
#'       \strong{\eqn{\theta}}\eqn{ = }(\strong{\eqn{\upsilon}}\eqn{_{1}},
#'       ...,\strong{\eqn{\upsilon}}\eqn{_{M}}, \eqn{\alpha_{1},...,\alpha_{M-1},}\strong{\eqn{\nu}}\eqn{)}, where
#'       \itemize{
#'         \item \strong{\eqn{\upsilon}}\eqn{_{m}} \eqn{ = (\phi_{m,0},}\strong{\eqn{\phi}}\eqn{_{m}}\eqn{,\sigma_{m})}
#'         \item \strong{\eqn{\phi}}\eqn{_{m}}\eqn{ = (vec(A_{m,1}),...,vec(A_{m,p})}
#'         \item and \eqn{\sigma_{m} = vech(\Omega_{m})}, m=1,...,M,
#'         \item \strong{\eqn{\nu}}\eqn{=(\nu_{M1+1},...,\nu_{M})}
#'         \item \eqn{M1} is the number of GMVAR type regimes.
#'       }
#'     }
#'     \item{\strong{For structural model:}}{
#'       Should have the form
#'       \strong{\eqn{\theta}}\eqn{ = (\phi_{1,0},...,\phi_{M,0},}\strong{\eqn{\phi}}\eqn{_{1},...,}\strong{\eqn{\phi}}\eqn{_{M},
#'       vec(W),}\strong{\eqn{\lambda}}\eqn{_{2},...,}\strong{\eqn{\lambda}}\eqn{_{M},\alpha_{1},...,\alpha_{M-1},}\strong{\eqn{\nu}}\eqn{)},
#'       where
#'       \itemize{
#'         \item\strong{\eqn{\lambda}}\eqn{_{m}=(\lambda_{m1},...,\lambda_{md})} contains the eigenvalues of the \eqn{m}th mixture component.
#'       }
#'     }
#'   }
#'   Above, \eqn{\phi_{m,0}} is the intercept parameter, \eqn{A_{m,i}} denotes the \eqn{i}th coefficient matrix of the \eqn{m}th
#'   mixture component, \eqn{\Omega_{m}} denotes the error term covariance matrix of the \eqn{m}:th mixture component, and
#'   \eqn{\alpha_{m}} is the mixing weight parameter. The \eqn{W} and \eqn{\lambda_{mi}} are structural parameters replacing the
#'   error term covariance matrices (see Virolainen, 2022). If \eqn{M=1}, \eqn{\alpha_{m}} and \eqn{\lambda_{mi}} are dropped.
#'   If \code{parametrization=="mean"}, just replace each \eqn{\phi_{m,0}} with regimewise mean \eqn{\mu_{m}}.
#'   \eqn{vec()} is vectorization operator that stacks columns of a given matrix into a vector. \eqn{vech()} stacks columns
#'   of a given matrix from the principal diagonal downwards (including elements on the diagonal) into a vector.
#'
#'   In the \strong{GMVAR model}, \eqn{M1=M} and \strong{\eqn{\nu}} is dropped from the parameter vector. In the \strong{StMVAR} model,
#'   \eqn{M1=0}. In the \strong{G-StMVAR} model, the first \code{M1} regimes are \emph{GMVAR type} and the rest \code{M2} regimes are
#'   \emph{StMVAR type}. In \strong{StMVAR} and \strong{G-StMVAR} models, the degrees of freedom parameters in \strong{\eqn{\nu}}
#'   # should be strictly larger than two.
#'
#'   The notation is similar to the cited literature.
#' @param all_boldA 3D array containing the \eqn{((dp)x(dp))} "bold A" matrices related to each mixture component VAR-process,
#'   obtained from \code{form_boldA}. Will be computed if not given.
#' @param tolerance Returns \code{FALSE} if modulus of any eigenvalue is larger or equal to \code{1-tolerance}.
#' @details If the model is constrained, remove the constraints first with the function \code{reform_constrained_pars}.
#' @return Returns \code{TRUE} if the model is stationary and \code{FALSE} if not. Based on the argument \code{tolerance},
#'   \code{is_stationary} may return \code{FALSE} when the parameter vector is in the stationarity region, but
#'   very close to the boundary (this is used to ensure numerical stability in estimation of the model parameters).
#' @section Warning:
#'  No argument checks!
#' @inherit in_paramspace_int references
#' @keywords internal

is_stationary <- function(p, M, d, params, all_boldA=NULL, structural_pars=NULL, tolerance=1e-3) {
  M <- sum(M)
  if(is.null(all_boldA)) {
    all_A <- pick_allA(p=p, M=M, d=d, params=params, structural_pars=structural_pars)
    all_boldA <- form_boldA(p=p, M=M, d=d, all_A=all_A)
  }
  for(m in 1:M) {
    if(any(abs(eigen(all_boldA[, , m], symmetric=FALSE, only.values=TRUE)$'values') >= 1 - tolerance)) {
      return(FALSE)
    }
  }
  TRUE
}


#' @title Determine whether the parameter vector lies in the parameter space
#'
#' @description \code{in_paramspace_int} checks whether the parameter vector lies in the parameter
#'   space.
#'
#' @inheritParams loglikelihood_int
#' @inheritParams is_stationary
#' @param alphas (Mx1) vector containing all mixing weight parameters, obtained from \code{pick_alphas}.
#' @param all_Omega 3D array containing all covariance matrices \eqn{\Omega_{m}}, obtained from \code{pick_Omegas}.
#' @param W_constraints set \code{NULL} for reduced form models. For structural models, this should be the
#'   constraint matrix \eqn{W} from the list of structural parameters.
#' @param stat_tol numerical tolerance for stationarity of the AR parameters: if the "bold A" matrix of any regime
#'   has eigenvalues larger that \code{1 - stat_tol} the model is classified as non-stationary. Note that if the
#'   tolerance is too small, numerical evaluation of the log-likelihood might fail and cause error.
#' @param posdef_tol numerical tolerance for positive definiteness of the error term covariance matrices: if
#'   the error term covariance matrix of any regime has eigenvalues smaller than this, the model is classified
#'   as not satisfying positive definiteness assumption. Note that if the tolerance is too small, numerical
#'   evaluation of the log-likelihood might fail and cause error.
#' @param df_tol the parameter vector is considered to be outside the parameter space if all degrees of
#'   freedom parameters are not larger than \code{2 + df_tol}.
#' @details The parameter vector in the argument \code{params} should be unconstrained and it is used for
#'   structural models only.
#' @return Returns \code{TRUE} if the given parameter values are in the parameter space and \code{FALSE} otherwise.
#'   This function does NOT consider the identifiability condition!
#' @references
#'  \itemize{
#'    \item Kalliovirta L., Meitz M. and Saikkonen P. 2016. Gaussian mixture vector autoregression.
#'          \emph{Journal of Econometrics}, \strong{192}, 485-498.
#'    \item Virolainen S. 2022. Structural Gaussian mixture vector autoregressive model with application to the asymmetric
#'      effects of monetary policy shocks. Unpublished working paper, available as arXiv:2007.04713.
#'    \item Virolainen S. 2022. Gaussian and Student's t mixture vector autoregressive model with application to the
#'      asymmetric effects of monetary policy shocks in the Euro area. Unpublished working
#'      paper, available as arXiv:2109.13648.
#'  }
#'  @keywords internal

in_paramspace_int <- function(p, M, d, params, model=c("GMVAR", "StMVAR", "G-StMVAR"), all_boldA, alphas, all_Omega,
                              W_constraints=NULL, stat_tol=1e-3, posdef_tol=1e-8, df_tol=1e-8) {

  model <- match.arg(model)
  if(model != "GMVAR") { # Check degrees of freedom parameters for StMVAR and G-StMVAR models
    all_df <- pick_df(M=M, params=params, model=model)
    if(any(all_df <= 2 + df_tol) || any(all_df > 1e+5)) {
      return(FALSE)
    }
  }
  M <- sum(M) # pick_W, pick_lambdas, and is_stationary all just use the number of mixture components, not M1 and M2 separately
  if(!is.null(W_constraints)) {
    W_pars <- pick_W(p=p, M=M, d=d, params=params, structural_pars=list(W=W_constraints))
    # No need to check zero constraints because the zeros are not parametrized
    if(any(W_pars[W_constraints < 0] > -1e-8, na.rm=TRUE)) {
      return(FALSE)
    } else if(any(W_pars[W_constraints > 0] < 1e-8, na.rm=TRUE)) {
      return(FALSE)
    }
    if(M > 1) {
      lambdas <- pick_lambdas(p=p, M=M, d=d, params=params, structural_pars=list(W=W_constraints))
      if(any(lambdas < 1e-8)) {
        return(FALSE)
      }
    }
  }
  if(M >= 2 & sum(alphas[-M]) >= 1) {
    return(FALSE)
  } else if(any(alphas <= 0)) {
    return(FALSE)
  } else if(!is_stationary(p=p, M=M, d=d, all_boldA=all_boldA, tolerance=stat_tol)) {
    return(FALSE)
  }
  for(m in 1:M) {
    if(any(eigen(all_Omega[, , m], symmetric=TRUE, only.values=TRUE)$values < posdef_tol)) {
      return(FALSE)
    }
  }
  TRUE
}


#' @title Determine whether the parameter vector lies in the parameter space
#'
#' @description \code{in_paramspace} checks whether the given parameter vector lies in
#'   the parameter space. Does NOT test the identification conditions!
#'
#' @inheritParams loglikelihood_int
#' @inheritParams is_stationary
#' @inheritParams in_paramspace_int
#' @return Returns \code{TRUE} if the given parameter vector lies in the parameter space
#'  and \code{FALSE} otherwise.
#' @inherit in_paramspace_int references
#' @examples
#' # GMVAR(1,1), d=2 model:
#' params11 <- c(1.07, 127.71, 0.99, 0.00, -0.01, 0.99, 4.05,
#'   2.22, 8.87)
#' in_paramspace(p=1, M=1, d=2, params=params11)
#'
#' # GMVAR(2,2), d=2 model:
#' params22 <- c(1.39, -0.77, 1.31, 0.14, 0.09, 1.29, -0.39,
#'  -0.07, -0.11, -0.28, 0.92, -0.03, 4.84, 1.01, 5.93, 1.25,
#'   0.08, -0.04, 1.27, -0.27, -0.07, 0.03, -0.31, 5.85, 3.57,
#'   9.84, 0.74)
#' in_paramspace(p=2, M=2, d=2, params=params22)
#'
#' # GMVAR(2,2), d=2 model with AR-parameters restricted to be
#' # the same for both regimes:
#' C_mat <- rbind(diag(2*2^2), diag(2*2^2))
#' params22c <- c(1.03, 2.36, 1.79, 3.00, 1.25, 0.06,0.04,
#'  1.34, -0.29, -0.08, -0.05, -0.36, 0.93, -0.15, 5.20,
#'  5.88, 3.56, 9.80, 0.37)
#' in_paramspace(p=2, M=2, d=2, params=params22c, constraints=C_mat)
#'
#' # Structural GMVAR(2, 2), d=2 model identified with sign-constraints:
#' params22s <- c(1.03, 2.36, 1.79, 3, 1.25, 0.06, 0.04, 1.34, -0.29,
#'  -0.08, -0.05, -0.36, 1.2, 0.05, 0.05, 1.3, -0.3, -0.1, -0.05, -0.4,
#'   0.89, 0.72, -0.37, 2.16, 7.16, 1.3, 0.37)
#' W_22 <- matrix(c(1, 1, -1, 1), nrow=2, byrow=FALSE)
#' in_paramspace(p=2, M=2, d=2, params=params22s,
#'   structural_pars=list(W=W_22))
#' @export

in_paramspace <- function(p, M, d, params, model=c("GMVAR", "StMVAR", "G-StMVAR"), constraints=NULL, same_means=NULL,
                          weight_constraints=NULL, structural_pars=NULL, stat_tol=1e-3, posdef_tol=1e-8, df_tol=1e-8) {
  model <- match.arg(model)
  check_pMd(p=p, M=M, d=d, model=model)
  check_constraints(p=p, M=M, d=d, constraints=constraints, same_means=same_means, weight_constraints=weight_constraints,
                    structural_pars=structural_pars)
  if(length(params) != n_params(p=p, M=M, d=d, model=model, constraints=constraints, same_means=same_means,
                                weight_constraints=weight_constraints, structural_pars=structural_pars)) {
    stop("The parameter vector has wrong length!")
  }
  W_constraints <- structural_pars$W
  params <- reform_constrained_pars(p=p, M=M, d=d, params=params, model=model, constraints=constraints,
                                    same_means=same_means, weight_constraints=weight_constraints,
                                    structural_pars=structural_pars)
  structural_pars <- get_unconstrained_structural_pars(structural_pars=structural_pars)
  all_A <- pick_allA(p=p, M=M, d=d, params=params, structural_pars=structural_pars)
  in_paramspace_int(p=p, M=M, d=d, params=params, model=model, all_boldA=form_boldA(p=p, M=M, d=d, all_A=all_A),
                    alphas=pick_alphas(p=p, M=M, d=d, params=params, model=model),
                    all_Omega=pick_Omegas(p=p, M=M, d=d, params=params, structural_pars=structural_pars),
                    W_constraints=W_constraints, stat_tol=stat_tol, posdef_tol=posdef_tol, df_tol=df_tol)
}


#' @title Check that the given parameter vector satisfies the model assumptions
#'
#' @description \code{check_parameters} checks whether the given parameter vector satisfies
#'   the model assumptions. Does NOT consider the identifiability condition!
#'
#' @inheritParams loglikelihood_int
#' @inheritParams is_stationary
#' @inheritParams in_paramspace_int
#' @return Throws an informative error if there is something wrong with the parameter vector.
#' @inherit in_paramspace references
#' @examples
#' \dontrun{
#' # These examples will cause an informative error
#'
#' # GMVAR(1, 1), d=2 model:
#' params11 <- c(1.07, 127.71, 0.99, 0.00, -0.01, 1.00, 4.05,
#'   2.22, 8.87)
#' check_parameters(p=1, M=1, d=2, params=params11)
#'
#' # GMVAR(2, 2), d=2 model:
#' params22 <- c(1.39, -0.77, 1.31, 0.14, 0.09, 1.29, -0.39,
#'  -0.07, -0.11, -0.28, 0.92, -0.03, 4.84, 1.01, 5.93, 1.25,
#'   0.08, -0.04, 1.27, -0.27, -0.07, 0.03, -0.31, 5.85, 10.57,
#'   9.84, 0.74)
#' check_parameters(p=2, M=2, d=2, params=params22)
#'
#' # GMVAR(2, 2), d=2 model with AR-parameters restricted to be
#' # the same for both regimes:
#' C_mat <- rbind(diag(2*2^2), diag(2*2^2))
#' params222c <- c(1.03, 2.36, 1.79, 3.00, 1.25, 0.06,0.04,
#'  1.34, -0.29, -0.08, -0.05, -0.36, 0.93, -0.15, 5.20,
#'  5.88, 3.56, 9.80, 1.37)
#' check_parameters(p=2, M=2, d=2, params=params22c, constraints=C_mat)
#'
#' # Structural GMVAR(2, 2), d=2 model identified with sign-constraints
#' # (no error):
#' params22s <- c(1.03, 2.36, 1.79, 3, 1.25, 0.06, 0.04, 1.34, -0.29,
#'  -0.08, -0.05, -0.36, 1.2, 0.05, 0.05, 1.3, -0.3, -0.1, -0.05, -0.4,
#'   0.89, 0.72, -0.37, 2.16, 7.16, 1.3, 0.37)
#' W_22 <- matrix(c(1, 1, -1, 1), nrow=2, byrow=FALSE)
#' check_parameters(p=2, M=2, d=2, params=params22s,
#'  structural_pars=list(W=W_22))
#' }
#' @export

check_parameters <- function(p, M, d, params, model=c("GMVAR", "StMVAR", "G-StMVAR"), parametrization=c("intercept", "mean"),
                             constraints=NULL, same_means=NULL, weight_constraints=NULL, structural_pars=NULL,
                             stat_tol=1e-3, posdef_tol=1e-8, df_tol=1e-8) {
  model <- match.arg(model)
  check_pMd(p=p, M=M, d=d, model=model)
  if(model != "GMVAR") { # Check degrees of freedom parameters for StMVAR and G-StMVAR models
    if(any(pick_df(M=M, params=params, model=model) <= 2 + df_tol)) {
      stop("The degrees of freedom parameters are not strictly larger than two (with large enough numerical tolerance)")
    } else if(any(pick_df(M=M, params=params, model=model) > 1e+5)) {
      stop("The degrees of freedom parameters should not be larger than 1e+5 (they cause numerical problems and are useless)")
    }
  }
  M_orig <- M
  M <- sum(M)
  parametrization <- match.arg(parametrization)
  check_same_means(parametrization=parametrization, same_means=same_means)
  check_constraints(p=p, M=M_orig, d=d, constraints=constraints, same_means=same_means, weight_constraints=weight_constraints,
                    structural_pars=structural_pars)
  if(length(params) != n_params(p=p, M=M_orig, d=d, model=model, constraints=constraints, same_means=same_means,
                                weight_constraints=weight_constraints, structural_pars=structural_pars)) {
    stop("The parameter vector has wrong dimension!")
  }
  params <- reform_constrained_pars(p=p, M=M_orig, d=d, params=params, model=model, constraints=constraints,
                                    same_means=same_means, weight_constraints=weight_constraints,
                                    structural_pars=structural_pars)
  W_constraints <- structural_pars$W
  structural_pars <- get_unconstrained_structural_pars(structural_pars=structural_pars)
  alphas <- pick_alphas(p=p, M=M_orig, d=d, params=params, model=model)

  if(!is.null(structural_pars)) {
    W_pars <- pick_W(p=p, M=M_orig, d=d, params=params, structural_pars=structural_pars)
    if(any(W_pars[W_constraints < 0] > -1e-8, na.rm=TRUE)) {
      stop("The W parameter does not satisfy the (strict) negative sign constraints (with large enough numerical tolerance)")
    } else if(any(W_pars[W_constraints > 0] < 1e-8, na.rm=TRUE)) {
      stop("The W parameter does not satisfy the (strict) positive sign constraints (with large enough numerical tolerance)")
    }
    if(M > 1) {
      lambdas <- pick_lambdas(p=p, M=M_orig, d=d, params=params, structural_pars=structural_pars)
      if(any(lambdas < 1e-8)) {
        stop("The lambda parameters are not strictly positive (with large enough numerical tolerance)")
      }
    }
  }

  if(M >= 2 & sum(alphas[-M]) >= 1) {
    stop("The mixing weight parameters don't sum to one")
  } else if(any(alphas <= 0)) {
    stop("The mixing weight parameters must be strictly positive")
  } else if(!is_stationary(p=p, M=M, d=d, params=params, structural_pars=structural_pars, tolerance=stat_tol)) {
    stop("The stationarity condition is not satisfied (with large enough numerical tolerance)")
  }
  all_Omega <- pick_Omegas(p=p, M=M, d=d, params=params, structural_pars=structural_pars)
  for(m in 1:M) {
    if(any(eigen(all_Omega[, , m], symmetric=TRUE, only.values=TRUE)$values < posdef_tol)) {
      stop(paste0("Error term covariance matrix of regime ", m, " is not (numerically enough) positive definite"))
    }
  }
}


#' @title Check the constraint matrix has the correct form
#'
#' @description \code{check_constraints} checks that the constraints are correctly set.
#'
#' @inheritParams loglikelihood_int
#' @inheritParams is_stationary
#' @return Checks the constraints and throws an error
#'   if something is wrong.
#' @keywords internal

check_constraints <- function(p, M, d, constraints=NULL, same_means=NULL, weight_constraints=NULL, structural_pars=NULL) {
  M <- sum(M)
  if(!is.null(constraints)) {
    if(!is.matrix(constraints) | !is.numeric(constraints)) {
      stop("The argument constraints should be a numeric matrix (or NULL if no constraints should be employed)")
    } else if(nrow(constraints) != M*p*d^2) {
      stop("The constraint matrix should have M*p*d^2 rows")
    } else if(ncol(constraints) > nrow(constraints)) {
      stop("The constraint matrix has more columns than rows! What are you doing??")
    } else if(qr(constraints)$rank != ncol(constraints)) {
      stop("The constraint matrix should have full column rank")
    }
  }
  if(!is.null(same_means)) {
    if(!is.list(same_means)) {
      stop("The argument same_means should a list (or null if mean parameters are not constrained)")
    } else if(length(same_means) == 0) {
      stop("The argument same_means should not of length zero")
    }
    for(i1 in 1:length(same_means)) {
      if(!is.numeric(same_means[[i1]]) || length(same_means[[i1]]) == 0) {
        stop("The elements of same_means should be numeric vectors with strictly positive length")
      }
    }
    tmp <- sort(unlist(same_means), decreasing=FALSE)
    if(length(tmp) != M || !all(tmp == 1:M)) {
      stop("The argument same_means should contains all regimes in some group exactly once")
    }
  }
  if(!is.null(weight_constraints)) {
    if(!is.numeric(weight_constraints)) {
      stop("The argument weight_constraints should be a numeric vector of length M - 1")
    } else if(length(weight_constraints) != M - 1) {
      stop("The argument weight_constraints should be a numeric vector of length M - 1")
    } else if(any(weight_constraints <= 0)) {
      stop("Each element of weight_constraints should be strictly larger than zero")
    } else if(sum(weight_constraints) >= 1) {
      stop("The elements of weight_constraints should sum to strictly less than one")
    }
  }
  if(!is.null(structural_pars)) {
    if(!is.list(structural_pars)) {
      stop("The argument structural_pars should be a list")
    } else if(is.null(structural_pars$W)) {
      stop("The list 'structural_pars' should contain an element 'W' imposing zero and/or sign constraints on the time-varying B-matrix")
    } else if(!is.matrix(structural_pars$W) || any(dim(structural_pars$W) != d)) {
      stop("The element 'W' in 'structural_pars' should be a (d x d) matrix")
    }
    n_zeros1 <- vapply(1:d, function(i1) sum(structural_pars$W[i1,] == 0, na.rm=TRUE), numeric(1))
    n_zeros2 <- vapply(1:d, function(i1) sum(structural_pars$W[,i1] == 0, na.rm=TRUE), numeric(1))
    if(any(n_zeros1 == d) || any(n_zeros2 == d)) {
      stop("The matrix 'W' is non-singular so you cannot constrain all elements in a row or in a column to be zero")
    }
    if(!is.null(structural_pars$C_lambda)) {
      if(M < 2) stop("There are not lambdas and thus not C_lambda when M < 2")
      C_lamb <- structural_pars$C_lambda
      if(nrow(C_lamb) != d*(M - 1)) {
        stop("The element 'C_lambda' in 'structural_pars' should have d*(M - 1) rows")
      } else if(ncol(C_lamb) > nrow(C_lamb)) {
        stop("The structural parameter constraint matrix 'C_lambda' has more columns than rows! What are you doing??")
      } else if(qr(C_lamb)$rank != ncol(C_lamb)) {
        stop("The structural parameter constraint matrix 'C_lambda' should have full column rank")
      } else if(any(C_lamb < 0)) {
        stop("Entries of the structural parameter constraint matrix 'C_lambda' should be positive or zero")
      }
    }
    if(!is.null(structural_pars$C_lambda) && !is.null(structural_pars$fixed_lambdas)) {
      stop("The constraints C_lambda and fixed_lambdas cannot be used at the same time, remove one of them from structural_pars")
    }
    if(!is.null(structural_pars$fixed_lambdas)) {
      if(!is.numeric(structural_pars$fixed_lambdas)) {
        stop("structural_pars$fixed_lambdas should be a numeric vector of length d*(M - 1)")
      } else if(length(structural_pars$fixed_lambdas) != d*(M - 1)) {
        stop("structural_pars$fixed_lambdas should be a numeric vector of length d*(M - 1)")
      } else if(any(structural_pars$fixed_lambdas <= 0)) {
        stop("Each element of structural_pars$fixed_lambdas should be strictly larger than zero.")
      }
    }
  }
}


#' @title Calculate the number of parameters in a GMVAR, StMVAR, or G-StMVAR model's parameter vector
#'
#' @description \code{n_params} calculates the number of parameters in the model.
#'
#' @inheritParams check_parameters
#' @return Returns the number of parameters in the parameter vector of the specified
#'  GMVAR, StMVAR, or G-StMVAR model.
#' @section Warning:
#'  No argument checks!
#' @inherit in_paramspace references
#' @keywords internal

n_params <- function(p, M, d, model=c("GMVAR", "StMVAR", "G-StMVAR"), constraints=NULL, same_means=NULL,
                     weight_constraints=NULL, structural_pars=NULL) {
  model <- match.arg(model)
  if(model == "GMVAR") {
    n_df <- 0
  } else if(model == "StMVAR") {
    n_df <- M
  } else { # model == "G-StMVAR"
    n_df <- M[2]
    M <- sum(M)
  }
  if(is.null(same_means)) {
    less_pars <- 0 # Number of parameters less compared to models without same intercept constraints
  } else {
    g <- length(same_means) # Number groups with the same mean parameters
    less_pars <- d*(M - g) # Number of parameters less compared to models without same intercept constraints
  }
  if(!is.null(weight_constraints)) {
    less_pars <- less_pars + M - 1
  }
  if(!is.null(structural_pars$fixed_lambdas)) {
    less_pars <- less_pars + d*(M - 1)
  }
  if(is.null(structural_pars)) {
    ret <- ifelse(is.null(constraints), M*(d^2*p + d + d*(d+1)/2 + 1) - 1, M*(d + d*(d + 1)/2 + 1) + ncol(constraints) - 1) - less_pars
  } else {
    q <- ifelse(is.null(constraints), M*p*d^2, ncol(constraints))
    n_Wpars <- length(Wvec(structural_pars$W))
    r <- ifelse(is.null(structural_pars$C_lambda), d*(M - 1), ncol(structural_pars$C_lambda))
    ret <- M*d + q + n_Wpars + r + M - 1 - less_pars
  }
  ret + n_df
}


#' @title Check the data is in the correct form
#'
#' @description \code{check_data} checks the data.
#'
#' @inheritParams loglikelihood_int
#' @return Checks the data and tries to correct it. Throws an error if something is wrong and
#'   returns the corrected data otherwise.
#' @keywords internal

check_data <- function(data, p) {
  if(is.data.frame(data)) {
    data <- as.matrix(data)
  }
  if(!is.matrix(data)) {
    stop("The data must be numeric matrix (possibly a class 'ts' object)!")
  } else {
    if(anyNA(data)) stop("The data contains NA values!")
    if(!is.numeric(data)) stop("The data must be numeric!")
    if(ncol(data) < 2) stop("The data matrix must contain at least two columns! For univariate analysis use the package 'uGMAR'.")
    if(nrow(data) < p + 1) stop("The data must contain at least p+1 observations!")
  }
  data
}


#' @title Check whether all arguments are positive integers
#'
#' @description \code{all_pos_ints} checks whether all the elements in a vector
#'   are positive integers.
#'
#' @param x a vector containing the elements to be tested.
#' @return Returns \code{TRUE} or \code{FALSE} accordingly.
#' @keywords internal

all_pos_ints <- function(x) {
  all(vapply(x, function(x1) x1 %% 1 == 0 && length(x1) == 1 && x1 >= 1, logical(1)))
}


#' @title Check that p, M, and d are correctly set
#'
#' @description \code{check_pMd} checks the arguments p, M, and d.
#'
#' @inheritParams is_stationary
#' @return Throws an error if something is wrong.
#' @keywords internal

check_pMd <- function(p, M, d, model=c("GMVAR", "StMVAR", "G-StMVAR")) {
  model <- match.arg(model)
  if(model == "G-StMVAR") {
    if(length(M) != 2 || !all_pos_ints(M)) {
      stop("For G-StMVAR model, the argument M must a length two vector with positive integer entries")
    }
  } else {
    if(!all_pos_ints(M) || length(M) != 1) {
      stop("For GMVAR and StMVAR models, the argument M must be a positive integer")
    }
  }
  if(!all_pos_ints(p) || length(p) != 1) {
    stop("The argument p must be a positive integer!")
  }
  if(!missing(d)) {
    if(d < 2 | d%%1 != 0) {
      stop("Argument d, the number of columns in the data matrix, has to be a positive integer larger than one!
           For univariate analysis, use the package 'uGMAR'")
    }
  }
}


#' @title Checks whether the given object has class attribute 'gsmvar'
#'
#' @description \code{check_gsmvar} checks that the object has class attribute 'gsmvar'.
#'
#' @param object S3 object to be tested
#' @param object_name what is the name of the object that should of class 'gsmvar'?
#' @return Throws an error if the object doesn't have the class attribute 'gsmvar'.
#' @keywords internal

check_gsmvar <- function(object, object_name) {
  if(missing(object_name)) object_name <- "gsmvar"
  if(!any(class(object) == "gsmvar")) {
    stop(paste("The object", object_name, "has to be of class 'gsmvar', typically created with the function 'GSMVAR' or 'fitGSMVAR'"))
  }
}


#' @title Checks whether the given object contains data
#'
#' @description \code{check_null_data} checks that the gsmvar object has data.
#'
#' @inheritParams quantile_residual_tests
#' @return Throws an error if is.null(gsmvar$data).
#' @keywords internal

check_null_data <- function(gsmvar) {
  if(is.null(gsmvar$data)) {
    stop("The model has to contain data! Data can be added without parameter estimation with the function 'add_data'")
  }
}

#' @title Check whether the parametrization is correct for usage of same means restrictions
#'
#' @description \code{check_same_means} checks whether the parametrization is correct for
#'  usage of same means restrictions
#'
#' @inheritParams loglikelihood_int
#' @return Throws an error if parametrization type is not "mean" and means are constrained
#' @keywords internal

check_same_means <- function(parametrization, same_means) {
  if(parametrization == "intercept" && !is.null(same_means)) {
    stop("The same_means constraints are available for models with mean-parametrization only (parametrization='mean')")
  }
}


#' @title Warn about near-unit-roots in some regimes
#'
#' @description \code{warn_ar_roots} warns if the model contains near-unit-roots in some regimes
#'
#' @inheritParams quantile_residual_tests
#' @details Warns if, for some regime, some moduli of "bold A" eigenvalues are larger than \code{1 - stat_tol} or
#'  some eigenvalue of the error term covariance matrix is smaller than \code{podef_tol}.
#' @return Doesn't return anything.
#' @keywords internal

warn_eigens <- function(gsmvar, stat_tol=0.0015, posdef_tol=0.0002) {
  boldA_eigens <- get_boldA_eigens(gsmvar)
  omega_eigens <- get_omega_eigens(gsmvar)
  near_nonstat <- vapply(1:sum(gsmvar$model$M), function(i1) any(abs(boldA_eigens[,i1]) > 1 - stat_tol), logical(1))
  near_singular <- vapply(1:sum(gsmvar$model$M), function(i1) any(abs(omega_eigens[,i1]) < posdef_tol), logical(1))
  if(any(near_nonstat)) {
    my_string1 <- ifelse(sum(near_nonstat) == 1,
                        paste("Regime", which(near_nonstat),"has near-unit-roots! "),
                        paste("Regimes", paste(which(near_nonstat), collapse=" and ") ,"have near-unit-roots! "))
  } else {
    my_string1 <- NULL
  }
  if(any(near_singular)) {
    my_string2 <- ifelse(sum(near_singular) == 1,
                         paste("Regime", which(near_singular),"has near-singular error term covariance matrix! "),
                         paste("Regimes", paste(which(near_singular), collapse=" and "),
                               "have near-singular error term covariance matrices! "))
  } else {
    my_string2 <- NULL
  }
  if(any(near_nonstat) || any(near_singular)) {
    warning(paste0(my_string1, my_string2, "Consider building a model from the next-largest",
                   " local maximum with the function 'alt_gsmvar' by adjusting its argument 'which_largest'."))
  }
}



#' @title Warn about large degrees of freedom parameter values
#'
#' @description \code{warn_df} warns if the model contains large degrees of freedom parameter values
#'
#' @inheritParams quantile_residual_tests
#' @inheritParams loglikelihood_int
#' @details Warns if, for some regime, the degrees of freedom parameter value is larger than 100.
#' @return Doesn't return anything.
#' @keywords internal

warn_df <- function(gsmvar, p, M, params, model=c("GMVAR", "StMVAR", "G-StMVAR")) {
  model <- match.arg(model)
  if(!missing(gsmvar)) { # If the model is given as a class 'gsmvar' object, extract the information from it
    M <- gsmvar$model$M
    params <- gsmvar$params
    model <- gsmvar$model$model
  }
  if(model == "StMVAR" || model == "G-StMVAR") { # Check whether there is large df parameter value in some regime
    all_df <- pick_df(M=M, params=params, model=model)
    if(any(all_df > 100)) {
      warning(paste0("The model contains overly large degrees of freedom parameters.",
                     "Consider switching to the appropriate G-StMVAR model by setting",
                     "the corresponding regimes to GMVAR type with the function 'stmvar_to_gstmvar'."))
    }
  }
}
