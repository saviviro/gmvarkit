#' @title Check the stationary condition of a given GMVAR model
#'
#' @description \code{is_stationary} checks the stationarity condition of a GMVAR model.
#'
#' @inheritParams loglikelihood_int
#' @param d the number of time series in the system.
#' @param params  a real valued vector specifying the parameter values.
#'   \describe{
#'     \item{\strong{For reduced form model:}}{
#'       Should be size
#'       \eqn{((M(pd^2+d+d(d+1)/2+1)-1)x1)} and have form \strong{\eqn{\theta}}\eqn{ = }(\strong{\eqn{\upsilon}}\eqn{_{1}},
#'       ...,\strong{\eqn{\upsilon}}\eqn{_{M}}, \eqn{\alpha_{1},...,\alpha_{M-1}}), where:
#'       \itemize{
#'         \item \strong{\eqn{\upsilon}}\eqn{_{m}} \eqn{ = (\phi_{m,0},}\strong{\eqn{\phi}}\eqn{_{m}}\eqn{,\sigma_{m})}
#'         \item \strong{\eqn{\phi}}\eqn{_{m}}\eqn{ = (vec(A_{m,1}),...,vec(A_{m,p})}
#'         \item and \eqn{\sigma_{m} = vech(\Omega_{m})}, m=1,...,M.
#'       }
#'     }
#'     \item{\strong{For structural GMVAR model:}}{
#'       Should have the form
#'       \strong{\eqn{\theta}}\eqn{ = (\phi_{1,0},...,\phi_{M,0},}\strong{\eqn{\phi}}\eqn{_{1},...,}\strong{\eqn{\phi}}\eqn{_{M},
#'       vec(W),}\strong{\eqn{\lambda}}\eqn{_{2},...,}\strong{\eqn{\lambda}}\eqn{_{M},\alpha_{1},...,\alpha_{M-1})}, where
#'       \itemize{
#'         \item\strong{\eqn{\lambda}}\eqn{_{m}=(\lambda_{m1},...,\lambda_{md})} contains the eigenvalues of the \eqn{m}th mixture component.
#'       }
#'     }
#'   }
#'   Above, \eqn{\phi_{m,0}} is the intercept parameter, \eqn{A_{m,i}} denotes the \eqn{i}:th coefficient matrix of
#'   the \eqn{m}:th mixture component, \eqn{\Omega_{m}} denotes the error term covariance matrix of the \eqn{m}:th
#'   mixture component, and \eqn{\alpha_{m}} is the mixing weight parameter. The \eqn{W} and \eqn{\lambda_{mi}} are
#'   structural parameters replacing the error term covariance matrices (see Virolainen, 2020). If \eqn{M=1}, \eqn{\alpha_{m}}
#'   and \eqn{\lambda_{mi}} are dropped.
#'
#'   If \code{parametrization=="mean"}, just replace each \eqn{\phi_{m,0}} with the regimewise mean \eqn{\mu_{m}}.
#'   \eqn{vec()} is vectorization operator that stacks columns of a given matrix into a vector. \eqn{vech()} stacks columns
#'   of a given matrix from the principal diagonal downwards (including elements on the diagonal) into a vector.
#'   The notation is in line with the cited article by KMS (2016) introducing the GMVAR model.
#' @param all_boldA 3D array containing the \eqn{((dp)x(dp))} "bold A" matrices related to each mixture component VAR-process,
#'   obtained from \code{form_boldA}. Will be computed if not given.
#' @param tolerance Returns \code{FALSE} if modulus of any eigenvalue is larger or equal to \code{1-tolerance}.
#' @details If the model is constrained, remove the constraints first with the function \code{reform_constrained_pars}.
#' @return Returns \code{TRUE} if the model is stationary and \code{FALSE} if not. Based on the argument \code{tolerance},
#'   \code{is_stationary} may return \code{FALSE} when the parameter vector is in the stationarity region, but
#'   very close to the boundary (this is used to ensure numerical stability in estimation of the model parameters).
#' @section Warning:
#'  No argument checks!
#' @inherit loglikelihood_int references

is_stationary <- function(p, M, d, params, all_boldA=NULL, structural_pars=NULL, tolerance=1e-3) {
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
#' @inheritParams is_stationary
#' @param alphas (Mx1) vector containing all mixing weight parameters, obtained from \code{pick_alphas}.
#' @param all_Omega 3D array containing all covariance matrices \eqn{\Omega_{m}}, obtained from \code{pick_Omegas}.
#' @param W_constraints set \code{NULL} for reduced form models. For structural models, this should be the
#'   constraint matrix \eqn{W} from the list of structural parameters.
#' @details The parameter vector in the argument \code{params} should be unconstrained and it is used for
#'   structural models only.
#' @return Returns \code{TRUE} if the given parameter values are in the parameter space and \code{FALSE} otherwise.
#'   This function does NOT consider the identifiability condition!
#' @references
#'  \itemize{
#'    \item Kalliovirta L., Meitz M. and Saikkonen P. 2016. Gaussian mixture vector autoregression.
#'          \emph{Journal of Econometrics}, \strong{192}, 485-498.
#'    \item Virolainen S. 2020. Structural Gaussian mixture vector autoregressive model. Unpublished working
#'      paper, available as arXiv:2007.04713.
#'  }

in_paramspace_int <- function(p, M, d, params, all_boldA, alphas, all_Omega, W_constraints=NULL) {

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
  } else if(!is_stationary(p=p, M=M, d=d, all_boldA=all_boldA)) {
    return(FALSE)
  }
  for(m in 1:M) {
    if(any(eigen(all_Omega[, , m], symmetric=TRUE, only.values=TRUE)$values < 1e-8)) {
      return(FALSE)
    }
  }
  TRUE
}


#' @title Determine whether the parameter vector lies in the parameter space
#'
#' @description \code{in_paramspace} checks whether the given parameter vector lies in
#'   the parameter space. Does NOT consider the identifiability condition!
#'
#' @inheritParams loglikelihood_int
#' @inheritParams is_stationary
#' @return Returns \code{TRUE} if the given parameter vector lies in the parameter space
#'  and \code{FALSE} otherwise.
#' @inherit in_paramspace_int references
#' @examples
#' # GMVAR(1,1), d=2 model:
#' params112 <- c(1.07, 127.71, 0.99, 0.00, -0.01, 0.99, 4.05,
#'   2.22, 8.87)
#' in_paramspace(p=1, M=1, d=2, params=params112)
#'
#' # GMVAR(2,2), d=2 model:
#' params222 <- c(1.39, -0.77, 1.31, 0.14, 0.09, 1.29, -0.39,
#'  -0.07, -0.11, -0.28, 0.92, -0.03, 4.84, 1.01, 5.93, 1.25,
#'   0.08, -0.04, 1.27, -0.27, -0.07, 0.03, -0.31, 5.85, 3.57,
#'   9.84, 0.74)
#' in_paramspace(p=2, M=2, d=2, params=params222)
#'
#' # GMVAR(2,2), d=2 model with AR-parameters restricted to be
#' # the same for both regimes:
#' C_mat <- rbind(diag(2*2^2), diag(2*2^2))
#' params222c <- c(1.03, 2.36, 1.79, 3.00, 1.25, 0.06,0.04,
#'  1.34, -0.29, -0.08, -0.05, -0.36, 0.93, -0.15, 5.20,
#'  5.88, 3.56, 9.80, 0.37)
#' in_paramspace(p=2, M=2, d=2, params=params222c, constraints=C_mat)
#'
#' # Structural GMVAR(2, 2), d=2 model identified with sign-constraints:
#' params222s <- c(1.03, 2.36, 1.79, 3, 1.25, 0.06, 0.04, 1.34, -0.29,
#'  -0.08, -0.05, -0.36, 1.2, 0.05, 0.05, 1.3, -0.3, -0.1, -0.05, -0.4,
#'   0.89, 0.72, -0.37, 2.16, 7.16, 1.3, 0.37)
#' W_222 <- matrix(c(1, NA, -1, 1), nrow=2, byrow=FALSE)
#' in_paramspace(p=2, M=2, d=2, params=params222s,
#'   structural_pars=list(W=W_222))
#' @export

in_paramspace <- function(p, M, d, params, constraints=NULL, structural_pars=NULL) {
  check_pMd(p=p, M=M, d=d)
  check_constraints(p=p, M=M, d=d, constraints=constraints, structural_pars=structural_pars)
  if(length(params) != n_params(p=p, M=M, d=d, constraints=constraints, structural_pars=structural_pars)) {
    stop("The parameter vector has wrong length!")
  }
  W_constraints <- structural_pars$W
  params <- reform_constrained_pars(p=p, M=M, d=d, params=params, constraints=constraints, structural_pars=structural_pars)
  structural_pars <- get_unconstrained_structural_pars(structural_pars=structural_pars)
  all_A <- pick_allA(p=p, M=M, d=d, params=params, structural_pars=structural_pars)
  in_paramspace_int(p=p, M=M, d=d, params=params, all_boldA=form_boldA(p=p, M=M, d=d, all_A=all_A),
                    alphas=pick_alphas(p=p, M=M, d=d, params=params),
                    all_Omega=pick_Omegas(p=p, M=M, d=d, params=params, structural_pars=structural_pars),
                    W_constraints=W_constraints)
}


#' @title Check that the given parameter vector satisfies the model assumptions
#'
#' @description \code{check_parameters} checks whether the given parameter vector satisfies
#'   the model assumptions. Does NOT consider the identifiability condition!
#'
#' @inheritParams loglikelihood_int
#' @inheritParams is_stationary
#' @return Throws an informative error if there is something wrong with the parameter vector.
#' @inherit in_paramspace references
#' @examples
#' \dontrun{
#' # These examples will cause an informative error
#'
#' # GMVAR(1, 1), d=2 model:
#' params112 <- c(1.07, 127.71, 0.99, 0.00, -0.01, 1.00, 4.05,
#'   2.22, 8.87)
#' check_parameters(p=1, M=1, d=2, params=params11)
#'
#' # GMVAR(2, 2), d=2 model:
#' params222 <- c(1.39, -0.77, 1.31, 0.14, 0.09, 1.29, -0.39,
#'  -0.07, -0.11, -0.28, 0.92, -0.03, 4.84, 1.01, 5.93, 1.25,
#'   0.08, -0.04, 1.27, -0.27, -0.07, 0.03, -0.31, 5.85, 10.57,
#'   9.84, 0.74)
#' check_parameters(p=2, M=2, d=2, params=params222)
#'
#' # GMVAR(2, 2), d=2 model with AR-parameters restricted to be
#' # the same for both regimes:
#' C_mat <- rbind(diag(2*2^2), diag(2*2^2))
#' params222c <- c(1.03, 2.36, 1.79, 3.00, 1.25, 0.06,0.04,
#'  1.34, -0.29, -0.08, -0.05, -0.36, 0.93, -0.15, 5.20,
#'  5.88, 3.56, 9.80, 1.37)
#' check_parameters(p=2, M=2, d=2, params=params222c, constraints=C_mat)
#'
#' # Structural GMVAR(2, 2), d=2 model identified with sign-constraints:
#' params222s <- c(1.03, 2.36, 1.79, 3, 1.25, 0.06, 0.04, 1.34, -0.29,
#'  -0.08, -0.05, -0.36, 1.2, 0.05, 0.05, 1.3, -0.3, -0.1, -0.05, -0.4,
#'   0.89, 0.72, -0.37, 2.16, 7.16, 1.3, 0.37)
#' W_222 <- matrix(c(1, NA, -1, 1), nrow=2, byrow=FALSE)
#' check_parameters(p=2, M=2, d=2, params=params222s,
#'  structural_pars=list(W=W_222))
#' }
#' @export

check_parameters <- function(p, M, d, params, constraints=NULL, structural_pars=NULL) {

  check_pMd(p=p, M=M, d=d)
  check_constraints(p=p, M=M, d=d, constraints=constraints, structural_pars=structural_pars)
  if(length(params) != n_params(p=p, M=M, d=d, constraints=constraints, structural_pars=structural_pars)) {
    stop("The parameter vector has wrong dimension!")
  }
  params <- reform_constrained_pars(p=p, M=M, d=d, params=params, constraints=constraints, structural_pars=structural_pars)
  W_constraints <- structural_pars$W
  structural_pars <- get_unconstrained_structural_pars(structural_pars=structural_pars)
  alphas <- pick_alphas(p=p, M=M, d=d, params=params)

  if(!is.null(structural_pars)) {
    W_pars <- pick_W(p=p, M=M, d=d, params=params, structural_pars=structural_pars)
    if(any(W_pars[W_constraints < 0] > -1e-8, na.rm=TRUE)) {
      stop("The W parameter does not satisfy the (strict) negative sign constraints (with large enough numerical tolerance)")
    } else if(any(W_pars[W_constraints > 0] < 1e-8, na.rm=TRUE)) {
      stop("The W parameter does not satisfy the (strict) positive sign constraints (with large enough numerical tolerance)")
    }
    if(M > 1) {
      lambdas <- pick_lambdas(p=p, M=M, d=d, params=params, structural_pars=structural_pars)
      if(any(lambdas < 1e-8)) {
        stop("The lambda parameters are not strictly positive (with large enough numerical tolerance)")
      }
    }
  }

  if(M >= 2 & sum(alphas[-M]) >= 1) {
    stop("The mixing weight parameters don't sum to one")
  } else if(any(alphas <= 0)) {
    stop("The mixing weight parameters must be strictly positive")
  } else if(!is_stationary(p=p, M=M, d=d, params=params, structural_pars=structural_pars)) {
    stop("The stationarity condition is not satisfied (with large enough numerical tolerance)")
  }
  all_Omega <- pick_Omegas(p=p, M=M, d=d, params=params, structural_pars=structural_pars)
  for(m in 1:M) {
    if(any(eigen(all_Omega[, , m], symmetric=TRUE, only.values=TRUE)$values < 1e-8)) {
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
#' @return Checks the constraint matrix \strong{C} and throws an error
#'   if something is wrong.
#' @details If \code{is.null(constraints)}, then this function doesn't do anything.

check_constraints <- function(p, M, d, constraints=NULL, structural_pars=NULL) {
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
  }
}


#' @title Calculate the number of parameters in GMVAR model parameter vector
#'
#' @description \code{n_params} calculates the number of parameters in the model.
#'
#' @inheritParams check_parameters
#' @return Returns the number of parameters in parameter vector of the specified GMVAR model.
#' @section Warning:
#'  No argument checks!
#' @inherit in_paramspace references

n_params <- function(p, M, d, constraints=NULL, structural_pars=NULL) {
  if(is.null(structural_pars)) {
    return(ifelse(is.null(constraints), M*(d^2*p + d + d*(d+1)/2 + 1) - 1, M*(d + d*(d + 1)/2 + 1) + ncol(constraints) - 1))
  } else {
    q <- ifelse(is.null(constraints), M*p*d^2, ncol(constraints))
    n_Wpars <- length(Wvec(structural_pars$W))
    r <- ifelse(is.null(structural_pars$C_lambda), d*(M - 1), ncol(structural_pars$C_lambda))
    return(M*d + q + n_Wpars + r + M - 1)
  }
}


#' @title Check the data is in the correct form
#'
#' @description \code{check_data} checks the data.
#'
#' @inheritParams loglikelihood_int
#' @return Checks the data and tries to correct it. Throws an error if something is wrong and
#'   returns the corrected data otherwise.

check_data <- function(data, p) {
  if(is.data.frame(data)) {
    data <- as.matrix(data)
  }
  if(!is.matrix(data)) {
    stop("The data must be numeric matrix!")
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

all_pos_ints <- function(x) {
  all(vapply(x, function(x1) x1 %% 1 == 0 && length(x1) == 1 && x1 >= 1, logical(1)))
}


#' @title Check that p, M, and d are correctly set
#'
#' @description \code{check_pMd} checks the arguments p, M, and d.
#'
#' @inheritParams is_stationary
#' @return Throws an error if something is wrong.

check_pMd <- function(p, M, d) {
  if(!all_pos_ints(c(p, M))) {
    stop("Arguments p and M have to be positive integers!")
  }
  if(!missing(d)) {
    if(d < 2 | d%%1 != 0) {
      stop("Argument d, number of columns in the data matrix, has to be positive integer larger than one!
           For univariate analysis, use the package 'uGMAR'")
    }
  }
}


#' @title Checks whether the given object has class attribute 'gmvar'
#'
#' @description \code{check_gmvar} checks that the object has class attribute 'gmvar'.
#'
#' @param object S3 object to be tested
#' @param object_name what is the name of the object that should of class 'gmvar'?
#' @return Throws an error if the object doesn't have the class attribute 'gmvar'.

check_gmvar <- function(object, object_name) {
  if(missing(object_name)) object_name <- "gmvar"
  if(!any(class(object) == "gmvar")) {
    stop(paste("The object", object_name, "has to be of class 'gmvar', typically created with the function 'GMVAR' or 'fitGMVAR'"))
  }
}


#' @title Checks whether the given object contains data
#'
#' @description \code{check_null_data} checks that the gmvar object has data.
#'
#' @inheritParams simulateGMVAR
#' @return Throws an error if is.null(gmvar$data).

check_null_data <- function(gmvar) {
  if(is.null(gmvar$data)) {
    stop("The model has to contain data! Data can be added without parameter estimation with the function 'add_data'")
  }
}

