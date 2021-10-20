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
#' @keywords internal

reform_data <- function(data, p) {
  d <- ncol(data)
  n_obs <- nrow(data)
  T_obs <- n_obs - p
  matrix(vapply(1:p, function(i1) as.vector(data[(p - i1 + 1):(T_obs + p - i1 + 1),]), numeric((n_obs - p + 1)*d)), nrow=n_obs - p + 1, byrow=FALSE)
}


#' @title Reform constrained parameter vector into the "standard" form
#'
#' @description \code{reform_constrained_pars} reforms constrained parameter vector
#'   into the form that corresponds to unconstrained parameter vectors.
#'
#' @inheritParams loglikelihood_int
#' @inheritParams is_stationary
#' @param change_na change NA parameter values of constrained models to -9.999?
#' @return Returns "regular model" parameter vector corresponding to the constraints.
#' @section Warning:
#'  No argument checks!
#' @inherit in_paramspace_int references
#' @keywords internal

reform_constrained_pars <- function(p, M, d, params, model=c("GMVAR", "StMVAR", "G-StMVAR"),
                                    constraints=NULL, same_means=NULL, structural_pars=NULL,
                                    change_na=FALSE) {
  model <- match.arg(model)
  if(is.null(constraints) && is.null(structural_pars) && is.null(same_means)) {
    return(params)
  } else if(is.null(constraints) && is.null(same_means) && !is.null(structural_pars) && !any(structural_pars$W == 0, na.rm=TRUE) && is.null(structural_pars$C_lambda)) {
    return(params)
  }
  all_df <- pick_df(M=M, params=params, model=model)
  M2 <- length(all_df)
  M <- sum(M)

  if(is.null(same_means)) {
    less_pars <- 0 # Number of parameters less compared to models without same mean constraints
  } else {
    g <- length(same_means) # Number groups with the same mean parameters
    less_pars <- d*(M - g) # Number of parameters less compared to models without same mean constraints
  }

  # Obtain the AR coefficients from the constraints
  if(is.null(constraints)) { # For structural model with constrained structural parameters but no AR constraints
    q <- M*p*d^2
    psi_expanded <- params[(d*M + 1 - less_pars):(d*M + d^2*p*M - less_pars)] # AR coefficients (without constraints)
    psiNA <- FALSE
  } else {
    q <- ncol(constraints)
    psi <- params[(M*d + 1 - less_pars):(M*d + q - less_pars)]
    if(change_na) {
      if(anyNA(psi)) {
        warning("Replaced some NA values with -9.999")
        psiNA <- TRUE
      } else {
        psiNA <- FALSE
      }
      psi[is.na(psi)] <- -9.999
    }
    psi_expanded <- constraints%*%psi
  }

  # Obtain the mean parameters from the constrained parameter vector
  if(is.null(same_means)) {
    if(is.null(structural_pars)) { # There is computationally more efficient way for structural models.
      all_phi0 <- matrix(params[1:(d*M)], nrow=d, ncol=M) # params[((m - 1)*d + 1):(m*d)]
    }
  } else {
    group_phi0 <- matrix(params[1:(d*g)], nrow=d, ncol=g) # Column for each group
    all_phi0 <- matrix(NA, nrow=d, ncol=M) # Storage for all phi0 (=mean parameters in this case)
    for(i1 in 1:g) {
      all_phi0[,same_means[[i1]]] <- group_phi0[,i1]
    }
  }

  if(is.null(structural_pars)) { # Reduced form model
    pars <- as.vector(vapply(1:M, function(m) c(all_phi0[,m], psi_expanded[((m - 1)*p*d^2 + 1):(m*p*d^2)],
                                                params[(M*d + q + (m - 1)*d*(d + 1)/2 + 1 - less_pars):(M*d + q + m*d*(d + 1)/2 - less_pars)]),
                             numeric(p*d^2 + d + d*(d + 1)/2)))
  } else { # Structural model
    W <- structural_pars$W # Obtain the indices with zero constraints (the zeros don't exist in params)
    n_zeros <- sum(W == 0, na.rm=TRUE)
    new_W <- numeric(d^2)
    W_pars <- params[(d*M + q + 1 - less_pars):(d*M + q + d^2 - n_zeros - less_pars)]
    new_W[W != 0 | is.na(W)] <- W_pars

    if(M > 1) {
      if(!is.null(structural_pars$C_lambda)) {
        r <- ncol(structural_pars$C_lambda)
        gamma <- params[(d*M + q + d^2 - n_zeros + 1 - less_pars):(d*M + q + d^2 - n_zeros + r - less_pars)]
        if(change_na) {
          if(anyNA(gamma) && !psiNA) warning("Replaced some NA values with -9.999")
          gamma[is.na(gamma)] <- -9.999
        }
        lambdas <- structural_pars$C_lambda%*%gamma
      } else {
        lambdas <- params[(d*M + q + d^2 - n_zeros + 1 - less_pars):(d*M + q + d^2 - n_zeros + d*(M - 1) - less_pars)]
      }
    } else {
      lambdas <- numeric(0)
    }
    if(is.null(same_means)) {
      pars <- c(params[1:(M*d)], psi_expanded, vec(new_W), lambdas)
    } else {
      pars <- c(as.vector(all_phi0), psi_expanded, vec(new_W), lambdas)
    }
  }
  if(M == 1) {
    return(c(pars, all_df))
  } else {
    return(c(pars, params[(length(params) - M - M2 + 2):length(params)]))
  }
}


#' @title Reform structural parameter vector into the "standard" form
#'
#' @description \code{reform_structural_pars} reforms (unconstrained) structural
#'   parameter vector into the form that corresponds to reduced form parameter vectors.
#'
#' @inheritParams loglikelihood_int
#' @inheritParams is_stationary
#' @details If the structural parameter vector is a constrained one, use
#'   \code{reform_constrained_pars} first to remove the constraints.
#' @return Returns (unconstrained) "reduced form model" parameter vector.
#' @section Warning:
#'  No argument checks!
#' @inherit in_paramspace_int references
#' @keywords internal

reform_structural_pars <- function(p, M, d, params, model=c("GMVAR", "StMVAR", "G-StMVAR"), structural_pars=NULL) {
  if(is.null(structural_pars)) {
    return(params)
  }
  model <- match.arg(model)
  Omegas <- pick_Omegas(p=p, M=M, d=d, params=params, structural_pars=structural_pars)
  allA <- pick_allA(p=p, M=M, d=d, params=params, structural_pars=structural_pars)
  phi0 <- pick_phi0(p=p, M=M, d=d, params=params, structural_pars=structural_pars)
  all_df <- pick_df(M=M, params=params, model=model)
  M_orig <- M
  M <- sum(M)
  make_upsilon <- function(phi0, Am, Omega) c(phi0, vec(Am), vech(Omega))
  n_upsilon_pars <- d^2*p + d + d*(d + 1)/2
  pars <- numeric(M*n_upsilon_pars) # no alphas
  for(m in 1:M) {
    pars[((m - 1)*n_upsilon_pars + 1):(m*n_upsilon_pars)] <- c(phi0[, m], as.vector(allA[, , , m]), vech(Omegas[, , m]))
  }
  if(M == 1) {
    ret <- pars
  } else {
    alphas <- pick_alphas(p=p, M=M_orig, d=d, params=params, model=model)
    ret <- c(pars, alphas[-M]) # + alphas
  }
  c(ret, all_df)
}


#' @title Get structural parameters that indicate there are no constraints
#'
#' @description \code{get_unconstrained_struct_pars} return structural parameters that indicate there are no constraints
#'  (except possibly sign constraints).
#'
#' @inheritParams loglikelihood_int
#' @return Returns a list with \code{$W} being \eqn{(d x d)} matrix of ones and \code{$C_lambda} being \code{NULL}. If the
#'   supplied argument is \code{NULL}, returns \code{NULL}.
#' @details Intended to be called after calling the function \code{reform_constrained_pars} to avoid remove the constraints
#'   again in any further function calls as this will create bugs. Sign constraints are irrelevant in this context.
#' @section Warning:
#'  No argument checks!
#' @keywords internal

get_unconstrained_structural_pars <- function(structural_pars=NULL) {
  if(is.null(structural_pars)) {
    return(NULL)
  } else {
    d <- nrow(structural_pars$W)
    new_W <- matrix(rep(1, d^2), nrow=d)
    return(list(W=new_W))
  }
}


#' @title Form the \eqn{((dp)x(dp))} "bold A" matrices related to the VAR processes
#'
#' @description \code{form_boldA} creates the "bold A" coefficient matrices related to
#'   VAR processes.
#'
#' @inheritParams pick_allA
#' @param all_A 4D array containing all coefficient matrices \eqn{A_{m,i}}, obtained from \code{pick_allA}.
#' @return Returns 3D array containing the \eqn{((dp)x(dp))} "bold A" matrices related to each component VAR-process.
#'  The matrix \strong{\eqn{A_{m}}} can be obtained by choosing \code{[, , m]}.
#' @section Warning:
#'  No argument checks!
#' @inherit is_stationary references
#' @keywords internal

form_boldA <- function(p, M, d, all_A) {
  M <- sum(M)
  I_all <- diag(nrow=d*(p - 1))
  ZER_all <- matrix(0, nrow=d*(p - 1), ncol=d)
  array(vapply(1:M, function(m) rbind(matrix(all_A[, , 1:p, m], nrow=d, byrow=FALSE), cbind(I_all, ZER_all)), numeric((d*p)^2)), dim=c(d*p, d*p, M))
}



#' @title Sort components in parameter vector according to mixing weights into a decreasing order
#'
#' @description \code{sort_components} sorts mixture components in the parameter vector according
#'   to mixing weights into a decreasing order.
#'
#' @inheritParams is_stationary
#' @details Constrained parameter vectors are not supported (expect for constraints in W but including
#'   constraining some mean parameters to be the same among different regimes)!
#'   For structural models, sorting the regimes in a decreasing order requires re-parametrizing the
#'   decomposition of the covariance matrices if the first regime changes. As a result, the sorted
#'   parameter vector will differ from the given one not only by the ordering of the elements but
#'   also by some of the parameter values.
#' @return Returns sorted parameter vector of the form described for the argument \code{params},
#'   with the mixture components sorted so that \eqn{\alpha_{1}>...>\alpha_{M}} for GMVAR and StMVAR
#'   models, and \eqn{\alpha_{1}>...>\alpha_{M1}} and \eqn{\alpha_{M1+1}>...>\alpha_{M}} for G-StMVAR models.
#' @section Warning:
#'  No argument checks!
#' @inherit in_paramspace_int references
#' @keywords internal

sort_components <- function(p, M, d, params, model=c("GMVAR", "StMVAR", "G-StMVAR"), structural_pars=NULL) {
  model <- match.arg(model)
  M_orig <- M
  M <- sum(M)
  if(M == 1) return(params)
  alphas <- pick_alphas(p=p, M=M_orig, d=d, params=params, model=model)
  all_df <- pick_df(M=M_orig, params=params, model=model)
  if(model != "G-StMVAR") {
    M2 <- ifelse(model == "StMVAR", M, 0)
    ord <- order(alphas, decreasing=TRUE, method="radix")
    if(all(ord == 1:M)) return(params)
    if(model == "StMVAR") all_df <- all_df[ord]
  } else {  # G-StMVAR model: sort the M1 GMVAR type and M2 StMVAR type components separately
    M1 <- M_orig[1]
    M2 <- M_orig[2]
    if(M1 == 1 && M2 == 1) return(params) # Only one component of each type - nothing to sort
    ord_M1 <- order(alphas[1:M1], decreasing=TRUE, method="radix")
    ord_M2 <- order(alphas[(M1 + 1):M], decreasing=TRUE, method="radix")
    ord <- c(ord_M1, M1 + ord_M2) # Overall ordering
    if(all(ord == 1:M)) return(params) # Already in the correct order
    all_df <- all_df[ord_M2] # Df ordered
  }
  alphas <- alphas[ord][-M]
  if(is.null(structural_pars)) { # Reduced form model
    pars <- matrix(params[1:(length(params) - M + 1 - M2)], ncol=M)[,ord]
    return(c(pars, alphas, all_df))
  } else { # Structural model
    n_zeros <- sum(structural_pars$W == 0, na.rm=TRUE)
    phi0 <- pick_phi0(p=p, M=M, d=d, params=params, structural_pars=structural_pars)
    A <- matrix(params[(d*M + 1):(d*M + d^2*p*M)], nrow=d^2*p, byrow=FALSE) # [, m]
    lambdas <- matrix(params[(d*M + d^2*p*M + d^2 - n_zeros + 1):(d*M + d^2*p*M + d^2 - n_zeros + d*(M - 1))],
                      nrow=d, byrow=FALSE)
    W_const <- as.vector(structural_pars$W)
    W_const[is.na(W_const)] <- 1 # Insert arbitrary non-NA and non-zero constraint where was NA
    old_W <- rep(0, times=d^2) # Include non-parametrized zeros here
    old_W[W_const != 0] <- params[(d*M + d^2*p*M + 1):(d*M + d^2*p*M + d^2 - n_zeros)] # Zeros where there are zero constaints
    return(c(phi0[, ord], # sorted phi0/mu parameters
             A[, ord], # sorted AR parameters
             Wvec(redecompose_Omegas(M=M, d=d, W=old_W, lambdas=lambdas, perm=ord)), # Sorted and possibly recomposed the covariance matrices
             alphas, all_df)) # sorted alphas, excluding the M:th one.
  }
}


#' @title Sort the columns of W matrix by sorting the lambda parameters of the second regime to increasing order
#'
#' @description \code{sort_W_and_lambdas} sorts the columns of W matrix by sorting the lambda parameters of
#'  the second regime to increasing order.
#'
#' @inheritParams is_stationary
#' @details Only structural models are supported (but there is no need to provide
#'  structural_pars).
#'  \strong{This function does not sort the constraints of the W matrix but just sorts
#'  the columns of the W matrix and the lambda parameters.} It is mainly used in the genetic
#'  algorithm to assist estimation with better identification when the constraints are not
#'  itself strong for identification of the parameters (but are invariant to different orderings
#'  of the columns of the W matrix).
#'
#'  Before using this function, make sure the parameter vector is sortable: the constraints on
#'  the W matrix is invariant to different orderings of the columns, there are no zero restrictions,
#'  and there are no constraints on the lambda parameters.
#' @return Returns the sorted parameter vector (that implies the same reduced form model).
#' @section Warning:
#'  No argument checks!
#' @references
#'  \itemize{
#'    \item Virolainen S. 2020. Structural Gaussian mixture vector autoregressive model. Unpublished working
#'      paper, available as arXiv:2007.04713.
#'  }
#' @keywords internal

sort_W_and_lambdas <- function(p, M, d, params, model=c("GMVAR", "StMVAR", "G-StMVAR")) {
  model <- match.arg(model)
  if(model == "GMVAR") {
    M2 <- 0
  } else if(model == "StMVAR") {
    M2 <- M
  } else { # model == "G-StMVAR"
    M2 <- M[2]
    M <- sum(M)
  }
  if(M == 1) return(params)
  # The last M2 parameters the degrees of freedom, the next M - 1 parameters are alphas, after that,
  # the last d^2 + d*(M - 1) params are the W and lambda parameters.
  W_and_lambda_inds <- (length(params) - (d^2 + d*(M - 1)) - M2 - M + 2):(length(params) - M2 - M + 1)
  W_and_lambdas <- params[W_and_lambda_inds]
  W <- matrix(W_and_lambdas[1:(d^2)], nrow=d, ncol=d, byrow=FALSE)
  lambdas <- matrix(W_and_lambdas[(d^2 + 1):length(W_and_lambdas)], nrow=d, ncol=M - 1, byrow=FALSE)
  new_order <- order(lambdas[,1], decreasing=FALSE) # Increasing order for second regime lambdas
  W <- W[,new_order] # Columns to the new order
  lambdas <- lambdas[new_order,] # Rows to the new order
  params[W_and_lambda_inds] <- c(W, lambdas)
  params
}


#' @title Change parametrization of a parameter vector
#'
#' @description \code{change_parametrization} changes the parametrization of the given parameter
#'   vector to \code{change_to}.
#'
#' @inheritParams loglikelihood_int
#' @inheritParams form_boldA
#' @param change_to either "intercept" or "mean" specifying to which parametrization it should be switched to.
#'   If set to \code{"intercept"}, it's assumed that \code{params} is mean-parametrized, and if set to \code{"mean"}
#'   it's assumed that \code{params} is intercept-parametrized.
#' @return Returns parameter vector described in \code{params}, but with parametrization changed from intercept to mean
#'   (when \code{change_to==mean}) or from mean to intercept (when \code{change_to==intercept}).
#' @details Parametrization cannot be changed for models with same_means constraints!
#' @section Warning:
#'  No argument checks!
#' @inherit is_stationary references
#' @keywords internal

change_parametrization <- function(p, M, d, params, model=c("GMVAR", "StMVAR", "G-StMVAR"), constraints=NULL, same_means=NULL,
                                   structural_pars=NULL, change_to=c("intercept", "mean")) {
  stopifnot(is.null(same_means))
  model <- match.arg(model)
  change_to <- match.arg(change_to)
  re_params <- params
  params <- reform_constrained_pars(p=p, M=M, d=d, params=params, model=model, constraints=constraints, same_means=same_means,
                                    structural_pars=structural_pars) # Parameters in "regular" form
  Id <- diag(nrow=d)
  all_A <- pick_allA(p=p, M=M, d=d, params=params, structural_pars=structural_pars)
  all_phi0_or_mu <- pick_phi0(p=p, M=M, d=d, params=params, structural_pars=structural_pars)
  M <- sum(M)

  if(is.null(constraints) && is.null(structural_pars)) {
    qm1 <- (1:M - 1)*(d + p*d^2 + d*(d + 1)/2)
  } else {
    qm1 <- (1:M - 1)*d
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
#'   of the given parameter vector
#'
#' @description \code{change_regime} changes the regime parameters (excluding mixing weights parameter)
#'   of the pointed regime to the new given parameters.
#'
#' @inheritParams pick_regime
#' @param regime_pars
#'   \describe{
#'     \item{For reduced form models:}{a size \eqn{((pd^2+d+d(d+1)/2)x1)} vector
#'       (\strong{\eqn{\upsilon_{m}}}\eqn{,\nu_m})\eqn{ = (\phi_{m,0},}\strong{\eqn{\phi_{m}}}\eqn{,\sigma_{m},\nu_m)}.}
#'     \item{For structural models:}{a length \eqn{pd^2 + d} vector \eqn{(\phi_{m,0},}\strong{\eqn{\phi_{m}}}\eqn{,\nu_m)}.}
#'     In the case of a GMVAR type regime, \eqn{\nu_m} is omitted.
#'   }
#' @return Returns parameter vector with \code{m}:th regime changed to \code{regime_pars}.
#' @details Does not currently support models with AR, mean, or lambda parameter constraints.
#' @section Warning:
#'  No argument checks!
#' @inherit in_paramspace_int references
#' @keywords internal

change_regime <- function(p, M, d, params, model=c("GMVAR", "StMVAR", "G-StMVAR"), m, regime_pars, structural_pars=NULL) {
  model <- match.arg(model)
  all_df <- pick_df(M=M, params=params, model=model)
  if(model == "StMVAR" || (model == "G-StMVAR" && m > M[1])) {
    M <- sum(M)
    df <- regime_pars[length(regime_pars)] # Extract df
    regime_pars <- regime_pars[-length(regime_pars)] # Remove df
    params[length(params) - M + m] <- df # Change df in the parameter vector (the index is the same for G-StMVAR model)
  } else if(model == "G-StMVAR") { # GMVAR type regime is changes
    M <- sum(M)
  }
  if(is.null(structural_pars)) {
    qm1 <- (m - 1)*(d + p*d^2 + d*(d + 1)/2)
    params[(qm1 + 1):(qm1 + d + p*d^2 + d*(d + 1)/2)] <- regime_pars
  } else {
    n_zeros <- sum(structural_pars$W == 0, na.rm=TRUE)
    params[(d*(m - 1) + 1):(d*m)] <- regime_pars[1:d] # phi0
    params[(d*M + d^2*p*(m - 1) + 1):(d*M + d^2*p*(m - 1) + d^2*p)] <- regime_pars[(d + 1):(d + d^2*p)] # AR coefs
  }
  params
}


#' @title Calculate "distance" between two (scaled) regimes
#'  \strong{\eqn{\upsilon_{m}}}\eqn{ = (\phi_{m,0},}\strong{\eqn{\phi_{m}}}\eqn{,\sigma_{m})}
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
#' @keywords internal

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
  c(sqrt(crossprod(regime_pars1/scales1 - regime_pars2/scales2)))
}


#' @title Sort mixing weight parameters in a decreasing order and standardize them
#'  to sum to one. For G-StMVAR models, the mixing weight parameters are sorted
#'  to a decreasing order for GMVAR and StMVAR type regimes separately.
#'
#' @description \code{sort_and_standardize_alphas} sorts mixing weight parameters in a decreasing
#'  order and standardizes them to sum to one. For G-StMVAR models, the mixing weight parameters are sorted
#'  to a decreasing order for GMVAR and StMVAR type regimes separately. Does not sort if AR constraints,
#'  lambda constraints, or same means are employed.
#'
#' @inheritParams loglikelihood_int
#' @param alphas mixing weights parameters alphas, \strong{INCLUDING} the one for the M:th regime (that is
#'  not parametrized in the model). Don't need to be standardized to sum to one.
#' @return Returns the given alphas in a (M x 1) vector sorted in decreasing order and their sum standardized to one.
#'  If AR constraints, lambda constraints, or same means are employed, does not sort but standardizes the alphas
#'  to sum to one.
#' @section Warning:
#'  No argument checks!
#' @keywords internal

sort_and_standardize_alphas <- function(alphas, M, model=c("GMVAR", "StMVAR", "G-StMVAR"),
                                        constraints=NULL, same_means=NULL, structural_pars=NULL) {
  model <- match.arg(model)
  if(is.null(constraints) && is.null(structural_pars$C_lambda) && is.null(same_means)) {
    if(model != "G-StMVAR") {
      alphas <- alphas[order(alphas, decreasing=TRUE, method="radix")]
    } else { # model == "G-StMVAR"
      alphas <- c(alphas[order(alphas[1:M[1]], decreasing=TRUE, method="radix")], # Sort GMVAR type alphas
                  alphas[M[1] + order(alphas[(M[1] + 1):sum(M)], decreasing=TRUE, method="radix")]) # Sort StMVAR type alphas
    }
  }
  alphas/sum(alphas)
}



#' @title Transform a StMVAR (or G-StMVAR) model parameter vector to the corresponding G-StMVAR model parameter vector
#'  with the large df parameters removed.
#'
#' @description \code{stmvarpars_to_gstmvar} transforms a StMVAR (or G-StMVAR) model parameter vector to the corresponding
#'   G-StMVAR model parameter vector with the large df parameters removed.
#'
#' @inheritParams loglikelihood_int
#' @param maxdf regimes with degrees of freedom parameter value larger than this will be turned into
#'  GMVAR type.
#' @return Returns a list with four elements: \code{$params} contains the corresponding G-StMVAR model
#'  parameter vector, \code{$reg_order} contains the permutation that was applied to the regimes
#'  (GMVAR type regimes first, and decreasing ordering by mixing weight parameters),
#'  \code{$M} a vector of length two containing the number of GMVAR type regimes in the first element
#'  and the number of StMVAR type regimes in the second, \code{$same_means} constains the \code{same_mean}
#'  constraints of the new model.
#' @section Warning:
#'  No argument checks!
#' @keywords internal

stmvarpars_to_gstmvar <- function(p, M, d, params, model=c("GMVAR", "StMVAR", "G-StMVAR"),
                                  constraints=NULL, same_means=NULL, structural_pars=NULL,
                                  maxdf=100) {
  model <- match.arg(model)
  stopifnot(model %in% c("StMVAR", "G-StMVAR"))
  constraints_orig <- constraints
  if(model == "StMVAR") {
    M1 <- 0
    M2 <- M
  } else { # model == "G-StMVAR"
    M1 <- M[1]
    M2 <- M[2]
  }
  M_orig <- M
  M <- sum(M)
  all_df <- pick_df(M=M_orig, params=params, model=model)
  if(!any(all_df > maxdf)) { # Nothing is changed in the model
    warning("No degrees of freedom parameter is larger than 'maxdf'. The original model is returned.")
    if(model == "StMVAR") {
      ret_M <- c(0, M)
    } else { # model == "G-StMVAR"
      ret_M <- M_orig
    }
    return(list(params=params,
                reg_order=1:M,
                M=ret_M,
                same_means=same_means))
  }

  # Which regimes will be GMVAR type in the new model
  regs_to_change <- which(all_df > maxdf) + M1 # Which regimes to change to GMVAR type
  if(length(regs_to_change) == M2) message("All regimes are changed to GMVAR type. The result is therefore a GMVAR model and not a G-StMVAR model.")
  if(model == "StMVAR") {
    gmvar_regs <- regs_to_change
  } else { # model == "G-StMVAR"
    gmvar_regs <- c(1:M1, regs_to_change)
  }
  new_M1 <- length(gmvar_regs) # M1 in the new model
  new_M2 <- M - new_M1
  ret_M <- c(new_M1, new_M2)

  # Obtain the new parameter vector, but excluding the mixing weight and degrees of freedom parameters
  params_std <- reform_constrained_pars(p=p, M=M_orig, d=d, params=params, model=model, constraints=constraints, same_means=same_means,
                                        structural_pars=structural_pars)
  alphas <- pick_alphas(p=p, M=M_orig, d=d, params=params_std, model=model)
  if(M > 1) {
    reg_order <- c(gmvar_regs[order(alphas[gmvar_regs], decreasing=TRUE)], # GMVAR type regimes
                   (1:M)[-gmvar_regs][order(alphas[-gmvar_regs], decreasing=TRUE)]) # StMVAR type regimes
    new_alphas <- alphas[reg_order][-M] # New mixing weight parameters
  } else { # Only one regime
    reg_order <- 1
    new_alphas <- numeric(0)
    new_df <- numeric(0)
  }
  if(new_M2 > 0) {
    new_df <-  all_df[-c(regs_to_change - M1)][order(reg_order[(new_M1 + 1):M])] # The new degrees of freedom parameter vector
  } else {
    new_df <- numeric(0)
  }

  # Construct the new parameter vectors, exxluding mixing weights and degrees of feedom parameters
  if(all(reg_order == 1:M)) {
    # The regime order is the same in the new model, so only the degrees of freedom parameters change (some are removed)
    new_params <- params[1:(length(params) - length(all_df) - (M - 1))] # Remove the mixing weights and degrees of freedom parameters
  } else {
    # Here the ordering of the regimes changes

    # If constraints are employed on the autoregressive parameters, the regimes cannot be sorted due to the
    # specific form of the constraints. In this case, all constraints are removed and we then do the switch
    # with the unconstrained model.
    strpars_unc <- get_unconstrained_structural_pars(structural_pars=structural_pars)
    if(!is.null(constraints)) {
      # We remove also same_mean constraints, if any
      params <- params_std
      constraints <- NULL
      same_means <- NULL
      if(is.null(structural_pars)) {
        message("All constraints were removed from the model when switching to the G-StMVAR model")
      } else {
        message("All constraints were removed from the model, except the ones imposed on the W-matrix, when switching to the G-StMVAR model")
      }
    } else if(!is.null(same_means)) {
      # If only mean parameters are constrained, the regimes can be sorted as long as the mean-parameter
      # constraints reformed accordingly.
      new_same_means <- vector(mode="list", length=length(same_means))
      for(i1 in 1:length(same_means)) { # Go though the same means groups
        new_same_means[[i1]] <- which(reg_order %in% same_means[[i1]]) # Create the new groups according to the new ordering
      }
      cat("The 'same_means' constraints were reconstructed when reordering the regimes. The new same_means is:\n")
      same_means <- new_same_means
      print(same_means)
    }
    if(!is.null(structural_pars$C_lambda) && is.null(constraints_orig)) {
      # In this case, removal of the (lambda) constraints has not yet been messaged
      # We will remove lambda constraints from the structural models where regimes are re-ordered.
      message("Lambda constraints were removed from the model when switching to the G-StMVAR model")
    }

    # Obtain the new parameters, form the new regimes, sort them, and form the new parameter vector
    all_phi0 <- pick_phi0(p=p, M=M_orig, d=d, params=params_std, structural_pars=strpars_unc) # [, m]
    all_A <- pick_allA(p=p, M=M_orig, d=d, params=params_std, structural_pars=strpars_unc) # [, , p, m]

    new_all_phi0 <- all_phi0[, reg_order] # [, m] with the new order
    new_allA <- matrix(all_A, nrow=p*d^2, ncol=M)[, reg_order] # [, m] with the new order

    # Go through the special cases, case by case
    if(is.null(structural_pars)) { ## Reduced form model
      all_Omega <- pick_Omegas(p=p, M=M_orig, d=d, params=params_std, structural_pars=structural_pars) # [, , m]
      new_all_Omega <- vapply(1:M, function(m) vech(all_Omega[, , m]), numeric(d*(d + 1)/2))[, reg_order] # [, m] with the new order
      if(!is.null(constraints) || (is.null(constraints) && is.null(same_means))) {
        # AR (and possibly mean constraints) employed and regime-order changes OR no constraints employed at all.
        # -> In this case, all constraints are removed and we work with the "standard form" parameter vector.
        new_params <- as.vector(rbind(new_all_phi0, new_allA, new_all_Omega)) # The new parameter vector, excluding alphas and df
      } else {
        # Here AR constraints are not employed but mean constraints are.
        # -> We proceed with mean constrained parameter vector.

        # The mean parameters should in the order of the groups: i:th mean is for group i.
        # Because the ordering of the groups did not change, although the regimes were re-ordered,
        # we do not re-order the mean-parameters.
        new_mean_params <- params[1:(d*length(same_means))] # The old mean params are also the new mean params
        new_params <- c(new_mean_params, new_allA, new_all_Omega)
      }
    } else { ## Structural model
      W_pars <- pick_W(p=p, M=M_orig, d=d, params=params_std, structural_pars=structural_pars)
      lambdas <- pick_lambdas(p=p, M=M_orig, d=d, params=params_std, structural_pars=structural_pars)
      new_W_and_lambdas <- redecompose_Omegas(M=M_orig, d=d, W=W_pars, lambdas=lambdas, perm=reg_order)
      new_W <- new_W_and_lambdas[1:(d^2)]
      new_W <- new_W[!structural_pars$W == 0] # Remove zeros when the corresponding elements of W are constrained to zero and hence not parametrized (note: the new W has zeros in the same entries as the old W)
      new_lambdas <- new_W_and_lambdas[(d^2 + 1):length(new_W_and_lambdas)]
      if(!is.null(same_means)) {
        new_mean_params <- params[1:(d*length(same_means))] # The old mean params are also the new mean params (see the reduced form case)
      } else {
        new_mean_params <- new_all_phi0
      }
      # In the structural models, constrained and unconstrainwd models have the same form for the parameter vector.
      # Note that if AR or lambda constraints are employed, they are removed.
      new_params <- c(new_mean_params, new_allA, new_W, new_lambdas)
    }
  }

  list(params=c(new_params, new_alphas, new_df), # Add new mixing weight and degrees of freedom parameter vectors
       reg_order=reg_order,
       M=ret_M,
       same_means=same_means)
}
