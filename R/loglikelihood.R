#' @import stats
#'
#' @title Compute log-likelihood of a Gaussian mixture vector autoregressive model
#'
#' @description \code{loglikelihood_int} computes log-likelihood of a GMVAR model.
#'
#' @inheritParams in_paramspace_int
#' @param data a matrix or class \code{'ts'} object with \code{d>1} columns. Each column is taken to represent
#'  a single time series. \code{NA} values are not supported.
#' @param p a positive integer specifying the autoregressive order of the model.
#' @param M a positive integer specifying the number of mixture components.
#' @param params a real valued vector specifying the parameter values.
#'   \describe{
#'     \item{\strong{For unconstrained models:}}{
#'       Should be size \eqn{((M(pd^2+d+d(d+1)/2+1)-1)x1)} and have the form
#'       \strong{\eqn{\theta}}\eqn{ = }(\strong{\eqn{\upsilon}}\eqn{_{1}},
#'       ...,\strong{\eqn{\upsilon}}\eqn{_{M}}, \eqn{\alpha_{1},...,\alpha_{M-1}}), where
#'       \itemize{
#'         \item \strong{\eqn{\upsilon}}\eqn{_{m}} \eqn{ = (\phi_{m,0},}\strong{\eqn{\phi}}\eqn{_{m}}\eqn{,\sigma_{m})}
#'         \item \strong{\eqn{\phi}}\eqn{_{m}}\eqn{ = (vec(A_{m,1}),...,vec(A_{m,p})}
#'         \item and \eqn{\sigma_{m} = vech(\Omega_{m})}, m=1,...,M.
#'       }
#'     }
#'     \item{\strong{For constrained models:}}{
#'       Should be size \eqn{((M(d+d(d+1)/2+1)+q-1)x1)} and have the form
#'       \strong{\eqn{\theta}}\eqn{ = (\phi_{1,0},...,\phi_{M,0},}\strong{\eqn{\psi},}
#'       \eqn{\sigma_{1},...,\sigma_{M},\alpha_{1},...,\alpha_{M-1})}, where
#'       \itemize{
#'         \item \strong{\eqn{\psi}} \eqn{(qx1)} satisfies (\strong{\eqn{\phi}}\eqn{_{1}}\eqn{,...,}
#'         \strong{\eqn{\phi}}\eqn{_{M}) =} \strong{\eqn{C \psi}} where \strong{\eqn{C}} is \eqn{(Mpd^2xq)}
#'         constraint matrix.
#'       }
#'     }
#'     \item{\strong{For same_means models:}}{
#'       Should have the form
#'       \strong{\eqn{\theta}}\eqn{ = (}\strong{\eqn{\mu},}\strong{\eqn{\psi},}
#'       \eqn{\sigma_{1},...,\sigma_{M},\alpha_{1},...,\alpha_{M-1})}, where
#'       \itemize{
#'         \item \strong{\eqn{\mu}}\eqn{= (\mu_{1},...,\mu_{g})} where
#'           \eqn{\mu_{i}} is the mean parameter for group \eqn{i} and
#'           \eqn{g} is the number of groups.
#'         \item If AR constraints are employed, \strong{\eqn{\psi}} is as for constrained
#'           models, and if AR constraints are not employed, \strong{\eqn{\psi}}\eqn{ = }
#'           (\strong{\eqn{\phi}}\eqn{_{1}}\eqn{,...,}\strong{\eqn{\phi}}\eqn{_{M})}.
#'       }
#'     }
#'     \item{\strong{For structural GMVAR model:}}{
#'       Should have the form
#'       \strong{\eqn{\theta}}\eqn{ = (\phi_{1,0},...,\phi_{M,0},}\strong{\eqn{\phi}}\eqn{_{1},...,}\strong{\eqn{\phi}}\eqn{_{M},
#'       vec(W),}\strong{\eqn{\lambda}}\eqn{_{2},...,}\strong{\eqn{\lambda}}\eqn{_{M},\alpha_{1},...,\alpha_{M-1})}, where
#'       \itemize{
#'         \item\strong{\eqn{\lambda}}\eqn{_{m}=(\lambda_{m1},...,\lambda_{md})} contains the eigenvalues of the \eqn{m}th mixture component.
#'       }
#'       \describe{
#'         \item{\strong{If AR parameters are constrained: }}{Replace \strong{\eqn{\phi}}\eqn{_{1}}\eqn{,...,}
#'         \strong{\eqn{\phi}}\eqn{_{M}} with \strong{\eqn{\psi}} \eqn{(qx1)} that satisfies (\strong{\eqn{\phi}}\eqn{_{1}}\eqn{,...,}
#'         \strong{\eqn{\phi}}\eqn{_{M}) =} \strong{\eqn{C \psi}}, as above.}
#'         \item{\strong{If same_means: }}{Replace \eqn{(\phi_{1,0},...,\phi_{M,0})} with \eqn{(\mu_{1},...,\mu_{g})},
#'           as above.}
#'         \item{\strong{If \eqn{W} is constrained:}}{Remove the zeros from \eqn{vec(W)} and make sure the other entries satisfy
#'          the sign constraints.}
#'         \item{\strong{If \eqn{\lambda_{mi}} are constrained:}}{Replace \strong{\eqn{\lambda}}\eqn{_{2},...,}\strong{\eqn{\lambda}}\eqn{_{M}}
#'          with \strong{\eqn{\gamma}} \eqn{(rx1)} that satisfies (\strong{\eqn{\lambda}}\eqn{_{2}}\eqn{,...,}
#'         \strong{\eqn{\lambda}}\eqn{_{M}) =} \strong{\eqn{C_{\lambda} \gamma}} where \eqn{C_{\lambda}} is a \eqn{(d(M-1) x r)}
#'          constraint matrix.}
#'       }
#'     }
#'   }
#'   Above, \eqn{\phi_{m,0}} is the intercept parameter, \eqn{A_{m,i}} denotes the \eqn{i}th coefficient matrix of the \eqn{m}th
#'   mixture component, \eqn{\Omega_{m}} denotes the error term covariance matrix of the \eqn{m}:th mixture component, and
#'   \eqn{\alpha_{m}} is the mixing weight parameter. The \eqn{W} and \eqn{\lambda_{mi}} are structural parameters replacing the
#'   error term covariance matrices (see Virolainen, 2020). If \eqn{M=1}, \eqn{\alpha_{m}} and \eqn{\lambda_{mi}} are dropped.
#'   If \code{parametrization=="mean"}, just replace each \eqn{\phi_{m,0}} with regimewise mean \eqn{\mu_{m}}.
#'   \eqn{vec()} is vectorization operator that stacks columns of a given matrix into a vector. \eqn{vech()} stacks columns
#'   of a given matrix from the principal diagonal downwards (including elements on the diagonal) into a vector.
#'   The notation is in line with the cited article by \emph{Kalliovirta, Meitz and Saikkonen (2016)} introducing the GMVAR model.
#' @param conditional a logical argument specifying whether the conditional or exact log-likelihood function
#'  should be used.
#' @param parametrization \code{"intercept"} or \code{"mean"} determining whether the model is parametrized with intercept
#'   parameters \eqn{\phi_{m,0}} or regime means \eqn{\mu_{m}}, m=1,...,M.
#' @param constraints a size \eqn{(Mpd^2 x q)} constraint matrix \strong{\eqn{C}} specifying general linear constraints
#'   to the autoregressive parameters. We consider constraints of form
#'   (\strong{\eqn{\phi}}\eqn{_{1}}\eqn{,...,}\strong{\eqn{\phi}}\eqn{_{M}) = }\strong{\eqn{C \psi}},
#'   where \strong{\eqn{\phi}}\eqn{_{m}}\eqn{ = (vec(A_{m,1}),...,vec(A_{m,p}) (pd^2 x 1), m=1,...,M},
#'   contains the coefficient matrices and \strong{\eqn{\psi}} \eqn{(q x 1)} contains the related parameters.
#'   For example, to restrict the AR-parameters to be the same for all regimes, set \strong{\eqn{C}}=
#'   [\code{I:...:I}]\strong{'} \eqn{(Mpd^2 x pd^2)} where \code{I = diag(p*d^2)}.
#'   Ignore (or set to \code{NULL}) if linear constraints should \strong{not} be employed.
#' @param same_means Restrict the mean parameters of some regimes to be the same? Provide a list of numeric vectors
#'   such that each numeric vector contains the regimes that should share the common mean parameters. For instance, if
#'   \code{M=3}, the argument \code{list(1, 2:3)} restricts the mean parameters of the second and third regime to be
#'   the same but the first regime has freely estimated (unconditional) mean. Ignore or set to \code{NULL} if mean parameters
#'   should not be restricted to be the same among any regimes. \strong{This constraint is available only for mean parametrized models;
#'   that is, when \code{parametrization="mean"}.}
#' @param structural_pars If \code{NULL} a reduced form model is considered. For structural model, should be a list containing
#'   the following elements:
#'   \itemize{
#'     \item \code{W} - a \eqn{(dxd)} matrix with its entries imposing constraints on \eqn{W}: \code{NA} indicating that the element is
#'       unconstrained, a positive value indicating strict positive sign constraint, a negative value indicating strict
#'       negative sign constraint, and zero indicating that the element is constrained to zero.
#'     \item \code{C_lambda} - a \eqn{(d(M-1) x r)} constraint matrix that satisfies (\strong{\eqn{\lambda}}\eqn{_{2}}\eqn{,...,}
#'       \strong{\eqn{\lambda}}\eqn{_{M}) =} \strong{\eqn{C_{\lambda} \gamma}} where \strong{\eqn{\gamma}} is the new \eqn{(r x 1)}
#'       parameter subject to which the model is estimated (similarly to AR parameter constraints). The entries of \code{C_lambda}
#'       must be either \strong{positive} or \strong{zero}. Ignore (or set to \code{NULL}) if the eigenvalues \eqn{\lambda_{mi}}
#'       should not be constrained.
#'   }
#'   See Virolainen (2020) for the conditions required to identify the shocks and for the B-matrix as well (it is \eqn{W} times
#'   a time-varying diagonal matrix with positive diagonal entries).
#' @param check_params should it be checked that the parameter vector satisfies the model assumptions? Can be skipped to save
#'   computation time if it does for sure.
#' @param minval the value that will be returned if the parameter vector does not lie in the parameter space
#'   (excluding the identification condition).
#' @param to_return should the returned object be the log-likelihood value, mixing weights, mixing weights including
#'   value for \eqn{alpha_{m,T+1}}, a list containing log-likelihood value and mixing weights, or
#'   the terms \eqn{l_{t}: t=1,..,T} in the log-likelihood function (see \emph{KMS 2016, eq.(9)})? Or should
#'   the regimewise conditional means, total conditional means, or total conditional covariance matrices
#'   be returned? Default is the log-likelihood value (\code{"loglik"}).
#' @details \code{loglikelihood_int} takes use of the function \code{dmvn} from the package \code{mvnfast}.
#' @return
#'  \describe{
#'   \item{By default:}{log-likelihood value of the specified GMVAR model,}
#'   \item{If \code{to_return=="mw"}:}{a size ((n_obs-p)xM) matrix containing the mixing weights: for m:th component in m:th column.}
#'   \item{If \code{to_return=="mw_tplus1"}:}{a size ((n_obs-p+1)xM) matrix containing the mixing weights: for m:th component in m:th column.
#'     The last row is for \eqn{\alpha_{m,T+1}}}.
#'   \item{If \code{to_return=="terms"}:}{a size ((n_obs-p)x1) numeric vector containing the terms \eqn{l_{t}}.}
#'   \item{if \code{to_return=="loglik_and_mw"}:}{a list of two elements. The first element contains the log-likelihood value and the
#'     second element contains the mixing weights.}
#'   \item{If \code{to_return=="regime_cmeans"}:}{an \code{[T-p, d, M]} array containing the regimewise conditional means
#'    (the first p values are used as the initial values).}
#'   \item{If \code{to_return=="total_cmeans"}:}{a \code{[T-p, d]} matrix containing the conditional means of the process
#'    (the first p values are used as the initial values).}
#'   \item{If \code{to_return=="total_ccov"}:}{an \code{[d, d, T-p]} array containing the conditional covariance matrices of the process
#'    (the first p values are used as the initial values).}
#'  }
#' @references
#'  \itemize{
#'    \item Kalliovirta L., Meitz M. and Saikkonen P. 2016. Gaussian mixture vector autoregression.
#'          \emph{Journal of Econometrics}, \strong{192}, 485-498.
#'    \item Lütkepohl H. 2005. New Introduction to Multiple Time Series Analysis,
#'          \emph{Springer}.
#'    \item McElroy T. 2017. Computation of vector ARMA autocovariances.
#'          \emph{Statistics and Probability Letters}, \strong{124}, 92-96.
#'    \item Virolainen S. 2020. Structural Gaussian mixture vector autoregressive model. Unpublished working
#'          paper, available as arXiv:2007.04713.
#'  }
#' @keywords internal

loglikelihood_int <- function(data, p, M, params, conditional=TRUE, parametrization=c("intercept", "mean"), constraints=NULL,
                              same_means=NULL, structural_pars=NULL,
                              to_return=c("loglik", "mw", "mw_tplus1", "loglik_and_mw", "terms", "regime_cmeans", "total_cmeans", "total_ccovs"),
                              check_params=TRUE, minval=NULL, stat_tol=1e-3, posdef_tol=1e-8) {

  # Compute required values
  epsilon <- round(log(.Machine$double.xmin) + 10) # Logarithm of the smallest value that can be handled normally
  d <- ncol(data)
  n_obs <- nrow(data)
  T_obs <- n_obs - p
  to_return <- match.arg(to_return)

  # Collect parameter values
  parametrization <- match.arg(parametrization)
  check_same_means(parametrization=parametrization, same_means=same_means)
  params <- reform_constrained_pars(p=p, M=M, d=d, params=params, constraints=constraints, same_means=same_means, structural_pars=structural_pars)
  W_constraints <- structural_pars$W
  structural_pars <- get_unconstrained_structural_pars(structural_pars=structural_pars)
  if(parametrization == "intercept") {
    all_phi0 <- pick_phi0(p=p, M=M, d=d, params=params, structural_pars=structural_pars)
  } else {
    mu <- pick_phi0(p=p, M=M, d=d, params=params, structural_pars=structural_pars) # mean parameters instead of phi0
  }
  all_A <- pick_allA(p=p, M=M, d=d, params=params, structural_pars=structural_pars) # A_{m,i}, m=1,...,M, i=1,..,p
  all_Omega <- pick_Omegas(p=p, M=M, d=d, params=params, structural_pars=structural_pars) # Omega_m
  all_boldA <- form_boldA(p=p, M=M, d=d, all_A=all_A) # The 'bold A' for each m=1,..,M, Lütkepohl 2005, eq.(2.1.8)
  alphas <- pick_alphas(p=p, M=M, d=d, params=params) # Mixing weight parameters

  # Check that the parameter vector lies in the parameter space (excluding indentifiability)
  if(check_params) {
    if(!in_paramspace_int(p=p, M=M, d=d, params=params, all_boldA=all_boldA, alphas=alphas, all_Omega=all_Omega,
                          W_constraints=W_constraints, stat_tol=stat_tol, posdef_tol=posdef_tol)) {
      return(minval)
    }
  }

  # An i:th row denotes the vector \bold{y_{i-1}} = (y_{i-1}',...,y_{i-p}') (dpx1),
  # assuming the observed data is y_{-p+1},...,y_0,y_1,...,y_{T}
  Y <- reform_data(data, p)

  # Calculate expected values (column per component) or phi0-parameters if using mean-parametrization
  Id <- diag(nrow=d)
  if(parametrization == "intercept") {
    mu <- vapply(1:M, function(m) solve(Id - rowSums(all_A[, , , m, drop=FALSE], dims=2), all_phi0[,m]), numeric(d)) # rowSums: sum over dims+1=3
  } else {
    all_phi0 <- vapply(1:M, function(m) (Id - rowSums(all_A[, , , m, drop=FALSE], dims=2))%*%mu[,m], numeric(d))
  }

  # Calculate the covariance matrices Sigma_{m,p} (Lutkepohl 2005, eq. (2.1.39) or the algorithm proposed by McElroy 2017)
  Sigmas <- get_Sigmas(p=p, M=M, d=d, all_A=all_A, all_boldA=all_boldA, all_Omega=all_Omega) # Store the (dpxdp) covariance matrices
  chol_Sigmas <- array(NA, dim=c(d*p, d*p, M))
  for(m in 1:M) {
    chol_Sigmas[, , m] <- chol(Sigmas[, , m]) # Take Cholesky here to avoid unnecessary warnings from mvnfast::dmvn
  }

  # Calculate the dp-dimensional multinormal densities (KMS 2016, eq.(6)), i:th row for index i-1 etc, m:th column for m:th component
  # Calculated in logarithm because same values may be too close to zero for machine accuracy
  log_mvnvalues <- vapply(1:M, function(m) mvnfast::dmvn(X=Y, mu=rep(mu[,m], p), sigma=chol_Sigmas[, , m], log=TRUE, ncores=1, isChol=TRUE), numeric(T_obs + 1))

  ## Calculate the mixing weights alpha_{m,t} (KMS 2016, eq.(7))
  if(to_return != "mw_tplus1") {
    log_mvnvalues <- log_mvnvalues[1:T_obs, , drop=FALSE] # alpha_mt uses y_{t-1} so the last row is not needed
  }
  alpha_mt_and_l_0 <- get_alpha_mt(M=M, log_mvnvalues=log_mvnvalues, alphas=alphas,
                                   epsilon=epsilon, conditional=conditional, also_l_0=TRUE)
  alpha_mt <- alpha_mt_and_l_0$alpha_mt
  l_0 <- alpha_mt_and_l_0$l_0 # The first term in the exact log-likelihood function (=0 for conditional)

  if(to_return == "mw" | to_return == "mw_tplus1") {
    return(alpha_mt)
  }

  # Calculate the conditional means mu_{m,t} (KMS 2016, Condition 1 (a)).
  # The dimensions of mu_mt will be: [t, p, m]
  all_A2 <- array(all_A, dim=c(d, d*p, M)) # cbind coefficient matrices of each component: m:th component is obtained at [, , m]
  Y2 <- Y[1:T_obs,] # Last row is not needed because mu_mt uses lagged values
  mu_mt <- array(vapply(1:M, function(m) t(all_phi0[, m] + tcrossprod(all_A2[, , m], Y2)), numeric(d*T_obs)), dim=c(T_obs, d, M)) # [, , m]

  if(to_return == "regime_cmeans") {
    return(mu_mt)
  } else if(to_return == "total_cmeans") { # KMS 2016, eq.(3)
    return(matrix(rowSums(vapply(1:M, function(m) alpha_mt[,m]*mu_mt[, , m], numeric(d*T_obs))), nrow=T_obs, ncol=d, byrow=FALSE))
  } else if(to_return == "total_ccovs") { # KMS 2016, eq.(4)
    first_term <- array(rowSums(vapply(1:M, function(m) rep(alpha_mt[, m], each=d*d)*as.vector(all_Omega[, , m]), numeric(d*d*T_obs))), dim=c(d, d, T_obs))
    sum_alpha_mu <- matrix(rowSums(vapply(1:M, function(m) alpha_mt[, m]*mu_mt[, , m], numeric(d*T_obs))), nrow=T_obs, ncol=d, byrow=FALSE)
    second_term <- array(rowSums(vapply(1:M, function(m) rep(alpha_mt[, m], each=d*d)*as.vector(vapply(1:nrow(alpha_mt), function(i1) tcrossprod((mu_mt[, , m] - sum_alpha_mu)[i1,]),
                                                                                                       numeric(d*d))), numeric(d*d*T_obs))), dim=c(d, d, T_obs))
    return(first_term + second_term)
  }

  # Calculate the second term of the log-likelihood (KMS 2016 eq.(10))
  dat <- data[(p + 1):n_obs,] # Initial values are not used here
  mvn_vals <- vapply(1:M, function(m) mvnfast::dmvn(X=dat - mu_mt[, , m], mu=rep(0, times=d), sigma=all_Omega[, , m], log=FALSE, ncores=1, isChol=FALSE), numeric(T_obs))
  weighted_mvn <- rowSums(alpha_mt*mvn_vals)
  weighted_mvn[weighted_mvn == 0] <- exp(epsilon)
  l_t <- log(weighted_mvn)

  if(to_return == "terms") {
     return(l_t)
  } else if(to_return == "loglik_and_mw") {
    return(list(loglik=l_0 + sum(l_t), mw=alpha_mt))
  } else {
    return(l_0 + sum(l_t))
  }
}


#' @title Get mixing weights alpha_mt (this function is for internal use)
#'
#' @description \code{get_alpha_mt} computes the mixing weights based on
#'   the logarithm of the multivariate normal densities in the definition of
#'   the mixing weights.
#'
#' @inheritParams loglikelihood_int
#' @param log_mvnvalues \eqn{T x M} matrix containing the log multivariate normal densities.
#' @param alphas \eqn{M x 1} vector containing the mixing weight pa
#' @param epsilon the smallest number such that its exponent is wont classified as numerically zero
#'   (around \code{-698} is used).
#' @param also_l_0 return also l_0 (the first term in the exact log-likelihood function)?
#' @details Note that we index the time series as \eqn{-p+1,...,0,1,...,T} as in Kalliovirta et al. (2016).
#' @return Returns the mixing weights a matrix of the same dimension as \code{log_mvnvalues} so
#'   that the t:th row is for the time point t and m:th column is for the regime m.
#' @inherit in_paramspace_int references
#' @seealso \code{\link{loglikelihood_int}}
#' @keywords internal

get_alpha_mt <- function(M, log_mvnvalues, alphas, epsilon, conditional, also_l_0=FALSE) {
  if(M == 1) {
    if(!is.matrix(log_mvnvalues)) log_mvnvalues <- as.matrix(log_mvnvalues) # Possibly many time points but only one regime
    alpha_mt <- as.matrix(rep(1, nrow(log_mvnvalues)))
  } else {
    if(!is.matrix(log_mvnvalues)) log_mvnvalues <- t(as.matrix(log_mvnvalues)) # Only one time point but multiple regimes

    log_mvnvalues_orig <- log_mvnvalues
    small_logmvns <- log_mvnvalues < epsilon
    if(any(small_logmvns)) {
      # If too small or large non-log-density values are present (i.e., that would yield -Inf or Inf),
      # we replace them with ones that are not too small or large but imply the same mixing weights
      # up to negligible numerical tolerance.
      which_change <- rowSums(small_logmvns) > 0 # Which rows contain too small  values
      to_change <- log_mvnvalues[which_change, , drop=FALSE]
      largest_vals <- do.call(pmax, split(to_change, f=rep(1:ncol(to_change), each=nrow(to_change)))) # The largest values of those rows
      diff_to_largest <- to_change - largest_vals # Differences to the largest value of the row

      # For each element in each row, check the (negative) distance from the largest value of the row. If the difference
      # is smaller than epsilon, replace the with epsilon. The results are then the new log_mvn values.
      diff_to_largest[diff_to_largest < epsilon] <- epsilon

      # Replace the old log_mvnvalues with the new ones
      log_mvnvalues[which_change,] <- diff_to_largest
    }

    mvnvalues <- exp(log_mvnvalues)
    denominator <- as.vector(mvnvalues%*%alphas)
    alpha_mt <- (mvnvalues/denominator)%*%diag(alphas)
  }
  if(!also_l_0) {
    return(alpha_mt)
  } else {
    # First term of the exact log-likelihood (Kalliovirta et al. 2016, eq.(9))
    l_0 <- 0
    if(M == 1 && conditional == FALSE) {
      l_0 <- log_mvnvalues[1,]
    } else if(M > 1 && conditional == FALSE) {
      if(any(log_mvnvalues_orig[1,] < epsilon)) { # Need to use Brobdingnag
        l_0 <- log(Reduce("+", lapply(1:M, function(i1) alphas[i1]*exp(Brobdingnag::as.brob(log_mvnvalues_orig[1, i1])))))
      } else {
        l_0 <- log(sum(alphas*mvnvalues[1,]))
      }
    }
    return(list(alpha_mt=alpha_mt,
                l_0=l_0))
  }
}


#' @title Compute log-likelihood of a GMVAR model using parameter vector
#'
#' @description \code{loglikelihood} computes log-likelihood of a GMVAR model using parameter vector
#'   instead of an object of class 'gmvar'. Exists for convenience if one wants to for example
#'   employ other estimation algorithms than the ones used in \code{fitGMVAR}. Use \code{minval} to
#'   control what happens when the parameter vector is outside the parameter space.
#'
#' @inheritParams loglikelihood_int
#' @return Returns log-likelihood if \code{params} is in the parameters space and \code{minval} if not.
#' @inherit loglikelihood_int details references
#' @seealso \code{\link{fitGMVAR}}, \code{\link{GMVAR}}, \code{\link{calc_gradient}}
#' @examples
#' # GMVAR(2, 2), d=2 model;
#' params22 <- c(0.36, 0.121, 0.223, 0.059, -0.151, 0.395, 0.406, -0.005,
#'  0.083, 0.299, 0.215, 0.002, 0.03, 0.484, 0.072, 0.218, 0.02, -0.119,
#'  0.722, 0.093, 0.032, 0.044, 0.191, 1.101, -0.004, 0.105, 0.58)
#' loglikelihood(data=gdpdef, p=2, M=2, params=params22)
#'
#' # Structural GMVAR(2, 2), d=2 model identified with sign-constraints:
#' params22s <- c(0.36, 0.121, 0.484, 0.072, 0.223, 0.059, -0.151, 0.395,
#'  0.406, -0.005, 0.083, 0.299, 0.218, 0.02, -0.119, 0.722, 0.093, 0.032,
#'  0.044, 0.191, 0.057, 0.172, -0.46, 0.016, 3.518, 5.154, 0.58)
#' W_22 <- matrix(c(1, 1, -1, 1), nrow=2, byrow=FALSE)
#' loglikelihood(data=gdpdef, p=2, M=2, params=params22s, structural_pars=list(W=W_22))
#' @export

loglikelihood <- function(data, p, M, params, conditional=TRUE, parametrization=c("intercept", "mean"), constraints=NULL,
                          same_means=NULL, structural_pars=NULL, minval=NA, stat_tol=1e-3, posdef_tol=1e-8) {
  if(!all_pos_ints(c(p, M))) stop("Arguments p and M must be positive integers")
  parametrization <- match.arg(parametrization)
  check_same_means(parametrization=parametrization, same_means=same_means)
  data <- check_data(data, p)
  d <- ncol(data)
  check_constraints(p=p, M=M, d=d, constraints=constraints, same_means=same_means, structural_pars=structural_pars)
  if(length(params) != n_params(p=p, M=M, d=d, constraints=constraints, same_means=same_means, structural_pars=structural_pars)) {
    stop("Parameter vector has wrong dimension")
  }
  loglikelihood_int(data=data, p=p, M=M, params=params, conditional=conditional, parametrization=parametrization,
                    constraints=constraints, same_means=same_means, structural_pars=structural_pars, to_return="loglik",
                    check_params=TRUE, minval=minval, stat_tol=stat_tol, posdef_tol=posdef_tol)
}


#' @title Compute conditional moments of a GMVAR model
#'
#' @description \code{loglikelihood} compute conditional regimewise means, conditional means, and conditional covariance matrices
#'  of a GMVAR model.
#'
#' @inheritParams loglikelihood_int
#' @param to_return should the regimewise conditional means, total conditional means, or total conditional covariance matrices
#'   be returned?
#' @details The first p values are used as the initial values, and by conditional we mean conditioning on the past. Formulas
#'   for the conditional means and covariance matrices are given in equations (3) and (4) of KMS (2016).
#' @return
#'  \describe{
#'   \item{If \code{to_return=="regime_cmeans"}:}{an \code{[T-p, d, M]} array containing the regimewise conditional means
#'    (the first p values are used as the initial values).}
#'   \item{If \code{to_return=="total_cmeans"}:}{a \code{[T-p, d]} matrix containing the conditional means of the process
#'    (the first p values are used as the initial values).}
#'   \item{If \code{to_return=="total_ccov"}:}{an \code{[d, d, T-p]} array containing the conditional covariance matrices of the process
#'    (the first p values are used as the initial values).}
#'  }
#' @inherit loglikelihood_int references
#' @family moment functions
#' @examples
#' # GMVAR(2, 2), d=2 model;
#' params22 <- c(0.36, 0.121, 0.223, 0.059, -0.151, 0.395, 0.406, -0.005,
#'  0.083, 0.299, 0.215, 0.002, 0.03, 0.484, 0.072, 0.218, 0.02, -0.119,
#'  0.722, 0.093, 0.032, 0.044, 0.191, 1.101, -0.004, 0.105, 0.58)
#' cond_moments(data=gdpdef, p=2, M=2, params=params22, to_return="regime_cmeans")
#' cond_moments(data=gdpdef, p=2, M=2, params=params22, to_return="total_cmeans")
#' cond_moments(data=gdpdef, p=2, M=2, params=params22, to_return="total_ccovs")
#' @export

cond_moments <- function(data, p, M, params, parametrization=c("intercept", "mean"), constraints=NULL, same_means=NULL,
                         structural_pars=NULL, to_return=c("regime_cmeans", "total_cmeans", "total_ccovs"),
                         stat_tol=1e-3, posdef_tol=1e-8) {
  if(!all_pos_ints(c(p, M))) stop("Arguments p and M must be positive integers")
  parametrization <- match.arg(parametrization)
  check_same_means(parametrization=parametrization, same_means=same_means)
  to_return <- match.arg(to_return)
  data <- check_data(data, p)
  d <- ncol(data)
  check_constraints(p=p, M=M, d=d, constraints=constraints, structural_pars=structural_pars)
  if(length(params) != n_params(p=p, M=M, d=d, constraints=constraints, same_means=same_means, structural_pars=structural_pars)) {
    stop("Parameter vector has wrong dimension")
  }
  loglikelihood_int(data=data, p=p, M=M, params=params, conditional=TRUE, parametrization=parametrization,
                    constraints=constraints, same_means=same_means, structural_pars=structural_pars, to_return=to_return, check_params=TRUE,
                    minval=NA, stat_tol=stat_tol, posdef_tol=posdef_tol)
}


#' @title Calculate AIC, HQIC, and BIC
#'
#' @description \code{get_IC} calculates the information criteria values AIC, HQIC, and BIC.
#'
#' @param loglik log-likelihood value
#' @param npars number of (freely estimated) parameters in the model
#' @param obs numbers of observations with starting values excluded for conditional models.
#' @details Note that for conditional models with different autoregressive order p the
#'  information criteria values are \strong{NOT} comparable.
#' @return Returns a data frame containing the information criteria values.
#' @keywords internal

get_IC <- function(loglik, npars, obs) {
  AIC <- -2*loglik + 2*npars
  HQIC <- -2*loglik + 2*npars*log(log(obs))
  BIC <- -2*loglik + npars*log(obs)
  data.frame(AIC=AIC, HQIC=HQIC, BIC=BIC)
}

