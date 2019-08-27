#' @import stats
#'
#' @title Compute the log-likelihood of Gaussian Mixture Vector Autoregressive model
#'
#' @description \code{loglikelihood_int} computes the log-likelihood of GMVAR model.
#'
#' @param data a matrix or class \code{'ts'} object with \code{d>1} columns. Each column is taken to represent
#'  a single time series. \code{NA} values are not supported.
#' @param p a positive integer specifying the autoregressive order of the model.
#' @param M a positive integer specifying the number of mixture components.
#' @param params a real valued vector specifying the parameter values.
#'   \describe{
#'     \item{\strong{For regular models:}}{
#'       Should be size \eqn{((M(pd^2+d+d(d+1)/2+1)-1)x1)} and have form
#'       \strong{\eqn{\theta}}\eqn{ = }(\strong{\eqn{\upsilon}}\eqn{_{1}},
#'       ...,\strong{\eqn{\upsilon}}\eqn{_{M}}, \eqn{\alpha_{1},...,\alpha_{M-1}}), where:
#'       \itemize{
#'         \item \strong{\eqn{\upsilon}}\eqn{_{m}} \eqn{ = (\phi_{m,0},}\strong{\eqn{\phi}}\eqn{_{m}}\eqn{,\sigma_{m})}
#'         \item \strong{\eqn{\phi}}\eqn{_{m}}\eqn{ = (vec(A_{m,1}),...,vec(A_{m,p})}
#'         \item and \eqn{\sigma_{m} = vech(\Omega_{m})}, m=1,...,M.
#'       }
#'     }
#'     \item{\strong{For constrained models:}}{
#'       Should be size \eqn{((M(d+d(d+1)/2+1)+q-1)x1)} and have form
#'       \strong{\eqn{\theta}}\eqn{ = (\phi_{1,0},...,\phi_{M,0},}\strong{\eqn{\psi}}
#'       \eqn{,\sigma_{1},...,\sigma_{M},\alpha_{1},...,\alpha_{M-1})}, where:
#'       \itemize{
#'         \item \strong{\eqn{\psi}} \eqn{(qx1)} satisfies (\strong{\eqn{\phi}}\eqn{_{1}}\eqn{,...,}
#'         \strong{\eqn{\phi}}\eqn{_{M}) =} \strong{\eqn{C \psi}}. Here \strong{\eqn{C}} is \eqn{(Mpd^2xq)}
#'         constraint matrix.
#'       }
#'     }
#'   }
#'   Above \eqn{\phi_{m,0}} is the intercept parameter, \eqn{A_{m,i}} denotes the \eqn{i}:th coefficient matrix of the \eqn{m}:th
#'   mixture component, \eqn{\Omega_{m}} denotes the error term covariance matrix of the \eqn{m}:th mixture component and
#'   \eqn{\alpha_{m}} is the mixing weight parameter.
#'   If \code{parametrization=="mean"}, just replace each \eqn{\phi_{m,0}} with regimewise mean \eqn{\mu_{m}}.
#'   \eqn{vec()} is vectorization operator that stacks columns of a given matrix into a vector. \eqn{vech()} stacks columns
#'   of a given matrix from the principal diagonal downwards (including elements on the diagonal) into a vector.
#'   The notations are in line with the cited article by \emph{Kalliovirta, Meitz and Saikkonen (2016)}.
#' @param conditional a logical argument specifying whether the conditional or exact log-likelihood function
#'  should be used. Default is \code{TRUE}.
#' @param parametrization \code{"mean"} or \code{"intercept"} determining whether the model is parametrized with regime means \eqn{\mu_{m}} or
#'   intercept parameters \eqn{\phi_{m,0}}, m=1,...,M. Default is \code{"intercept"}.
#' @param constraints a size \eqn{(Mpd^2 x q)} constraint matrix \strong{\eqn{C}} specifying general linear constraints
#'   to the autoregressive parameters. We consider constraints of form
#'   (\strong{\eqn{\phi}}\eqn{_{1}}\eqn{,...,}\strong{\eqn{\phi}}\eqn{_{M}) = }\strong{\eqn{C \psi}},
#'   where \strong{\eqn{\phi}}\eqn{_{m}}\eqn{ = (vec(A_{m,1}),...,vec(A_{m,p}) (pd^2 x 1), m=1,...,M}
#'   contains the coefficient matrices and \strong{\eqn{\psi}} \eqn{(q x 1)} contains the constrained parameters.
#'   For example, to restrict the AR-parameters to be the same for all regimes, set \strong{\eqn{C}}=
#'   [\code{I:...:I}]\strong{'} \eqn{(Mpd^2 x pd^2)} where \code{I = diag(p*d^2)}.
#'   Ignore (or set to \code{NULL}) if linear constraints should \strong{not} be employed.
#' @param check_params should it be checked that the parameter vector satisfies the model assumptions? Can be skipped to save
#'   computation time if it does for sure. Default \code{TRUE}.
#' @param minval value that will be returned if the parameter vector does not lie in the parameter space
#'   (excluding the identification condition).
#' @param to_return should the returned object be log-likelihood value, mixing weights, mixing weights including
#'   value for \eqn{alpha_{m,T+1}}, a list containing log-likelihood value and mixing weights or
#'   the terms \eqn{l_{t}: t=1,..,T} in the log-likelihood function (see \emph{KMS 2016, eq.(9)})? Or should
#'   the regimewise conditional means, total conditional means, or total conditional covariance matrices
#'   be returned? Default is the log-likelihood value (\code{"loglik"}).
#' @details Takes use of the function \code{dmvn} from the package \code{mvnfast} to cut down computation time.
#'   Values extremely close to zero are handled with the package \code{Brobdingnag}.
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
#'            \emph{Journal of Econometrics}, \strong{192}, 485-498.
#'    \item Lutkepohl H. 2005. New Introduction to Multiple Time Series Analysis,
#'            \emph{Springer}.
#'  }

loglikelihood_int <- function(data, p, M, params, conditional=TRUE, parametrization=c("intercept", "mean"), constraints=NULL,
                              to_return=c("loglik", "mw", "mw_tplus1", "loglik_and_mw", "terms", "regime_cmeans", "total_cmeans", "total_ccovs"),
                              check_params=TRUE, minval=NULL) {

  # Compute required values
  epsilon <- round(log(.Machine$double.xmin) + 10) # Logarithm of the smallest value that can be handled normally
  d <- ncol(data)
  n_obs <- nrow(data)
  T_obs <- n_obs - p
  to_return <- match.arg(to_return)

  # Collect parameter values
  parametrization <- match.arg(parametrization)
  params <- reform_constrained_pars(p=p, M=M, d=d, params=params, constraints=constraints)
  if(parametrization=="intercept") {
    all_phi0 <- pick_phi0(p=p, M=M, d=d, params=params)
  } else {
    mu <- pick_phi0(p=p, M=M, d=d, params=params) # mean parameters instead of phi0
  }
  all_A <- pick_allA(p=p, M=M, d=d, params=params)
  all_Omega <- pick_Omegas(p=p, M=M, d=d, params=params)
  all_boldA <- form_boldA(p=p, M=M, d=d, all_A=all_A)
  alphas <- pick_alphas(p=p, M=M, d=d, params=params)

  # Check that the parameter vector lies in the parameter space (excluding indentifiability)
  if(check_params) {
    if(!in_paramspace_int(p=p, M=M, d=d, all_boldA=all_boldA, alphas=alphas, all_Omega=all_Omega)) {
      return(minval)
    }
  }

  # An i:th row denotes the vector \bold{y_{i-1}} = (y_{i-1}',...,y_{i-p}') (dpx1), assuming the observed data is y_{-p+1},...,y_0,y_1,...,y_{T}
  Y <- reform_data(data, p)

  # Calculate expected values (column per component) or phi0-parameters if using mean-parametrization
  Id <- diag(nrow=d)
  if(parametrization == "intercept") {
    mu <- vapply(1:M, function(m) solve(Id - rowSums(all_A[, , , m, drop=FALSE], dims=2), all_phi0[,m]), numeric(d)) # rowSums: sum over dims+1=3
  } else {
    all_phi0 <- vapply(1:M, function(m) (Id - rowSums(all_A[, , , m, drop=FALSE], dims=2))%*%mu[,m], numeric(d))
  }

  # Calculate the covariance matrices Sigma_{m,p} (Lutkepohl 2005, eq. (2.1.39))
  I_dp2 <- diag(nrow=(d*p)^2)
  ZER_lower <- matrix(0, nrow=d*(p-1), ncol=d*p)
  ZER_right <- matrix(0, nrow=d, ncol=d*(p-1))
  Sigmas <- array(NA, dim=c(d*p, d*p, M)) # Store the (dpxdp) covariance matrices
  for(m in 1:M) {
    kronmat <- I_dp2 - kronecker(all_boldA[, , m], all_boldA[, , m])
    sigma_epsm <- rbind(cbind(all_Omega[, , m], ZER_right), ZER_lower)
    Sigma_m <- solve(kronmat, vec(sigma_epsm))
    Sigmas[, , m] <- Sigma_m
  }

  # Calculate the dp-dimensional multinormal densities (Kalliovirta ym. 2016, eq.(6)), i:th row for index i-1 etc, m:th column for m:th component
  # Calculated in logarithm because same values may be too close to zero for machine accuracy
  log_mvnvalues <- vapply(1:M, function(m) mvnfast::dmvn(X=Y, mu=rep(mu[,m], p), sigma=Sigmas[, , m], log=TRUE, ncores=1, isChol=FALSE), numeric(T_obs+1))


  # Calculate the mixing weights alpha_{m,t} (Kalliovirta et al. 2016, eq.(7))
  if(to_return != "mw_tplus1") {
    log_mvnvalues <- log_mvnvalues[1:T_obs, , drop=FALSE] # alpha_mt uses y_{t-1} so the last row is not needed
  }

  l_0 <- 0 # First term of the exact log-likelihood (Kalliovirta et al. 2016, eq.(9))
  if(M == 1) {
    alpha_mt <- as.matrix(rep(1, nrow(log_mvnvalues)))
    if(conditional == FALSE) {
      l_0 <- log_mvnvalues[1,]
    }
  } else if(any(log_mvnvalues < epsilon)) { # If some values are too close to zero use the package Brobdingnag
    numerators <- lapply(1:M, function(m) alphas[m]*exp(Brobdingnag::as.brob(log_mvnvalues[,m]))) # lapply(1:M, function(m) alphas[m]*Brobdingnag::as.brob(exp(1))^log_mvnvalues[,m])
    denominator <- Reduce('+', numerators)
    alpha_mt <- vapply(1:M, function(m) as.numeric(numerators[[m]]/denominator), numeric(nrow(log_mvnvalues)))

    if(conditional == FALSE) {
      l_0 <- log(Reduce('+', lapply(1:M, function(m) numerators[[m]][1])))
    }
  } else {
    mvnvalues <- exp(log_mvnvalues)
    denominator <- as.vector(mvnvalues%*%alphas)
    alpha_mt <- (mvnvalues/denominator)%*%diag(alphas)

    if(conditional == FALSE) {
      l_0 <- log(sum(alphas*mvnvalues[1,]))
    }
  }
  if(to_return == "mw" | to_return == "mw_tplus1") {
    return(alpha_mt)
  }

  # Calculate the conditional means mu_{m,t} (Kalliovirta et al. 2016, Condition 1 (a)).
  # Map to dimensions of mu_mt: [t, p, m]
  all_A2 <- array(all_A, dim=c(d, d*p, M)) # cbind coefficient matrices of each component: m:th component is obtained at [, , m]
  Y2 <- Y[1:T_obs,] # Last row is not needed because mu_mt uses lagged values
  mu_mt <- array(vapply(1:M, function(m) t(all_phi0[, m] + tcrossprod(all_A2[, , m], Y2)), numeric(d*T_obs)), dim=c(T_obs, d, M)) # [, , m]

  if(to_return == "regime_cmeans") {
    return(mu_mt)
  } else if(to_return == "total_cmeans") {
    return(matrix(rowSums(vapply(1:M, function(m) alpha_mt[,m]*mu_mt[, , m], numeric(d*T_obs))), nrow=T_obs, ncol=d, byrow=FALSE))
  } else if(to_return == "total_ccovs") {
    # array(vapply(1:nrow(alpha_mt), function(i1) rowSums(vapply(1:M, function(m) alpha_mt[i1, m]*all_Omega[, , m], numeric(d*d))),  numeric(d*d)), dim=c(d, d, T_obs))
    first_term <- array(rowSums(vapply(1:M, function(m) rep(alpha_mt[, m], each=d*d)*as.vector(all_Omega[, , m]), numeric(d*d*T_obs))), dim=c(d, d, T_obs))
    sum_alpha_mu <- matrix(rowSums(vapply(1:M, function(m) alpha_mt[, m]*mu_mt[, , m], numeric(d*T_obs))), nrow=T_obs, ncol=d, byrow=FALSE)
    # array(vapply(1:nrow(alpha_mt), function(i1) rowSums(vapply(1:M, function(m) alpha_mt[i1, m]*tcrossprod((mu_mt[, , m] - sum_alpha_mu)[i1,]), numeric(d*d))), numeric(d*d)), dim=c(d, d, T_obs))
    second_term <- array(rowSums(vapply(1:M, function(m) rep(alpha_mt[, m], each=d*d)*as.vector(vapply(1:nrow(alpha_mt), function(i1) tcrossprod((mu_mt[, , m] - sum_alpha_mu)[i1,]),
                                                                                                       numeric(d*d))), numeric(d*d*T_obs))), dim=c(d, d, T_obs))
    return(first_term + second_term)
  }

  # Calculate the second term of the log-likelihood (Kalliovirta et al. 2016 eq.(10))
  dat <- data[(p+1):n_obs,] # Initial values are not used here
  mvn_vals <- vapply(1:M, function(m) mvnfast::dmvn(X=dat-mu_mt[, , m], mu=rep(0, times=d), sigma=all_Omega[, , m], log=FALSE, ncores=1, isChol=FALSE), numeric(T_obs))
  l_t <- log(rowSums(alpha_mt*mvn_vals))

  if(to_return == "terms") {
     return(l_t)
  } else if(to_return == "loglik_and_mw") {
    return(list(loglik=l_0 + sum(l_t), mw=alpha_mt))
  } else {
    return(l_0 + sum(l_t))
  }
}



#' @title Compute log-likelihood of GMVAR model using parameter vector
#'
#' @description \code{loglikelihood} computes log-likelihood of GMVAR model by using parameter vector
#'   instead of object of class 'gmvar'. Exists for convenience if one wants to for example
#'   plot profile log-likelihoods or employ other estimation algorithms than used in \code{fitGMVAR}.
#'   Use \code{minval} to control what happens when the parameter vector is outside the parameter space.
#'
#' @inheritParams loglikelihood_int
#' @return Returns log-likelihood if \code{params} is in the parameters space and \code{minval} if not.
#' @inherit loglikelihood_int details references
#' @seealso \code{\link{fitGMVAR}}, \code{\link{GMVAR}}, \code{\link{calc_gradient}}
#' @examples
#' data <- cbind(10*eurusd[,1], 100*eurusd[,2])
#' params222 <- c(-11.904, 154.684, 1.314, 0.145, 0.094, 1.292, -0.389,
#'  -0.070, -0.109, -0.281, 0.920, -0.025, 4.839, 11.633, 124.983, 1.248,
#'   0.077, -0.040, 1.266, -0.272, -0.074, 0.034, -0.313, 5.855, 3.570,
#'   9.838, 0.740)
#' loglikelihood(data=data, p=2, M=2, params=params222, parametrization="mean")
#' @export

loglikelihood <- function(data, p, M, params, conditional=TRUE, parametrization=c("intercept", "mean"), constraints=NULL, minval=NA) {
  if(!all_pos_ints(c(p, M))) stop("Arguments p and M must be positive integers")
  parametrization <- match.arg(parametrization)
  data <- check_data(data, p)
  d <- ncol(data)
  check_constraints(p=p, M=M, d=d, constraints=constraints)
  if(length(params) != n_params(p=p, M=M, d=d, constraints=constraints)) stop("Parameter vector has wrong dimension")
  loglikelihood_int(data, p, M, params, conditional=conditional, parametrization=parametrization,
                    constraints=constraints, to_return="loglik", check_params=TRUE, minval=minval)
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
#'   for the conditional means and covariance matrices are given in equations (3) and (4) of Kalliovirta et al. (2016).
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
#' data <- cbind(10*eurusd[,1], 100*eurusd[,2])
#' params222 <- c(-11.904, 154.684, 1.314, 0.145, 0.094, 1.292, -0.389,
#'  -0.070, -0.109, -0.281, 0.920, -0.025, 4.839, 11.633, 124.983, 1.248,
#'   0.077, -0.040, 1.266, -0.272, -0.074, 0.034, -0.313, 5.855, 3.570,
#'   9.838, 0.740)
#' cond_moments(data=data, p=2, M=2, params=params222, parametrization="mean",
#'   to_return="regime_cmeans")
#' cond_moments(data=data, p=2, M=2, params=params222, parametrization="mean",
#'   to_return="total_cmeans")
#' cond_moments(data=data, p=2, M=2, params=params222, parametrization="mean",
#'   to_return="total_ccovs")
#' @export

cond_moments <- function(data, p, M, params, parametrization=c("intercept", "mean"), constraints=NULL,
                         to_return=c("regime_cmeans", "total_cmeans", "total_ccovs")) {
  if(!all_pos_ints(c(p, M))) stop("Arguments p and M must be positive integers")
  parametrization <- match.arg(parametrization)
  to_return <- match.arg(to_return)
  data <- check_data(data, p)
  d <- ncol(data)
  check_constraints(p=p, M=M, d=d, constraints=constraints)
  if(length(params) != n_params(p=p, M=M, d=d, constraints=constraints)) stop("Parameter vector has wrong dimension")
  loglikelihood_int(data, p, M, params, conditional=TRUE, parametrization=parametrization,
                    constraints=constraints, to_return=to_return, check_params=TRUE, minval=NA)
}


#' @title Calculate AIC, HQIC and BIC
#'
#' @description \code{get_IC} calculates information criteria values AIC, HQIC and BIC.
#'
#' @param loglik log-likelihood value
#' @param npars number of (freely estimated) parameters in the model
#' @param obs numbers of observations with starting values excluded for conditional models.
#' @return Returns a data frame containing the information criteria values.

get_IC <- function(loglik, npars, obs) {
  AIC <- -2*loglik + 2*npars
  HQIC <- -2*loglik + 2*npars*log(log(obs))
  BIC <- -2*loglik + npars*log(obs)
  data.frame(AIC=AIC, HQIC=HQIC, BIC=BIC)
}

