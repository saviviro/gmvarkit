#' @title Create random mean-parametrized parameter vector of a GMVAR, StMVAR, or G-StMVAR model that may not be stationary
#'
#' @description \code{random_ind} generates random mean-parametrized parameter vector that may not be stationary.
#'
#' @inheritParams GAfit
#' @inheritParams is_stationary
#' @return Returns random mean-parametrized parameter vector that has the same form as the argument \code{params}
#'   in the other functions, for instance, in the function \code{loglikelihood}.
#' @inherit in_paramspace references
#' @keywords internal

random_ind <- function(p, M, d, model=c("GMVAR", "StMVAR", "G-StMVAR"), constraints=NULL, same_means=NULL,
                       weight_constraints=NULL, structural_pars=NULL,
                       mu_scale, mu_scale2, omega_scale, W_scale, lambda_scale, ar_scale2=1) {
  model <- match.arg(model)
  M_orig <- M
  M <- sum(M)
  scale_A <- ar_scale2*ifelse(is.null(constraints),
                              1 + log(2*mean(c((p - 0.2)^(1.25), d))),
                              1 + (sum(constraints)/(M*d^2))^0.85)
  g <- ifelse(is.null(same_means), M, length(same_means)) # Number of groups of regimes with the same mean parameters
  # AR and covmat pars
  if(is.null(constraints)) {
    if(is.null(structural_pars)) {
      if(is.null(same_means)) { # No AR constraints, reduced form, no same_means
        x <- as.vector(vapply(1:M, function(m) c(rnorm(d, mean=mu_scale, sd=mu_scale2),
                                                 random_coefmats(d=d, how_many=p, scale=scale_A),
                                                 random_covmat(d=d, omega_scale=omega_scale)), numeric(p*d^2 + d + d*(d+1)/2)))
      } else { # No AR constraints, reduced form, same_means
        x <- c(rnorm(d*g, mean=mu_scale, sd=mu_scale2),
               replicate(n=M, expr=random_coefmats(d=d, how_many=p, scale=scale_A)),
               replicate(n=M, expr=random_covmat(d=d, omega_scale=omega_scale)))
      }
    } else { # No AR constraints, structural model, possibly with same_means
      x <- c(rnorm(d*g, mean=mu_scale, sd=mu_scale2),
             replicate(n=M, random_coefmats(d=d, how_many=p, scale=scale_A)),
             random_covmat(d=d, M=M, W_scale=W_scale, lambda_scale=lambda_scale, structural_pars=structural_pars))
    }
  } else { # AR constraints employed
    q <- ncol(constraints)
    psi <- rnorm(q, mean=0, sd=0.5/scale_A) # random psi
    all_phi0 <- rnorm(d*g, mean=mu_scale, sd=mu_scale2)
    if(is.null(structural_pars)) {
      x <- c(all_phi0, psi, replicate(n=M, random_covmat(d=d, omega_scale=omega_scale)))
    } else {
      x <- c(all_phi0, psi, random_covmat(d=d, M=M, W_scale=W_scale, lambda_scale=lambda_scale, structural_pars=structural_pars))
    }
  }
  if(M > 1 && is.null(weight_constraints)) {
    alphas <- runif(n=M)
    alphas <- sort_and_standardize_alphas(alphas=alphas, constraints=constraints, same_means=same_means,
                                          structural_pars=structural_pars)
    ret <- c(x, alphas[-M])
  } else { # No alpha params
    ret <- x
  }
  c(ret, random_df(M=M_orig, model=model))
}


#' @title Create random parameter vector of a GMVAR, StMVAR, or G-StMVAR model fairly close to a given
#'   parameter vector
#'
#' @description \code{smart_ind} creates random mean-parametrized parameter vector of a GMVAR, StMVAR, or G-StMVAR
#'   model fairly close to a given parameter vector. The result may not be stationary.
#'
#' @inheritParams is_stationary
#' @inheritParams GAfit
#' @inheritParams random_coefmats2
#' @param accuracy a positive real number adjusting how close to the given parameter vector the returned individual should be.
#'   Larger number means larger accuracy. Read the source code for details.
#' @param which_random a vector with length between 1 and M specifying the mixture components that should be random instead of
#'   close to the given parameter vector. This does not consider constrained AR or lambda parameters.
#' @section Warning:
#'   No argument checks!
#' @inherit random_ind return references
#' @keywords internal

smart_ind <- function(p, M, d, params, model=c("GMVAR", "StMVAR", "G-StMVAR"), constraints=NULL, same_means=NULL,
                      weight_constraints=NULL, structural_pars=NULL, accuracy=1, which_random=numeric(0),
                      mu_scale, mu_scale2, omega_scale, ar_scale=1, ar_scale2=1, W_scale, lambda_scale) {
  model <- match.arg(model)
  M_orig <- M
  M <- sum(M)
  scale_A <- ar_scale2*(1 + log(2*mean(c((p - 0.2)^(1.25), d))))
  params_std <- reform_constrained_pars(p=p, M=M_orig, d=d, params=params, model=model,
                                        constraints=constraints, same_means=same_means,
                                        weight_constraints=weight_constraints,
                                        structural_pars=structural_pars)
  unc_structural_pars <- get_unconstrained_structural_pars(structural_pars=structural_pars)
  alphas <- pick_alphas(p=p, M=M_orig, d=d, params=params_std, model=model)
  all_df <- pick_df(M=M_orig, params=params_std, model=model)
  if(is.null(structural_pars)) {
    all_Omega <- pick_Omegas(p=p, M=M_orig, d=d, params=params_std)
  }
  if(is.null(constraints) && is.null(structural_pars) && is.null(same_means)) {
    # No AR constraints, reduced form model, no same_means
    all_phi0_A <- pick_all_phi0_A(p=p, M=M_orig, d=d, params=params_std) # all_mu if called from GA
    pars <- vapply(1:M, function(m) {
      if(any(which_random == m)) {
        if(runif(1) > 0.5) { # Use algorithm to force stationarity of coefficient matrices
          coefmats <- random_coefmats2(p=p, d=d, ar_scale=ar_scale)
        } else {
          coefmats <- random_coefmats(d=d, how_many=p, scale=scale_A)
        }
        c(rnorm(d, mean=mu_scale, sd=mu_scale2),
          coefmats,
          random_covmat(d=d, omega_scale=omega_scale))
      } else {
        c(rnorm(length(all_phi0_A[,m]), mean=all_phi0_A[,m], sd=pmax(0.2, abs(all_phi0_A[,m]))/accuracy),
          smart_covmat(d=d, Omega=all_Omega[, , m], accuracy=accuracy))
      }
    }, numeric(d + p*d^2 + d*(d + 1)/2))
  } else { # AR constraints, structural model, or same_means
    g <- ifelse(is.null(same_means), M, length(same_means)) # Number of groups of regimes with the same mean parameters
    if(length(which_random) == 0) {
      smart_regs <- 1:M
    } else {
      smart_regs <- (1:M)[-which_random]
    }
    less_pars <- ifelse(is.null(same_means), 0, d*(M - g)) # Number of (mean) parameters less in same_means models
    all_phi0 <- matrix(params[1:(d*g)], nrow=d, ncol=g, byrow=FALSE) # Always mean parameters (not intercept) when called from GA
    phi0_pars <- as.vector(vapply(1:g, function(m) {
      which_reg <- ifelse(is.null(same_means), m, same_means[[m]]) # Can be many if same_means used
      if(any(which_reg %in% smart_regs)) { # Smart parameters
        rnorm(d, mean=all_phi0[,m], sd=abs(all_phi0[,m]/accuracy))
      } else { # Random parameters
        rnorm(d, mean=mu_scale, sd=mu_scale2)
      }
    }, numeric(d)))
    if(is.null(constraints)) { # Structural model with AR parameters not constrained, possibly with same_means
      all_A <- pick_allA(p=p, M=M, d=d, params=params_std, structural_pars=unc_structural_pars)
      AR_pars <- as.vector(vapply(1:M, function(m) {
        if(any(which_random == m)) {
          if(runif(1) > 0.5) { # Use algorithm to force stationarity of coefficient matrices
            random_coefmats2(p=p, d=d, ar_scale=ar_scale)
          } else {
            random_coefmats(d=d, how_many=p, scale=scale_A)
          }
        } else {
          all_Am <- as.vector(all_A[, , , m])
          rnorm(n=length(all_Am), mean=all_Am, sd=pmax(0.2, abs(all_Am))/accuracy)
        }
      }, numeric(p*d^2)))
    } else { # Structural or reduced form model with AR parameters constrained, possibly with same_means
      q <- ncol(constraints)
      psi <- params[(M*d + 1 - less_pars):(M*d + q - less_pars)]
      AR_pars <- rnorm(q, mean=psi, sd=pmax(0.2, abs(psi))/accuracy)
    }
    if(is.null(structural_pars)) { # Reduced form model, possibly with same_means
      covmat_pars <- as.vector(vapply(1:M, function(m) {
        if(any(which_random == m)) {
          random_covmat(d=d, omega_scale=omega_scale)
        } else {
          smart_covmat(d=d, Omega=all_Omega[, , m], accuracy=accuracy)
        }
      }, numeric(d*(d + 1)/2)))
    } else { # Structural model, possibly with same_means
      # For the covarince matrix parameters, there is the problem that the parameters of the first regime (W)
      # affect the covariance matrices of the other regimes as well. So using random W messes up the covariance
      # matrices of the other regimes too even if the lambda-parameters are "smart mutated". For this reason,
      # if the first regime is "redundant" and randomly mutated, with probability 0.5 all the covariance matrix parameters
      # will be random, and with probability 0.5 the W parameters will be smart but lambda parameters of redundant
      # regimes are still random.
      if(any(which_random == 1)  && runif(1) < 0.5) {
        # If first regime is random, then W must be random so the lambdas may as well be random too.
        covmat_pars <- random_covmat(d=d, M=M, W_scale=W_scale, lambda_scale=lambda_scale, structural_pars=structural_pars)
      } else { # First regime is smart
        W_pars <- Wvec(pick_W(p=p, M=M, d=d, params=params_std, structural_pars=unc_structural_pars))
        if(M > 1 && is.null(structural_pars$fixed_lambdas)) {
          n_lambs <- ifelse(is.null(structural_pars$C_lambda), d*(M - 1), ncol(structural_pars$C_lambda))
          W_and_lambdas <- c(W_pars, params[(length(params) - (M - 1) - n_lambs + 1):(length(params) - (M - 1))])
        } else {
          W_and_lambdas <- W_pars # No lambdas when M == 1 or fixed_lambdas are used
        }
        covmat_pars <- smart_covmat(d=d, M=M, W_and_lambdas=W_and_lambdas, accuracy=accuracy, structural_pars=structural_pars)
      }
      if(is.null(structural_pars$C_lambda) && is.null(structural_pars$fixed_lambdas) && M > 1) {
        # If lambdas are not constrained, we can replace smart lambdas of some regimes with random lambdas
        for(m in 2:M) {
          if(any(which_random == m)) {
            covmat_pars[(length(covmat_pars) - d*(M - 1) + d*(m - 2) +
                           1):(length(covmat_pars) - d*(M - 1) + d*(m - 2) + d)] <- abs(rnorm(n=d, mean=0, sd=lambda_scale[m - 1]))
          }
        }
      }
    }
    pars <- c(phi0_pars, AR_pars, covmat_pars)
  }
  if(M > 1 && is.null(weight_constraints)) {
    alphas <- abs(rnorm(M, mean=alphas, sd=0.1))
    ret <- c(pars, (alphas/sum(alphas))[-M])
  } else { # No alpha params
    ret <- pars
  }
  c(ret, smart_df(M=M_orig, df=all_df, accuracy=accuracy, which_random=which_random, model=model))
}


#' @title Create somewhat random parameter vector of a GMVAR, StMVAR, or G-StMVAR model that is always stationary
#'
#' @description \code{random_ind2} generates random mean-parametrized parameter vector
#'  that is always stationary.
#'
#' @inheritParams GAfit
#' @inheritParams is_stationary
#' @details The coefficient matrices are generated using the algorithm proposed by Ansley
#'   and Kohn (1986) which forces stationarity. It's not clear in detail how \code{ar_scale}
#'   exactly affects the coefficient matrices but larger \code{ar_scale} seems to result in larger
#'   AR coefficients. Read the cited article by Ansley and Kohn (1986) and the source code
#'   for more information.
#'
#'   The covariance matrices are generated from (scaled) Wishart distribution.
#'
#'   Models with AR parameters constrained are not supported!
#' @inherit random_ind return
#' @references
#'  \itemize{
#'    \item Ansley C.F., Kohn R. 1986. A note on reparameterizing a vector autoregressive
#'       moving average model to enforce stationarity.
#'       \emph{Journal of statistical computation and simulation}, \strong{24}:2, 99-106.
#'    \item Kalliovirta L., Meitz M. and Saikkonen P. 2016. Gaussian mixture vector autoregression.
#'          \emph{Journal of Econometrics}, \strong{192}, 485-498.
#'    \item Virolainen S. (forthcoming). A statistically identified structural vector autoregression with endogenously
#'           switching volatility regime. \emph{Journal of Business & Economic Statistics}.
#'    \item Virolainen S. 2022. Gaussian and Student's t mixture vector autoregressive model with application to the
#'      asymmetric effects of monetary policy shocks in the Euro area. Unpublished working
#'      paper, available as arXiv:2109.13648.
#'  }
#'  @keywords internal

random_ind2 <- function(p, M, d, model=c("GMVAR", "StMVAR", "G-StMVAR"), same_means=NULL,
                        weight_constraints=NULL, structural_pars=NULL,
                        mu_scale, mu_scale2, omega_scale, ar_scale=1, W_scale, lambda_scale) {
  model <- match.arg(model)
  M_orig <- M
  M <- sum(M)
  if(is.null(structural_pars) && is.null(same_means)) { # Reduced form, no same_means
      x <- as.vector(vapply(1:M, function(m) c(rnorm(d, mean=mu_scale, sd=mu_scale2),
                                               random_coefmats2(p=p, d=d, ar_scale=ar_scale),
                                               random_covmat(d=d, omega_scale=omega_scale, structural_pars=structural_pars)),
                            numeric(p*d^2 + d + d*(d + 1)/2)))
  } else { # Structural model or same_means
   g <- ifelse(is.null(same_means), M, length(same_means)) # Number of groups of regimes with the same mean parameters
   n_covmats <- ifelse(is.null(structural_pars), M, 1) # Only one covmat for structural models (that includes all lambdas)
   x <- c(rnorm(d*g, mean=mu_scale, sd=mu_scale2),
          replicate(n=M, random_coefmats2(p=p, d=d, ar_scale=ar_scale)),
          replicate(n=n_covmats,
                    random_covmat(d=d, M=M, omega_scale=omega_scale, W_scale=W_scale,
                                  lambda_scale=lambda_scale, structural_pars=structural_pars)))
  }
  if(M > 1 && is.null(weight_constraints)) {
    alphas <- runif(n=M)
    alphas <- sort_and_standardize_alphas(alphas=alphas, constraints=NULL, same_means=same_means,
                                          weight_constraints=weight_constraints,
                                          structural_pars=structural_pars)
    ret <- c(x, alphas[-M])
  } else {
    ret <- x
  }
  c(ret, random_df(M=M_orig, model=model))
}



#' @title Create random VAR-model \eqn{(dxd)} coefficient matrices \eqn{A}.
#'
#' @description \code{random_coefmats} generates random VAR model coefficient matrices.
#'
#' @inheritParams is_stationary
#' @param how_many how many \eqn{(d\times d)} coefficient matrices \eqn{A} should be drawn?
#' @param scale non-diagonal elements will be drawn from mean zero normal distribution
#'   with \code{sd=0.3/scale} and diagonal elements from one with \code{sd=0.6/scale}.
#'   Larger scale will hence more likely result stationary coefficient matrices, but
#'   will explore smaller area of the parameter space. Can be for example
#'   \code{1 + log(2*mean(c((p-0.2)^(1.25), d)))}.
#' @return Returns \eqn{((how_many*d^2)\times 1)} vector containing vectorized coefficient
#'  matrices \eqn{(vec(A_{1}),...,vec(A_{how_many}))}. Note that if \code{how_many==p},
#'  then the returned vector equals \strong{\eqn{\phi_{m}}}.
#' @keywords internal

random_coefmats <- function(d, how_many, scale) {
  as.vector(vapply(1:how_many, function(i1) {
    x <- rnorm(d*d, mean=0, sd=0.3/scale)
    x[1 + 0:(d - 1) * (d + 1)] <- rnorm(d, mean=0, sd=0.6/scale)
    x
  }, numeric(d*d)))
}


#' @title Create random stationary VAR model \eqn{(dxd)} coefficient matrices \eqn{A}.
#'
#' @description \code{random_coefmats2} generates random VAR model coefficient matrices.
#'
#' @inheritParams is_stationary
#' @param ar_scale a positive real number between zero and one. Larger values will typically result
#'   larger AR coefficients.
#' @details The coefficient matrices are generated using the algorithm proposed by Ansley
#'   and Kohn (1986) which forces stationarity. It's not clear in detail how \code{ar_scale}
#'   affects the coefficient matrices. Read the cited article by Ansley and Kohn (1986) and
#'   the source code for more information.
#'
#'   Note that when using large \code{ar_scale} with large \code{p} or \code{d}, numerical
#'   inaccuracies caused by the imprecision of the float-point presentation may result in errors
#'   or nonstationary AR-matrices. Using smaller \code{ar_scale} facilitates the usage of larger
#'   \code{p} or \code{d}.
#' @return Returns \eqn{((pd^2)\times 1)} vector containing stationary vectorized coefficient
#'  matrices \eqn{(vec(A_{1}),...,vec(A_{p})}.
#' @references
#'  \itemize{
#'    \item Ansley C.F., Kohn R. 1986. A note on reparameterizing a vector autoregressive
#'       moving average model to enforce stationarity.
#'       \emph{Journal of statistical computation and simulation}, \strong{24}:2, 99-106.
#'  }
#' @keywords internal

random_coefmats2 <- function(p, d, ar_scale=1) {
  # First generate matrices P_1,..,P_p with singular values less than one
  stopifnot(ar_scale > 0 && ar_scale <= 1)
  Id <- diag(x=1, nrow=d)
  all_P <- array(dim=c(d, d, p))
  for(i1 in 1:p) {
    A <- matrix(rnorm(d*d, sd=ar_scale), nrow=d)
    B <- t(chol(Id + tcrossprod(A, A)))
    all_P[, , i1] <- solve(B, A)
  }

  all_phi <- array(dim=c(d, d, p, p)) # [ , , i, j] for phi_{i, j}
  all_phi_star <- array(dim=c(d, d, p, p)) # [ , , i, j] for phi_{i, j}*

  # Set initial values
  L <- L_star <- Sigma <- Sigma_star <- Gamma <- Id

  # Recursion algorithm (Ansley and Kohn 1986, lemma 2.1)
  for(s in 0:(p - 1)) {
    all_phi[, , s+1, s+1] <- L%*%all_P[, , s+1]%*%solve(L_star)
    all_phi_star[, , s+1, s+1] <- tcrossprod(L_star, all_P[, , s+1])%*%solve(L)

    if(s >= 1) {
      for(k in 1:s) {
        all_phi[, , s+1, k] <- all_phi[, , s, k] - all_phi[, , s+1, s+1]%*%all_phi_star[, , s, s-k+1]
        all_phi_star[, , s+1, k] <- all_phi_star[, , s, k] - all_phi_star[, , s+1, s+1]%*%all_phi[, , s, s-k+1]
      }
    }

    if(s < p - 1) { # These are not needed in the last round because only coefficient matrices will be returned.
      Sigma_next <- Sigma - all_phi[, , s+1, s+1]%*%tcrossprod(Sigma_star, all_phi[, , s+1, s+1])
      if(s < p + 1) {
        Sigma_star <- Sigma_star - all_phi_star[, , s+1, s+1]%*%tcrossprod(Sigma, all_phi_star[, , s+1, s+1])
        L_star <- t(chol(Sigma_star))
      }
      Sigma <- Sigma_next
      L <- t(chol(Sigma))
    }
  }
  #all_A <- all_phi[, , p, 1:p]
  #as.vector(all_A)
  as.vector(all_phi[, , p, 1:p])
}


#' @title Create random VAR model error term covariance matrix
#'
#' @description \code{random_covmat} generates random VAR model \eqn{(d\times d)} error term covariance matrix \eqn{\Omega}
#'   from (scaled) Wishart distribution for reduced form models and the parameters \eqn{W},\eqn{\lambda_1,...,\lambda_M}
#'   for structural models (from normal distributions).
#'
#' @inheritParams is_stationary
#' @inheritParams GAfit
#' @details Note that for StMVAR type regimes, the error term covariance matrix is consists of an ARCH type scalar that
#'   multiplies a constant covariance matrix. This function generates the constant covariance matrix part of the
#'   error term covariance matrix.
#' @return
#'   \describe{
#'     \item{For \strong{reduced form models}:}{Returns a \eqn{(d(d+1)/2\times 1)} vector containing vech-vectorized covariance matrix
#'       \eqn{\Omega}.}
#'     \item{For \strong{structural models}:}{Returns a length \eqn{d^2 - n_zeros - d*(M - 1)} vector of the form
#'       \eqn{(Wvec(W),\lambda_2,...,\lambda_M)} where \eqn{\lambda_m=(\lambda_{m1},...,\lambda_{md})}
#'       contains the eigenvalue parameters of the \eqn{m}th regime \eqn{(m>1)} and \eqn{n_zeros} is the number of zero constraints
#'       in \eqn{W}. If lambdas are \code{C_lambda} constrained, replace \eqn{d*(M - 1)} in the length with \eqn{r} and
#'       \eqn{\lambda_2,...,\lambda_M)} with \strong{\eqn{\gamma}}. If \code{fixed_lambdas} are used, the \eqn{\lambda_{mi}} parameters
#'        are not included. The operator \eqn{Wvec()} vectorizes a matrix and removes zeros.}
#'   }
#' @keywords internal

random_covmat <- function(d, M, omega_scale, W_scale, lambda_scale, structural_pars=NULL) {
  if(is.null(structural_pars)) {
    return(smart_covmat(d=d, Omega=diag(x=omega_scale), accuracy=1))
  } else {
    M <- sum(M)
    W <- structural_pars$W
    n_zeros <- vapply(1:d, function(i1) sum(W[i1,] == 0, na.rm=TRUE), numeric(1)) # Number of zeros in each row
    std_devs <- sqrt(W_scale/(d - n_zeros))
    new_W <- matrix(rnorm(n=d^2, mean=0, sd=std_devs), nrow=d, byrow=FALSE) # The standard deviations are recycled
    new_W[W == 0 & !is.na(W)] <- 0
    new_W[W > 0 & !is.na(W)] <- abs(new_W[W > 0 & !is.na(W)])
    new_W[W < 0 & !is.na(W)] <- -abs(new_W[W < 0 & !is.na(W)])
    if(M > 1 && is.null(structural_pars$fixed_lambdas)) {
      if(is.null(structural_pars$C_lambda)) {
        lambdas <- as.vector(abs(vapply(1:(M - 1), function(i1) rnorm(n=d, mean=0, sd=lambda_scale[i1]), numeric(d))))
      } else {
        lambdas <- abs(rnorm(n=ncol(structural_pars$C_lambda), mean=0, sd=lambda_scale)) # gammas
      }
      return(c(Wvec(new_W), lambdas))
    } else { # No lambda params
      return(Wvec(new_W))
    }
  }
}


#' @title Create random VAR-model \eqn{(d\times d)} error term covariance matrix \eqn{\Omega}
#'   fairly close to a given \strong{positive definite} covariance matrix using (scaled)
#'   Wishart distribution
#'
#' @description \code{random_covmat} generates random VAR model \eqn{(dxd)} error term covariance matrix \eqn{\Omega}
#'   from (scaled) Wishart distribution that is fairly close to the given matrix.
#'
#' @inheritParams is_stationary
#' @param Omega a symmetric positive definite \eqn{(d \times d)} covariance matrix specifying
#'   expected value of the matrix to be generated.
#' @param W_and_lambdas the mean of the normal distribution the new parameters are generated
#'   from.
#'   \describe{
#'     \item{If lambdas are \strong{not constrained}:}{a size \eqn{(d^2 - n_zeros + d*(M - 1))} vector
#'       \eqn{(Wvec(W),\lambda_{2},...,\lambda{M})}, where \eqn{n_zeros} is the number of zero constraints
#'       in \eqn{W} and \eqn{\lambda_m=(\lambda_{m1},...,\lambda_{md})}.}
#'     \item{If lambdas are \strong{constrained}:}{a size \eqn{(d^2 - n_zeros + r)} vector
#'       \eqn{(Wvec(W),\gamma)}, where \eqn{C_{\lambda}\gamma =(\lambda_2,....,\lambda_M)}, \eqn{\gamma}
#'       is of the size \eqn{(r x 1)}, and \eqn{C_{\lambda}} of the size \eqn{(d*(M - 1) x r}).}
#'   }
#' @param accuracy a positive real number adjusting how close to the given covariance matrix
#'   the returned individual should be.
#'
#'   For \strong{reduced form models} standard deviation of each diagonal element is for reduced form
#'   models
#'   \itemize{
#'     \item \eqn{\omega_{i,i}/}\code{accuracy} when \code{accuracy > d/2}
#'     \item and \code{sqrt(2/d)*}\eqn{\omega_{i,i}} when \code{accuracy <= d/2}.
#'   }
#'   Wishart distribution is used for reduced form models, but for more details read the source code.
#'
#'   For \strong{structural models}, the parameters are generated from normal distribution with mean given
#'   by the argument \code{W_and_lambdas} and the standard deviation is \code{sqrt(abs(W_and_lambdas)/(d + accuracy))}.
#' @inherit random_covmat return
#' @keywords internal

smart_covmat <- function(d, M, Omega, W_and_lambdas, accuracy, structural_pars=NULL) {
  if(is.null(structural_pars)) {
    if(accuracy <= d/2) {
      covmat <- rWishart(n=1, df=d, Sigma=Omega/d)[, , 1]
    } else {
      covmat <- (1 - d/(2*accuracy))*Omega + rWishart(n=1, df=(d^2)/2, Sigma=Omega/(d*accuracy))[, , 1]
    }
    return(vech(covmat))
  } else {
    M <- sum(M)
    pars <- rnorm(n=length(W_and_lambdas), mean=W_and_lambdas, sd=sqrt(abs(W_and_lambdas)/(d + accuracy)))
    W_const <- structural_pars$W[structural_pars$W != 0] # Non-zero constrnts/free values, first length(W_const) values in pars are W params
    # We enforce W to satisfy the sign constraints
    pars[1:length(W_const)][W_const > 0 & !is.na(W_const)] <- abs(pars[1:length(W_const)][W_const > 0 & !is.na(W_const)])
    pars[1:length(W_const)][W_const < 0 & !is.na(W_const)] <- -abs(pars[1:length(W_const)][W_const < 0 & !is.na(W_const)])
    if(M > 1 && is.null(structural_pars$fixed_lambdas)) {
      n_lambdas <- ifelse(is.null(structural_pars$C_lambda), d*(M - 1), ncol(structural_pars$C_lambda))
      lambdas <- abs(pars[(length(pars) - n_lambdas + 1):length(pars)]) # Make lambdas positive
      pars[(length(pars) - n_lambdas + 1):length(pars)] <- lambdas
    }
    return(pars)
  }
}



#' @title Create random degrees of freedom parameter values
#'
#' @description \code{random_df} generates random \code{M2} degrees of freedom parameter values, where
#' \code{M2} is number of StMVAR type regimes in the model-
#'
#' @inheritParams loglikelihood_int
#' @return
#'   \describe{
#'     \item{\strong{GMVAR models}:}{a numeric vector of length zero.}
#'     \item{\strong{StMVAR models}:}{a numeric vector of length \code{M} with random entries strictly larger than two.}
#'     \item{\strong{G-StMVAR models}:}{a numeric vector of length \code{M2} with random entries strictly larger than two.}
#'   }
#' @keywords internal

random_df <- function(M, model=c("GMVAR", "StMVAR", "G-StMVAR")) {
  model <- match.arg(model)
  if(model == "GMVAR") {
    return(numeric(0))
  } else if(model == "StMVAR") {
    M2 <- M
  } else { # model == "G-StMVAR"
    M2 <- M[2]
  }
  2.000001 + rgamma(M2, shape=0.3, rate=0.007)
}


#' @title Create random degrees of freedom parameter values close to given values
#'
#' @description \code{random_df} generates random \code{M2} degrees of freedom parameter values
#'   close to given values, where \code{M2} is number of StMVAR type regimes in the model.
#'
#' @inheritParams loglikelihood_int
#' @param df the old degrees of freedom parameters (of all regimes)
#' @param which_random a vector with length between 1 and M specifying the mixture components that should be random instead of
#'   close to the given degrees of freedom.
#' @param accuracy a positive real number adjusting how close to the given degrees of freedom parameters
#'   the returned df should be.
#' @param which_random a vector with length between 1 and M specifying the mixture components that should be random instead of
#'   close to the given degrees of freedom.
#' @inherit random_df return
#' @keywords internal

smart_df <- function(M, df, accuracy, which_random=numeric(0), model=c("GMVAR", "StMVAR", "G-StMVAR")) {
  model <- match.arg(model)
  if(model == "GMVAR") {
    return(numeric(0))
  }
  M1 <- ifelse(model == "StMVAR", 0, M[1])
  new_df <- rnorm(length(df), mean=df, sd=pmax(0.2, abs(df))/accuracy) # smart df
  if(model == "G-StMVAR") which_random <- which_random[which_random > M1]
  if(length(which_random) > 0) { # Some df should be random instead of smart?
    rand_df <- 2 + rgamma(length(which_random), shape=0.3, rate=0.007)
    new_df[which_random - M1] <- rand_df
  }
  pmax(2.01, new_df) # Make sure all df are above the strict lower bound 2
}
