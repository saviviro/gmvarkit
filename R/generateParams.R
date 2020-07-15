
#' @title Create random mean-parametrized parameter vector of a GMVAR model that may not be stationary
#'
#' @description \code{random_ind} generates random mean-parametrized parameter vector that may not be stationary.
#'
#' @inheritParams GAfit
#' @inheritParams is_stationary
#' @return Returns random mean-parametrized parameter vector that has the same form as the argument \code{params}
#'   in the other functions, for instance, in the function \code{loglikelihood}.
#' @inherit in_paramspace references

random_ind <- function(p, M, d, constraints=NULL, mu_scale, mu_scale2, omega_scale, W_scale, lambda_scale, structural_pars=NULL) {
  scale_A <- ifelse(is.null(constraints),
                    1 + log(2*mean(c((p - 0.2)^(1.25), d))),
                    1 + (sum(constraints)/(M*d^2))^0.85)
  if(is.null(constraints)) {
    if(is.null(structural_pars)) {
      x <- as.vector(vapply(1:M, function(m) c(rnorm(d, mean=mu_scale, sd=mu_scale2),
                                               random_coefmats(d=d, how_many=p, scale=scale_A),
                                               random_covmat(d=d, omega_scale=omega_scale)), numeric(p*d^2 + d + d*(d+1)/2)))
    } else {
      x <- c(rnorm(d*M, mean=mu_scale, sd=mu_scale2), as.vector(replicate(n=M, random_coefmats(d=d, how_many=p, scale=scale_A))),
             random_covmat(d=d, M=M, W_scale=W_scale, lambda_scale=lambda_scale, structural_pars=structural_pars))
    }
  } else {
    q <- ncol(constraints)
    psi <- rnorm(q, mean=0, sd=0.5/scale_A) # random psi
    all_phi0 <- rnorm(d*M, mean=mu_scale, sd=mu_scale2) # as.vector(replicate(n=M, rnorm(d, mean=mu_scale, sd=mu_scale2)))
    if(is.null(structural_pars)) {
      x <- c(all_phi0, psi, as.vector(replicate(n=M, random_covmat(d=d, omega_scale=omega_scale))))
    } else {
      x <- c(all_phi0, psi, random_covmat(d=d, M=M, W_scale=W_scale, lambda_scale=lambda_scale, structural_pars=structural_pars))
    }
  }
  if(M > 1) {
    alphas <- runif(n=M)
    return(c(x, (alphas[order(alphas, decreasing=TRUE, method="radix")]/sum(alphas))[-M]))
  } else {
    return(x)
  }
}


#' @title Create random parameter vector of a GMVAR model fairly close to a given
#'   parameter vector
#'
#' @description \code{smart_ind} creates random mean-parametrized parameter vector of a GMVAR model fairly
#'  close to a given parameter vector. The result may not be stationary.
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

smart_ind <- function(p, M, d, params, constraints=NULL, accuracy=1, which_random=numeric(0), mu_scale, mu_scale2, omega_scale,
                      ar_scale=1, W_scale, lambda_scale, structural_pars=NULL) {
    scale_A <- 1 + log(2*mean(c((p - 0.2)^(1.25), d)))
    params_std <- reform_constrained_pars(p=p, M=M, d=d, params=params, constraints=constraints, structural_pars=structural_pars)
    unc_structural_pars <- get_unconstrained_structural_pars(structural_pars=structural_pars)
    alphas <- pick_alphas(p=p, M=M, d=d, params=params_std)
    if(is.null(structural_pars)) {
      all_Omega <- pick_Omegas(p=p, M=M, d=d, params=params_std)
    }
    if(is.null(constraints) && is.null(structural_pars)) { # If there are no AR constraints and a reduced form model is considered
      all_phi0_A <- pick_all_phi0_A(p=p, M=M, d=d, params=params_std) # all_mu if called from GA
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
    } else { # If AR parameters are constrained or a structural model is considered
      all_phi0 <- pick_phi0(p=p, M=M, d=d, params=params_std, structural_pars=unc_structural_pars) # all_mu if called from GA
      phi0_pars <- as.vector(vapply(1:M, function(m) {
        if(any(which_random == m)) {
          rnorm(d, mean=mu_scale, sd=mu_scale2)
        } else {
          rnorm(d, mean=all_phi0[,m], sd=abs(all_phi0[,m]/accuracy))
        }
      }, numeric(d)))

      if(is.null(constraints)) { # Structural model with AR parameters not constrained
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
      } else { # Structural or reduced form model with AR parameters constrained
        q <- ncol(constraints)
        psi <- params[(M*d + 1):(M*d + q)]
        AR_pars <- rnorm(q, mean=psi, sd=pmax(0.2, abs(psi))/accuracy)
      }
      if(is.null(structural_pars)) { # Reduced form model
        covmat_pars <- as.vector(vapply(1:M, function(m) {
          if(any(which_random == m)) {
            random_covmat(d=d, omega_scale=omega_scale)
          } else {
            smart_covmat(d=d, Omega=all_Omega[, , m], accuracy=accuracy)
          }
        }, numeric(d*(d + 1)/2)))
      } else { # Structural model
        if(any(which_random == 1)) {
          # If first regime is random, then W must be random so the lambdas may as well be random too.
          covmat_pars <- random_covmat(d=d, M=M, W_scale=W_scale, lambda_scale=lambda_scale, structural_pars=structural_pars)

        } else { # First regime is smart
          W_pars <- Wvec(pick_W(p=p, M=M, d=d, params=params_std, structural_pars=unc_structural_pars))
          if(M > 1) {
            n_lambs <- ifelse(is.null(structural_pars$C_lambda), d*(M - 1), ncol(structural_pars$C_lambda))
            W_and_lambdas <- c(W_pars, params[(length(params) - (M - 1) - n_lambs + 1):(length(params) - (M - 1))])
          } else {
            W_and_lambdas <- W_pars # No lambdas when M == 1
          }
          covmat_pars <- smart_covmat(d=d, M=M, W_and_lambdas=W_and_lambdas, accuracy=accuracy, structural_pars=structural_pars)
        }
        if(is.null(structural_pars$C_lambda) && M > 1) {
          # If lambdas are not constrained, we can replace smart lambdas of some regimes with random lambdas
          for(m in 2:M) {
            if(any(which_random == m)) {
              covmat_pars[(length(covmat_pars) - d*(M - 1) + d*(m - 2) + 1):(length(covmat_pars) - d*(M - 1) + d*(m - 2) + d)] <- rnorm(n=d, mean=0, sd=lambda_scale[m - 1])
            }
          }
        }
      }
      pars <- c(phi0_pars, AR_pars, covmat_pars)
    }
    if(M > 1) {
      alphas2 <- abs(rnorm(M, mean=alphas, sd=0.2))
      return(c(pars, (alphas2[order(alphas2, decreasing=TRUE, method="radix")]/sum(alphas2))[-M]))
    } else {
      return(pars)
    }
  }


#' @title Create somewhat random parameter vector of a GMVAR model that is always stationary
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
#'    \item Virolainen S. 2020. Structural Gaussian mixture vector autoregressive model. Unpublished working
#'      paper, available as arXiv:2007.04713.
#'  }

random_ind2 <- function(p, M, d, mu_scale, mu_scale2, omega_scale, ar_scale=1, W_scale, lambda_scale, structural_pars=NULL) {
  if(is.null(structural_pars)) {
    x <- as.vector(vapply(1:M, function(m) c(rnorm(d, mean=mu_scale, sd=mu_scale2),
                                             random_coefmats2(p=p, d=d, ar_scale=ar_scale),
                                             random_covmat(d=d, omega_scale=omega_scale)), numeric(p*d^2 + d + d*(d + 1)/2)))
  } else {
   x <- c(rnorm(d*M, mean=mu_scale, sd=mu_scale2), as.vector(replicate(n=M, random_coefmats2(p=p, d=d, ar_scale=ar_scale))),
          random_covmat(d=d, M=M, W_scale=W_scale, lambda_scale=lambda_scale, structural_pars=structural_pars))
  }
  if(M > 1) {
    alphas <- runif(n=M)
    return(c(x, (alphas[order(alphas, decreasing=TRUE, method="radix")]/sum(alphas))[-M]))
  } else {
    return(x)
  }
}



#' @title Create random VAR-model \eqn{(dxd)} coefficient matrices \eqn{A}.
#'
#' @description \code{random_coefmats} generates random VAR model coefficient matrices.
#'
#' @inheritParams is_stationary
#' @param how_many how many \eqn{(dxd)} coefficient matrices \eqn{A} should be drawn?
#' @param scale non-diagonal elements will be drawn from mean zero normal distribution
#'   with \code{sd=0.3/scale} and diagonal elements from one with \code{sd=0.6/scale}.
#'   Larger scale will hence more likely result stationary coefficient matrices, but
#'   will explore smaller area of the parameter space. Can be for example
#'   \code{1 + log(2*mean(c((p-0.2)^(1.25), d)))}.
#' @return Returns \eqn{((how_many*d^2)x1)} vector containing vectorized coefficient
#'  matrices \eqn{(vec(A_{1}),...,vec(A_{how_many}))}. Note that if \code{how_many==p},
#'  then the returned vector equals \strong{\eqn{\phi_{m}}}.

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
#' @param ar_scale a positive real number. Larger values will typically result larger AR coefficients.
#' @details The coefficient matrices are generated using the algorithm proposed by Ansley
#'   and Kohn (1986) which forces stationarity. It's not clear in detail how \code{ar_scale}
#'   affects the coefficient matrices. Read the cited article by Ansley and Kohn (1986) and
#'   the source code for more information.
#' @return Returns \eqn{((pd^2)x1)} vector containing stationary vectorized coefficient
#'  matrices \eqn{(vec(A_{1}),...,vec(A_{p})}.
#' @references
#'  \itemize{
#'    \item Ansley C.F., Kohn R. 1986. A note on reparameterizing a vector autoregressive
#'       moving average model to enforce stationarity.
#'       \emph{Journal of statistical computation and simulation}, \strong{24}:2, 99-106.
#'  }

random_coefmats2 <- function(p, d, ar_scale=1) {
  # First generate matrices P_1,..,P_p with singular values less than one
  stopifnot(ar_scale > 0)
  Id <- diag(nrow=d)
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
  all_A <- all_phi[, , p, 1:p]
  as.vector(all_A)
}


#' @title Create random VAR model error term covariance matrix
#'
#' @description \code{random_covmat} generates random VAR model \eqn{(dxd)} error term covariance matrix \eqn{\Omega}
#'   from (scaled) Wishart distribution for reduced form models and the parameters \eqn{W},\eqn{\lambda_1,...,\lambda_M}
#'   for structural models (from normal distributions).
#'
#' @inheritParams is_stationary
#' @inheritParams GAfit
#' @return
#'   \describe{
#'     \item{For \strong{reduced form models}:}{Returns a \eqn{(d(d+1)/2x1)} vector containing vech-vectorized covariance matrix
#'       \eqn{\Omega}.}
#'     \item{For \strong{structural models}:}{Returns a length \eqn{d^2 - n_zeros - d*(M - 1)} vector of the form
#'       \eqn{(Wvec(W),\lambda_2,...,\lambda_M)} where \eqn{\lambda_m=(\lambda_{m1},...,\lambda_{md})}
#'       contains the eigenvalue parameters of the \eqn{m}th regime \eqn{(m>1)} and \eqn{n_zeros} is the number of zero constraints
#'       in \eqn{W}. If lambdas are constrained, replacce \eqn{d*(M - 1)} in the length with \eqn{r} and
#'       \eqn{\lambda_2,...,\lambda_M)} with \strong{\eqn{\gamma}}. The operator \eqn{Wvec()} vectorizes a matrix and removes zeros.}
#'   }

random_covmat <- function(d, M, omega_scale, W_scale, lambda_scale, structural_pars=NULL) {
  if(is.null(structural_pars)) {
    return(smart_covmat(d=d, Omega=diag(x=omega_scale), accuracy=1))
  } else {
    W <- structural_pars$W
    n_zeros <- vapply(1:d, function(i1) sum(W[i1,] == 0, na.rm=TRUE), numeric(1))
    std_devs <- sqrt(W_scale/(d - n_zeros))
    new_W <- matrix(rnorm(n=d^2, mean=0, sd=std_devs), nrow=d, byrow=FALSE) # The standard deviations are recycled
    new_W[W == 0 & !is.na(W)] <- 0
    new_W[W > 0 & !is.na(W)] <- abs(new_W[W > 0 & !is.na(W)])
    new_W[W < 0 & !is.na(W)] <- -abs(new_W[W < 0 & !is.na(W)])
    if(M > 1) {
      if(is.null(structural_pars$C_lambda)) {
        lambdas <- as.vector(abs(vapply(1:(M - 1), function(i1) rnorm(n=d, mean=0, sd=lambda_scale[i1]), numeric(d))))
      } else {
        lambdas <- abs(rnorm(n=ncol(structural_pars$C_lambda), mean=0, sd=lambda_scale)) # gammas
      }
      lambdas[lambdas == Inf | lambdas == -Inf] <- 1 # If the df is very close to zero, Inf values may appear
      # lambdas <- sort(lambdas, decreasing=TRUE)

      return(c(Wvec(new_W), lambdas))
    } else {
      return(Wvec(new_W))
    }
  }
}


#' @title Create random VAR-model \eqn{(dxd)} error term covariance matrix \eqn{\Omega}
#'   fairly close to a given \strong{positive definite} covariance matrix using (scaled)
#'   Wishart distribution
#'
#' @description \code{random_covmat} generates random VAR model \eqn{(dxd)} error term covariance matrix \eqn{\Omega}
#'   from (scaled) Wishart distribution that is fairly close to the given matrix.
#'
#' @inheritParams is_stationary
#' @param Omega a symmetric positive definite \eqn{(dxd)} covariance matrix specifying
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

smart_covmat <- function(d, M, Omega, W_and_lambdas, accuracy, structural_pars=NULL) {
  if(is.null(structural_pars)) {
    if(accuracy <= d/2) {
      covmat <- rWishart(n=1, df=d, Sigma=Omega/d)[, , 1]
    } else {
      covmat <- (1 - d/(2*accuracy))*Omega + rWishart(n=1, df=(d^2)/2, Sigma=Omega/(d*accuracy))[, , 1]
    }
    return(vech(covmat))
  } else {
    pars <- rnorm(n=length(W_and_lambdas), mean=W_and_lambdas, sd=sqrt(abs(W_and_lambdas)/(d + accuracy)))
    W_const <- structural_pars$W
    pars[1:(d^2)][W_const > 0 & !is.na(W_const)] <- abs(pars[1:(d^2)][W_const > 0 & !is.na(W_const)]) # We enforce W to satisfy the sign constraints
    pars[1:(d^2)][W_const < 0 & !is.na(W_const)] <- -abs(pars[1:(d^2)][W_const < 0 & !is.na(W_const)])
    if(M > 1) {
      n_lambdas <- ifelse(is.null(structural_pars$C_lambda), d*(M - 1), ncol(structural_pars$C_lambda))
      lambdas <- abs(pars[(length(pars) - n_lambdas + 1):length(pars)]) # Make lambdas positive
      # lambdas <- sort(lambdas, decreasing=TRUE)
      pars[(length(pars) - n_lambdas + 1):length(pars)] <- lambdas
    }
    return(pars)
  }
}
