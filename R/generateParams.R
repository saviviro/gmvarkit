
#' @title Create random mean-parametrized parameter vector of a GMVAR model that may not be stationary
#'
#' @description \code{random_ind} generates random mean-parametrized parameter vector that may not be stationary.
#'
#' @inheritParams GAfit
#' @inheritParams is_stationary
#' @return Returns random mean-parametrized parameter vector that has form \strong{\eqn{\theta}}\eqn{ = }(\strong{\eqn{\upsilon_{1}}},
#'  ...,\strong{\eqn{\upsilon_{M}}}, \eqn{\alpha_{1},...,\alpha_{M-1}}), where:
#'  \itemize{
#'    \item \strong{\eqn{\upsilon_{m}}} \eqn{ = (\mu_{m},}\strong{\eqn{\phi_{m}}}\eqn{,\sigma_{m})}
#'    \item \strong{\eqn{\phi_{m}}}\eqn{ = (vec(A_{m,1}),...,vec(A_{m,1})}
#'    \item and \eqn{\sigma_{m} = vech(\Omega_{m})}, m=1,...,M.
#'  }
#' @inherit in_paramspace references

random_ind <- function(p, M, d, constraints=NULL, mu_scale, mu_scale2, omega_scale) {
  scale_A <- ifelse(is.null(constraints),
                    1 + log(2*mean(c((p - 0.2)^(1.25), d))),
                    1 + (sum(constraints)/(M*d^2))^0.85)
  if(is.null(constraints)) {
    x <- as.vector(vapply(1:M, function(m) c(rnorm(d, mean=mu_scale, sd=mu_scale2),
                                             random_coefmats(d=d, how_many=p, scale=scale_A),
                                             random_covmat(d=d, omega_scale=omega_scale)), numeric(p*d^2 + d + d*(d+1)/2)))
  } else {
    q <- ncol(constraints)
    x <- c(as.vector(replicate(n=M, rnorm(d, mean=mu_scale, sd=mu_scale2))),
           rnorm(q, mean=0, sd=0.5/scale_A), # random psi
           as.vector(replicate(n=M, random_covmat(d=d, omega_scale=omega_scale))))
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
#'   close to the given parameter vector. If constraints are employed, then this does not consider AR coefficients. Default is \code{NULL}.
#' @section Warning:
#'   No argument checks!
#' @inherit random_ind return references

smart_ind <- function(p, M, d, params, constraints=NULL, accuracy=1, which_random=numeric(0), mu_scale, mu_scale2, omega_scale, ar_scale=1) {
    scale_A <- 1 + log(2*mean(c((p - 0.2)^(1.25), d)))
    params_std <- reform_constrained_pars(p=p, M=M, d=d, params=params, constraints=constraints) # parameters in the standard form
    alphas <- pick_alphas(p=p, M=M, d=d, params=params_std)
    all_Omega <- pick_Omegas(p=p, M=M, d=d, params=params_std)

    if(is.null(constraints)) { # If there are no contraints
      all_phi0_A <- pick_all_phi0_A(p=p, M=M, d=d, params=params_std) # or all_mu
      pars <- vapply(1:M, function(m) {
        if(any(which_random == m)) {
          if(runif(1) > 0.5) { # Use algorithm to force stationarity for coefficient matrices
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
    } else { # If linear constraints are employed
      all_phi0 <- pick_phi0(p=p, M=M, d=d, params=params_std) # or all_mu
      q <- ncol(constraints)
      psi <- params[(M*d + 1):(M*d + q)]
      pars_phi0 <- as.vector(vapply(1:M, function(m) {
        if(any(which_random == m)) {
          rnorm(d, mean=mu_scale, sd=mu_scale2)
        } else {
          rnorm(d, mean=all_phi0[,m], sd=abs(all_phi0[,m]/accuracy))
        }
      }, numeric(d)))
      pars_psi <- rnorm(q, mean=psi, sd=pmax(0.2, abs(psi))/accuracy)
      pars_sigma <- as.vector(vapply(1:M, function(m) {
        if(any(which_random == m)) {
          random_covmat(d=d, omega_scale=omega_scale)
        } else {
          smart_covmat(d=d, Omega=all_Omega[, , m], accuracy=accuracy)
        }
      }, numeric(d*(d + 1)/2)))
      pars <- c(pars_phi0, pars_psi, pars_sigma)
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
#'   affects the coefficient matrices but larger \code{ar_scale} seems to result in larger
#'   AR coefficients. Read the cited article by Ansley and Kohn (1986) and the source code
#'   for more information.
#'
#'   The covariance matrices are generated from (scaled) Wishart distribution.
#'
#'   Constrained models are not supported!
#' @inherit random_ind return
#' @references
#'  \itemize{
#'    \item Ansley C.F., Kohn R. 1986. A note on reparameterizing a vector autoregressive
#'       moving average model to enforce stationarity.
#'       \emph{Journal of statistical computation and simulation}, \strong{24}:2, 99-106.
#'    \item Kalliovirta L., Meitz M. and Saikkonen P. 2016. Gaussian mixture vector autoregression.
#'          \emph{Journal of Econometrics}, \strong{192}, 485-498.
#'  }

random_ind2 <- function(p, M, d, mu_scale, mu_scale2, omega_scale, ar_scale=1) {
  x <- as.vector(vapply(1:M, function(m) c(rnorm(d, mean=mu_scale, sd=mu_scale2),
                                           random_coefmats2(p=p, d=d, ar_scale=ar_scale),
                                           random_covmat(d=d, omega_scale=omega_scale)), numeric(p*d^2 + d +d*(d+1)/2)))
  if(M > 1) {
    alphas <- runif(n=M)
    return(c(x, (alphas[order(alphas, decreasing=TRUE, method="radix")]/sum(alphas))[-M]))
  } else {
    return(x)
  }
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


#' @title Create random VAR model error term covariance matrix
#'
#' @description \code{random_covmat} generates random VAR model \eqn{(dxd)} error term covariance matrix \eqn{\Omega}
#'   from (scaled) Wishart distribution.
#'
#' @inheritParams is_stationary
#' @inheritParams GAfit
#' @return Returns \eqn{(d(d+1)/2x1)} vector containing vech-vectorized covariance matrix \eqn{\Omega}.

random_covmat <- function(d, omega_scale) {
  smart_covmat(d=d, Omega=diag(x=omega_scale), accuracy=1)
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
#' @param accuracy a positive real number adjusting how close to the given covariance matrix
#'   the returned individual should be. Standard deviation of each diagonal element is
#'   \itemize{
#'    \item \eqn{\omega_{i,i}/}\code{accuracy} when \code{accuracy > d/2}
#'    \item and \code{sqrt(2/d)*}\eqn{\omega_{i,i}} when \code{accuracy <= d/2}.
#'  }
#'  Wishart distribution is used, but for more details read the source code.
#' @return Returns \eqn{(d(d+1)/2x1)} vector containing vech-vectorized covariance matrix \eqn{\Omega}.

smart_covmat <- function(d, Omega, accuracy) {
  if(accuracy <= d/2) {
    covmat <- rWishart(n=1, df=d, Sigma=Omega/d)[, , 1]
  } else {
    covmat <- (1-d/(2*accuracy))*Omega + rWishart(n=1, df=(d^2)/2, Sigma=Omega/(d*accuracy))[, , 1]
  }
  vech(covmat)
}

