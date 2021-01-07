#' @title Calculate regime means \eqn{\mu_{m}}
#'
#' @description \code{get_regime_means} calculates regime means \eqn{\mu_{m} = (I - \sum A)^(-1))}
#'   from the given parameter vector.
#'
#' @inheritParams loglikelihood_int
#' @inheritParams is_stationary
#' @return Returns a \eqn{(dxM)} matrix containing regime mean \eqn{\mu_{m}} in the m:th column, \eqn{m=1,..,M}.
#' @section Warning:
#'  No argument checks!
#' @inherit is_stationary references

get_regime_means_int <- function(p, M, d, params, parametrization=c("intercept", "mean"),
                                 constraints=NULL, same_means=NULL, structural_pars=NULL) {
  parametrization <- match.arg(parametrization)
  params <- reform_constrained_pars(p=p, M=M, d=d, params=params, constraints=constraints,
                                    same_means=same_means, structural_pars=structural_pars)
  structural_pars <- get_unconstrained_structural_pars(structural_pars=structural_pars)
  if(parametrization == "intercept") {
    params <- change_parametrization(p=p, M=M, d=d, params=params, constraints=NULL, structural_pars=structural_pars, change_to="mean")
  }
  pick_phi0(p=p, M=M, d=d, params=params, structural_pars=structural_pars)
}


#' @title Calculate regime means \eqn{\mu_{m}}
#'
#' @description \code{get_regime_means} calculates regime means \eqn{\mu_{m} = (I - \sum A_{m,i})^(-1))}
#'   for the given GMVAR model.
#'
#' @inheritParams simulateGMVAR
#' @return Returns a \eqn{(dxM)} matrix containing regime mean \eqn{\mu_{m}} in the m:th column, \eqn{m=1,..,M}.
#' @family moment functions
#' @seealso \code{\link{uncond_moments}}, \code{\link{get_regime_autocovs}}, \code{\link{cond_moments}}
#' @inherit is_stationary references
#' @examples
#' # These examples use the data 'eurusd' which comes with the
#' # package, but in a scaled form.
#' data <- cbind(10*eurusd[,1], 100*eurusd[,2])
#' colnames(data) <- colnames(eurusd)
#'
#' # GMVAR(1,2), d=2 model:
#' params122 <- c(0.623, -0.129, 0.959, 0.089, -0.006, 1.006, 1.746,
#'  0.804, 5.804, 3.245, 7.913, 0.952, -0.037, -0.019, 0.943, 6.926,
#'  3.982, 12.135, 0.789)
#' mod122 <- GMVAR(data, p=1, M=2, params=params122)
#' mod122
#' get_regime_means(mod122)
#'
#' # GMVAR(2,2), d=2 model with mean-parametrization:
#' params222 <- c(-11.904, 154.684, 1.314, 0.145, 0.094, 1.292, -0.389,
#'  -0.070, -0.109, -0.281, 0.920, -0.025, 4.839, 11.633, 124.983, 1.248,
#'   0.077, -0.040, 1.266, -0.272, -0.074, 0.034, -0.313, 5.855, 3.570,
#'   9.838, 0.740)
#' mod222 <- GMVAR(data, p=2, M=2, params=params222, parametrization="mean")
#' mod222
#' get_regime_means(mod222)
#'
#' # Structural GMVAR(2, 2), d=2 model identified with sign-constraints:
#' params222s <- c(1.03, 2.36, 1.79, 3, 1.25, 0.06, 0.04, 1.34, -0.29,
#'  -0.08, -0.05, -0.36, 1.2, 0.05, 0.05, 1.3, -0.3, -0.1, -0.05, -0.4,
#'   0.89, 0.72, -0.37, 2.16, 7.16, 1.3, 0.37)
#' W_222 <- matrix(c(1, 1, -1, 1), nrow=2, byrow=FALSE)
#' mod222s <- GMVAR(data, p=2, M=2, params=params222s, structural_pars=list(W=W_222))
#' mod222s
#' get_regime_means(mod222s)
#' @export

get_regime_means <- function(gmvar) {
  check_gmvar(gmvar)
  get_regime_means_int(p=gmvar$model$p, M=gmvar$model$M, d=gmvar$model$d, params=gmvar$params,
                       parametrization=gmvar$model$parametrization, constraints=gmvar$model$constraints,
                       same_means=gmvar$model$same_means, structural_pars=gmvar$model$structural_pars)
}



#' @title Calculate regimewise autocovariance matrices
#'
#' @description \code{get_regime_autocovs_int} calculates the regimewise autocovariance matrices \eqn{\Gamma_{m}(j)}
#'  \eqn{j=0,1,...,p} for the given GMVAR model.
#'
#' @inheritParams loglikelihood_int
#' @inheritParams reform_constrained_pars
#' @return Returns an \eqn{(d x d x p+1 x M)} array containing the first p regimewise autocovariance matrices.
#'   The subset \code{[, , j, m]} contains the j-1:th lag autocovariance matrix of the m:th regime.
#' @inherit loglikelihood_int references

get_regime_autocovs_int <- function(p, M, d, params, constraints=NULL, same_means=NULL, structural_pars=NULL) {

  params <- reform_constrained_pars(p=p, M=M, d=d, params=params, constraints=constraints,
                                    same_means=same_means, structural_pars=structural_pars)
  structural_pars <- get_unconstrained_structural_pars(structural_pars=structural_pars)
  all_A <- pick_allA(p=p, M=M, d=d, params=params, structural_pars=structural_pars)
  all_Omega <- pick_Omegas(p=p, M=M, d=d, params=params, structural_pars=structural_pars)
  all_boldA <- form_boldA(p=p, M=M, d=d, all_A=all_A)
  I_dp2 <- diag(nrow=(d*p)^2)
  ZER_lower <- matrix(0, nrow=d*(p-1), ncol=d*p)
  ZER_right <- matrix(0, nrow=d, ncol=d*(p - 1))
  all_Gammas <- array(NA, dim=c(d, d, p + 1, M)) # For each m=1,..,M, store the (dxd) covariance matrices Gamma_{y,m}(0),...,Gamma{y,m}(p-1),,Gamma{y,m}(p)
  for(m in 1:M) {
    # Calculate the (dpxdp) Gamma_{Y,m}(0) covariance matrix (Lütkepohl 2005, eq. (2.1.39))
    kronmat <- I_dp2 - kronecker(all_boldA[, , m], all_boldA[, , m])
    sigma_epsm <- rbind(cbind(all_Omega[, , m], ZER_right), ZER_lower)
    Gamma_m <- matrix(solve(kronmat, vec(sigma_epsm)), nrow=d*p, ncol=d*p, byrow=FALSE)

    # Obtain the Gamma_{y,m}(0),...,Gamma_{y,m}(p-1) covariance matrices from Gamma_{Y,m}(0)
    all_Gammas[, , , m] <- c(as.vector(Gamma_m[1:d,]), rep(NA, d*d))

    # Calculate the Gamma{y,m}(p) recursively from Gamma_{y,m}(0),...,Gamma_{y,m}(p-1) (Lütkepohl 2005, eq. (2.1.37))
    all_Gammas[, , p + 1, m] <- rowSums(vapply(1:p, function(i1) all_A[, ,i1 , m]%*%all_Gammas[, , p + 1 - i1, m], numeric(d*d)))
  }
  all_Gammas
}


#' @title Calculate regimewise autocovariance matrices
#'
#' @description \code{get_regime_autocovs} calculates the first p regimewise autocovariance
#'  matrices \eqn{\Gamma_{m}(j)} for the given GMVAR model.
#'
#' @inheritParams simulateGMVAR
#' @family moment functions
#' @inherit get_regime_autocovs_int return
#' @inherit loglikelihood_int references
#' @examples
#' # GMVAR(1,2), d=2 model:
#' params122 <- c(0.623, -0.129, 0.959, 0.089, -0.006, 1.006, 1.746,
#'  0.804, 5.804, 3.245, 7.913, 0.952, -0.037, -0.019, 0.943, 6.926,
#'  3.982, 12.135, 0.789)
#' mod122 <- GMVAR(p=1, M=2, d=2, params=params122)
#' get_regime_autocovs(mod122)
#'
#' # GMVAR(2,2), d=2 model with AR-parameters restricted to be
#' # the same for both regimes:
#' C_mat <- rbind(diag(2*2^2), diag(2*2^2))
#' params222c <- c(1.031, 2.356, 1.786, 3.000, 1.250, 0.060, 0.036,
#'  1.335, -0.290, -0.083, -0.047, -0.356, 0.934, -0.152, 5.201, 5.883,
#'  3.560, 9.799, 0.368)
#' mod222c <- GMVAR(p=2, M=2, d=2, params=params222c, constraints=C_mat)
#' get_regime_autocovs(mod222c)
#'
#' # Structural GMVAR(2, 2), d=2 model identified with sign-constraints:
#' params222s <- c(1.03, 2.36, 1.79, 3, 1.25, 0.06, 0.04, 1.34, -0.29,
#'  -0.08, -0.05, -0.36, 1.2, 0.05, 0.05, 1.3, -0.3, -0.1, -0.05, -0.4,
#'   0.89, 0.72, -0.37, 2.16, 7.16, 1.3, 0.37)
#' W_222 <- matrix(c(1, 1, -1, 1), nrow=2, byrow=FALSE)
#' mod222s <- GMVAR(p=2, M=2, d=2, params=params222s, structural_pars=list(W=W_222))
#' get_regime_autocovs(mod222s)
#' @export

get_regime_autocovs <- function(gmvar) {
  check_gmvar(gmvar)
  get_regime_autocovs_int(p=gmvar$model$p, M=gmvar$model$M, d=gmvar$model$d, params=gmvar$params,
                          constraints=gmvar$model$constraints, same_means=gmvar$model$same_means,
                          structural_pars=gmvar$model$structural_pars)
}


#' @title Calculate the unconditional mean, variance, the first p autocovariances, and the first p autocorrelations
#'  of a GMVAR process
#'
#' @description \code{uncond_moments_int} calculates the unconditional mean, variance, the first p autocovariances,
#'  and the first p autocorrelations of the specified GMVAR process.
#'
#' @inheritParams loglikelihood_int
#' @inheritParams reform_constrained_pars
#' @details The unconditional moments are based on the stationary distribution of the process.
#' @return Returns a list with three components:
#'   \describe{
#'     \item{\code{$uncond_mean}}{a length d vector containing the unconditional mean of the process.}
#'     \item{\code{$autocovs}}{an \eqn{(d x d x p+1)} array containing the lag 0,1,...,p autocovariances of
#'       the process. The subset \code{[, , j]} contains the lag \code{j-1} autocovariance matrix (lag zero for the variance).}
#'     \item{\code{$autocors}}{the autocovariance matrices scaled to autocorrelation matrices.}
#'   }
#' @inherit loglikelihood_int references

uncond_moments_int <- function(p, M, d, params, parametrization=c("intercept", "mean"), constraints=NULL,
                               same_means=NULL, structural_pars=NULL) {
  parametrization <- match.arg(parametrization)
  params <- reform_constrained_pars(p=p, M=M, d=d, params=params, constraints=constraints, same_means=same_means,
                                    structural_pars=structural_pars) # Remove any constraints
  structural_pars <- get_unconstrained_structural_pars(structural_pars=structural_pars)
  alphas <- pick_alphas(p=p, M=M, d=d, params=params)
  reg_means <- get_regime_means_int(p=p, M=M, d=d, params=params, parametrization=parametrization, constraints=NULL,
                                    same_means=NULL, structural_pars=structural_pars)
  uncond_mean <- colSums(alphas*t(reg_means))
  tmp <- rowSums(vapply(1:M, function(m) alphas[m]*tcrossprod(reg_means[,m] - uncond_mean), numeric(d*d))) # Vectorized matrix
  reg_autocovs <- get_regime_autocovs_int(p=p, M=M, d=d, params=params, constraints=NULL, structural_pars=structural_pars)
  autocovs <- array(rowSums(vapply(1:M, function(m) alphas[m]*reg_autocovs[, , , m], numeric(d*d*(p + 1)))) + tmp, dim=c(d, d, p + 1))
  ind_vars <- diag(autocovs[, , 1])
  autocors <- array(vapply(1:(p + 1), function(i1) {
    acor_mat <- matrix(NA, nrow=d, ncol=d)
    for(i2 in 1:d) {
      for(i3 in 1:d) {
        acor_mat[i2, i3] <- autocovs[i2, i3, i1]/sqrt(ind_vars[i2]*ind_vars[i3])
      }
    }
    acor_mat
    }, numeric(d*d)), dim=c(d, d, p + 1))

  list(uncond_mean=uncond_mean,
       autocovs=autocovs,
       autocors=autocors)
}


#' @title Calculate the unconditional mean, variance, the first p autocovariances, and the first p autocorrelations
#'  of a GMVAR process
#'
#' @description \code{uncond_moments} calculates the unconditional mean, variance, the first p autocovariances,
#'  and the first p autocorrelations of the given GMVAR process.
#'
#' @inheritParams simulateGMVAR
#' @family moment functions
#' @inherit uncond_moments_int return
#' @inherit uncond_moments_int details references
#' @examples
#' # GMVAR(1,2), d=2 model:
#' params122 <- c(0.623, -0.129, 0.959, 0.089, -0.006, 1.006, 1.746,
#'  0.804, 5.804, 3.245, 7.913, 0.952, -0.037, -0.019, 0.943, 6.926,
#'  3.982, 12.135, 0.789)
#' mod122 <- GMVAR(p=1, M=2, d=2, params=params122)
#' uncond_moments(mod122)
#'
#' # GMVAR(2,2), d=2 model with AR-parameters restricted to be
#' # the same for both regimes:
#' C_mat <- rbind(diag(2*2^2), diag(2*2^2))
#' params222c <- c(1.031, 2.356, 1.786, 3.000, 1.250, 0.060, 0.036,
#'  1.335, -0.290, -0.083, -0.047, -0.356, 0.934, -0.152, 5.201, 5.883,
#'  3.560, 9.799, 0.368)
#' mod222c <- GMVAR(p=2, M=2, d=2, params=params222c, constraints=C_mat)
#' uncond_moments(mod222c)
#'
#' # Structural GMVAR(2, 2), d=2 model identified with sign-constraints:
#' params222s <- c(1.03, 2.36, 1.79, 3, 1.25, 0.06, 0.04, 1.34, -0.29,
#'  -0.08, -0.05, -0.36, 1.2, 0.05, 0.05, 1.3, -0.3, -0.1, -0.05, -0.4,
#'   0.89, 0.72, -0.37, 2.16, 7.16, 1.3, 0.37)
#' W_222 <- matrix(c(1, 1, -1, 1), nrow=2, byrow=FALSE)
#' mod222s <- GMVAR(p=2, M=2, d=2, params=params222s, structural_pars=list(W=W_222))
#' uncond_moments(mod222s)
#' @export

uncond_moments <- function(gmvar) {
  check_gmvar(gmvar)
  uncond_moments_int(p=gmvar$model$p, M=gmvar$model$M, d=gmvar$model$d, params=gmvar$params,
                     parametrization=gmvar$model$parametrization, constraints=gmvar$model$constraints,
                     same_means=gmvar$model$same_means, structural_pars=gmvar$model$structural_pars)
}




#' @title Calculate the initial (lag zero)  dp-dimensional autocovariance matrix of p consecutive
#'  observations of VAR process
#'
#' @description \code{VAR_init_acov} calculate the initial (lag zero) dp-dimensional autocovariance
#'  matrix of a VAR process (that is, for p consecutive observations) by taking advantage of the algorithm
#'  proposed by McElroy (2017).
#'
#' @inheritParams loglikelihood_int
#' @param all_Am HERE PARAMETERS
#' @param Omega_m the dxd error term covarian
#' @details Copyright (2015) Tucker McElroy
#'
#'	This program is free software; you can redistribute it and/or
#'	modify it under the terms of the GNU General Public License
#'	as published by the Free Software Foundation; either version 2
#'	of the License, or (at your option) any later version.
#'
#'	This program is distributed in the hope that it will be useful,
#'	but WITHOUT ANY WARRANTY; without even the implied warranty of
#'	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#'	GNU General Public License for more details.
#'
#'	You should have received a copy of the GNU General Public License
#'	along with this program; if not, write to the Free Software
#'	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#' @references HERE THE REFERENCE

VAR_init_acov <- function(p, d, all_Am, Omega_m) {
  # all_Am = all_A[, , , m]
  # sigma = all_Omega[, , m] = Omega_m

  # The K commutator matrix
  Kcommut <- function(vect, d, n) matrix(t(matrix(vect, nrow=d, ncol=n)), ncol=1)
  Kmat <- apply(diag(d^2), MARGIN=1, FUN=Kcommut, d, d)

  # Step 1: vectorized error term covariance matrix for lag zero
  gamMAvec <- matrix(Omega_m, nrow=d^2, ncol=1)

  # Step 2: error term versus y_t covariance matrix for lag zero
  Amat <- array(0, dim=c(d^2, p + 1, d^2, 2*p + 1))
  Arow <- diag(d^2)
  for(i1 in 1:p) { # TÄTÄ VOI NOPEUTTAA ALUSTAMALLA Arow:n CBINDAAMISEN SIJAAN
    Arow <- cbind(-1*diag(d)%x%all_Am[, , i1], Arow) # Kronecker product
  }
  for(i1 in 1:(p + 1)) {
    Amat[, i1, , i1:(i1 + p)] <- Arow
  }
  newA <- array(Amat[, 1:(p + 1), , 1:p], dim=c(d^2, p + 1, d^2, p))
  for(i1 in 1:(p + 1)) {
    for(i2 in 1:p) {
      newA[, i1, , i2] <- newA[, i1, , i2]%*%Kmat
    }
  }
  Amat <- cbind(matrix(Amat[, , , p + 1], nrow=d^2*(p + 1), ncol=d^2),
                matrix(Amat[, , , (p + 2):(2*p + 1)], nrow=d^2*(p + 1), ncol=d^2*(p)) + matrix(newA[, , , p:1], nrow=d^2*(p + 1), ncol=d^2*(p)))

  Bmat <- array(0, dim=c(d^2, 1, d^2, p + 1))
  Brow <- diag(d^2)
  for(i1 in 1:p) { # TÄTÄ VOI NOPEUTTAA ALUSTAMALLA Brow:n CBINDAAMISEN SIJAAN
    Brow <- cbind(Brow, -1*all_Am[, , i1] %x% diag(d))
  }
  Bmat[, 1, , 1:(1 + p)] <- Brow
  Bmat <- matrix(Bmat[, , , 1], nrow=d^2, ncol=d^2)
  #Bmat <- matrix(Bmat, nrow=d^2, ncol=d^2)
#  Binv <- solve(Bmat)
  gamMix <- solve(Bmat, gamMAvec) #Binv%*%gamMAvec

  # Step 3: we pad Gamma_WX(0) with zeros
  gamMixTemp <- c(gamMix, rep(0, times=p*d^2))

  # Step 4: compute the Gamma_XX(h) autocovariances for lags h=0,...,p
  gamARMA <- solve(Amat, gamMixTemp)
  gamMix <- array(gamMix, dim=c(d, d, 1))
  gamARMA <- array(gamARMA, dim=c(d, d, p + 1))
  gamARMA <- gamARMA[, , 1:p]

  # After obtaining the autocovariances for lags h=0,...,p-1,
  # we construct the dp-dimensional covariance matrix for p consecutive
  # observations of a VAR process.

  # FILL IN

  gamARMA
}
