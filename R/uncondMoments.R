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
#' @keywords internal

get_regime_means_int <- function(p, M, d, params, model=c("GMVAR", "StMVAR", "G-StMVAR"), parametrization=c("intercept", "mean"),
                                 constraints=NULL, same_means=NULL, weight_constraints=NULL, structural_pars=NULL) {
  parametrization <- match.arg(parametrization)
  model <- match.arg(model)
  params <- reform_constrained_pars(p=p, M=M, d=d, params=params, model=model, constraints=constraints,
                                    same_means=same_means, weight_constraints=weight_constraints,
                                    structural_pars=structural_pars)
  structural_pars <- get_unconstrained_structural_pars(structural_pars=structural_pars)
  if(parametrization == "intercept") {
    params <- change_parametrization(p=p, M=M, d=d, params=params, model=model, constraints=NULL,
                                     weight_constraints=NULL, structural_pars=structural_pars, change_to="mean")
  }
  pick_phi0(p=p, M=M, d=d, params=params, structural_pars=structural_pars)
}


#' @title Calculate regime means \eqn{\mu_{m}}
#'
#' @description \code{get_regime_means} calculates regime means \eqn{\mu_{m} = (I - \sum A_{m,i})^(-1))}
#'   for the given GMVAR, StMVAR, or G-StMVAR model.
#'
#' @inheritParams quantile_residual_tests
#' @return Returns a \eqn{(dxM)} matrix containing regime mean \eqn{\mu_{m}} in the m:th column, \eqn{m=1,..,M}.
#' @family moment functions
#' @seealso \code{\link{uncond_moments}}, \code{\link{get_regime_autocovs}}, \code{\link{cond_moments}}
#' @inherit is_stationary references
#' @examples
#' # GMVAR(1,2), d=2 model:
#' params12 <- c(0.55, 0.112, 0.344, 0.055, -0.009, 0.718, 0.319, 0.005,
#'  0.03, 0.619, 0.173, 0.255, 0.017, -0.136, 0.858, 1.185, -0.012,
#'  0.136, 0.674)
#' mod12 <- GSMVAR(gdpdef, p=1, M=2, params=params12)
#' mod12
#' get_regime_means(mod12)
#'
#' # Structural GMVAR(2, 2), d=2 model identified with sign-constraints:
#' params22s <- c(0.36, 0.121, 0.484, 0.072, 0.223, 0.059, -0.151, 0.395,
#'  0.406, -0.005, 0.083, 0.299, 0.218, 0.02, -0.119, 0.722, 0.093, 0.032,
#'  0.044, 0.191, 0.057, 0.172, -0.46, 0.016, 3.518, 5.154, 0.58)
#' W_22 <- matrix(c(1, 1, -1, 1), nrow=2, byrow=FALSE)
#' mod22s <- GSMVAR(gdpdef, p=2, M=2, params=params22s, structural_pars=list(W=W_22))
#' mod22s
#' get_regime_means(mod22s)
#' @export

get_regime_means <- function(gsmvar) {
  gsmvar <- gmvar_to_gsmvar(gsmvar) # Backward compatibility
  check_gsmvar(gsmvar)
  get_regime_means_int(p=gsmvar$model$p, M=gsmvar$model$M, d=gsmvar$model$d, params=gsmvar$params,
                       model=gsmvar$model$model, parametrization=gsmvar$model$parametrization,
                       constraints=gsmvar$model$constraints, same_means=gsmvar$model$same_means,
                       weight_constraints=gsmvar$model$weight_constraints,
                       structural_pars=gsmvar$model$structural_pars)
}



#' @title Calculate regimewise autocovariance matrices
#'
#' @description \code{get_regime_autocovs_int} calculates the regimewise autocovariance matrices \eqn{\Gamma_{m}(j)}
#'  \eqn{j=0,1,...,p} for the given GSMVAR model.
#'
#' @inheritParams loglikelihood_int
#' @inheritParams reform_constrained_pars
#' @return Returns an \eqn{(d x d x p+1 x M)} array containing the first p regimewise autocovariance matrices.
#'   The subset \code{[, , j, m]} contains the j-1:th lag autocovariance matrix of the m:th regime.
#' @inherit loglikelihood_int references
#' @keywords internal

get_regime_autocovs_int <- function(p, M, d, params, model=c("GMVAR", "StMVAR", "G-StMVAR"), constraints=NULL, same_means=NULL,
                                    weight_constraints=NULL, structural_pars=NULL) {
  model <- match.arg(model)
  params <- reform_constrained_pars(p=p, M=M, d=d, params=params, model=model, constraints=constraints,
                                    same_means=same_means, weight_constraints=weight_constraints, structural_pars=structural_pars)
  structural_pars <- get_unconstrained_structural_pars(structural_pars=structural_pars)
  all_A <- pick_allA(p=p, M=M, d=d, params=params, structural_pars=structural_pars)
  all_Omega <- pick_Omegas(p=p, M=M, d=d, params=params, structural_pars=structural_pars)
  all_boldA <- form_boldA(p=p, M=M, d=d, all_A=all_A)
  M <- sum(M)
  I_dp2 <- diag(nrow=(d*p)^2)
  ZER_lower <- matrix(0, nrow=d*(p-1), ncol=d*p)
  ZER_right <- matrix(0, nrow=d, ncol=d*(p - 1))
  # For each m=1,..,M, store the (dxd) covariance matrices Gamma_{y,m}(0),...,Gamma{y,m}(p-1),,Gamma{y,m}(p):
  all_Gammas <- array(NA, dim=c(d, d, p + 1, M))
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
#'  matrices \eqn{\Gamma_{m}(j)} for the given GMVAR, StMVAR, or G-StMVAR model.
#'
#' @inheritParams quantile_residual_tests
#' @family moment functions
#' @inherit get_regime_autocovs_int return
#' @inherit loglikelihood_int references
#' @examples
#' # GMVAR(1,2), d=2 model:
#' params12 <- c(0.55, 0.112, 0.344, 0.055, -0.009, 0.718, 0.319, 0.005,
#'  0.03, 0.619, 0.173, 0.255, 0.017, -0.136, 0.858, 1.185, -0.012,
#'  0.136, 0.674)
#' mod12 <- GSMVAR(gdpdef, p=1, M=2, params=params12)
#' get_regime_autocovs(mod12)
#'
#' # Structural GMVAR(2, 2), d=2 model identified with sign-constraints:
#' params22s <- c(0.36, 0.121, 0.484, 0.072, 0.223, 0.059, -0.151, 0.395,
#'  0.406, -0.005, 0.083, 0.299, 0.218, 0.02, -0.119, 0.722, 0.093, 0.032,
#'  0.044, 0.191, 0.057, 0.172, -0.46, 0.016, 3.518, 5.154, 0.58)
#' W_22 <- matrix(c(1, 1, -1, 1), nrow=2, byrow=FALSE)
#' mod22s <- GSMVAR(gdpdef, p=2, M=2, params=params22s, structural_pars=list(W=W_22))
#' mod22s
#' get_regime_autocovs(mod22s)
#' @export

get_regime_autocovs <- function(gsmvar) {
  gsmvar <- gmvar_to_gsmvar(gsmvar) # Backward compatibility
  check_gsmvar(gsmvar)
  get_regime_autocovs_int(p=gsmvar$model$p, M=gsmvar$model$M, d=gsmvar$model$d, params=gsmvar$params,
                          model=gsmvar$model$model, constraints=gsmvar$model$constraints,
                          same_means=gsmvar$model$same_means, weight_constraints=gsmvar$model$weight_constraints,
                          structural_pars=gsmvar$model$structural_pars)
}


#' @title Calculate the unconditional mean, variance, the first p autocovariances, and the first p autocorrelations
#'  of a GMVAR, StMVAR, or G-StMVAR process
#'
#' @description \code{uncond_moments_int} calculates the unconditional mean, variance, the first p autocovariances,
#'  and the first p autocorrelations of the specified GMVAR, StMVAR, or G-StMVAR process.
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
#' @keywords internal

uncond_moments_int <- function(p, M, d, params, model=c("GMVAR", "StMVAR", "G-StMVAR"), parametrization=c("intercept", "mean"),
                               constraints=NULL, same_means=NULL, weight_constraints=NULL, structural_pars=NULL) {
  parametrization <- match.arg(parametrization)
  model <- match.arg(model)
  params <- reform_constrained_pars(p=p, M=M, d=d, params=params, model=model, constraints=constraints, same_means=same_means,
                                    weight_constraints=weight_constraints, structural_pars=structural_pars) # Remove any constraints
  structural_pars <- get_unconstrained_structural_pars(structural_pars=structural_pars)
  alphas <- pick_alphas(p=p, M=M, d=d, params=params, model=model)
  reg_means <- get_regime_means_int(p=p, M=M, d=d, params=params, model=model, parametrization=parametrization,
                                    constraints=NULL, same_means=NULL, weight_constraints=weight_constraints,
                                    structural_pars=structural_pars)
  uncond_mean <- colSums(alphas*t(reg_means))
  tmp <- rowSums(vapply(1:sum(M), function(m) alphas[m]*tcrossprod(reg_means[,m] - uncond_mean), numeric(d*d))) # Vectorized matrix
  reg_autocovs <- get_regime_autocovs_int(p=p, M=M, d=d, params=params, model=model, constraints=NULL, structural_pars=structural_pars)
  autocovs <- array(rowSums(vapply(1:sum(M), function(m) alphas[m]*reg_autocovs[, , , m], numeric(d*d*(p + 1)))) + tmp, dim=c(d, d, p + 1))
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
#'  of a GMVAR, StMVAR, or G-StMVAR process
#'
#' @description \code{uncond_moments} calculates the unconditional mean, variance, the first p autocovariances,
#'  and the first p autocorrelations of the given GMVAR, StMVAR, or G-StMVAR process.
#'
#' @inheritParams quantile_residual_tests
#' @family moment functions
#' @inherit uncond_moments_int return
#' @inherit uncond_moments_int details references
#' @examples
#' # GMVAR(1,2), d=2 model:
#' params12 <- c(0.55, 0.112, 0.344, 0.055, -0.009, 0.718, 0.319, 0.005,
#'  0.03, 0.619, 0.173, 0.255, 0.017, -0.136, 0.858, 1.185, -0.012,
#'  0.136, 0.674)
#' mod12 <- GSMVAR(gdpdef, p=1, M=2, params=params12)
#' uncond_moments(mod12)
#'
#' # Structural GMVAR(2, 2), d=2 model identified with sign-constraints:
#' params22s <- c(0.36, 0.121, 0.484, 0.072, 0.223, 0.059, -0.151, 0.395,
#'  0.406, -0.005, 0.083, 0.299, 0.218, 0.02, -0.119, 0.722, 0.093, 0.032,
#'  0.044, 0.191, 0.057, 0.172, -0.46, 0.016, 3.518, 5.154, 0.58)
#' W_22 <- matrix(c(1, 1, -1, 1), nrow=2, byrow=FALSE)
#' mod22s <- GSMVAR(gdpdef, p=2, M=2, params=params22s, structural_pars=list(W=W_22))
#' mod22s
#' uncond_moments(mod22s)
#' @export

uncond_moments <- function(gsmvar) {
  gsmvar <- gmvar_to_gsmvar(gsmvar) # Backward compatibility
  check_gsmvar(gsmvar)
  uncond_moments_int(p=gsmvar$model$p, M=gsmvar$model$M, d=gsmvar$model$d, params=gsmvar$params,
                     model=gsmvar$model$model, parametrization=gsmvar$model$parametrization,
                     constraints=gsmvar$model$constraints, same_means=gsmvar$model$same_means,
                     weight_constraints=gsmvar$model$weight_constraints,
                     structural_pars=gsmvar$model$structural_pars)
}


#' @title Calculate the dp-dimensional covariance matrix of p consecutive
#'  observations of a VAR process
#'
#' @description \code{VAR_pcovmat} calculate the dp-dimensional covariance matrix of p consecutive
#'  observations of a VAR process with the algorithm proposed by McElroy (2017).
#'
#' @inheritParams loglikelihood_int
#' @inheritParams is_stationary
#' @param all_Am \code{[d, d, p]} array containing the AR coefficient matrices
#' @param Omega_m the dxd error term covariance matrix
#' @details
#'  Most of the code in this function is adapted from the one provided in the
#'  supplementary material of McElroy (2017). Reproduced under GNU General
#'  Public License, Copyright (2015) Tucker McElroy.
#' @return Returns the (dp x dp) covariance matrix.
#' @references
#'  \itemize{
#'    \item McElroy T. 2017. Computation of vector ARMA autocovariances.
#'          \emph{Statistics and Probability Letters}, \strong{124}, 92-96.
#'  }
#' @keywords internal

VAR_pcovmat <- function(p, d, all_Am, Omega_m) {
  # all_Am = all_A[, , , m]
  # Omega_m = all_Omega[, , m]

  # The K commutation matrix
  Kcommut <- function(vect) matrix(t(matrix(vect, nrow=d, ncol=d)), ncol=1)
  Kmat <- apply(diag(d^2), MARGIN=1, FUN=Kcommut)

  # Step 1: vectorized error term covariance matrix for lag zero
  gamMAvec <- matrix(Omega_m, nrow=d^2, ncol=1)

  # Step 2: error term versus y_t covariance matrix for lag zero
  Amat <- array(0, dim=c(d^2, p + 1, d^2, 2*p + 1))
  Arow <-  matrix(nrow=d^2, ncol=d^2*(p + 1))
  start_inds <- seq(from=ncol(Arow) - d^2 + 1, to=1, by=-d^2)
  end_inds <- seq(from=ncol(Arow), to=d^2, by=-d^2)
  Arow[,start_inds[1]:end_inds[1]] <- diag(d^2)
  for(i1 in 1:p) {
    Arow[,start_inds[i1 + 1]:end_inds[i1 + 1]] <- -1*diag(d)%x%all_Am[, , i1]
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
                matrix(Amat[, , , (p + 2):(2*p + 1)], nrow=d^2*(p + 1), ncol=d^2*(p)) + matrix(newA[, , , p:1], nrow=d^2*(p + 1), ncol=d^2*p))

  Bmat <- array(0, dim=c(d^2, 1, d^2, p + 1))
  Brow <-  matrix(nrow=d^2, ncol=d^2*(p + 1))
  start_inds <- rev(start_inds)
  end_inds <- rev(end_inds)
  Brow[,start_inds[1]:end_inds[1]] <- diag(d^2)
  for(i1 in 1:p) {
    Brow[,start_inds[i1 + 1]:end_inds[i1 + 1]] <- -1*diag(d)%x%all_Am[, , i1]
  }
  Bmat[, 1, , 1:(1 + p)] <- Brow
  gamMix <- solve(matrix(Bmat[, , , 1], nrow=d^2, ncol=d^2), gamMAvec)

  # Step 3: we pad out Gamma_WX(0) with zeros
  gamMixTemp <- c(gamMix, rep(0, times=p*d^2))

  # Step 4: compute the Gamma_XX(h) autocovariances for lags h=0,...,p
  gamVAR <- array(array(solve(Amat, gamMixTemp), dim=c(d, d, p + 1))[, , 1:p], dim=c(d, d, p))

  # After obtaining the autocovariances for lags h=0,...,p-1,
  # we construct the dp-dimensional covariance matrix for p consecutive
  # observations of a VAR process.

  # gamVAR contains lag=0 autocovariance in [, , 1], and lag=i in [, , i + 1].
  # Moreover, we use Gamma_Y(-h) = t(Gamma_Y(h)) and store the transposes
  # (as taking transpose multiple times uses more computation time):
  tgamVAR <- array(dim=c(d, d, p))
  for(i1 in 1:p) {
    tgamVAR[, , i1] <- t(gamVAR[, , i1])
  }

  # Finally, we fill in the covariance matrix for p consecutive observations
  # of the VAR process:
  Sigma_m <- matrix(nrow=d*p, ncol=d*p)
  start_inds <- seq(from=1, to=d*(p - 1) + 1, by=d)
  end_inds <- seq(from=d, to=d*p, by=d)
  for(i1 in 1:p) { # Go through row blocks
    for(i2 in 1:p) { # Go through column blocks
      # If i1 end_ind is larger than i2 end_ind, we consider blocks below block
      # diagonal (use transpose), and if i1 end_ind is smaller than i2 end_ind,
      # we consider blocks above block diagonal. If i1 end_ind == i2 end_ind,
      # we are at the block diagonal.
      if(end_inds[i1] > end_inds[i2]) { # Below diagonal, use transpose
        Sigma_m[start_inds[i1]:end_inds[i1], start_inds[i2]:end_inds[i2]] <- tgamVAR[, , abs(i1 - i2) + 1]
      } else {
        Sigma_m[start_inds[i1]:end_inds[i1], start_inds[i2]:end_inds[i2]] <- gamVAR[, , abs(i1 - i2) + 1]
      }
    }
  }
  Sigma_m
}


#' @title Calculate the dp-dimensional covariance matrices \eqn{\Sigma_{m,p}} in the mixing weights
#'  of the GMVAR, StMVAR, or G-StMVAR model.
#'
#' @description \code{get_Sigmas} calculates the dp-dimensional covariance matrices \eqn{\Sigma_{m,p}}
#'  in the mixing weights of the GMVAR, StMVAR, or G-StMVAR model so that the algorithm proposed by McElroy (2017) employed
#'  whenever it reduces the computation time.
#'
#' @inheritParams is_stationary
#' @inheritParams form_boldA
#' @param all_Omega a \code{[d, d, M]} array containing the covariance matrix Omegas
#' @details
#'  Calculates the dp-dimensional covariance matrix using the formula (2.1.39) in Lütkepohl (2005) when
#'  \code{d*p < 12} and using the algorithm proposed by McElroy (2017) otherwise.
#'
#'  The code in the implementation of the McElroy's (2017) algorithm (in the function \code{VAR_pcovmat}) is
#'  adapted from the one provided in the supplementary material of McElroy (2017). Reproduced under GNU General
#'  Public License, Copyright (2015) Tucker McElroy.
#' @return Returns a \code{[dp, dp, M]} array containing the dp-dimensional covariance matrices for each regime.
#' @references
#'  \itemize{
#'    \item Kalliovirta L., Meitz M. and Saikkonen P. 2016. Gaussian mixture vector autoregression.
#'            \emph{Journal of Econometrics}, \strong{192}, 485-498.
#'    \item Lütkepohl H. 2005. New Introduction to Multiple Time Series Analysis,
#'            \emph{Springer}.
#'    \item McElroy T. 2017. Computation of vector ARMA autocovariances.
#'          \emph{Statistics and Probability Letters}, \strong{124}, 92-96.
#'    \item Virolainen S. 2025. A statistically identified structural vector autoregression with endogenously
#'           switching volatility regime. \emph{Journal of Business & Economic Statistics}. \strong{43}:1, 44-54.
#'    \item Virolainen S. in press. A Gaussian and Student’s mixture vector autoregressive model with an application
#'     to monetary policy shocks. Econometrics and Statistics.
#'  }
#' @keywords internal

get_Sigmas <- function(p, M, d, all_A, all_boldA, all_Omega) {
  M <- sum(M)
  Sigmas <- array(NA, dim=c(d*p, d*p, M)) # Store the (dpxdp) covariance matrices
  if(d*p < 12) { # d*p < 12
    # Calculate the covariance matrices Sigma_{m,p} using the equation (2.1.39) in Lütkepohl (2005).
    I_dp2 <- diag(nrow=(d*p)^2)
    ZER_lower <- matrix(0, nrow=d*(p - 1), ncol=d*p)
    ZER_right <- matrix(0, nrow=d, ncol=d*(p - 1))
    for(m in 1:M) {
      kronmat <- I_dp2 - kronecker(all_boldA[, , m], all_boldA[, , m])
      sigma_epsm <- rbind(cbind(all_Omega[, , m], ZER_right), ZER_lower)
      Sigma_m <- solve(kronmat, vec(sigma_epsm))
      Sigmas[, , m] <- Sigma_m
    }
 } else { # d*p >= 12
   # Calculate the covariance matrices Sigma_{m,p} using the algorithm proposed my McElroy (2017).
   for(m in 1:M) {
     Sigmas[, , m] <- VAR_pcovmat(p=p, d=d, all_Am=all_A[, , , m], Omega_m=all_Omega[, , m])
   }
 }
  Sigmas
}
