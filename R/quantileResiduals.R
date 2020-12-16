
#' @title Calculate multivariate quantile residuals of a GMVAR model
#'
#' @description \code{quantile_residuals} calculates multivariate quantile residuals
#'  (described by \emph{Kalliovirta and Saikkonen 2010}) for a GMVAR model.
#'
#' @inheritParams simulateGMVAR
#' @return Returns \eqn{((n_obs-p) x d)} matrix containing the multivariate quantile residuals,
#'   \eqn{j}:th column corresponds to the time series in the \eqn{j}:th column of the data. The multivariate
#'   quantile residuals are calculated so that the first column quantile residuals are the "unconditioned ones"
#'   and the rest condition on all the previous ones in numerical order. Read the cited article by
#'   \emph{Kalliovirta and Saikkonen 2010} for details.
#' @inherit GMVAR references
#' @seealso \code{\link{fitGMVAR}}, \code{\link{GMVAR}}, \code{\link{quantile_residual_tests}},
#'   \code{\link{diagnostic_plot}}, \code{\link{predict.gmvar}}, \code{\link{profile_logliks}}
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
#' quantile_residuals(mod122)
#'
#' # GMVAR(2,2), d=2 model with mean-parametrization:
#' params222 <- c(-11.904, 154.684, 1.314, 0.145, 0.094, 1.292, -0.389,
#'  -0.070, -0.109, -0.281, 0.920, -0.025, 4.839, 11.633, 124.983, 1.248,
#'   0.077, -0.040, 1.266, -0.272, -0.074, 0.034, -0.313, 5.855, 3.570,
#'   9.838, 0.740)
#' mod222 <- GMVAR(data, p=2, M=2, params=params222, parametrization="mean")
#' quantile_residuals(mod222)
#'
#' # Structural GMVAR(2, 2), d=2 model identified with sign-constraints:
#' params222s <- c(1.03, 2.36, 1.79, 3, 1.25, 0.06, 0.04, 1.34, -0.29,
#'  -0.08, -0.05, -0.36, 1.2, 0.05, 0.05, 1.3, -0.3, -0.1, -0.05, -0.4,
#'   0.89, 0.72, -0.37, 2.16, 7.16, 1.3, 0.37)
#' W_222 <- matrix(c(1, 1, -1, 1), nrow=2, byrow=FALSE)
#' mod222s <- GMVAR(data, p=2, M=2, params=params222s, structural_pars=list(W=W_222))
#' quantile_residuals(mod222s)
#'
#' # GMVAR(2,2), d=2 model with AR-parameters restricted to be
#' # the same for both regimes:
#' C_mat <- rbind(diag(2*2^2), diag(2*2^2))
#' params222c <- c(1.031, 2.356, 1.786, 3.000, 1.250, 0.060, 0.036,
#'  1.335, -0.290, -0.083, -0.047, -0.356, 0.934, -0.152, 5.201, 5.883,
#'  3.560, 9.799, 0.368)
#' mod222c <- GMVAR(data, p=2, M=2, params=params222c, constraints=C_mat)
#' quantile_residuals(mod222c)
#'
#' # GMVAR(2,2), d=2 model with AR-parameters restricted to be
#' # the same for both regimes and the non-diagonal elements of
#' # the coefficient matrices constrained to zero.
#' tmp <- matrix(c(1, rep(0, 10), 1, rep(0, 8), 1, rep(0, 10), 1),
#'  nrow=2*2^2, byrow=FALSE)
#' C_mat2 <- rbind(tmp, tmp)
#' params222c2 <- c(0.355, 3.193, -0.114, 2.829, 1.263, 1.338, -0.292,
#'  -0.362, 5.597, 3.456, 9.622, 0.982, -0.327, 5.236, 0.650)
#' mod222c2 <- GMVAR(data, p=2, M=2, params=params222c2,
#'   constraints=C_mat2)
#' quantile_residuals(mod222c)
#' @export

quantile_residuals <- function(gmvar) {

  # Collect preliminary values
  check_gmvar(gmvar)
  epsilon <- round(log(.Machine$double.xmin) + 10)
  p <- gmvar$model$p
  M <- gmvar$model$M
  d <- gmvar$model$d
  data <- gmvar$data
  constraints <- gmvar$model$constraints
  structural_pars <- gmvar$model$structural_pars
  n_obs <- nrow(data)
  T_obs <- n_obs - p

  # Collect parameter values
  params <- gmvar$params
  if(gmvar$model$parametrization == "mean") {
    params <- change_parametrization(p=p, M=M, d=d, params=params, constraints=constraints,
                                     structural_pars=structural_pars, change_to="intercept")
  }
  params <- reform_constrained_pars(p=p, M=M, d=d, params=params, constraints=constraints,
                                    structural_pars=structural_pars)
  structural_pars <- get_unconstrained_structural_pars(structural_pars=structural_pars)

  all_mu <- get_regime_means(gmvar)
  all_phi0 <- pick_phi0(p=p, M=M, d=d, params=params, structural_pars=structural_pars)
  all_A <- pick_allA(p=p, M=M, d=d, params=params, structural_pars=structural_pars)
  all_Omega <- pick_Omegas(p=p, M=M, d=d, params=params, structural_pars=structural_pars)
  all_boldA <- form_boldA(p=p, M=M, d=d, all_A=all_A)
  alphas <- pick_alphas(p=p, M=M, d=d, params=params)

  # An i:th row denotes the vector \bold{y_{i-1}} = (y_{i-1}',...,y_{i-p}') (dpx1),
  # assuming the observed data is y_{-p+1},...,y_0,y_1,...,y_{T}
  Y <- reform_data(data, p)

  # Mixing weights
  alpha_mt <- gmvar$mixing_weights

  # Calculate the conditional means mu_{m,t} (KMS 2016, Condition 1 (a)).
  # Map to dimensions of mu_mt: [t, p, m]
  all_A2 <- array(all_A, dim=c(d, d*p, M)) # cbind coefficient matrices of each component: m:th component is obtained at [, , m]
  Y2 <- Y[1:T_obs,] # Last row is not needed because mu_mt uses lagged values
  mu_mt <- array(vapply(1:M, function(m) t(all_phi0[, m] + tcrossprod(all_A2[, , m], Y2)), numeric(d*T_obs)), dim=c(T_obs, d, M)) # [, , m]

  ## Start computing the multivariate quantile residuals (Kalliovirta and Saikkonen 2010, eq.(4))
  # using properties of marginal and conditional distributions of multinormal random variables and
  # applying them to the mixture components at each time point t.

  # Compute partitions and matrix products of partioned covariance matrices Omega that will be used multiple times.
  upleft_jjmat <- function(mat, j) mat[1:j, 1:j, drop=FALSE]

  # Storage for variances and conditional means (obtained from multinormal)
  variances <- matrix(nrow=d, ncol=M) # Column per mixture component and row per component
  mu_mtj <- array(dim=c(T_obs, d, M)) # [t, j, m]

  # Calculate variances and means (both conditional for j=2,...,d) for t=1,...,T, m=1,...,M, j=1,...,d
  dat <- data[(p + 1):nrow(data),] # Remove the initial values
  for(m in 1:M) {
    variances[1, m] <- all_Omega[1, 1, m] # "Marginal" variance
    mu_dif <- dat - mu_mt[, , m] # mu_mt - y_t, t = 1,...,T
    mu_mtj[, 1, m] <- mu_mt[, 1, m] # "Marginal" means

    for(j in 2:d) {
      Omega_mj <- upleft_jjmat(all_Omega[, , m], j)
      up_left <- upleft_jjmat(Omega_mj, j - 1)
      up_right <- Omega_mj[1:(j - 1), j, drop=FALSE] # (j-1 x 1)
      low_left <- t(up_right)
      low_right <- Omega_mj[j, j]

      matprod <- low_left%*%solve(up_left) # (1 x j-1)
      variances[j, m] <- low_right - matprod%*%up_right # Conditional variance conditioned to components 1,..,j-1
      mu_mtj[, j, m] <- mu_mt[, j, m] + t(tcrossprod(matprod, mu_dif[, 1:(j - 1), drop=FALSE])) #t(matprod%*%t(mu_dif[,1:(j-1)]))
    }
  }

  # Calculate beta_{m,t,j} for j=2,...,d (Virolainen 2018, eq. (1.12) unpublished work paper)
  beta_mtj <- array(dim=c(T_obs, M, d)) # [t, m, j] j=1,...,d
  beta_mtj[, , 1] <- alpha_mt
  for(j in 2:d) {
    log_mvnvalues <- vapply(1:M, function(m) dlogmultinorm(y=dat[,1:(j - 1), drop=FALSE],
                                                           mu=as.matrix(mu_mt[, 1:(j - 1), m]),
                                                           Omega=upleft_jjmat(all_Omega[, , m], j - 1)),
                                             numeric(T_obs))

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

    numerators <- as.matrix(alpha_mt*exp(log_mvnvalues))
    denominator <- rowSums(numerators)
    beta_mtj[, , j] <- numerators/denominator

  }

  # Then calculate (y_{i_j,t} - mu_mtj)/sqrt(variance_mtj) for m=1,...,M, t=1,...,T, j=1,...,d
  points_mtj <- array(vapply(1:M, function(m) t(t(dat - mu_mtj[, , m])/sqrt(variances[, m])),
                      numeric(d*T_obs)), dim=c(T_obs, d, M))

  F_values <- vapply(1:d, function(j) rowSums(as.matrix(beta_mtj[, , j]*pnorm(points_mtj[, j, ]))), numeric(T_obs))

  # Values too close to 0 or 1 will be scaled so that qnorm doesn't return inf-values
  F_values[F_values >= 1 - .Machine$double.eps/2] <- 1 - .Machine$double.eps/2
  F_values[F_values <= .Machine$double.eps/2] <- .Machine$double.eps/2

  qnorm(F_values)
}


#' @title Calculate multivariate quantile residuals of GMVAR model
#'
#' @description \code{quantile_residuals_int} is a wrapper for \code{quantile_residuals} to compute
#'   quantile residuals using parameter values instead of class \code{gmvar} object.
#'
#' @inheritParams loglikelihood_int
#' @section Warning:
#'   No argument checks!
#' @inherit quantile_residuals return references

quantile_residuals_int <- function(data, p, M, params, conditional, parametrization, constraints=NULL,
                                   structural_pars=NULL, stat_tol=1e-3, posdef_tol=1e-8) {
  lok_and_mw <- loglikelihood_int(data=data, p=p, M=M, params=params, conditional=conditional,
                                  parametrization=parametrization, constraints=constraints,
                                  structural_pars=structural_pars, to_return="loglik_and_mw",
                                  check_params=TRUE, minval=NA, stat_tol=stat_tol, posdef_tol=posdef_tol)
  d <- ncol(data)
  npars <- n_params(p=p, M=M, d=d, constraints=constraints)
  mod <- structure(list(data=data,
                        model=list(p=p,
                                   M=M,
                                   d=ncol(data),
                                   conditional=conditional,
                                   parametrization=parametrization,
                                   constraints=constraints,
                                   structural_pars=structural_pars),
                        params=params,
                        std_errors=rep(NA, npars),
                        mixing_weights=lok_and_mw$mw,
                        quantile_residuals=NA,
                        loglik=structure(lok_and_mw$loglik,
                                         class="logLik",
                                         df=npars),
                        IC=NA,
                        all_estimates=NULL,
                        all_logliks=NULL,
                        which_converged=NULL,
                        num_tols=list(stat_tol=stat_tol,
                                      posdef_tol=posdef_tol)),
                   class="gmvar")
  quantile_residuals(mod)
}


#' @title Calculate logarithms of multiple multivariate normal densities with varying
#'  mean and constant covariance matrix
#'
#' @description \code{dlogmultinorm} calculates logarithms of multiple multivariate normal
#'   densities with varying mean and constant covariance matrix.
#'
#' @param y dimension \eqn{(T x k)} matrix where each row is a k-dimensional random vector
#' @param mu dimension \eqn{(T x k)} matrix where each row is the mean of the k-dimensional
#'   random vector in corresponding row of \code{y}.
#' @param Omega the \eqn{(k x k)} covariance matrix Omega.
#' @return Returns a size \eqn{(T x 1)} vector containing the multinormal densities in logarithm.

dlogmultinorm <- function(y, mu, Omega) {
  tmp <- -0.5*ncol(y)*log(2*pi) - 0.5*log(det(Omega))
  ymu <- y - mu
  tmp - 0.5*rowSums(ymu%*%solve(Omega)*ymu)
}
