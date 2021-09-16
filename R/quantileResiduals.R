
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
#' # GMVAR(1,2), d=2 model:
#' params12 <- c(0.55, 0.112, 0.344, 0.055, -0.009, 0.718, 0.319, 0.005, 0.03,
#'  0.619, 0.173, 0.255, 0.017, -0.136, 0.858, 1.185, -0.012, 0.136, 0.674)
#' mod12 <- GMVAR(gdpdef, p=1, M=2, params=params12)
#' quantile_residuals(mod12)
#'
#' # GMVAR(2,2), d=2 model with mean-parametrization:
#' params22 <- c(0.869, 0.549, 0.223, 0.059, -0.151, 0.395, 0.406, -0.005,
#'  0.083, 0.299, 0.215, 0.002, 0.03, 0.576, 1.168, 0.218, 0.02, -0.119,
#'  0.722, 0.093, 0.032, 0.044, 0.191, 1.101, -0.004, 0.105, 0.58)
#' mod22 <- GMVAR(gdpdef, p=2, M=2, params=params22, parametrization="mean")
#' quantile_residuals(mod22)
#'
#' # Structural GMVAR(2, 2), d=2 model identified with sign-constraints:
#' params22s <- c(0.36, 0.121, 0.484, 0.072, 0.223, 0.059, -0.151, 0.395,
#'  0.406, -0.005, 0.083, 0.299, 0.218, 0.02, -0.119, 0.722, 0.093, 0.032,
#'  0.044, 0.191, 0.057, 0.172, -0.46, 0.016, 3.518, 5.154, 0.58)
#' W_22 <- matrix(c(1, 1, -1, 1), nrow=2, byrow=FALSE)
#' mod22s <- GMVAR(gdpdef, p=2, M=2, params=params22s, structural_pars=list(W=W_22))
#' quantile_residuals(mod22s)
#' @export

quantile_residuals <- function(gmvar) {

  # Collect preliminary values
  check_gmvar(gmvar)
  epsilon <- round(log(.Machine$double.xmin) + 10)
  p <- gmvar$model$p
  M <- gmvar$model$M
  d <- gmvar$model$d
  data <- gmvar$data
  structural_pars <- gmvar$model$structural_pars
  n_obs <- nrow(data)
  T_obs <- n_obs - p

  # Collect parameter values
  params <- gmvar$params
  params <- reform_constrained_pars(p=p, M=M, d=d, params=params, constraints=gmvar$model$constraints,
                                    same_means=gmvar$model$same_means, structural_pars=structural_pars)
  structural_pars <- get_unconstrained_structural_pars(structural_pars=structural_pars)
  if(gmvar$model$parametrization == "mean") {
    params <- change_parametrization(p=p, M=M, d=d, params=params, constraints=NULL,
                                     structural_pars=structural_pars, change_to="intercept")
  }
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
#' @keywords internal

quantile_residuals_int <- function(data, p, M, params, conditional, parametrization, constraints=NULL,
                                   same_means=NULL, structural_pars=NULL, stat_tol=1e-3, posdef_tol=1e-8) {
  lok_and_mw <- loglikelihood_int(data=data, p=p, M=M, params=params, conditional=conditional,
                                  parametrization=parametrization, constraints=constraints,
                                  same_means=same_means, structural_pars=structural_pars,
                                  to_return="loglik_and_mw", check_params=TRUE, minval=NA,
                                  stat_tol=stat_tol, posdef_tol=posdef_tol)
  d <- ncol(data)
  npars <- n_params(p=p, M=M, d=d, constraints=constraints)
  mod <- structure(list(data=data,
                        model=list(p=p,
                                   M=M,
                                   d=ncol(data),
                                   conditional=conditional,
                                   parametrization=parametrization,
                                   constraints=constraints,
                                   same_means=same_means,
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
#' @keywords internal

dlogmultinorm <- function(y, mu, Omega) {
  tmp <- -0.5*ncol(y)*log(2*pi) - 0.5*log(det(Omega))
  ymu <- y - mu
  tmp - 0.5*rowSums(ymu%*%solve(Omega)*ymu)
}
