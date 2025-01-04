#' @title Calculate multivariate quantile residuals of a GMVAR, StMVAR, or G-StMVAR model
#'
#' @description \code{quantile_residuals} calculates multivariate quantile residuals
#'  (proposed by \emph{Kalliovirta and Saikkonen 2010}) for a GMVAR, StMVAR, or G-StMVAR model.
#'
#' @inheritParams quantile_residual_tests
#' @return Returns \eqn{((n_obs-p) x d)} matrix containing the multivariate quantile residuals,
#'   \eqn{j}:th column corresponds to the time series in the \eqn{j}:th column of the data. The multivariate
#'   quantile residuals are calculated so that the first column quantile residuals are the "unconditioned ones"
#'   and the rest condition on all the previous ones in numerical order. Read the cited article by
#'   \emph{Kalliovirta and Saikkonen 2010} for details.
#' @inherit GSMVAR references
#' @seealso \code{\link{fitGSMVAR}}, \code{\link{GSMVAR}}, \code{\link{quantile_residual_tests}},
#'   \code{\link{diagnostic_plot}}, \code{\link{predict.gsmvar}}, \code{\link{profile_logliks}}
#' @examples
#' # GMVAR(1,2), d=2 model:
#' params12 <- c(0.55, 0.112, 0.344, 0.055, -0.009, 0.718, 0.319, 0.005, 0.03,
#'  0.619, 0.173, 0.255, 0.017, -0.136, 0.858, 1.185, -0.012, 0.136, 0.674)
#' mod12 <- GSMVAR(gdpdef, p=1, M=2, params=params12)
#' quantile_residuals(mod12)
#'
#' # GMVAR(2,2), d=2 model with mean-parametrization:
#' params22 <- c(0.869, 0.549, 0.223, 0.059, -0.151, 0.395, 0.406, -0.005,
#'  0.083, 0.299, 0.215, 0.002, 0.03, 0.576, 1.168, 0.218, 0.02, -0.119,
#'  0.722, 0.093, 0.032, 0.044, 0.191, 1.101, -0.004, 0.105, 0.58)
#' mod22 <- GSMVAR(gdpdef, p=2, M=2, params=params22, parametrization="mean")
#' quantile_residuals(mod22)
#'
#' # Structural GMVAR(2, 2), d=2 model identified with sign-constraints:
#' params22s <- c(0.36, 0.121, 0.484, 0.072, 0.223, 0.059, -0.151, 0.395,
#'  0.406, -0.005, 0.083, 0.299, 0.218, 0.02, -0.119, 0.722, 0.093, 0.032,
#'  0.044, 0.191, 0.057, 0.172, -0.46, 0.016, 3.518, 5.154, 0.58)
#' W_22 <- matrix(c(1, 1, -1, 1), nrow=2, byrow=FALSE)
#' mod22s <- GSMVAR(gdpdef, p=2, M=2, params=params22s, structural_pars=list(W=W_22))
#' quantile_residuals(mod22s)
#' @export

quantile_residuals <- function(gsmvar) {
  gsmvar <- gmvar_to_gsmvar(gsmvar) # Backward compatibility

  # Collect preliminary values
  check_gsmvar(gsmvar)
  epsilon <- round(log(.Machine$double.xmin) + 10)
  p <- gsmvar$model$p
  M <- gsmvar$model$M
  d <- gsmvar$model$d
  model <- gsmvar$model$model
  data <- gsmvar$data
  structural_pars <- gsmvar$model$structural_pars
  n_obs <- nrow(data)
  T_obs <- n_obs - p

  # Collect parameter values
  params <- reform_constrained_pars(p=p, M=M, d=d, params=gsmvar$params, model=model, constraints=gsmvar$model$constraints,
                                    same_means=gsmvar$model$same_means, weight_constraints=gsmvar$model$weight_constraints,
                                    structural_pars=structural_pars)
  structural_pars <- get_unconstrained_structural_pars(structural_pars=structural_pars)
  if(gsmvar$model$parametrization == "mean") {
    params <- change_parametrization(p=p, M=M, d=d, params=params, model=model, constraints=NULL, weight_constraints=NULL,
                                     structural_pars=structural_pars, change_to="intercept")
  }
  all_mu <- get_regime_means(gsmvar)
  all_phi0 <- pick_phi0(p=p, M=M, d=d, params=params, structural_pars=structural_pars)
  all_A <- pick_allA(p=p, M=M, d=d, params=params, structural_pars=structural_pars)
  all_Omega <- pick_Omegas(p=p, M=M, d=d, params=params, structural_pars=structural_pars)
  all_boldA <- form_boldA(p=p, M=M, d=d, all_A=all_A)
  alphas <- pick_alphas(p=p, M=M, d=d, params=params, model=model)
  all_df <- pick_df(M=M, params=params, model=model)

  if(model == "GMVAR") {
    M1 <- M
    M2 <- 0
  } else if(model == "StMVAR") {
    M1 <- 0
    M2 <- M
  } else { # model == "G-StMVAR"
    M1 <- M[1]
    M2 <- M[2]
  }
  M <- sum(M) # The total number of mixture components

  # An i:th row denotes the vector \bold{y_{i-1}} = (y_{i-1}',...,y_{i-p}') (dpx1),
  # assuming the observed data is y_{-p+1},...,y_0,y_1,...,y_{T}
  Y <- reform_data(data, p)

  # Mixing weights
  alpha_mt <- gsmvar$mixing_weights

  # Calculate the conditional means mu_{m,t} (KMS 2016, Condition 1 (a)).
  # Map to dimensions of mu_mt: [t, p, m]
  all_A2 <- array(all_A, dim=c(d, d*p, M)) # cbind coefficient matrices of each component: m:th component is obtained at [, , m]
  Y2 <- Y[1:T_obs,] # Last row is not needed because mu_mt uses lagged values
  mu_mt <- array(vapply(1:M, function(m) t(all_phi0[, m] + tcrossprod(all_A2[, , m], Y2)), numeric(d*T_obs)), dim=c(T_obs, d, M)) # [, , m]

  # The arch scalars multiplying the error covariance matrices (ones for GMVAR type regimes)
  arch_scalars <- gsmvar$arch_scalars

  ## Start computing the multivariate quantile residuals (Kalliovirta and Saikkonen 2010, eq.(4))
  # using properties of marginal and conditional distributions of multinormal random variables and
  # applying them to the mixture components at each time point t.

  # Compute partitions and matrix products of partitioned covariance matrices Omega_m that will be used multiple times.
  upleft_jjmat <-  function(mat, j) mat[1:j, 1:j , drop=FALSE] # function(arr, j) arr[1:j, 1:j, , drop=FALSE]
  # Above returns upper-left (j x j) block matrix from each slice of the 3d array "arr"

  # Storage for variances and conditional means
  # (obtained from properties of multinormal/multistudent for 1-dimensional conditional distribution)
  Omega_mtj <- array(dim=c(T_obs, d, M)) # [t, j, m]; time-varying for StMVAR type regimes, conditional on F_{t-1} and A_{j-1}
  mu_mtj <- array(dim=c(T_obs, d, M)) # [t, j, m], conditional on F_{t-1} and A_{j-1}

  # Inverses of the upper-left (j-1 x j-1) blocks of Omega; and logarithms of the determinants of the non-inverted block matrices
  all_inv_Omega_j_minus_1 <- array(dim=c(d-1, d-1, d-1, M)) # [d, d, j-1, m], inverses to be filled in the upper-left (j-1 x j-1) block
  log_det_Omega_j_minus_1 <- matrix(nrow=d-1, ncol=M) # [j-1, M]

  # Calculate 1-dimensional conditional variances and means, conditional on F_{t-1} and A_{j-1}, for t=1,...,T, m=1,...,M, and j=1,...,d
  dat <- data[(p + 1):nrow(data),] # Remove the initial values
  for(m in 1:M) {
    # j = 1; conditional on previous observations only
    Omega_mtj[, 1, m] <- arch_scalars[,m]*all_Omega[1, 1, m]
    mu_mtj[, 1, m] <- mu_mt[, 1, m] # Marginal conditional means
    mu_dif <- dat - mu_mt[, , m] # mu_mt - y_t, t = 1,...,T


    # j=2,...,d; conditional on previous observations and y_{1,t},...,y_{j-1,t}
    for(j in 2:d) {
      Omega_mj <- upleft_jjmat(all_Omega[, , m], j) # (j x j)
      up_left <- upleft_jjmat(Omega_mj, j - 1) # (j-1 x j-1)
      up_right <- Omega_mj[1:(j - 1), j, drop=FALSE] # (j-1 x 1)
      low_right <- Omega_mj[j, j] # (1 x 1)

      chol_up_left <- chol(up_left)
      all_inv_Omega_j_minus_1[1:(j-1), 1:(j-1), j-1,
                              m] <- inv_Omega_j_minus_1 <- chol2inv(chol_up_left) # Inv of upper left (j-1 x j-1) Omega_m block matrix
      log_det_Omega_j_minus_1[j-1, m] <- 2*log(prod(diag(chol_up_left))) # Log of the det of inv of upper left (j-1 x j-1) Omega_m block matrix
      matprod <- crossprod(up_right, inv_Omega_j_minus_1) # (1 x j-1)

      # Common formula for GMVAR and StMVAR type regimes
      mu_dif_j_minus_1 <- mu_dif[, 1:(j - 1), drop=FALSE]
      mu_mtj[, j, m] <- mu_mt[, j, m] + t(tcrossprod(matprod, mu_dif_j_minus_1))

      if(m <= M1) { # Constant conditional variance for GMVAR type regimes Omega_mtj [t, j, m];
        Omega_mtj[, j, m] <- low_right - matprod%*%up_right
      } else { # Time-varying conditional variance for StMVAR type regimes
        arch_scalars2 <- (all_df[m - M1] + d*p +
                            rowSums(mu_dif_j_minus_1%*%inv_Omega_j_minus_1*mu_dif_j_minus_1)/arch_scalars[,m])/(all_df[m - M1] + d*p + j - 3)
        Omega_mtj[, j, m] <- arch_scalars2*arch_scalars[,m]*c(low_right - matprod%*%up_right)
      }
    }
  }

  # Calculate beta_{m,t,j} for j=2,...,d (Virolainen 2022, eq. (2.6),(2.7) unpublished work paper)
  beta_mtj <- array(dim=c(T_obs, M, d)) # [t, m, j] j=1,...,d
  beta_mtj[, , 1] <- alpha_mt
  for(j in 2:d) {
    log_mvdvalues <- matrix(nrow=T_obs, ncol=M)
    for(m in 1:M) {
      if(m <= M1) { # GMVAR type regime
        log_mvdvalues[,m] <- dlogmultinorm(y=dat[,1:(j - 1), drop=FALSE],
                                           mu=as.matrix(mu_mt[, 1:(j - 1), m]),
                                           inv_Omega=as.matrix(all_inv_Omega_j_minus_1[1:(j-1), 1:(j-1), j-1, m]),
                                           log_det_Omega=log_det_Omega_j_minus_1[j-1, m])
      } else { # StMVAR type regime
      log_mvdvalues[,m] <- dlogmultistudent(y=dat[,1:(j - 1), drop=FALSE],
                                            mu=as.matrix(mu_mt[, 1:(j - 1), m]),
                                            inv_Omega=as.matrix(all_inv_Omega_j_minus_1[1:(j-1), 1:(j-1), j-1, m]),
                                            log_det_Omega=log_det_Omega_j_minus_1[j-1, m],
                                            arch_scalars=arch_scalars[,m],
                                            df=all_df[m - M1] + d*p)
      }
    }

    small_logmvds <- log_mvdvalues < epsilon
    if(any(small_logmvds)) {
      # If too small or large non-log-density values are present (i.e., that would yield -Inf or Inf),
      # we replace them with ones that are not too small or large but imply the same mixing weights
      # up to negligible numerical tolerance.
      which_change <- rowSums(small_logmvds) > 0 # Which rows contain too small  values
      to_change <- log_mvdvalues[which_change, , drop=FALSE]
      largest_vals <- do.call(pmax, split(to_change, f=rep(1:ncol(to_change), each=nrow(to_change)))) # The largest values of those rows
      diff_to_largest <- to_change - largest_vals # Differences to the largest value of the row

      # For each element in each row, check the (negative) distance from the largest value of the row. If the difference
      # is smaller than epsilon, replace the with epsilon. The results are then the new log_mvd values.
      diff_to_largest[diff_to_largest < epsilon] <- epsilon

      # Replace the old log_mvdvalues with the new ones
      log_mvdvalues[which_change,] <- diff_to_largest
    }

    numerators <- as.matrix(alpha_mt*exp(log_mvdvalues))
    denominator <- rowSums(numerators)
    beta_mtj[, , j] <- numerators/denominator

  }

  # Then, we calculate the conditional cumulative distribution function values, for all m=1,...,M, t=1,...,T, and j=1,...,d.
  F_values_regime <- array(dim=c(T_obs, d, M)) # [t, d, m]
  for(m in 1:M) {
    ydiff <- dat - mu_mtj[, , m] # [t, j]
    if(m <= M1) { # GMVAR type regimes
      F_values_regime[, , m] <- pnorm(ydiff/sqrt(Omega_mtj[, , m])) # Omega_mtj, mu_mtj = [t, j, m]
    } else { # StMVAR type regimes
      # Function for numerical integration of the Student's t pdf (when closed form expression if not available)
      my_integral <- function(t) { # Takes the observation index t (1,...,T) as formal argument and the rest from parent frame
        tryCatch(integrate(function(y_t) { # The conditional density function to be integrated numerically;
          C0/sqrt(Omega_mtj[t, j, m])*(1 + ((y_t - mu_mtj[t, j, m])^2)/(Omega_mtj[t, j, m]*(df - 2)))^(-0.5*(1 + df))
        }, lower=-Inf, upper=dat[t, j])$value,
        error=function(e) {
          warning("Couldn't analytically nor numerically integrate all quantile residuals:")
          warning(e)
          return(NA)
        })
      }

      # Go through dimensions
      for(j in 1:d) {
        df <- all_df[m - M1] + d*p + j - 1
        which_def <- which(abs(ydiff[,j]) < sqrt(Omega_mtj[, j, m]*(df - 2)))
        if(length(which_def) > 0) {
          which_not_def <- (1:nrow(ydiff))[-which_def]
        } else {
          which_not_def <- 1:nrow(ydiff)
        }
        C0 <- exp(lgamma(0.5*(1 + df)) - 0.5*log(base::pi) - 0.5*log(df - 2) - lgamma(0.5*df))

        # Calculate CDF values using hypergeometric function whenever it is defined
        if(length(which_def) > 0) {
          ydiff0 <- ydiff[which_def, j]
          F_values_regime[which_def, j, m] <- 0.5 +
            C0/sqrt(Omega_mtj[which_def, j, m])*ydiff0*gsl::hyperg_2F1(a=0.5, b=0.5*(1 + df), c=1.5,
                                                                       x=-(ydiff0^2)/(Omega_mtj[which_def, j, m]*(df - 2)),
                                                                       give=FALSE, strict=TRUE)
        }

        # Calculate CDF values by numerically integrating the t-densities whenever hypergeometric function is not defined
        if(length(which_not_def) > 0) {
          for(t in which_not_def) {
            F_values_regime[t, j, m] <- my_integral(t=t)
          }
        }
      }
    }
  }

  # Then calculate the conditional (marginal) CDFs of the process
  F_values <- vapply(1:d, function(j) rowSums(as.matrix(beta_mtj[, , j]*F_values_regime[, j, ])), numeric(T_obs))

  # Values too close to 0 or 1 will be scaled so that qnorm doesn't return inf-values
  F_values[F_values >= 1 - .Machine$double.eps/2] <- 1 - .Machine$double.eps/2
  F_values[F_values <= .Machine$double.eps/2] <- .Machine$double.eps/2

  qnorm(F_values)
}


#' @title Calculate multivariate quantile residuals of GMVAR, StMVAR, or G-StMVAR model
#'
#' @description \code{quantile_residuals_int} is a wrapper for \code{quantile_residuals} to compute
#'   quantile residuals using parameter values instead of class \code{gsmvar} object.
#'
#' @inheritParams loglikelihood_int
#' @section Warning:
#'   No argument checks!
#' @inherit quantile_residuals return references
#' @keywords internal

quantile_residuals_int <- function(data, p, M, params, model=c("GMVAR", "StMVAR", "G-StMVAR"), conditional, parametrization, constraints=NULL,
                                   same_means=NULL, weight_constraints=NULL, structural_pars=NULL,
                                   stat_tol=1e-3, posdef_tol=1e-8, df_tol=1e-8) {
  model <- match.arg(model)
  loglik_mw_archscalars <- loglikelihood_int(data=data, p=p, M=M, params=params, model=model, conditional=conditional,
                                             parametrization=parametrization, constraints=constraints,
                                             weight_constraints=weight_constraints,
                                             same_means=same_means, structural_pars=structural_pars,
                                             to_return="loglik_mw_archscalars", check_params=TRUE, minval=NA,
                                             stat_tol=stat_tol, posdef_tol=posdef_tol, df_tol=df_tol)
  d <- ncol(data)
  npars <- n_params(p=p, M=M, d=d, model=model, constraints=constraints)

  mod <- structure(list(data=data,
                        model=list(p=p,
                                   M=M,
                                   d=ncol(data),
                                   model=model,
                                   conditional=conditional,
                                   parametrization=parametrization,
                                   constraints=constraints,
                                   same_means=same_means,
                                   weight_constraints=weight_constraints,
                                   structural_pars=structural_pars),
                        params=params,
                        std_errors=rep(NA, npars),
                        mixing_weights=loglik_mw_archscalars$mw,
                        arch_scalars=loglik_mw_archscalars$arch_scalars,
                        quantile_residuals=NA,
                        loglik=structure(loglik_mw_archscalars$loglik,
                                         class="logLik",
                                         df=npars),
                        IC=NA,
                        all_estimates=NULL,
                        all_logliks=NULL,
                        which_converged=NULL,
                        num_tols=list(stat_tol=stat_tol,
                                      posdef_tol=posdef_tol,
                                      df_tol=df_tol)),
                   class="gsmvar")

  quantile_residuals(mod)
}


#' @title Calculate logarithms of multiple multivariate normal densities with varying
#'  mean and constant covariance matrix
#'
#' @description \code{dlogmultinorm} calculates logarithms of multiple multivariate normal
#'   densities with varying mean and constant covariance matrix.
#'
#' @param y dimension \eqn{(T x k)} matrix where each row is a k-dimensional vector
#' @param mu dimension \eqn{(T x k)} matrix where each row is the mean of the k-dimensional
#'   vector in corresponding row of \code{y}.
#' @param inv_Omega inverse of the \eqn{(k x k)} covariance matrix Omega.
#' @param log_det_Omega logarithm of the determinant of the covariance matrix Omega.
#' @return Returns a size \eqn{(T x 1)} vector containing the multinormal densities in logarithm.
#' @keywords internal

dlogmultinorm <- function(y, mu, inv_Omega, log_det_Omega) {
  tmp <- -0.5*ncol(y)*log(2*base::pi) - 0.5*log_det_Omega
  ymu <- y - mu
  tmp - 0.5*rowSums(ymu%*%inv_Omega*ymu)
}


#' @title Calculate logarithms of multiple multivariate Student's t densities with varying
#'  mean and covariance matrix of specific structure, but constant degrees of freedom.
#'
#' @description \code{dlogmultistudent} calculates logarithms of multiple multivariate
#'  Student's t densities with varying mean and covaraince matrix of specific structure.
#'
#' @inheritParams dlogmultinorm
#' @param arch_scalars length \eqn{T} numeric vector containing the coefficients that multiply
#'   the covariance matrix \code{Omega}.
#' @param df the degrees of freedom parameter that is common for all \eqn{t=1,...,T}.
#' @return Returns a size \eqn{(T x 1)} vector containing the  multivariate Student's t
#'   densities in logarithm.
#' @keywords internal

dlogmultistudent <- function(y, mu, inv_Omega, log_det_Omega, arch_scalars, df) {
  logC <- lgamma(0.5*(ncol(y) + df)) - 0.5*ncol(y)*log(base::pi) - 0.5*ncol(y)*log(df - 2) - lgamma(0.5*df)
  tmp <- -0.5*ncol(y)*log(arch_scalars) - 0.5*log_det_Omega
  ymu <- y - mu
  logC + tmp - 0.5*(ncol(y) + df)*log(1 + rowSums(ymu%*%inv_Omega*ymu)/(arch_scalars*(df - 2)))
}
