#' @title Simulate method for class 'gsmvar' objects
#'
#' @description \code{simulate.gsmvar} is a simulate method for class 'gsmvar' objects.
#'   It allows to simulate observations from a GMVAR, StMVAR, or G-StMVAR process.
#'
#' @param object an object of class \code{'gsmvar'}, typically created with \code{fitGSMVAR} or \code{GSMVAR}.
#' @param nsim number of observations to be simulated.
#' @param seed set seed for the random number generator?
#' @param ... currently not in use.
#' @param init_values a size \eqn{(pxd)} matrix specifying the initial values, where d is the number
#'   of time series in the system. The \strong{last} row will be used as initial values for the first lag,
#'   the second last row for second lag etc. If not specified, initial values will be drawn according to
#'   mixture distribution specifed by the argument \code{init_regimes}.
#' @param init_regimes a numeric vector of length at most \eqn{M} and elements
#'   in \eqn{1,...,M} specifying the regimes from which the initial values
#'   should be generated from. The initial values will be generated from a
#'   mixture distribution with the mixture components being the stationary
#'   distributions of the specific regimes and the (proportional) mixing weights
#'   given by the mixing weight parameters of those regimes. Note that if
#'   \code{init_regimes=1:M}, the initial values are generated from the
#'   stationary distribution of the process and if \code{init_regimes=m}, the
#'   initial value are generated from the stationary distribution of the
#'   \eqn{m}th regime. Ignored if the argument \code{init_values} is specified.
#' @param ntimes how many sets of simulations should be performed?
#' @param drop if \code{TRUE} (default) then the components of the returned list are coerced to lower dimension if \code{ntimes==1}, i.e.,
#'   \code{$sample} and \code{$mixing_weights} will be matrices, and \code{$component} will be vector.
#' @param girf_pars This argument is used internally in the estimation of generalized impulse response functions (see \code{?GIRF}). You
#'   should ignore it.
#' @details The argument \code{ntimes} is intended for forecasting: a GMVAR, StMVAR, or G-StMVAR process can be forecasted by simulating its possible future values.
#'  One can easily perform a large number simulations and calculate the sample quantiles from the simulated values to obtain prediction
#'  intervals (see the forecasting example).
#' @return If \code{drop==TRUE} and \code{ntimes==1} (default): \code{$sample}, \code{$component}, and \code{$mixing_weights} are matrices.
#'   Otherwise, returns a list with...
#'   \describe{
#'     \item{\code{$sample}}{a size (\code{nsim}\eqn{ x d x }\code{ntimes}) array containing the samples: the dimension \code{[t, , ]} is
#'      the time index, the dimension \code{[, d, ]} indicates the marginal time series, and the dimension \code{[, , i]} indicates
#'      the i:th set of simulations.}
#'     \item{\code{$component}}{a size (\code{nsim}\eqn{ x }\code{ntimes}) matrix containing the information from which mixture component each
#'      value was generated from.}
#'     \item{\code{$mixing_weights}}{a size (\code{nsim}\eqn{ x M x }\code{ntimes}) array containing the mixing weights corresponding to the
#'      sample: the dimension \code{[t, , ]} is the time index, the dimension \code{[, m, ]} indicates the regime, and the dimension
#'      \code{[, , i]} indicates the i:th set of simulations.}
#'   }
#' @seealso \code{\link{fitGSMVAR}}, \code{\link{GSMVAR}}, \code{\link{diagnostic_plot}}, \code{\link{predict.gsmvar}},
#'  \code{\link{profile_logliks}}, \code{\link{quantile_residual_tests}}, \code{\link{GIRF}}, \code{\link{GFEVD}}
#' @inherit loglikelihood_int references
#' @examples
#'  # GMVAR(1,2), d=2 process, initial values from the stationary
#'  # distribution
#'  params12 <- c(0.55, 0.112, 0.344, 0.055, -0.009, 0.718, 0.319, 0.005,
#'   0.03, 0.619, 0.173, 0.255, 0.017, -0.136, 0.858, 1.185, -0.012, 0.136,
#'   0.674)
#'  mod12 <- GSMVAR(p=1, M=2, d=2, params=params12)
#'  set.seed(1)
#'  sim12 <- simulate(mod12, nsim=500)
#'  plot.ts(sim12$sample)
#'  ts.plot(sim12$mixing_weights, col=c("blue", "red"), lty=2)
#'  plot(sim12$component, type="l")
#'
#'  # StMVAR(2, 2), d=2 model
#'  params22t <- c(0.554, 0.033, 0.184, 0.005, -0.186, 0.683, 0.256, 0.031,
#'   0.026, 0.204, 0.583, -0.002, 0.048, 0.697, 0.154, 0.049, 0.374, 0.476,
#'    0.318, -0.645, -0.302, -0.222, 0.193, 0.042, -0.013, 0.048, 0.818,
#'     4.334, 20)
#'  mod22t <- GSMVAR(gdpdef, p=2, M=2, params=params22t, model="StMVAR")
#'  sim22t <- simulate(mod22t, nsim=100)
#'  plot.ts(sim22t$mixing_weights)
#'
#'  ## FORECASTING EXAMPLE ##
#'  # Forecast 5-steps-ahead, 500 sets of simulations with initial
#'  # values from the data:
#'  # GMVAR(2,2), d=2 model
#'  params22 <- c(0.36, 0.121, 0.223, 0.059, -0.151, 0.395, 0.406, -0.005,
#'   0.083, 0.299, 0.215, 0.002, 0.03, 0.484, 0.072, 0.218, 0.02, -0.119,
#'    0.722, 0.093, 0.032, 0.044, 0.191, 1.101, -0.004, 0.105, 0.58)
#'  mod22 <- GSMVAR(gdpdef, p=2, M=2, params=params22)
#'  sim22 <- simulate(mod22, nsim=5, ntimes=500)
#'
#'  # Point forecast + 95% prediction intervals:
#'  apply(sim22$sample, MARGIN=1:2, FUN=quantile, probs=c(0.025, 0.5, 0.972))
#'
#'  # Similar forecast for the mixing weights:
#'  apply(sim22$mixing_weights, MARGIN=1:2, FUN=quantile,
#'        probs=c(0.025, 0.5, 0.972))
#' @export

simulate.gsmvar <- function(object, nsim=1, seed=NULL, ..., init_values=NULL, init_regimes=1:sum(gsmvar$model$M),
                            ntimes=1, drop=TRUE, girf_pars=NULL) {
  # girf_pars$shock_numb - which shock?
  # girf_pars$shock_size - size of the structural shock?
  # girf_pars$include_mixweights - should GIRFs be estimated for the mixing weights as well? TRUE or FALSE
  # If !is.null(girf_pars) and girf_pars$include_mixweights == TRUE, returns a size (N+1 x d+M) vector containing
  #                                                                  the estimated GIRFs for the variables and
  #                                                                  and the mixing weights (column d+m for the m:th regime).
  # If !is.null(girf_pars) and girf_pars$include_mixweights == FALSE, returns a size (N+1 x d) vector containing
  #                                                                   the estimated GIRFs for the variables only.
  # The first row for response at impact
  gsmvar <- object
  if(!is.null(seed)) set.seed(seed)
  epsilon <- round(log(.Machine$double.xmin) + 10)
  p <- gsmvar$model$p
  M_orig <- gsmvar$model$M
  M <- sum(M_orig)
  d <- gsmvar$model$d
  model <- gsmvar$model$model
  constraints <- gsmvar$model$constraints
  same_means <- gsmvar$model$same_means
  structural_pars <- gsmvar$model$structural_pars
  if(model == "GMVAR") { # The number of GMVAR type regimes
    M1 <- M
  } else if(model == "StMVAR") {
    M1 <- 0
  } else { # model == "G-StMVAR"
    M1 <- M_orig[1]
  }

  if(!all_pos_ints(c(nsim, ntimes))) stop("Arguments n and ntimes must be positive integers")
  if(!is.null(init_values)) {
    if(!is.matrix(init_values)) stop("init_values must be a numeric matrix")
    if(anyNA(init_values)) stop("init_values contains NA values")
    if(ncol(init_values) != d | nrow(init_values) < p) stop("init_values must contain d columns and at least p rows")
  }

  # Collect parameter values
  params <- gsmvar$params
  params <- reform_constrained_pars(p=p, M=M_orig, d=d, params=params, model=model,
                                    constraints=gsmvar$model$constraints,
                                    same_means=gsmvar$model$same_means,
                                    structural_pars=structural_pars)
  structural_pars <- get_unconstrained_structural_pars(structural_pars=structural_pars)
  if(gsmvar$model$parametrization == "mean") {
    params <- change_parametrization(p=p, M=M, d=d, params=params, constraints=NULL,
                                     structural_pars=structural_pars, change_to="intercept")
  }

  all_mu <- get_regime_means(gsmvar)
  all_phi0 <- pick_phi0(p=p, M=M_orig, d=d, params=params, structural_pars=structural_pars)
  all_A <- pick_allA(p=p, M=M_orig, d=d, params=params, structural_pars=structural_pars)
  all_Omega <- pick_Omegas(p=p, M=M_orig, d=d, params=params, structural_pars=structural_pars)
  all_boldA <- form_boldA(p=p, M=M_orig, d=d, all_A=all_A)
  alphas <- pick_alphas(p=p, M=M_orig, d=d, params=params, model=model)
  all_df <- pick_df(M=M_orig, params=params, model=model)

  all_Bm <- array(dim=c(d, d, M)) # Matrices such that all_Bm[, , m]%*%rnorm(d) follows N(0, \Omega_m)
  if(!is.null(structural_pars)) { # Calculate  W%*%Lambda_m^{1/2} where Lambda_1 = I_d (not to be confused with the time-varying B-matrix)
    W <- pick_W(p=p, M=M_orig, d=d, params=params, structural_pars=structural_pars)
    all_Bm[, , 1] <- W
    if(M > 1) {
      lambdas <- matrix(pick_lambdas(p=p, M=M_orig, d=d, params=params, structural_pars=structural_pars), nrow=d, byrow=FALSE)
      for(m in 2:M) {
        Lambda_m <- diag(lambdas[, m - 1])
        all_Bm[, , m] <- W%*%sqrt(Lambda_m)
      }
    }
  } else { # Reduced form model
    for(m in 1:M) { # t(chol(all_Omega[, , m]))
      eig <- eigen(all_Omega[, , m], symmetric=TRUE) # Orthogonal eigenvalue decomposition of the error term covariance matrix
      all_Bm[, , m] <- eig$vectors%*%tcrossprod(diag(sqrt(eig$values)), eig$vectors) # Symmetric square root matrix of the error term covariance matrix
    }
  }

  # Calculate the covariance matrices Sigma_{m,p} (Lutkepohl 2005, eq. (2.1.39) or the algorithm proposed by McElroy 2017)
  Sigmas <- get_Sigmas(p=p, M=M_orig, d=d, all_A=all_A, all_boldA=all_boldA, all_Omega=all_Omega) # Store the (dpxdp) covariance matrices
  inv_Sigmas <- array(NA, dim=c(d*p, d*p, M)) # Store inverses of the (dpxdp) covariance matrices
  det_Sigmas <- numeric(M) # Store determinants of the (dpxdp) covariance matrices
  chol_Sigmas <- array(dim=c(d*p, d*p, M)) # Cholesky decompositions of the  (dpxdp) covariance matrices
  for(m in 1:M) {
    chol_Sigmas[, , m] <- chol(Sigmas[, , m]) # Upper triangle
    inv_Sigmas[, , m] <- chol2inv(chol_Sigmas[, , m]) # Faster inverse
    det_Sigmas[m] <- prod(diag(chol_Sigmas[, , m]))^2 # Faster determinant
  }

  if(M1 < M) { # Calculate statistics used for StMVAR type regimes
    all_logC <- lgamma(0.5*(d*p + all_df)) - 0.5*d*p*log(base::pi) - 0.5*d*p*log(all_df - 2) - lgamma(0.5*all_df)
  }


  if(!is.null(girf_pars)) {
    R1 <- girf_pars$R1
  }

  # Relative mixing probabilities of the initial regimes
  reg_probs <- alphas[init_regimes]

  # Set/generate initial values
  if(is.null(init_values)) {
    # Generate initial values from the stationary distribution of the process or of the given regime
    m <- ifelse(length(init_regimes) == 1, init_regimes, sample(x=init_regimes, size=1, replace=TRUE, prob=reg_probs))  # From which mixture component the initial values are drawn from?
    mu <- rep(all_mu[, m], p)
    if(m <= M1) { # Draw initial values from a GMVAR type regime
      L <- t(chol_Sigmas[, , m]) # Lower triangle
      init_values <- matrix(mu + L%*%rnorm(d*p), nrow=p, ncol=d, byrow=TRUE) # i:th row for the i:th length d random vector
    } else { # m > M1, draw initial values from a StMVAR type regime
      # Note that we use parametrization: mean, covariance (=v/(v - 2)*Scale), df
      nu <- all_df[m - M1] # df
      L <- t(chol((nu - 2)/nu*Sigmas[, , m])) # Lower triangle
      Z <- mu + L%*%rnorm(d*p) # Sample from N(0, ((v - 2)/v)*Sigmas[, , m]
      X <- rchisq(1, df=nu) # Sample from Chisq(v)
      init_values <- matrix(mu + Z*sqrt(nu/X), nrow=p, ncol=d, byrow=TRUE) # i:th row for the i:th length d random vector
    }
  } else {
    # Fixed initial values given as argument
    init_values <- init_values[(nrow(init_values) - p + 1):nrow(init_values), , drop=FALSE]
  }

  # Container for the simulated values and initial values. First row row initial values vector, and t:th row for (y_{i-1},...,y_{i-p})
  Y <- matrix(nrow=nsim + 1, ncol=d*p)
  Y[1,] <- reform_data(init_values, p=p)
  if(!is.null(girf_pars)) Y2 <- Y # Storage for the second sample path in the GIRF algorithm

  # Initialize data structures
  sample <- array(dim=c(nsim, d, ntimes))
  component <- matrix(nrow=nsim, ncol=ntimes)
  mixing_weights <- array(dim=c(nsim, M, ntimes))
  if(!is.null(girf_pars)) {
    sample2 <- array(dim=c(nsim, d, ntimes))
    mixing_weights2 <- array(dim=c(nsim, M, ntimes))
  }

  # Some functions to be used
  get_matprods <- function(Y) vapply(1:M, function(m) crossprod(Y[i1,] - rep(all_mu[, m], p), inv_Sigmas[, , m])%*%(Y[i1,] - rep(all_mu[, m], p)), numeric(1))
  get_logmvdvalues <- function(matprods, arch_scalars) vapply(1:M, function(m) {
    if(m <= M1) { # GMVAR type regime
      return(-0.5*d*p*log(2*pi) - 0.5*log(det_Sigmas[m]) - 0.5*matprods[m]) # Log multivariate normal density
    } else { # StMVAR type regime
      return(all_logC[m - M1] - 0.5*log(det_Sigmas[m]) - 0.5*(d*p + all_df[m - M1])*log(1 + matprods[m]/(all_df[m - M1] - 2)))
    }
  }, numeric(1))

  for(j1 in seq_len(ntimes)) {
    for(i1 in seq_len(nsim)) {
      # Calculate the dp-dimensional multinormal densities (KMS 2016, eq.(6)).
      # Calculated in logarithm because same values may be too close to zero for machine accuracy.
      matprods <- get_matprods(Y) # Not that get_matprods obtains i1 as non-formal argument
      log_mvdvalues <- get_logmvdvalues(matprods)
      alpha_mt <- get_alpha_mt(M=M, log_mvdvalues=log_mvdvalues, alphas=alphas,
                               epsilon=epsilon, also_l_0=FALSE)

      # Calculate the time varying arch-scalars for StMVAR type regimes
      if(M1 < M) { # StMVAR and G-StMVAR models
        arch_scalars <- (all_df - 2 + matprods[(M1 + 1):M])/(all_df - 2 + d*p) # [m - M1] indexing
      }

      # Draw a mixture component and store the values
      m <- sample.int(n=M, size=1, replace=TRUE, prob=alpha_mt)
      component[i1, j1] <- m
      mixing_weights[i1, , j1] <- alpha_mt

      # Calculate conditional mean for regime m
      A2 <- matrix(all_A[, , , m], nrow=d, byrow=FALSE) # (A_1:...:A_p)
      mu_mt <- all_phi0[, m] + A2%*%Y[i1,]

      # Draw the sample and store it
      eps_t <- rnorm(d) # We use the same std normal shocks to create the Student't variables as well to control random variation across the sample paths
      if(M1 < M) { # StMVAR and G-StMVAR models
        # We generate a chi^2 variable with df nu + dp for all StMVAR type regimes and when computing GIRF use the same ones in both sample paths
        all_chisq_rv <- vapply(all_df + d*p, function(i1) rchisq(n=1, df=i1), numeric(1)) # [m - M1] indexed
      }

      if(m <= M1) { # GMVAR type regime
        sample[i1, , j1] <- mu_mt + all_Bm[, , m]%*%eps_t # Conditional mean + reduced form shock
      } else { # StMVAR type regime
        # Here, we obtain a (reduced form) multivariate Student's t shock with time-varying conditional variance
        df_to_use <- all_df[m - M1] + d*p
        Z <- sqrt(arch_scalars[m - M1])*(df_to_use - 2)/df_to_use*all_Bm[, , m]%*%eps_t # Sample from N(0, arch_scalar*(df - 2)/df*Omega_m))
        sample[i1, , j1] <- mu_mt + Z*sqrt(df_to_use/all_chisq_rv[m - M1]) # Sample from t_d(mu_mt, arch_scalar*Omega_m, all_df[m - M1] + d*p)
      }

      # Update storage matrix Y (overwrites when ntimes > 1)
      if(p == 1) {
        Y[i1 + 1,] <- sample[i1, , j1]
      } else {
        Y[i1 + 1,] <- c(sample[i1, , j1], Y[i1 , 1:(d*p - d)])
      }

      # For the second sample path in GIRF (with a specific structural shock occurring)
      if(!is.null(girf_pars)) {
        matprods2 <- get_matprods(Y2)
        log_mvdvalues2 <- get_logmvdvalues(matprods2)
        alpha_mt2 <- get_alpha_mt(M=M, log_mvdvalues=log_mvdvalues2, alphas=alphas,
                                  epsilon=epsilon, also_l_0=FALSE)
        mixing_weights2[i1, , j1] <- alpha_mt2

        # Calculate the time varying arch-scalars for StMVAR type regimes
        if(M1 < M) { # StMVAR and G-StMVAR models
          arch_scalars2 <- (all_df - 2 + matprods2[(M1 + 1):M])/(all_df - 2 + d*p) # [m - M1] indexing
        }
        if(i1 == 1) {
          m2 <- m # Common regime at impact (the mixing weights are the same so m is drawn from the same distribution)
        } else {
          m2 <- sample.int(n=M, size=1, replace=TRUE, prob=alpha_mt2) # Draw the regime m
        }
        A22 <- matrix(all_A[, , , m2], nrow=d, byrow=FALSE) # (A_1:...:A_p)
        mu_mt2 <- all_phi0[, m2] + A22%*%Y2[i1,]

        # Compute reduced form shock from the same NID(0,I) (and Chisq) shock as the other sample path
        if(m2 <= M1) { # GMVAR type regime
          u_t <- all_Bm[, , m2]%*%eps_t
        } else {# StMVAR type regime
          df_to_use <- all_df[m2 - M1] + d*p
          Z <- sqrt(arch_scalars[m2 - M1])*(df_to_use - 2)/df_to_use*all_Bm[, , m2]%*%eps_t # Sample from N(0, arch_scalar*(v - 2)/v*Omega_m))
          u_t <- Z*sqrt(df_to_use/all_chisq_rv[m2 - M1]) # Sample from t_d(0, arch_scalar*Omega_m, all_df[m - M1] + d*p)
        }

        if(i1 == 1) { # At impact, obtain reduced form shock from the specific structural shock

          # Calculate the time-varying B-matrix
          if(M == 1) {
            B_t <- W
          } else {
            tmp <- array(dim=c(d, d, M))
            if(model == "StMVAR") { # The first regime is StMVAR type
              multiplier <- arch_scalars[1]*alpha_mt2[1]
            } else {
              multiplier <-  alpha_mt2[1]
            }
            tmp[, , 1] <- multiplier*diag(d)

            for(m in 2:M) {
              if(m <= M1) { # GMVAR type regime
                multiplier <- alpha_mt2[m]
              } else { # StMVAR type regime
                multiplier <- arch_scalars[m - M1]*alpha_mt2[m]
              }
              tmp[, , m] <- multiplier*diag(lambdas[, m - 1])
            }
            B_t <- W%*%sqrt(apply(tmp, MARGIN=1:2, FUN=sum))
          }
          e_t <- solve(B_t, u_t) # Structural shock
          e_t[girf_pars$shock_numb] <- girf_pars$shock_size # Impose the size of a shock
          u_t <- B_t%*%e_t # The reduced form shock corresponding to the specific sized structural shock in the j:th element
        }

        sample2[i1, , j1] <- mu_mt2 + u_t

        if(p == 1) {
          Y2[i1 + 1,] <- sample2[i1, , j1]
        } else {
          Y2[i1 + 1,] <- c(sample2[i1, , j1], Y2[i1 , 1:(d*p - d)])
        }
      }
    }
  }

  # Calculate a single GIRF for the given structural shock: (N + 1 x d) matrix
  if(!is.null(girf_pars)) {
    one_girf <- apply(X=sample2 - sample, MARGIN=1:2, FUN=mean)
    if(!is.null(gsmvar$data)) {
      colnames(one_girf) <- colnames(gsmvar$data)
    } else {
      colnames(one_girf) <- paste("shock", 1:d)
    }
    if(girf_pars$include_mixweights) {
      mix_girf <- apply(X=mixing_weights2 - mixing_weights, MARGIN=1:2, FUN=mean)
      colnames(mix_girf) <- paste("mw reg.", 1:M)
      return(cbind(one_girf, mix_girf))
    } else {
      return(one_girf)
    }
  }

  if(ntimes == 1 & drop) {
    sample <- matrix(sample, nrow=nsim, ncol=d)
    component <- as.vector(component)
    mixing_weights <- matrix(mixing_weights, nrow=nsim, ncol=M)
  }

  list(sample=sample,
       component=component,
       mixing_weights=mixing_weights)
}
