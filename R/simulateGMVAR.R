

#' @title Simulate from GMVAR process
#'
#' @description \code{simulateGMVAR} simulates observations from a GMVAR process.
#'
#' @param gmvar an object of class \code{'gmvar'} created with \code{fitGMVAR} or \code{GMVAR}.
#' @param nsimu number of observations to be simulated.
#' @param init_values a size \eqn{(pxd)} matrix specifying the initial values to be used in the simulation, where
#'   d is the number of time series in the system.
#'   The \strong{last} row will be used as initial values for the first lag, the second last row for second lag etc. If not
#'   specified, initial values will be drawn from the stationary distribution of the process.
#' @param ntimes how many sets of simulations should be performed?
#' @param drop if \code{TRUE} (default) then the components of the returned list are coerced to lower dimension if \code{ntimes==1}, i.e.,
#'   \code{$sample} and \code{$mixing_weights} will be matrices, and \code{$component} will be vector.
#' @param seed set seed for the random number generator?
#' @param girf_pars This argument is used internally in the estimation of generalized impulse response functions (see \code{?GIRF}). You
#'   should ignore it.
#' @details The argument \code{ntimes} is intended for forecasting: a GMVAR process can be forecasted by simulating its possible future values.
#'  One can easily perform a large number simulations and calculate the sample quantiles from the simulated values to obtain prediction
#'  intervals (see the forecasting example).
#' @return If \code{drop==TRUE} and \code{ntimes==1} (default): \code{$sample}, \code{$component}, and \code{$mixing_weights} are matrices.
#'   Otherwise, returns a list with...
#'   \describe{
#'     \item{\code{$sample}}{a size (\code{nsimu}\eqn{ x d x }\code{ntimes}) array containing the samples: the dimension \code{[t, , ]} is
#'      the time index, the dimension \code{[, d, ]} indicates the marginal time series, and the dimension \code{[, , i]} indicates
#'      the i:th set of simulations.}
#'     \item{\code{$component}}{a size (\code{nsimu}\eqn{ x }\code{ntimes}) matrix containing the information from which mixture component each
#'      value was generated from.}
#'     \item{\code{$mixing_weights}}{a size (\code{nsimu}\eqn{ x M x }\code{ntimes}) array containing the mixing weights corresponding to the
#'      sample: the dimension \code{[t, , ]} is the time index, the dimension \code{[, m, ]} indicates the regime, and the dimension
#'      \code{[, , i]} indicates the i:th set of simulations.}
#'   }
#' @seealso \code{\link{fitGMVAR}}, \code{\link{GMVAR}}, \code{\link{diagnostic_plot}}, \code{\link{predict.gmvar}},
#'  \code{\link{profile_logliks}}, \code{\link{quantile_residual_tests}}, \code{\link{GIRF}}
#' @inherit loglikelihood_int references
#' @examples
#'  # These examples use the data 'eurusd' which comes with the
#'  # package, but in a scaled form.
#'  data <- cbind(10*eurusd[,1], 100*eurusd[,2])
#'  colnames(data) <- colnames(eurusd)
#'
#'  # GMVAR(1,2), d=2 process, initial values from the stationary
#'  # distribution
#'  params122 <- c(0.623, -0.129, 0.959, 0.089, -0.006, 1.006, 1.746,
#'   0.804, 5.804, 3.245, 7.913, 0.952, -0.037, -0.019, 0.943, 6.926,
#'   3.982, 12.135, 0.789)
#'  mod122 <- GMVAR(p=1, M=2, d=2, params=params122)
#'  set.seed(1)
#'  sim122 <- simulateGMVAR(mod122, nsimu=500)
#'  plot.ts(sim122$sample)
#'  ts.plot(sim122$mixing_weights, col=c("blue", "red"), lty=2)
#'  plot(sim122$component, type="l")
#'
#'  # Structural GMVAR(2, 2), d=2 model identified with sign-constraints:
#'  params222s <- c(-11.964, 155.024, 11.636, 124.988, 1.314, 0.145, 0.094, 1.292,
#'    -0.389, -0.07, -0.109, -0.281, 1.248, 0.077, -0.04, 1.266, -0.272, -0.074,
#'    0.034, -0.313, 0.903, 0.718, -0.324, 2.079, 7.00, 1.44, 0.742)
#'  W_222 <- matrix(c(1, 1, -1, 1), nrow=2, byrow=FALSE)
#'  mod222s <- GMVAR(data, p=2, M=2, params=params222s, parametrization="mean",
#'    structural_pars=list(W=W_222))
#'  sim222s <- simulateGMVAR(mod222s, nsimu=100)
#'  plot.ts(sim222s$sample)
#'
#'  ## FORECASTING EXAMPLE ##
#'  # Forecast 5-steps-ahead, 500 sets of simulations with initial
#'  # values from the data:
#'  # GMVAR(2,2), d=2 model with mean-parametrization:
#'  params222 <- c(-11.904, 154.684, 1.314, 0.145, 0.094, 1.292, -0.389,
#'   -0.070, -0.109, -0.281, 0.920, -0.025, 4.839, 11.633, 124.983, 1.248,
#'    0.077, -0.040, 1.266, -0.272, -0.074, 0.034, -0.313, 5.855, 3.570,
#'    9.838, 0.740)
#'  mod222 <- GMVAR(data, p=2, M=2, params=params222, parametrization="mean")
#'  sim222 <- simulateGMVAR(mod222, nsimu=5, ntimes=500)
#'
#'  # Point forecast + 95% prediction intervals:
#'  apply(sim222$sample, MARGIN=1:2, FUN=quantile, probs=c(0.025, 0.5, 0.972))
#'
#'  # Similar forecast for the mixing weights:
#'  apply(sim222$mixing_weights, MARGIN=1:2, FUN=quantile,
#'        probs=c(0.025, 0.5, 0.972))
#'
#'
#'  # GMVAR(2,2), d=2 model with AR parameters restricted to be
#'  # the same for both regimes, custom inital values:
#'  C_mat <- rbind(diag(2*2^2), diag(2*2^2))
#'  params222c <- c(1.031, 2.356, 1.786, 3.000, 1.250, 0.060, 0.036,
#'   1.335, -0.290, -0.083, -0.047, -0.356, 0.934, -0.152, 5.201, 5.883,
#'   3.560, 9.799, 0.368)
#'  mod222c <- GMVAR(data, p=2, M=2, params=params222c, constraints=C_mat)
#'  sim222c <- simulateGMVAR(mod222c, nsimu=100,
#'               init_values=matrix(c(30, 30, 80, 80), nrow=2))
#'  plot.ts(sim222c$sample)
#'  ts.plot(sim222c$mixing_weights, col=c("blue", "red"), lty=2)
#'  plot(sim222c$component, type="l")
#' @export

simulateGMVAR <- function(gmvar, nsimu, init_values=NULL, ntimes=1, drop=TRUE, seed=NULL, girf_pars=NULL) {
  # girf_pars$shock_numb - which shock?
  # girf_pars$shock_size - size of the structural shock?
  # girf_pars$init_regimes - init values generated from which regimes? Ignored if !is.null(init_values)
  # girf_pars$include_mixweights - should GIRFs be estimated for the mixing weights as well? TRUE or FALSE
  # If !is.null(girf_pars) and girf_pars$include_mixweights == TRUE, returns a size (N+1 x d+M) vector containing
  #                                                                  the estimated GIRFs for the variables and
  #                                                                  and the mixing weights (column d+m for the m:th regime).
  # If !is.null(girf_pars) and girf_pars$include_mixweights == FALSE, returns a size (N+1 x d) vector containing
  #                                                                   the estimated GIRFs for the variables only.
  # The first row for response at impact
  if(!is.null(seed)) set.seed(seed)

  check_gmvar(gmvar)
  epsilon <- round(log(.Machine$double.xmin) + 10)
  p <- gmvar$model$p
  M <- gmvar$model$M
  d <- gmvar$model$d
  constraints <- gmvar$model$constraints
  same_means <- gmvar$model$same_means
  structural_pars <- gmvar$model$structural_pars
  if(!all_pos_ints(c(nsimu, ntimes))) stop("Arguments n and ntimes must be positive integers")
  if(!is.null(init_values)) {
    if(!is.matrix(init_values)) stop("init_values must be numeric matrix")
    if(anyNA(init_values)) stop("init_values contains NA values")
    if(ncol(init_values) != d | nrow(init_values) < p) stop("init_values must contain d columns and at least p rows")
  }

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

  all_Bm <- array(dim=c(d, d, M)) # Matrices such that all_Bm[, , m]%*%rnorm(d) follows N(0, \Omega_m)
  if(!is.null(structural_pars)) { # Calculate  W%*%Lambda_m^{1/2} where Lambda_1 = I_d (not to be confused with the time-varying B-matrix)
    W <- pick_W(p=p, M=M, d=d, params=params, structural_pars=structural_pars)
    all_Bm[, , 1] <- W
    if(M > 1) {
      lambdas <- matrix(pick_lambdas(p=p, M=M, d=d, params=params, structural_pars=structural_pars), nrow=d, byrow=FALSE)
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
  Sigmas <- get_Sigmas(p=p, M=M, d=d, all_A=all_A, all_boldA=all_boldA, all_Omega=all_Omega) # Store the (dpxdp) covariance matrices
  inv_Sigmas <- array(NA, dim=c(d*p, d*p, M)) # Store inverses of the (dpxdp) covariance matrices
  det_Sigmas <- numeric(M) # Store determinants of the (dpxdp) covariance matrices
  for(m in 1:M) {
    inv_Sigmas[, , m] <- solve(Sigmas[, , m])
    det_Sigmas[m] <- det(Sigmas[, , m])
  }

  if(is.null(girf_pars)) {
    init_regimes <- 1:M
    reg_probs <- alphas
  } else {
    R1 <- girf_pars$R1
    init_regimes <- girf_pars$init_regimes
    reg_probs <- alphas[init_regimes]
  }

  # Set/generate initial values
  if(is.null(init_values)) {
    m <- ifelse(length(init_regimes) == 1 && init_regimes != 1, init_regimes, sample(x=init_regimes, size=1, replace=TRUE, prob=reg_probs)) # From which mixture component the initial values are drawn from?
    mu <- rep(all_mu[, m], p)
    L <- t(chol(Sigmas[, , m]))
    init_values <- matrix(mu + L%*%rnorm(d*p), nrow=p, ncol=d, byrow=TRUE)
  } else {
    init_values <- init_values[(nrow(init_values) - p + 1):nrow(init_values), , drop=FALSE]
  }

  # Container for the simulated values and initial values. First row row initial values vector, and t:th row for (y_{i-1}',...,y_{i-p}')
  Y <- matrix(nrow=nsimu + 1, ncol=d*p)
  Y[1,] <- reform_data(init_values, p=p)
  if(!is.null(girf_pars)) Y2 <- Y # Storage for the second sample path in the GIRF algorithm

  # Initialize data structures
  sample <- array(dim=c(nsimu, d, ntimes))
  component <- matrix(nrow=nsimu, ncol=ntimes)
  mixing_weights <- array(dim=c(nsimu, M, ntimes))
  if(!is.null(girf_pars)) {
    sample2 <- array(dim=c(nsimu, d, ntimes))
    mixing_weights2 <- array(dim=c(nsimu, M, ntimes))
  }

  # Some functions to be used
  get_matprods <- function(Y) vapply(1:M, function(m) crossprod(Y[i1,] - rep(all_mu[, m], p), inv_Sigmas[, , m])%*%(Y[i1,] - rep(all_mu[, m], p)), numeric(1))
  get_logmvnvalues <- function(matprods) vapply(1:M, function(m) -0.5*d*p*log(2*pi) - 0.5*log(det_Sigmas[m]) - 0.5*matprods[m], numeric(1))

  for(j1 in seq_len(ntimes)) {
    for(i1 in seq_len(nsimu)) {
      # Calculate the dp-dimensional multinormal densities (KMS 2016, eq.(6)).
      # Calculated in logarithm because same values may be too close to zero for machine accuracy.
      matprods <- get_matprods(Y)
      log_mvnvalues <- get_logmvnvalues(matprods)
      alpha_mt <- get_alpha_mt(M=M, log_mvnvalues=log_mvnvalues, alphas=alphas,
                               epsilon=epsilon, also_l_0=FALSE)

      # Draw a mixture component and store the values
      m <- sample.int(n=M, size=1, replace=TRUE, prob=alpha_mt)
      component[i1, j1] <- m
      mixing_weights[i1, , j1] <- alpha_mt

      # Calculate conditional mean for regime m
      A2 <- matrix(all_A[, , , m], nrow=d, byrow=FALSE) # (A_1:...:A_p)
      mu_mt <- all_phi0[, m] + A2%*%Y[i1,]

      # Draw the sample and store it
      eps_t <- rnorm(d)
      sample[i1, , j1] <- mu_mt + all_Bm[, , m]%*%eps_t

      # Update storage matrix Y (overwrites when ntimes > 1)
      if(p == 1) {
        Y[i1 + 1,] <- sample[i1, , j1]
      } else {
        Y[i1 + 1,] <- c(sample[i1, , j1], Y[i1 , 1:(d*p - d)])
      }

      # For the second sample path in GIRF (with a specific structural shock occurring)
      if(!is.null(girf_pars)) {
        matprods2 <- get_matprods(Y2)
        log_mvnvalues2 <- get_logmvnvalues(matprods2)
        alpha_mt2 <- get_alpha_mt(M=M, log_mvnvalues=log_mvnvalues2, alphas=alphas,
                                  epsilon=epsilon, also_l_0=FALSE)
        mixing_weights2[i1, , j1] <- alpha_mt2

        if(i1 == 1) {
          m2 <- m # Common regime at impact (the mixing weights are the same)
        } else {
          m2 <- sample.int(n=M, size=1, replace=TRUE, prob=alpha_mt2)
        }

        A22 <- matrix(all_A[, , , m2], nrow=d, byrow=FALSE) # (A_1:...:A_p)
        mu_mt2 <- all_phi0[, m2] + A22%*%Y2[i1,]

        u_t <- all_Bm[, , m2]%*%eps_t # Reduced form shock from the same NID(0,I) shock

        if(i1 == 1) { # At impact, obtain reduced form shock from the specific structural shock

          # Calculate the time-varying B-matrix
          if(M == 1) {
            B_t <- W
          } else {
            tmp <- array(dim=c(d, d, M))
            tmp[, , 1] <- alpha_mt2[1]*diag(d)
            for(m in 2:M) {
              tmp[, , m] <- alpha_mt2[m]*diag(lambdas[, m - 1])
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
    if(!is.null(gmvar$data)) {
      colnames(one_girf) <- colnames(gmvar$data)
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
    sample <- matrix(sample, nrow=nsimu, ncol=d)
    component <- as.vector(component)
    mixing_weights <- matrix(mixing_weights, nrow=nsimu, ncol=M)
  }

  list(sample=sample,
       component=component,
       mixing_weights=mixing_weights)
}
