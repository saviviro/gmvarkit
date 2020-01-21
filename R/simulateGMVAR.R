

#' @title Simulate from GMVAR process
#'
#' @description \code{simulateGMVAR} simulates observations from a GMVAR process.
#'
#' @param gmvar an object of class \code{'gmvar'} created with \code{fitGMVAR} or \code{GMVAR}.
#' @param nsimu number of observations to be simulated.
#' @param init_values a size \eqn{(pxd)} matrix specifying the initial values to be used in the simulation, where
#'   d is the number of time series in the system.
#'   The \strong{last} row will be used as initial values for the first lag, the second last row for second lag etc. If not
#'   specified, initial values will be drawn from the stationary distribution.
#' @param ntimes how many sets of simulations should be performed?
#' @param drop if \code{TRUE} (default) then the components of the returned list are coerced to lower dimension if \code{ntimes==1}, i.e.,
#'   \code{$sample} and \code{$mixing_weights} will be matrices, and \code{$component} will be vector.
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
#' @inherit in_paramspace_int references
#' @examples
#'  \donttest{
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
#'  ## FORECASTING EXAMPLE
#'  # Forecast 5-steps-ahead, 10000 sets of simulations with initial
#'  # values from the data:
#'  # GMVAR(2,2), d=2 model with mean-parametrization:
#'  params222 <- c(-11.904, 154.684, 1.314, 0.145, 0.094, 1.292, -0.389,
#'   -0.070, -0.109, -0.281, 0.920, -0.025, 4.839, 11.633, 124.983, 1.248,
#'    0.077, -0.040, 1.266, -0.272, -0.074, 0.034, -0.313, 5.855, 3.570,
#'    9.838, 0.740)
#'  mod222 <- GMVAR(data, p=2, M=2, params=params222, parametrization="mean")
#'  sim222 <- simulateGMVAR(mod222, nsimu=5, ntimes=10000)
#'
#'  # Point forecast + 95% prediction intervals:
#'  apply(sim222$sample, 1:2, quantile, probs=c(0.025, 0.5, 0.972))
#'
#'  # Similar forecast for the mixing weights:
#'  apply(sim222$mixing_weights, 1:2, quantile, probs=c(0.025, 0.5, 0.972))
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
#'  }
#' @export

simulateGMVAR <- function(gmvar, nsimu, init_values=NULL, ntimes=1, drop=TRUE) {
  check_gmvar(gmvar)
  epsilon <- round(log(.Machine$double.xmin) + 10)
  p <- gmvar$model$p
  M <- gmvar$model$M
  d <- gmvar$model$d
  constraints <- gmvar$model$constraints
  if(!all_pos_ints(c(nsimu, ntimes))) stop("Arguments n and ntimes must be positive integers")
  if(!is.null(init_values)) {
    if(!is.matrix(init_values)) stop("init_values must be numeric matrix")
    if(anyNA(init_values)) stop("init_values contains NA values")
    if(ncol(init_values) != d | nrow(init_values) < p) stop("init_values must contain d columns and at least p rows")
  }

  # Collect parameter values
  params <- gmvar$params
  if(gmvar$model$parametrization == "mean") {
    params <- change_parametrization(p=p, M=M, d=d, params=params, constraints=constraints,
                                     change_to="intercept")
  }
  params <- reform_constrained_pars(p=p, M=M, d=d, params=params, constraints=constraints)

  all_mu <- get_regime_means(gmvar)
  all_phi0 <- pick_phi0(p=p, M=M, d=d, params=params)
  all_A <- pick_allA(p=p, M=M, d=d, params=params)
  all_Omega <- pick_Omegas(p=p, M=M, d=d, params=params)
  all_boldA <- form_boldA(p=p, M=M, d=d, all_A=all_A)
  alphas <- pick_alphas(p=p, M=M, d=d, params=params)

  # Calculate the covariance matrices Sigma_{m,p} (Lutkepohl 2005, eq. (2.1.39))
  I_dp2 <- diag(nrow=(d*p)^2)
  ZER_lower <- matrix(0, nrow=d*(p-1), ncol=d*p)
  ZER_right <- matrix(0, nrow=d, ncol=d*(p-1))
  Sigmas <- array(NA, dim=c(d*p, d*p, M)) # Store the (dpxdp) covariance matrices
  inv_Sigmas <- array(NA, dim=c(d*p, d*p, M)) # Store inverses of the (dpxdp) covariance matrices
  det_Sigmas <- numeric(M) # Store determinants of the (dpxdp) covariance matrices
  for(m in 1:M) {
    kronmat <- I_dp2 - kronecker(all_boldA[, , m], all_boldA[, , m])
    sigma_epsm <- rbind(cbind(all_Omega[, , m], ZER_right), ZER_lower)
    Sigma_m <- unvec(d=d*p, a=solve(kronmat, vec(sigma_epsm)))
    Sigmas[, , m] <- Sigma_m
    inv_Sigmas[, , m] <- solve(Sigma_m)
    det_Sigmas[m] <- det(Sigma_m)
  }

  # Set/generate initial values
  if(is.null(init_values)) {
    m <- sample.int(n=M, size=1, replace=TRUE, prob=alphas) # From which mixture component the initial values are drawn from?
    mu <- rep(all_mu[, m], p)
    L <- t(chol(Sigmas[, , m]))
    init_values <- matrix(mu + L%*%rnorm(d*p), nrow=p, ncol=d, byrow=TRUE)
  } else {
    init_values <- init_values[(nrow(init_values) - p + 1):nrow(init_values), , drop=FALSE]
  }

  # Container for the simulated values and initial values. First row row initival values vector, and t:th row for (y_{i-1}',...,y_{i-p}')
  Y <- matrix(nrow=nsimu + 1, ncol=d*p)
  Y[1,] <- reform_data(init_values, p=p)

  # Initialize data structures
  sample <- array(dim=c(nsimu, d, ntimes))
  component <- matrix(nrow=nsimu, ncol=ntimes)
  mixing_weights <- array(dim=c(nsimu, M, ntimes))

  for(j1 in seq_len(ntimes)) {
    for(i1 in seq_len(nsimu)) {
      # Calculate the dp-dimensional multinormal densities (KMS 2016, eq.(6)).
      # Calculated in logarithm because same values may be too close to zero for machine accuracy.
      matprods <- vapply(1:M, function(m) crossprod(Y[i1,] - rep(all_mu[, m], p), inv_Sigmas[, , m])%*%(Y[i1,] - rep(all_mu[, m], p)), numeric(1))
      log_mvnvalues <- vapply(1:M, function(m) -0.5*d*p*log(2*pi) - 0.5*log(det_Sigmas[m]) - 0.5*matprods[m], numeric(1))

      # Calculate mixing weights
      if(M == 1) {
        alpha_mt <- 1
      } else if(any(log_mvnvalues < epsilon)) { # If some values are too close to zero use the package Brobdingnag
        numerators <- lapply(1:M, function(m) alphas[m]*exp(Brobdingnag::as.brob(log_mvnvalues[m])))
        denominator <- Reduce('+', numerators)
        alpha_mt <- vapply(1:M, function(m) as.numeric(numerators[[m]]/denominator), numeric(1))
      } else {
        mvnvalues <- exp(log_mvnvalues)
        denominator <- as.vector(mvnvalues%*%alphas)
        alpha_mt <- alphas*(mvnvalues/denominator)
      }

      # Draw a mixture component and store the values
      m <- sample.int(n=M, size=1, replace=TRUE, prob=alpha_mt)
      component[i1, j1] <- m
      mixing_weights[i1, , j1] <- alpha_mt

      # Calculate conditional mean for regime m
      A2 <- matrix(all_A[, , , m], nrow=d, byrow=FALSE) # (A_1:...:A_p)
      mu_mt <- all_phi0[, m] + A2%*%Y[i1,]

      # Draw the sample and store it
      L <- t(chol(all_Omega[, , m])) # Lower triangular cholesky of error term covariance matrix
      sample[i1, , j1] <- mu_mt + L%*%rnorm(d)

      # Update storage matrix Y (overwrites when ntimes > 1)
      if(p == 1) {
        Y[i1 + 1,] <- sample[i1, , j1]
      } else {
        Y[i1 + 1,] <- c(sample[i1, , j1], Y[i1 , 1:(d*p - d)])
      }
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
