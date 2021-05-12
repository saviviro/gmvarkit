
#' @title Estimate generalized impulse response function for a structural GMVAR model.
#'
#' @description \code{GIRF} estimates generalized impulse response function for
#'  a structural GMVAR model.
#'
#' @inheritParams simulateGMVAR
#' @inheritParams fitGMVAR
#' @param which_shocks a numeric vector of length at most \eqn{d} (\code{=ncol(data)})
#'   and elements in \eqn{1,...,d} specifying the structural shocks for which the GIRF
#'   should be estimated.
#' @param shock_size a scalar value specifying the common size for all scalar components of the structural shock.
#'   Note that the conditional covariance matrix of the structural shock is an identity matrix
#'   and that the (generalized) impulse responses may not be symmetric to the sign and size of
#'   the shock.
#' @param N a positive integer specifying the horizon how far ahead should the generalized
#'   impulse responses be calculated.
#' @param R1 the number of repetitions used to estimate GIRF for each initial value.
#' @param R2 the number of initial values to be drawn from a stationary distribution of the process
#'   or of a specific regime? The confidence bounds will be sample quantiles of the GIRFs based on
#'   different initial values. Ignored if the argument \code{init_value} is specified.
#' @param init_regimes a numeric vector of length at most \eqn{M} and elements in \eqn{1,...,M}
#'   specifying the regimes from which the initial values should be generated from. The initial
#'   values will be generated from a mixture distribution with the mixture components being the
#'   stationary distributions of the specific regimes and the (proportional) mixing weights given
#'   by the mixing weight parameters of those regimes. Note that if \code{init_regimes=1:M}, the
#'   initial values are generated from the stationary distribution of the process and if
#'   \code{init_regimes=m}, the initial value are generated from the stationary distribution
#'   of the \eqn{m}th regime. Ignored if \code{init_value} is specified.
#' @param init_values a matrix or a multivariate class \code{'ts'} object with \eqn{d} columns
#'   and at least \eqn{p} rows specifying an initial value for the GIRF. The last \eqn{p} rows
#'   are taken to be the initial value assuming that the \strong{last} row is the most recent observation.
#' @param which_cumulative a numeric vector with values in \eqn{1,...,d} (\code{d=ncol(data)}) specifying
#'   which the variables for which the impulse responses should be cumulative. Default is none.
#' @param ci a numeric vector with elements in \eqn{(0, 1)} specifying the confidence levels of the
#'   confidence intervals.
#' @param include_mixweights should the generalized impulse response be calculated for the mixing weights
#'   as well? \code{TRUE} or \code{FALSE}.
#' @param ncores the number CPU cores to be used in parallel computing. Only single core computing is
#'   supported if an initial value is specified (and the GIRF won't thus be estimated multiple times).
#' @param plot \code{TRUE} if the results should be plotted, \code{FALSE} if not.
#' @param seeds a length \code{R2} vector containing the random number generator seed for estimation
#'   of each GIRF. A single number of an initial value is specified.
#'  or \code{NULL} for not initializing the seed. Exists for creating reproducible results.
#' @details The model needs to be structural in order for this function to be applicable. A structural
#'   GMVAR model can be estimated by specifying the argument \code{structural_pars} in the function \code{fitGMVAR}.
#'
#'   The confidence bounds reflect uncertainty about the initial state (but currently not about the parameter
#'   estimates) if initial values are not specified. If initial values are specified, there won't currently
#'   be confidence intervals. See the cited paper by Virolainen (2020) for details about the algorithm.
#' @return Returns a class \code{'girf'} list with the GIRFs in the first element (\code{$girf_res}) and the used
#'   arguments the rest. The first element containing the GIRFs is a list with the \eqn{m}th element containing
#'   the point estimates for the GIRF in \code{$point_est} (the first element) and confidence intervals in
#'   \code{$conf_ints} (the second element). The first row is for the GIRF at impact \eqn{(n=0)}, the second for
#'   \eqn{n=1}, the third for \eqn{n=2}, and so on.
#' @seealso \code{\link{GFEVD}}, \code{\link{fitGMVAR}}, \code{\link{GMVAR}}, \code{\link{gmvar_to_sgmvar}}, \code{\link{reorder_W_columns}},
#'  \code{\link{swap_W_signs}}, \code{\link{simulateGMVAR}}, \code{\link{predict.gmvar}},
#'  \code{\link{profile_logliks}}, \code{\link{quantile_residual_tests}}, \code{\link{LR_test}}, \code{\link{Wald_test}}
#' @inherit in_paramspace_int references
#' @examples
#'  \donttest{
#'  # These are long-running examples that use parallel computing.
#'  # It takes approximately 30 seconds to run all the below examples.
#'
#'  data(eurusd, package="gmvarkit")
#'  data <- cbind(10*eurusd[,1], 100*eurusd[,2])
#'  colnames(data) <- colnames(eurusd)
#'
#'  # Structural GMVAR(2, 2), d=2 model identified with sign-constraints:
#'  params22s <- c(1.386, -0.766, 1.005, 5.928, 1.314, 0.145, 0.094, 1.292,
#'   -0.389, -0.07, -0.109, -0.281, 1.248, 0.077, -0.04, 1.266, -0.272, -0.074,
#'    0.034, -0.313, 0.903, 0.718, -0.324, 2.079, 7.001, 1.44, 0.741)
#'  W_22 <- matrix(c(1, 1, -1, 1), nrow=2, byrow=FALSE)
#'  mod22s <- GMVAR(data, p=2, M=2, params=params22s,
#'   structural_pars=list(W=W_22))
#'  mod22s
#'  # Alternatively, use:
#'  # fit22s <- fitGMVAR(data, p=2, M=2, structural_pars=list(W=W_22),
#'  #                    ncalls=40, seeds=1:40)
#'  # To obtain an estimated version of the same model.
#'
#'  # Estimating the GIRFs of both structural shocks with default arguments
#'  # (initial values are drawn from the stationary distribution of the process,
#'  # 30 periods ahead, confidence levels 0.95 and 0.8):
#'  girf1 <- GIRF(mod22s, N=12, R1=100, R2=100)
#'  girf1
#'  plot(girf1)
#'
#'  # Estimating the GIRF of the second shock only, 36 periods ahead
#'  # and shock size 1, initial values drawn from the stationary distribution
#'  # of the first regime, confidence level 0.9:
#'  girf2 <- GIRF(mod22s, which_shocks=2, shock_size=1, N=12, init_regimes=1,
#'                ci=0.9, R1=100, R2=100)
#'  plot(girf2)
#'
#'  # Estimating the GIRFs of both structural shocks, shock sizes 1 and 3, N=20
#'  # periods ahead, estimation based on 200 Monte Carlo simulations, and fixed
#'  # initial values given by the last p observations of the data:
#'  girf3 <- GIRF(mod22s, shock_size=c(1, 3), N=20, R1=200,
#'                init_values=mod22s$data)
#'  plot(girf3)
#'  }
#' @export

GIRF <- function(gmvar, which_shocks, shock_size=1, N=30, R1=250, R2=250, init_regimes=1:gmvar$model$M, init_values=NULL,
                 which_cumulative=numeric(0), ci=c(0.95, 0.80), include_mixweights=TRUE, ncores=2,
                 plot=TRUE, seeds=NULL) {
  p <- gmvar$model$p
  M <- gmvar$model$M
  d <- gmvar$model$d
  if(is.null(gmvar$model$structural_pars)) stop("Only structural models are supported")
  if(M == 1) include_mixweights <- FALSE
  if(missing(which_shocks)) {
    which_shocks <- 1:d
  } else {
    stopifnot(all(which_shocks %in% 1:d))
    which_shocks <- unique(which_shocks)
  }
  if(!is.null(init_values)) R2 <- 1
  if(!is.null(seeds) && length(seeds) != R2) stop("The argument 'seeds' needs be NULL or a vector of length 'R2'")
  stopifnot(all(init_regimes %in% 1:M))
  init_regimes <- unique(init_regimes)
  stopifnot(length(ci) > 0 && all(ci > 0 & ci < 1))
  if(length(shock_size) != 1) {
    warning("The argument shock_size should be a numeric scalar. Using the first value.")
    shock_size <- shock_size[1]
  }
  if(length(which_cumulative) > 0) {
    which_cumulative <- unique(which_cumulative)
    stopifnot(all(which_cumulative %in% 1:d))
  }

  # Function that estimates GIRF
  get_one_girf <- function(shock_numb, shock_size, seed) {
    simulateGMVAR(gmvar, nsimu=N + 1, init_values=init_values, ntimes=R1, seed=seed, girf_pars=list(shock_numb=shock_numb,
                                                                                                    shock_size=shock_size,
                                                                                                    init_regimes=init_regimes,
                                                                                                    include_mixweights=include_mixweights))
  }

  if(ncores > parallel::detectCores()) {
    ncores <- parallel::detectCores()
    message("ncores was set to be larger than the number of cores detected")
  }
  if(is.null(init_values)) {
    cat(paste("Using", ncores, "cores to estimate", R2,"GIRFs for", length(which_shocks), "structural shocks,", "each based on", R1, "Monte Carlo repetitions."), "\n")
  } else {
    cat(paste("Using", ncores, "cores to estimate one GIRF for", length(which_shocks), "structural shocks, each based on", R1, "Monte Carlo repetitions."), "\n")
  }

  ### Calculate the GIRFs ###
  cl <- parallel::makeCluster(ncores)
  on.exit(try(parallel::stopCluster(cl), silent=TRUE)) # Close the cluster on exit, if not already closed.
  parallel::clusterExport(cl, ls(environment(GIRF)), envir = environment(GIRF)) # assign all variables from package:gmvarkit
  parallel::clusterEvalQ(cl, c(library(mvnfast), library(pbapply)))

  GIRF_shocks <- vector("list", length=length(which_shocks))

  for(i1 in 1:length(which_shocks)) {
    cat(paste0("Estimating GIRFs for structural shock ", which_shocks[i1], "..."), "\n")
    GIRF_shocks[[i1]] <- pbapply::pblapply(1:R2, function(i2) get_one_girf(shock_numb=which_shocks[i1], shock_size=shock_size, seed=seeds[i2]), cl=cl)
  }
  parallel::stopCluster(cl=cl)

  GIRF_results <- vector("list", length=length(which_shocks))
  names(GIRF_results) <- paste0("shock", which_shocks)

  for(i1 in 1:length(which_shocks)) {
    res_in_array <- array(unlist(GIRF_shocks[[i1]]), dim=c(N + 1, d + ifelse(include_mixweights, M, 0), R2))
    if(length(which_cumulative) > 0) {
      for(i2 in which_cumulative) {
        res_in_array[, i2, ] <- apply(res_in_array[, i2, , drop=FALSE], MARGIN=3, FUN=cumsum) # Replace GIRF with cumulative GIRF
      }
    }
    colnames(res_in_array) <- colnames(GIRF_shocks[[1]][[1]])
    point_estimate <- apply(X=res_in_array, MARGIN=1:2, FUN=mean)
    lower <- (1 - ci)/2
    upper <- rev(1 - lower)
    q_tocalc <- c(lower, upper)
    q_tocalc <- sort(q_tocalc, decreasing=FALSE)
    conf_ints <- apply(res_in_array, MARGIN=1:2, FUN=quantile, probs=q_tocalc)
    conf_ints <- aperm(conf_ints, perm=c(2, 1, 3))
    rownames(point_estimate) <- rownames(conf_ints) <- 0:N
    GIRF_results[[i1]] <- list(point_est=point_estimate,
                               conf_ints=conf_ints)
  }

  cat("Finished!\n")
  ret <- structure(list(girf_res=GIRF_results,
                        shocks=which_shocks,
                        shock_size=shock_size,
                        N=N,
                        R1=R1,
                        R2=R2,
                        ci=ci,
                        which_cumulative=which_cumulative,
                        init_regimes=init_regimes,
                        init_values=init_values,
                        include_mixweights=include_mixweights,
                        seeds=seeds,
                        gmvar=gmvar),
                   class="girf")
  if(plot) plot(ret)
  ret
}



#' @title Estimate generalized forecast error variance decomposition for a structural GMVAR model.
#'
#' @description \code{GFEVD} estimates generalized generalized forecast error variance decomposition for
#'  a structural GMVAR model.
#'
#' @inheritParams GIRF
#' @param shock_size What shocks size should be used for all shocks? By the definition of the SGMVAR model,
#'   the conditional covariance matrix of the structural shock is identity matrix.
#' @param N a positive integer specifying the horizon how far ahead should the GFEVD be calculated.
#' @param initval_type What type initial values are used for estimating the GIRFs that the GFEVD is based on?
#'   \describe{
#'     \item{\code{"data"}:}{Estimate the GIRF for all the possible length \eqn{p} histories in the data.}
#'     \item{\code{"random"}:}{Estimate the GIRF for several random initial values generated from the stationary
#'        distribution of the process or from the stationary distribution of specific regime(s) chosen with the
#'        argument \code{init_regimes}. The number of initial values is set with the argument \code{R2}.}
#'     \item{\code{"fixed"}:}{Estimate the GIRF for a fixed initial value only, which is specified with the argument
#'        \code{init_values}.}
#'   }
#' @param R2 the number of initial values to be drawn if \code{initval_type="random"}.
#' @param include_mixweights should the GFEVD be estimated for the mixing weights as well? Note that this is
#'   ignored if \code{M=1} and if \code{M=2} the GFEVD will be the same for both regime's mixing weights.
#' @param seeds a numeric vector containing the random number generator seed for estimation
#'   of each GIRF. Should have the length...
#'   \itemize{
#'     \item ...\code{nrow(data) - p + 1} if \code{initval_type="data"}.
#'     \item ...\code{R2} if \code{initval_type="random"}.
#'     \item ...\code{1} if \code{initval_type="fixed."}.
#'   }
#'   Set to \code{NULL} for not initializing the seed. Exists for creating reproducible results.
#' @details The model needs to be structural in order for this function to be applicable. A structural
#'   GMVAR model can be estimated by specifying the argument \code{structural_pars} in the function \code{fitGMVAR}.
#'
#'   The GFEVD is a forecast error variance decomposition calculated with the generalized impulse response function (GIRF).
#'   Lanne and Nyberg (2016) for details. Note, however, that the related GIRFs are calculated using the algorithm given in
#'   Virolainen (2020).
#' @return Returns and object of class 'gfevd' containing the GFEVD for all the variables and if
#'   \code{include_mixweights=TRUE} also to the mixing weights. Note that the decomposition does not
#'   exist at horizon zero for mixing weights because the related GIRFs are always zero at impact.
#' @seealso \code{\link{GIRF}}, \code{\link{fitGMVAR}}, \code{\link{GMVAR}}, \code{\link{gmvar_to_sgmvar}},
#'  \code{\link{reorder_W_columns}}, \code{\link{swap_W_signs}}, \code{\link{simulateGMVAR}}
#' @references
#' @references
#'  \itemize{
#'    \item Lanne M. and Nyberg H. 2016. Generalized Forecast Error Variance Decomposition for Linear
#'      and Nonlineae Multivariate Models. \emph{Oxford Bulletin of Economics and Statistics}, \strong{78}, 4, 595-603.
#'    \item Kalliovirta L., Meitz M. and Saikkonen P. 2016. Gaussian mixture vector autoregression.
#'          \emph{Journal of Econometrics}, \strong{192}, 485-498.
#'    \item Virolainen S. 2020. Structural Gaussian mixture vector autoregressive model. Unpublished working
#'      paper, available as arXiv:2007.04713.
#'  }
#' @examples
#'  \donttest{
#'  # These are long-running examples that use parallel computing.
#'  # It takes approximately 30 seconds to run all the below examples.
#'
#'  data(eurusd, package="gmvarkit")
#'  data <- cbind(10*eurusd[,1], 100*eurusd[,2])
#'  colnames(data) <- colnames(eurusd)
#'
#'  # Structural GMVAR(2, 2), d=2 model identified with sign-constraints:
#'  params22s <- c(1.386, -0.766, 1.005, 5.928, 1.314, 0.145, 0.094, 1.292,
#'   -0.389, -0.07, -0.109, -0.281, 1.248, 0.077, -0.04, 1.266, -0.272, -0.074,
#'    0.034, -0.313, 0.903, 0.718, -0.324, 2.079, 7.001, 1.44, 0.741)
#'  W_22 <- matrix(c(1, 1, -1, 1), nrow=2, byrow=FALSE)
#'  mod22s <- GMVAR(data, p=2, M=2, params=params22s,
#'   structural_pars=list(W=W_22))
#'  mod22s
#'
#'  ## NOTE: Use larger R1 is empirical applications! Small R1 is used
#'  ## Below only to fasten the execution time of the examples.
#'
#'  # Estimating the GFEVD using all possible histories in the data as the
#'  # initial values:
#'  gfevd1 <- GFEVD(mod22s, N=24, R1=20, initval_type="data")
#'  gfevd1
#'  plot(gfevd1)
#'
#'  # Estimate GFEVD with the initial values generated from the stationary
#'  # distribution of the process:
#'  gfevd2 <- GFEVD(mod22s, N=24, R1=20, R2=100, initval_type="random")
#'  gfevd2
#'  plot(gfevd2)
#'
#'  # Estimate GFEVD with fixed hand specified initial values. We use the
#'  # unconditional mean of the process:
#'  myvals <- rbind(mod22s$uncond_moments$uncond_mean,
#'                  mod22s$uncond_moments$uncond_mean)
#'  gfevd3 <- GFEVD(mod22s, N=36, R1=50, initval_type="fixed",
#'   init_values=myvals, include_mixweights=TRUE)
#'  gfevd3
#'  plot(gfevd3)
#'  }
#' @export

GFEVD <- function(gmvar, shock_size=1, N=30, initval_type=c("data", "random", "fixed"), R1=250, R2=250,
                  init_regimes=NULL, init_values=NULL, which_cumulative=numeric(0), include_mixweights=FALSE,
                  ncores=2, seeds=NULL) {
  initval_type <- match.arg(initval_type)
  p <- gmvar$model$p
  M <- gmvar$model$M
  d <- gmvar$model$d
  if(is.null(gmvar$model$structural_pars)) stop("Only structural models are supported")
  if(M == 1) include_mixweights <- FALSE
  if(initval_type == "data") {
    if(is.null(gmvar$data)) stop("The model does not contain data! Add data with the function 'add_data' or select another 'initval_type'.")
    stopifnot(nrow(gmvar$data) >= p)
    R2 <- nrow(gmvar$data) - p + 1 # The number of length p histories in the data
    all_initvals <- array(vapply(1:R2, function(i1) gmvar$data[i1:(i1 + p - 1),], numeric(p*d)), dim=c(p, d, R2)) # [, , i1] for i1:th initval
  } else if(initval_type == "random") {
    if(is.null(init_regimes)) {
      message("Initial regimes were not specified, so the initial values are generated from the stationary distribution of the process.")
      init_regimes <- 1:M
    } else {
      stopifnot(all(init_regimes %in% 1:M))
      init_regimes <- unique(init_regimes)
    }
    all_initvals <- array(dim=c(1, 1, R2)) # This won't be used. NULL init_values will be used.
  } else if(initval_type == "fixed") {
    R2 <- 1
    if(is.null(init_values)) stop("Initial values were not specified! Specify the initial values with the argument 'init_values' or choose another 'initval_type'.")
    if(!is.matrix(init_values)) stop("init_values must be a numeric matrix")
    if(anyNA(init_values)) stop("init_values contains NA values")
    if(ncol(init_values) != d | nrow(init_values) < p) stop("init_values must contain d columns and at least p rows")
    init_values <- init_values[(nrow(init_values) - p + 1):nrow(init_values), , drop=FALSE]
    all_initvals <- array(init_values, dim=c(p, d, R2))
  }
  if(!is.null(seeds) && length(seeds) != R2) stop("The argument 'seeds' has wrong length!")
  stopifnot(length(shock_size) == 1)
  if(length(which_cumulative) > 0) {
    which_cumulative <- unique(which_cumulative)
    stopifnot(all(which_cumulative %in% 1:d))
  }

  # Function that estimates GIRF
  get_one_girf <- function(shock_numb, shock_size, seed, init_values_for_1girf) {
    if(initval_type == "random") init_values_for_1girf <- NULL
    simulateGMVAR(gmvar, nsimu=N + 1, init_values=init_values_for_1girf, ntimes=R1, seed=seed,
                  girf_pars=list(shock_numb=shock_numb,
                                 shock_size=shock_size,
                                 init_regimes=init_regimes,
                                 include_mixweights=include_mixweights))
  }

  if(ncores > parallel::detectCores()) {
    ncores <- parallel::detectCores()
    message("ncores was set to be larger than the number of cores detected")
  }
  if(initval_type != "fixed") {
    cat(paste("Using", ncores, "cores to estimate", R2,"GIRFs for", d, "structural shocks,", "each based on", R1, "Monte Carlo repetitions."), "\n")
  } else {
    cat(paste("Using", ncores, "cores to estimate one GIRF for", d, "structural shocks, each based on", R1, "Monte Carlo repetitions."), "\n")
  }

  ### Calculate the GIRFs ###
  cl <- parallel::makeCluster(ncores)
  on.exit(try(parallel::stopCluster(cl), silent=TRUE)) # Close the cluster on exit, if not already closed.
  parallel::clusterExport(cl, ls(environment(GFEVD)), envir = environment(GFEVD)) # assign all variables from package:gmvarkit
  parallel::clusterEvalQ(cl, c(library(mvnfast), library(pbapply)))

  GIRF_shocks <- vector("list", length=d)
  for(i1 in 1:d) {
    cat(paste0("Estimating GIRFs for structural shock ", i1, "..."), "\n")
    GIRF_shocks[[i1]] <- pbapply::pblapply(1:R2, function(i2) get_one_girf(init_values=matrix(all_initvals[, , i2], nrow=p, ncol=d),
                                                                           shock_numb=i1,
                                                                           shock_size=shock_size,
                                                                           seed=seeds[i2]), cl=cl)
  }
  parallel::stopCluster(cl=cl)

  if(is.null(colnames(gmvar$data))) {
    varnames <- paste0("Variable", 1:d)
  } else {
    varnames <- colnames(gmvar$data)
  }
  if(include_mixweights) {
    varnames <- c(varnames, paste0("mw", 1:M))
  }
  shocknames <- paste0("Shock", 1:d)

  GIRF_square_cumsum <- array(dim=c(N + 1, length(varnames), d),
                              dimnames=list(0:N, varnames, shocknames)) # [horizon, variable, shock]
  for(i1 in 1:d) { # Go through the shocks
    res_in_array <- array(unlist(GIRF_shocks[[i1]]), dim=c(N + 1, length(varnames), R2))
    if(length(which_cumulative) > 0) {
      for(i2 in which_cumulative) {
        res_in_array[, i2, ] <- apply(res_in_array[, i2, , drop=FALSE], MARGIN=3, FUN=cumsum) # Replace GIRF with cumulative GIRF
      }
    }
    GIRF_square_cumsum[, , i1] <-  apply(apply(X=res_in_array, MARGIN=1:2, FUN=mean)^2, MARGIN=2, FUN=cumsum)
  }

  GFEVD_results <- array(dim=c(N + 1, d, length(varnames)),
                         dimnames=list(0:N, shocknames, varnames)) # [horizon, shock, variable]
  for(i1 in 1:ncol(GIRF_square_cumsum)) { # Go through the variables and possibly mixing weights
    denominator <- rowSums(GIRF_square_cumsum[, i1, ]) # The denominators for h=0,1,...,N
    for(i2 in 1:d) { # Go through the shocks
      GFEVD_results[ , i2, i1] <- GIRF_square_cumsum[, i1, i2]/denominator
    }
  }

  cat("Finished!\n")
  ret <- structure(list(gfevd_res=GFEVD_results,
                        shock_size=shock_size,
                        N=N,
                        initval_type=initval_type,
                        R1=R1,
                        R2=R2,
                        which_cumulative=which_cumulative,
                        init_regimes=init_regimes,
                        init_values=init_values,
                        include_mixweights=include_mixweights,
                        seeds=seeds,
                        gmvar=gmvar),
                   class="gfevd")
  ret
}
