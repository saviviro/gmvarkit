
#' @title Estimate generalized impulse response function for a structural GMVAR model.
#'
#' @description \code{GIRF} estimate generalized impulse response function for
#'  a structural GMVAR model.
#'
#' @inheritParams simulateGMVAR
#' @param which_shocks a numeric vector of length at most \eqn{d} (\code{=ncol(data)})
#'   and elements in \eqn{1,...,d} specifying the structural shocks for which the GIRF
#'   should be estimated.
#' @param shock_size a vector with the same length as \code{which_shocks} specifying
#'   the size of each structural shock. Alternatively, is a scalar value that specifies a
#'   common shock size for all structural shocks. By default, the shock size is
#'   one, which is then amplified by the B-matrix according to the conditional standard deviation
#'   of the model.
#' @param N a positive integer specifying the horizon how far ahead should the generalized
#'   impulse responses be calculated?
#' @param R1 the number of repetitions used to estimate GIRF for each initial value?
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
#' @seealso \code{\link{fitGMVAR}}, \code{\link{GMVAR}}, \code{\link{gmvar_to_sgmvar}}, \code{\link{reorder_W_columns}},
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

GIRF <- function(gmvar, which_shocks, shock_size, N=30, R1=250, R2=250, init_regimes=1:gmvar$model$M, init_values=NULL, which_cumulative,
                 ci=c(0.95, 0.80), include_mixweights=TRUE, ncores=min(2, parallel::detectCores()), plot=TRUE, seeds=NULL) {
  on.exit(closeAllConnections())

  p <- gmvar$model$p
  M <- gmvar$model$M
  d <- gmvar$model$d
  if(M == 1) include_mixweights <- FALSE
  if(missing(which_shocks)) {
    which_shocks <- 1:d
  } else {
    stopifnot(length(which_shocks) <= d && all(which_shocks %in% 1:d) && length(unique(which_shocks)) == length(which_shocks))
  }
  if(!is.null(init_values)) R2 <- 1
  if(!is.null(seeds) && length(seeds) != R2) stop("The argument 'seeds' needs be NULL or a vector of length 'R2'")
  stopifnot(length(init_regimes) <= M && all(init_regimes %in% 1:M) && length(unique(init_regimes)) == length(init_regimes))
  if(is.null(gmvar$model$structural_pars)) stop("Only structural models are supported")
  stopifnot(length(ci) > 0 && all(ci > 0 & ci < 1))

  if(missing(shock_size)) {
    shock_size <- rep(1, times=length(which_shocks))
  } else {
    stopifnot(length(shock_size) == length(which_shocks) | length(shock_size) == 1)
    if(length(shock_size) == 1) shock_size <- rep(shock_size, times=length(which_shocks))
  }
  if(missing(which_cumulative)) {
    which_cumulative <- numeric(0)
  } else {
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
  parallel::clusterExport(cl, ls(environment(GIRF)), envir = environment(GIRF)) # assign all variables from package:gmvarkit
  parallel::clusterEvalQ(cl, c(library(mvnfast), library(pbapply)))

  GIRF_shocks <- vector("list", length=length(which_shocks))

  for(i1 in 1:length(which_shocks)) {
    cat(paste0("Estimating GIRFs for structural shock ", which_shocks[i1], "..."), "\n")
    GIRF_shocks[[i1]] <- pbapply::pblapply(1:R2, function(i2) get_one_girf(shock_numb=which_shocks[i1], shock_size=shock_size[i1], seed=seeds[i2]), cl=cl)
  }
  parallel::stopCluster(cl=cl)

  GIRF_results <- vector("list", length=length(which_shocks))
  if(!is.null(gmvar$data) && !is.null(colnames(gmvar$data))) {
    names(GIRF_results) <- colnames(gmvar$data)[which_shocks]
  } else {
    names(GIRF_results) <- paste("shock", which_shocks)
  }

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
                        include_mixweights=include_mixweights),
                   class="girf")
  if(plot) plot(ret)
  ret
}
