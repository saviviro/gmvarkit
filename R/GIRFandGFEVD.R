#' @title Estimate generalized impulse response function for a structural GMVAR,
#'   StMVAR, or G-StMVAR model.
#'
#' @description \code{GIRF} estimates generalized impulse response function for
#'   a structural GMVAR, StMVAR, or G-StMVAR model.
#'
#' @inheritParams quantile_residual_tests
#' @inheritParams simulate.gsmvar
#' @inheritParams fitGSMVAR
#' @param which_shocks a numeric vector of length at most \eqn{d}
#'   (\code{=ncol(data)}) and elements in \eqn{1,...,d} specifying the
#'   structural shocks for which the GIRF should be estimated.
#' @param shock_size a non-zero scalar value specifying the common size for all scalar
#'   components of the structural shock. Note that the conditional covariance
#'   matrix of the structural shock is an identity matrix and that the
#'   (generalized) impulse responses may not be symmetric to the sign and size
#'   of the shock.
#' @param N a positive integer specifying the horizon how far ahead should the
#'   generalized impulse responses be calculated.
#' @param R1 the number of repetitions used to estimate GIRF for each initial
#'   value.
#' @param R2 the number of initial values to be drawn from a stationary
#'   distribution of the process or of a specific regime? The confidence bounds
#'   will be sample quantiles of the GIRFs based on different initial values.
#'   Ignored if the argument \code{init_value} is specified.
#' @param which_cumulative a numeric vector with values in \eqn{1,...,d}
#'   (\code{d=ncol(data)}) specifying which the variables for which the impulse
#'   responses should be cumulative. Default is none.
#' @param scale should the GIRFs to some of the shocks be scaled so that they
#'   correspond to a specific magnitude of instantaneous or peak response
#'   of some specific variable (see the argument \code{scale_type})?
#'   Provide a length three vector where the shock of interest
#'   is given in the first element (an integer in \eqn{1,...,d}), the variable of
#'   interest is given in the second element (an integer in \eqn{1,...,d}), and
#'   the magnitude of its instantaneous or peak response in the third element
#'   (a non-zero real number). If the GIRFs of multiple shocks should be scaled, provide
#'   a matrix which has one column for each of the shocks with the columns being
#'   the length three vectors described above.
#' @param scale_type If argument \code{scale} is specified, should the GIRFs be
#'   scaled to match an instantaneous response (\code{"instant"}) or peak response
#'   (\code{"peak"}). If \code{"peak"}, the scale is based on the largest magnitude
#'   of peak response in absolute value. Ignored if \code{scale} is not specified.
#' @param scale_horizon If \code{scale_type == "peak"} what the maximum horizon up
#'   to which peak response is expected? Scaling won't based on values after this.
#' @param ci a numeric vector with elements in \eqn{(0, 1)} specifying the
#'   confidence levels of the confidence intervals.
#' @param include_mixweights should the generalized impulse response be
#'   calculated for the mixing weights as well? \code{TRUE} or \code{FALSE}.
#' @param ncores the number CPU cores to be used in parallel computing. Only
#'   single core computing is supported if an initial value is specified (and
#'   the GIRF won't thus be estimated multiple times).
#' @param plot_res \code{TRUE} if the results should be plotted, \code{FALSE} if
#'   not.
#' @param seeds a length \code{R2} vector containing the random number generator
#'   seed for estimation of each GIRF. A single number of an initial value is
#'   specified. or \code{NULL} for not initializing the seed. Exists for
#'   creating reproducible results.
#' @param ... parameters passed to the plot method \code{plot.girf} that plots
#'   the results.
#' @details The model needs to be structural in order for this function to be
#'   applicable. A structural GMVAR, StMVAR, or G-StMVAR model can be estimated
#'   by specifying the argument \code{structural_pars} in the function \code{fitGSMVAR}.
#'
#'   The confidence bounds reflect uncertainty about the initial state (but
#'   currently not about the parameter estimates) if initial values are not
#'   specified. If initial values are specified, there won't currently be
#'   confidence intervals. See the cited paper by Virolainen (2021) for details
#'   about the algorithm.
#'
#'   Note that if the argument \code{scale} is used, the scaled responses of
#'   the mixing weights might be more than one in absolute value.
#' @return Returns a class \code{'girf'} list with the GIRFs in the first
#'   element (\code{$girf_res}) and the used arguments the rest. The first
#'   element containing the GIRFs is a list with the \eqn{m}th element
#'   containing the point estimates for the GIRF in \code{$point_est} (the first
#'   element) and confidence intervals in \code{$conf_ints} (the second
#'   element). The first row is for the GIRF at impact \eqn{(n=0)}, the second
#'   for \eqn{n=1}, the third for \eqn{n=2}, and so on.
#'
#'   The element \code{$all_girfs} is a list containing results from all the individual GIRFs
#'   obtained from the MC repetitions. Each element is for one shock and results are in
#'   array of the form \code{[horizon, variables, MC-repetitions]}.
#' @seealso \code{\link{GFEVD}}, \code{\link{fitGSMVAR}}, \code{\link{GSMVAR}},
#'   \code{\link{gsmvar_to_sgsmvar}}, \code{\link{reorder_W_columns}},
#'   \code{\link{swap_W_signs}}, \code{\link{simulate.gsmvar}},
#'   \code{\link{predict.gsmvar}}, \code{\link{profile_logliks}},
#'   \code{\link{quantile_residual_tests}}, \code{\link{LR_test}},
#'   \code{\link{Wald_test}}
#' @inherit in_paramspace_int references
#' @examples
#'  \donttest{
#'  # These are long-running examples that use parallel computing.
#'  # It takes approximately 30 seconds to run all the below examples.
#'
#'  # Structural GMVAR(2, 2), d=2 model identified with sign-constraints:
#'  params22s <- c(0.36, 0.121, 0.484, 0.072, 0.223, 0.059, -0.151, 0.395,
#'   0.406, -0.005, 0.083, 0.299, 0.218, 0.02, -0.119, 0.722, 0.093, 0.032,
#'    0.044, 0.191, 0.057, 0.172, -0.46, 0.016, 3.518, 5.154, 0.58)
#'  W_22 <- matrix(c(1, 1, -1, 1), nrow=2, byrow=FALSE)
#'  mod22s <- GSMVAR(gdpdef, p=2, M=2, params=params22s,
#'   structural_pars=list(W=W_22))
#'  mod22s
#'  # Alternatively, use:
#'  #fit22s <- fitGSMVAR(gdpdef, p=2, M=2, structural_pars=list(W=W_22),
#'  #                   ncalls=20, seeds=1:20)
#'  # To obtain an estimated version of the same model.
#'
#'  # Estimating the GIRFs of both structural shocks with initial values
#'  # drawn from the stationary distribution of the process,
#'  # 12 periods ahead, confidence levels 0.95 and 0.8:
#'  girf1 <- GIRF(mod22s, N=12, R1=100, R2=100)
#'  girf1
#'  plot(girf1)
#'
#'  # Estimating the GIRF of the second shock only, 12 periods ahead
#'  # and shock size 1, initial values drawn from the stationary distribution
#'  # of the first regime, confidence level 0.9:
#'  girf2 <- GIRF(mod22s, which_shocks=2, shock_size=1, N=12, init_regimes=1,
#'                ci=0.9, R1=100, R2=100)
#'
#'  # Estimating the GIRFs of both structural shocks, negative one standard
#'  # error shock, N=20 periods ahead, estimation based on 200 Monte Carlo
#'  # simulations, and fixed initial values given by the last p observations
#'  # of the data:
#'  girf3 <- GIRF(mod22s, shock_size=-1, N=20, R1=200,
#'                init_values=mod22s$data)
#'  }
#' @export

GIRF <- function(gsmvar, which_shocks, shock_size=1, N=30, R1=250, R2=250, init_regimes=1:sum(gsmvar$model$M), init_values=NULL,
                 which_cumulative=numeric(0), scale=NULL, scale_type=c("instant", "peak"), scale_horizon=N,
                 ci=c(0.95, 0.80), include_mixweights=TRUE, ncores=2, plot_res=TRUE, seeds=NULL, ...) {
  scale_type <- match.arg(scale_type)
  gsmvar <- gmvar_to_gsmvar(gsmvar) # Backward compatibility
  p <- gsmvar$model$p
  M <- sum(gsmvar$model$M)
  d <- gsmvar$model$d

  stopifnot(N %% 1 == 0 && N > 0)
  stopifnot(scale_horizon %in% 0:N)
  if(is.null(gsmvar$model$structural_pars)) stop("Only structural models are supported")
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
  if(!is.null(scale)) {
    scale <- as.matrix(scale)
    stopifnot(all(scale[1,] %in% 1:d)) # All shocks in 1,...,d
    stopifnot(length(unique(scale[1,])) == length(scale[1,])) # No duplicate scales for the same shock
    stopifnot(all(scale[2,] %in% 1:d)) # All variables in 1,...,d
    stopifnot(all(scale[3,] != 0)) # No zero initial magnitudes

    # For the considered shocks, check that there are not zero-constraints for the variable whose initial response is scaled.
    for(i1 in 1:ncol(scale)) {
      if(!is.na(gsmvar$model$structural_pars$W[scale[2, i1], scale[1, i1]]) && gsmvar$model$structural_pars$W[scale[2, i1], scale[1, i1]] == 0) {
        stop(paste("Instantaneous response of the variable that has a zero constraint for the considered shock cannot be scaled"))
      }
    }
  }
  stopifnot(shock_size != 0)
  if(length(which_cumulative) > 0) {
    which_cumulative <- unique(which_cumulative)
    stopifnot(all(which_cumulative %in% 1:d))
  }

  # Function that estimates GIRF
  get_one_girf <- function(shock_numb, shock_size, seed) {
    simulate.gsmvar(gsmvar, nsim=N + 1, seed=seed, init_values=init_values, init_regimes=init_regimes, ntimes=R1,
                    girf_pars=list(shock_numb=shock_numb,
                                   shock_size=shock_size,
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
  all_GIRFS <- vector("list", length=length(which_shocks))
  names(all_GIRFS) <- paste0("shock", which_shocks)
  names(GIRF_results) <- paste0("shock", which_shocks)

  for(i1 in 1:length(which_shocks)) { # Go through shocks
    res_in_array <- array(unlist(GIRF_shocks[[i1]]), dim=c(N + 1, d + ifelse(include_mixweights, M, 0), R2)) # [horizon, variables, MC-repetitions]
    if(length(which_cumulative) > 0) {
      for(i2 in which_cumulative) {
        res_in_array[, i2, ] <- apply(res_in_array[, i2, , drop=FALSE], MARGIN=3, FUN=cumsum) # Replace GIRF with cumulative GIRF
      }
    }

    # Scale the GIRFs if specified
    if(which_shocks[i1] %in% scale[1,]) { # GIRF of this shock should be scaled
      which_col <- which(which_shocks[i1] == scale[1,]) # which column of the scale-matrix contains the argument for this specific shock
      which_var <- scale[2, which_col] # According to initial/peak response of which variable the GIRFs should be scaled
      magnitude <- scale[3, which_col] # What should be the magnitude of the initial/peak response of this variable

      my_comparison_fun <- function(vec1, scalar1) which(abs(vec1 - scalar1) < .Machine$double.eps)[1] # To avoid potential problems with using == to compare numerical values
      for(i2 in 1:R2) { # Go through the MC repetitions
        # The scaling scalar is different for each MC repetition, because the instantaneous/peak movement is generally
        # different with different starting values.
        if(scale_type == "instant") {  # Scale by initial response
          one_scale <- magnitude/res_in_array[1, which_var, i2]
        } else {  # scale_type == "peak", "peak_max" or "peak_min", scale by peak response
          inds <- 1:(scale_horizon + 1) # +1 for period 0
          one_scale <- magnitude/res_in_array[my_comparison_fun(vec1=abs(res_in_array[inds, which_var, i2]),
                                                                scalar1=max(abs(res_in_array[inds, which_var, i2]))), which_var, i2]
          #one_scale <- magnitude/res_in_array[which(abs(res_in_array[, which_var, i2]) == max(abs(res_in_array[, which_var, i2]))), which_var, i2]
        }
        res_in_array[, , i2] <- one_scale*res_in_array[, , i2]
      }
    }
    all_GIRFS[[i1]] <- res_in_array

    # Point estimates, confidence intervals
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
                        all_girfs=all_GIRFS,
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
                        gsmvar=gsmvar),
                   class="girf")
  if(plot_res) plot(ret, ...)
  ret
}



#' @title Estimate generalized forecast error variance decomposition for a structural
#'   GMVAR, StMVAR, or G-StMVAR model.
#'
#' @description \code{GFEVD} estimates generalized generalized forecast error variance decomposition for
#'  a structural GMVAR, StMVAR, or G-StMVAR model.
#'
#' @inheritParams GIRF
#' @param shock_size What shocks size should be used for all shocks? By the definition of the SGMVAR,
#'   SStMVAR, and SG-StMVAR model, the conditional covariance matrix of the structural shock is identity matrix.
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
#'   GMVAR, StMVAR, or G-StMVAR model can be estimated by specifying the argument \code{structural_pars} in the function \code{fitGSMVAR}.
#'
#'   The GFEVD is a forecast error variance decomposition calculated with the generalized impulse response function (GIRF).
#'   Lanne and Nyberg (2016) for details. Note, however, that the related GIRFs are calculated using the algorithm given in
#'   Virolainen (2021).
#' @return Returns and object of class 'gfevd' containing the GFEVD for all the variables and if
#'   \code{include_mixweights=TRUE} also to the mixing weights. Note that the decomposition does not
#'   exist at horizon zero for mixing weights because the related GIRFs are always zero at impact.
#' @seealso \code{\link{GIRF}}, \code{\link{fitGSMVAR}}, \code{\link{GSMVAR}}, \code{\link{gsmvar_to_sgsmvar}},
#'  \code{\link{reorder_W_columns}}, \code{\link{swap_W_signs}}, \code{\link{simulate.gsmvar}}
#' @references
#' @references
#'  \itemize{
#'    \item Lanne M. and Nyberg H. 2016. Generalized Forecast Error Variance Decomposition for Linear
#'      and Nonlineae Multivariate Models. \emph{Oxford Bulletin of Economics and Statistics}, \strong{78}, 4, 595-603.
#'    \item Kalliovirta L., Meitz M. and Saikkonen P. 2016. Gaussian mixture vector autoregression.
#'          \emph{Journal of Econometrics}, \strong{192}, 485-498.
#'    \item Virolainen S. 2021. Structural Gaussian mixture vector autoregressive model. Unpublished working
#'      paper, available as arXiv:2007.04713.
#'    \item Virolainen S. 2021. Gaussian and Student's t mixture vector autoregressive model. Unpublished working
#'      paper, available as arXiv:2109.13648.
#'  }
#' @examples
#'  \donttest{
#'  # These are long-running examples that use parallel computing.
#'  # It takes approximately 30 seconds to run all the below examples.
#'
#'  # Structural GMVAR(2, 2), d=2 model identified with sign-constraints:
#'  params22s <- c(0.36, 0.121, 0.484, 0.072, 0.223, 0.059, -0.151, 0.395,
#'   0.406, -0.005, 0.083, 0.299, 0.218, 0.02, -0.119, 0.722, 0.093, 0.032,
#'   0.044, 0.191, 0.057, 0.172, -0.46, 0.016, 3.518, 5.154, 0.58)
#'  W_22 <- matrix(c(1, 1, -1, 1), nrow=2, byrow=FALSE)
#'  mod22s <- GSMVAR(gdpdef, p=2, M=2, params=params22s,
#'   structural_pars=list(W=W_22))
#'  mod22s
#'  # Alternatively, use:
#'  #fit22s <- fitGSMVAR(gdpdef, p=2, M=2, structural_pars=list(W=W_22),
#'  #                   ncalls=20, seeds=1:20)
#'  # To obtain an estimated version of the same model.
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

GFEVD <- function(gsmvar, shock_size=1, N=30, initval_type=c("data", "random", "fixed"), R1=250, R2=250,
                  init_regimes=NULL, init_values=NULL, which_cumulative=numeric(0), include_mixweights=FALSE,
                  ncores=2, seeds=NULL) {
  gsmvar <- gmvar_to_gsmvar(gsmvar) # Backward compatibility
  initval_type <- match.arg(initval_type)
  p <- gsmvar$model$p
  M <- sum(gsmvar$model$M)
  d <- gsmvar$model$d
  if(is.null(gsmvar$model$structural_pars)) stop("Only structural models are supported")
  if(M == 1) include_mixweights <- FALSE
  if(initval_type == "data") {
    if(is.null(gsmvar$data)) stop("The model does not contain data! Add data with the function 'add_data' or select another 'initval_type'.")
    stopifnot(nrow(gsmvar$data) >= p)
    R2 <- nrow(gsmvar$data) - p + 1 # The number of length p histories in the data
    all_initvals <- array(vapply(1:R2, function(i1) gsmvar$data[i1:(i1 + p - 1),], numeric(p*d)), dim=c(p, d, R2)) # [, , i1] for i1:th initval
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
    simulate.gsmvar(gsmvar, nsim=N + 1, init_values=init_values_for_1girf, ntimes=R1, seed=seed, init_regimes=init_regimes,
                    girf_pars=list(shock_numb=shock_numb,
                                   shock_size=shock_size,
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

  if(is.null(colnames(gsmvar$data))) {
    varnames <- paste0("Variable", 1:d)
  } else {
    varnames <- colnames(gsmvar$data)
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
                        gsmvar=gsmvar),
                   class="gfevd")
  ret
}
