
#' @title Estimate generalized impulse response function for a structural GMVAR model.
#'
#' @description \code{GIRF} estimate generalized impulse response function for
#'  a structural GMVAR model.
#'
#' @inheritParams simulateGMVAR
#' @param variables a numeric vector of length at most \eqn{d} (\code{=ncol(data)})
#'   and elements in \eqn{1,...,d} specifying the variables for which the GIRF
#'   should be estimated.
#' @param shock_size a vector with the same length as \code{variables} specifying
#'   the size of the structural shock for each variable. By default, the size of one standard
#'   deviation is used, calculated as the weighted average of the component process error term
#'   standard deviations with weights given by the mixing weight parameters.
#' @param N a positive integer specifying the horizon how far ahead should the generalized
#'   impulse responses be calculated?
#' @param R1 the number of repetitions used to estimate GIRF for each initial value?
#' @param R2 the number of initial values to be drawn from a stationary distribution of the process
#'   or of a specific regime? The confidence bounds will be sample quantiles of the GIRFs based on
#'   different initial values. Ignored if the argument \code{init_value} is specified.
#' @param init_regimes a numeric vector of length at most \eqn{M} and elements in \eqn{1,...,M}
#'   specifying the regimes from which the initial values should be generated from. The initial
#'   values will be generated from a mixture distribution with the mixture components being the
#'   stationary distirbutions of the specific regimes and the (proportiional) mixing weights given
#'   by the mixing weight parameters of those regimes. Note that if \code{init_regimes=1:M}, the
#'   initial values are generated from the stationary distribution of the process and if
#'   \code{init_regimes=m}, the initial value are generated from the stationary distribution
#'   of the \eqn{m}th regime. Ignored if \code{init_value} is specified.
#' @param init_values a matrix or a multivariate class \code{'ts'} object with \eqn{d} columns
#'   and at least \eqn{p} rows specifying an initial value for the GIRF. The last \eqn{p} rows
#'   are taken to be the initial value assuming that the \strong{last} row is the most recent observation.
#' @param ci a numeric vector with elements in \eqn{(0, 1)} specifying the confidence levels of the
#'   confidence intervals.
#' @param include_mixweights should the generalized impulse response be calculated for the mixing weights
#'   as well? \code{TRUE} or \code{FALSE}.
#' @param ncores the number CPU cores to be used in parallel computing. Only single core computing is
#'   supported if an initial value is specified (and the GIRF won't thus be estimated multiple times).
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
#'  data(eurusd, package="gmvarkit")
#'  data <- cbind(10*eurusd[,1], 100*eurusd[,2])
#'  colnames(data) <- colnames(eurusd)
#'
#'  # Structural GMVAR(2, 2), d=2 model identified with sign-constraints:
#'  params222s <- c(1.428, -0.808, 1.029, 5.84, 1.314, 0.145, 0.094, 1.292,
#'    -0.389, -0.07, -0.109, -0.281, 1.248, 0.077, -0.04, 1.266, -0.272,
#'    -0.074, 0.034, -0.313, 0.903, 0.718, -0.324, 2.079, 7, 1.44, 0.742)
#'  W_222 <- matrix(c(1, NA, -1, 1), nrow=2, byrow=FALSE)
#'  mod222s <- GMVAR(data, p=2, M=2, params=params222s, structural_pars=list(W=W_222))
#'  mod222s
#'  # Alternatively, use:
#'  # fit222s <- fitGMVAR(data, p=2, M=2, structural_pars=list(W=W_222),
#'  #                     ncalls=20, seeds=1:20)
#'  # To obtain an estimated version of the same model.
#'
#'  # Estimating the GIRFs of both variables with default arguments
#'  # (initial values are drawn from the stationary distribution of the process,
#'  # 30 periods ahead, confidence levels 0.95 and 0.8):
#'  girf1 <- GIRF(mod222s)
#'  plot(girf1)
#'  girf1
#'
#'  # Estimating the GIRF of the second variable only, 36 periods ahead
#'  # and shock size 1, initial values drawn from the stationary distribution
#'  # of the first regime, confidence level 0.9:
#'  girf2 <- GIRF(mod222s, variables=2, shock_size=1, N=36, init_regimes=1, ci=0.9)
#'  plot(girf2)
#'
#'  # Estimating the GIRFs of both variables, shock sizes 1 and 3, N=50 periods ahead,
#'  # estimation based on 1000 Monte Carlo simulations, and fixed initial values given
#'  # by the last p observations of the data:
#'  girf3 <- GIRF(mod222s, shock_size=c(1, 3), N=50, R1=1000, init_values=mod222s$data)
#'  }
#' @export

GIRF <- function(gmvar, variables, shock_size, N=30, R1=250, R2=250, init_regimes=1:M, init_values=NULL,
                 ci=c(0.95, 0.80), include_mixweights=TRUE, ncores=min(2, parallel::detectCores()), seeds=NULL) {
  on.exit(closeAllConnections())

  p <- gmvar$model$p
  M <- gmvar$model$M
  d <- gmvar$model$d
  if(M == 1) include_mixweights <- FALSE
  if(missing(variables)) {
    variables <- 1:d
  } else {
    stopifnot(length(variables) <= d && all(variables %in% 1:d) && length(unique(variables)) == length(variables))
  }
  if(!is.null(init_values)) R2 <- 1
  if(!is.null(seeds) && length(seeds) != R2) stop("The argument 'seeds' needs be NULL or a vector of length 'R2'")
  stopifnot(length(init_regimes) <= M && all(init_regimes %in% 1:M) && length(unique(init_regimes)) == length(init_regimes))
  stopifnot(!is.null(gmvar$model$structural_pars))
  stopifnot(length(ci) > 0 && all(ci > 0 & ci < 1))

  # Function that estimates GIRF
  get_one_girf <- function(variable, shock_size, seed) {
    simulateGMVAR(gmvar, nsimu=N + 1, init_values=init_values, ntimes=R1, seed=seed, girf_pars=list(variable=variable,
                                                                                                    shock_size=shock_size,
                                                                                                    init_regimes=init_regimes,
                                                                                                    include_mixweights=include_mixweights))
  }

  # Calculate shock sizes if not specified
  if(missing(shock_size)) {
    params <- reform_constrained_pars(p=p, M=M, d=d, params=gmvar$params, constraints=gmvar$model$constraints,
                                      structural_pars=gmvar$model$structural_pars)
    structural_pars <- get_unconstrained_structural_pars(gmvar$model$structural_pars)
    all_Omega <- pick_Omegas(p=p, M=M, d=d, params=params, structural_pars=structural_pars)
    alphas <- pick_alphas(p=p, M=M, d=d, params=params)
    shock_size_tmp <- matrix(nrow=M, ncol=d)
    for(m in 1:M) {
      shock_size_tmp[m,] <- diag(all_Omega[, , m])*alphas[m]
    }
    shock_size <- colMeans(shock_size_tmp)
  } else {
    stopifnot(length(shock_size) == length(variables))
  }


  if(ncores > parallel::detectCores()) {
    ncores <- parallel::detectCores()
    message("ncores was set to be larger than the number of cores detected")
  }
  if(is.null(init_values)) {
    cat(paste("Using", ncores, "cores to estimate", R2,"GIRFs for", length(variables), "variables,", "each based on", R1, "Monte Carlo repetitions."), "\n")
  } else {
    cat(paste("Using", ncores, "cores to estimate one GIRF for", length(variables), "variables, each based on", R1, "Monte Carlo repetitions."), "\n")
  }

  ### Calculate the GIRFs ###
  cl <- parallel::makeCluster(ncores)
  parallel::clusterExport(cl, ls(environment(GIRF)), envir = environment(GIRF)) # assign all variables from package:gmvarkit
  parallel::clusterEvalQ(cl, c(library(Brobdingnag), library(mvnfast), library(pbapply)))

  GIRF_variables <- vector("list", length=length(variables))

  for(i1 in 1:length(variables)) {
    cat(paste0("Estimating GIRFs for variable ", variables[i1], "..."), "\n")
    GIRF_variables[[i1]] <- pbapply::pblapply(1:R2, function(i2) get_one_girf(variable=variables[i1], shock_size=shock_size[i1], seed=seeds[i2]), cl=cl)
  }
  parallel::stopCluster(cl=cl)

  GIRF_results <- vector("list", length=length(variables))
  if(!is.null(gmvar$data) && !is.null(colnames(gmvar$data))) {
    names(GIRF_results) <- colnames(gmvar$data)[variables]
  } else {
    names(GIRF_results) <- paste("variable", variables)
  }

  for(i1 in 1:length(variables)) {
    res_in_array <- array(unlist(GIRF_variables[[i1]]), dim=c(N + 1, d + ifelse(include_mixweights, M, 0), R2))
    colnames(res_in_array) <- colnames(GIRF_variables[[1]][[1]])
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
  structure(list(girf_res=GIRF_results,
                 variables=variables,
                 shock_size=shock_size,
                 N=N,
                 R1=R1,
                 R2=R2,
                 ci=ci,
                 init_regimes=init_regimes,
                 init_values=init_values,
                 include_mixweights=include_mixweights),
            class="girf")
}
