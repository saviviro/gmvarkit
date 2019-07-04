#' @title Two-phase maximum likelihood estimation of GMVAR model
#'
#' @description \code{fitGMVAR} estimates GMVAR model in two phases: in the first phase it uses genetic algorithm
#'   to find starting values for gradient based variable metric algorithm, which it then uses to finalize the estimation in the second
#'   phase. Parallel computing is used to perform multiple rounds of estimations in parallel.
#'
#' @inheritParams GAfit
#' @param ncalls number of estimation rounds that should be performed.
#' @param ncores number cores to be used in parallel computing.
#' @param maxit maximum number of iterations in the variable metric algorithm.
#' @param print_res should summaries of estimation results be printed?
#' @param ... additional settings passed to the function \code{GAfit} employing the genetic algorithm.
#' @details
#'    Because of complexity and multimodality of the log-likelihood function, it's \strong{not certain} that the estimation
#'    algorithms will end up in the global maximum point. It's expected that most of the estimation rounds will end up in some local maximum
#'    point instead. Therefore a number of estimation rounds is required for reliable results. Because of the nature of the model,
#'    the estimation may fail especially in the cases where the number of mixture components is chosen too large.
#'
#'    Overall the estimation process is computationally heavy and it might take considerably long time for large models with
#'    large number of observations. If the iteration limit \code{maxit} in the variable metric algorithm is reached, one can continue
#'    the estimation by iterating more with the function \code{iterate_more}.
#'
#'    The genetic algorithm is mostly based on the description by \emph{Dorsey and Mayer (1995)}, but it includes some extra
#'    functionality designed for this particular estimation problem. The genetic algorithm uses (slightly modified) individually
#'    adaptive crossover and mutation rates described by \emph{Patnaik and Srinivas (1994)} and employs (50\%) fitness
#'    inheritance discussed by \emph{Smith, Dike and Stegmann (1995)}.
#'
#'    The gradient based variable metric algorithm used in the second phase is implemented with function \code{optim}
#'    from the package \code{stats}.
#' @return Returns an object of class \code{'gmvar'} defining the estimated GMVAR model. Multivariate quantile residuals
#'   (Kalliovirta and Saikkonen 2010) are also computed and included in the returned object. In addition, the returned
#'   object contains the estimates and log-likelihood values from all the estimation rounds performed.
#'   The estimated parameter vector can be obtained at \code{gmvar$params} (and corresponding approximate standard errors
#'   at \code{gmvar$std_errors}) and it is...
#'   \describe{
#'     \item{\strong{Regular models:}}{
#'       a size \eqn{((M(pd^2+d+d(d+1)/2+1)-1)x1)} vector that has form
#'       \strong{\eqn{\theta}}\eqn{ = }(\strong{\eqn{\upsilon}}\eqn{_{1}},
#'       ...,\strong{\eqn{\upsilon}}\eqn{_{M}}, \eqn{\alpha_{1},...,\alpha_{M-1}}), where:
#'       \itemize{
#'         \item \strong{\eqn{\upsilon}}\eqn{_{m}} \eqn{ = (\phi_{m,0},}\strong{\eqn{\phi}}\eqn{_{m}}\eqn{,\sigma_{m})}
#'         \item \strong{\eqn{\phi}}\eqn{_{m}}\eqn{ = (vec(A_{m,1}),...,vec(A_{m,p})}
#'         \item and \eqn{\sigma_{m} = vech(\Omega_{m})}, m=1,...,M.
#'       }
#'     }
#'     \item{\strong{Constrained models:}}{
#'       a size \eqn{((M(d+d(d+1)/2+1)+q-1)x1)} vector that has form
#'       \strong{\eqn{\theta}}\eqn{ = (\phi_{1,0},...,\phi_{M,0},}\strong{\eqn{\psi}}
#'       \eqn{,\sigma_{1},...,\sigma_{M},\alpha_{1},...,\alpha_{M-1})}, where:
#'       \itemize{
#'         \item \strong{\eqn{\psi}} \eqn{(qx1)} satisfies (\strong{\eqn{\phi}}\eqn{_{1}}\eqn{,...,}
#'         \strong{\eqn{\phi}}\eqn{_{M}) =} \strong{\eqn{C \psi}}. Here \strong{\eqn{C}} is \eqn{(Mpd^2xq)}
#'         constraint matrix.
#'       }
#'     }
#'   }
#'   Above \eqn{\phi_{m,0}} is the intercept parameter, \eqn{A_{m,i}} denotes the \eqn{i}:th coefficient matrix of the \eqn{m}:th
#'   mixture component, \eqn{\Omega_{m}} denotes the error term covariance matrix of the \eqn{m}:th mixture component and
#'   \eqn{\alpha_{m}} is the mixing weight parameter.
#'   If \code{parametrization=="mean"}, just replace each \eqn{\phi_{m,0}} with regimewise mean \eqn{\mu_{m}}.
#'   \eqn{vec()} is vectorization operator that stacks columns of a given matrix into a vector. \eqn{vech()} stacks columns
#'   of a given matrix from the principal diagonal downwards (including elements on the diagonal) into a vector.
#'   The notations are in line with the cited article by \emph{Kalliovirta, Meitz and Saikkonen (2016)}.
#'
#'   Remark that the first autocovariance/correlation matrix in \code{$uncond_moments} is for the lag zero,
#'   the second one for the lag one, etc.
#' @section S3 methods:
#'   The following S3 methods are supported for class \code{'gmvar'}: \code{logLik}, \code{residuals}, \code{print}, \code{summary},
#'    \code{predict} and \code{plot}.
#' @seealso \code{\link{GMVAR}}, \code{\link{iterate_more}}, \code{\link{predict.gmvar}},
#'   \code{\link{simulateGMVAR}}, \code{\link{quantile_residual_tests}}, \code{\link{print_std_errors}},
#'   \code{\link{swap_parametrization}}, \code{\link{get_gradient}}
#' @references
#'  \itemize{
#'    \item Dorsey R. E. and Mayer W. J. 1995. Genetic algorithms for estimation problems with multiple optima,
#'          nondifferentiability, and other irregular features. \emph{Journal of Business & Economic Statistics},
#'         \strong{13}, 53-66.
#'    \item Kalliovirta L., Meitz M. and Saikkonen P. 2016. Gaussian mixture vector autoregression.
#'          \emph{Journal of Econometrics}, \strong{192}, 485-498.
#'    \item Kalliovirta L. and Saikkonen P. 2010. Reliable Residuals for Multivariate Nonlinear
#'          Time Series Models. \emph{Unpublished Revision of HECER Discussion Paper No. 247}.
#'    \item Patnaik L.M. and Srinivas M. 1994. Adaptive Probabilities of Crossover and Mutation in Genetic Algorithms.
#'          \emph{Transactions on Systems, Man and Cybernetics} \strong{24}, 656-667.
#'    \item Smith R.E., Dike B.A., Stegmann S.A. 1995. Fitness inheritance in genetic algorithms.
#'          \emph{Proceedings of the 1995 ACM Symposium on Applied Computing}, 345-350.
#'  }
#' @examples
#' \donttest{
#' ## These are long running examples that use parallel computing!
#'
#' # These examples use the data 'eurusd' which comes with the
#' # package, but in a scaled form (similar to Kalliovirta et al. 2016).
#' data(eurusd, package="gmvarkit")
#' data <- cbind(10*eurusd[,1], 100*eurusd[,2])
#' colnames(data) <- colnames(eurusd)
#'
#' # GMVAR(1,2) model with default settings
#' fit12 <- fitGMVAR(data, p=1, M=2)
#' fit12
#' plot(fit12)
#' summary(fit12)
#'
#' # GMVAR(2,2) model with mean parametrization
#' fit22 <- fitGMVAR(data, p=2, M=2, parametrization="mean")
#' fit22
#'
#' # GMVAR(2,2) model with autoregressive parameters restricted
#' # to be the same for both regimes
#' C_mat <- rbind(diag(2*2^2), diag(2*2^2))
#' fit22c <- fitGMVAR(data, p=2, M=2, constraints=C_mat)
#' fit22c
#'
#' # GMVAR(2,2) model with autoregressive parameters restricted
#' # to be the same for both regimes and non-diagonl elements
#' # the coefficient matrices constrained to zero. Estimation
#' # with only 10 estimation rounds.
#' tmp <- matrix(c(1, rep(0, 10), 1, rep(0, 8), 1, rep(0, 10), 1),
#'  nrow=2*2^2, byrow=FALSE)
#' C_mat2 <- rbind(tmp, tmp)
#' fit22c2 <- fitGMVAR(data, p=2, M=2, constraints=C_mat2, ncalls=10)
#' fit22c2
#' }
#' @export


fitGMVAR <- function(data, p, M, conditional=TRUE, parametrization=c("intercept", "mean"), constraints=NULL, ncalls=round(10 + 9*log(M)),
                     ncores=min(ncalls, parallel::detectCores()), maxit=300, print_res=TRUE, ...) {

  on.exit(closeAllConnections())
  parametrization <- match.arg(parametrization)
  if(!all_pos_ints(c(p, M, ncalls, ncores, maxit))) stop("Arguments p, M, ncalls, ncores and maxit must be positive integers")
  data <- check_data(data=data, p=p)
  d <- ncol(data)
  n_obs <- nrow(data)
  npars <- n_params(p=p, M=M, d=d, constraints=constraints)
  if(npars >= d*nrow(data)) stop("There are at least as many parameters in the model as there are observations in the data")
  dot_params <- list(...)
  minval <- ifelse(is.null(dot_params$minval), -(10^(ceiling(log10(n_obs)) + d) - 1), dot_params$minval)
  red_criteria <- ifelse(rep(is.null(dot_params$red_criteria), 2), c(0.05, 0.01), dot_params$red_criteria)

  if(ncores > parallel::detectCores()) {
    ncores <- parallel::detectCores()
    message("ncores was set to be larger than the number of cores detected")
  }
  if(ncores > ncalls) {
    ncores <- ncalls
    message("ncores was set to be larger than the number of estimation rounds")
  }
  cat(paste("Using", ncores, "cores for", ncalls, "estimations rounds..."), "\n")

  ### Genetic algorithm optimization ###
  cl <- parallel::makeCluster(ncores)
  parallel::clusterExport(cl, ls(environment(fitGMVAR)), envir = environment(fitGMVAR)) # assign all variables from package:gmvarkit
  parallel::clusterEvalQ(cl, c(library(Brobdingnag), library(mvnfast), library(pbapply)))

  cat("Optimizing with genetic algorithm...", "\n")
  GAresults <- pbapply::pblapply(1:ncalls, function(x) GAfit(data=data, p=p, M=M, conditional=conditional, parametrization=parametrization,
                                                              constraints=constraints, ...), cl=cl)
  parallel::stopCluster(cl=cl)

  loks <- vapply(1:ncalls, function(i1) loglikelihood_int(data, p, M, params=GAresults[[i1]], conditional=conditional,
                                                          parametrization=parametrization, constraints=constraints,
                                                          check_params=TRUE, to_return="loglik", minval=minval), numeric(1))

  if(print_res == TRUE) {
    cat("Results from genetic algorithm:", "\n")
    cat(paste("lowest value: ", round(min(loks), 3)), "\n")
    cat(paste("mean value:   ", round(mean(loks), 3)), "\n")
    cat(paste("largest value:", round(max(loks), 3)), "\n")
  }


  ### Variable metric algorithm optimization ###
  loglik_fn <- function(params) {
    tryCatch(loglikelihood_int(data, p, M, params=params, conditional=conditional, parametrization=parametrization,
                               constraints=constraints, check_params=TRUE, to_return="loglik", minval=minval), error=function(e) minval)
  }

  h <- 6e-6
  I <- diag(rep(1, npars))
  loglik_grad <- function(params) {
    vapply(1:npars, function(i1) (loglik_fn(params + I[i1,]*h) - loglik_fn(params - I[i1,]*h))/(2*h), numeric(1))
  }

  cl <- parallel::makeCluster(ncores)
  parallel::clusterExport(cl, ls(environment(fitGMVAR)), envir = environment(fitGMVAR)) # assign all variables from package:gmvarkit
  parallel::clusterEvalQ(cl, c(library(Brobdingnag), library(mvnfast), library(pbapply)))

  cat("Optimizing with variable metric algorithm...\n")
  NEWTONresults <- pbapply::pblapply(1:ncalls, function(i1) optim(par=GAresults[[i1]], fn=loglik_fn, gr=loglik_grad, method=c("BFGS"),
                                                                   control=list(fnscale=-1, maxit=maxit)), cl=cl)
  parallel::stopCluster(cl=cl)

  converged <- vapply(1:ncalls, function(i1) NEWTONresults[[i1]]$convergence == 0, logical(1))

  loks <- vapply(1:ncalls, function(i1) loglikelihood_int(data=data, p=p, M=M, params=NEWTONresults[[i1]]$par, conditional=conditional,
                                                          constraints=constraints, parametrization=parametrization, check_params=TRUE,
                                                          to_return="loglik", minval=minval), numeric(1))
  if(print_res == TRUE) {
    cat("Results from variable metric algorithm:\n")
    cat(paste("lowest value: ", round(min(loks), 3)), "\n")
    cat(paste("mean value:   ", round(mean(loks), 3)), "\n")
    cat(paste("largest value:", round(max(loks), 3)), "\n")
  }


  ### Obtain estimates and standard errors, calculate IC ###
  all_estimates <- lapply(NEWTONresults, function(x) x$par)
  best_fit <- NEWTONresults[[which(loks == max(loks))[1]]]
  params <- best_fit$par
  if(is.null(constraints)) {
    params <- sort_components(p=p, M=M, d=d, params=params)
    all_estimates <- lapply(all_estimates, function(pars) sort_components(p=p, M=M, d=d, params=pars))
  }
  if(best_fit$convergence == 1) {
    message("Iteration limit was reached when estimating the best fitting individual! Consider further estimations with the function 'iterate_more()'")
  }
  mixing_weights <- loglikelihood_int(data=data, p=p, M=M, params=params, conditional=conditional,
                                      parametrization=parametrization, constraints=constraints,
                                      to_return="mw", check_params=TRUE, minval=NULL)
  if(any(vapply(1:M, function(i1) sum(mixing_weights[,i1] > red_criteria[1]) < red_criteria[2]*n_obs, logical(1)))) {
    message("At least one of the mixture components in the estimated model seems to be wasted!")
  }


  ### Wrap up ###
  cat("Calculating approximate standard errors...\n")
  ret <- GMVAR(data=data, p=p, M=M, d=d, params=params, conditional=conditional, parametrization=parametrization,
               constraints=constraints, calc_std_errors=TRUE)
  ret$all_estimates <- all_estimates
  ret$all_logliks <- loks
  ret$which_converger <- converged

  cat("Finished!\n")
  ret
}


#' @title Maximum likelihood estimation of GMVAR model with preliminary estimates
#'
#' @description \code{iterate_more} uses variable metric algorithm to finalize maximum likelihood
#'  estimation of GMVAR model (object of class \code{'gmvar'}) which already has preliminary estimates.
#'
#' @inheritParams simulateGMVAR
#' @inheritParams fitGMVAR
#' @inheritParams GMVAR
#' @details The main purpose of \code{iterate_more} is to provide a simple and convenient tool to finalize
#'   the estimation when the maximum number of iterations is reached when estimating a GMVAR model with the
#'   main estimation function \code{fitGMVAR}. It's just a simple wrapper around function \code{optim}
#'   from the package \code{stats} and \code{GMVAR} from the package \code{gmvarkit}.
#' @return Returns an object of class \code{'gmvar'} defining the estimated GMVAR model. Can be used
#'   to work with other functions provided in \code{gmvarkit}.
#' @seealso \code{\link{fitGMVAR}}, \code{\link{GMVAR}}, \code{\link[stats]{optim}}
#' @inherit GMVAR references
#' @examples
#' \donttest{
#' ## These are long running examples that use parallel computing!
#'
#' # These examples use the data 'eurusd' which comes with the
#' # package, but in a scaled form.
#' data <- cbind(10*eurusd[,1], 100*eurusd[,2])
#' colnames(data) <- colnames(eurusd)
#'
#' # GMVAR(1,2) model, only 5 iterations of the variable metric
#' # algorithm
#' fit12 <- fitGMVAR(data, p=1, M=2, maxit=5)
#' fit12
#'
#' # Iterate more:
#' fit12_2 <- iterate_more(fit12)
#' fit12_2
#'
#'
#' # GMVAR(2,2) model with autoregressive parameters restricted
#' # to be the same for all regimes, only 10 iterations of the
#' # variable metric algorithm
#' C_mat <- rbind(diag(2*2^2), diag(2*2^2))
#' fit22c <- fitGMVAR(data, p=2, M=2, constraints=C_mat, maxit=10)
#' fit22c
#'
#' # Iterate more:
#' fit22c_2 <- iterate_more(fit22c)
#' fit22c_2
#'
#' # GMVAR(3,2) model, only 10 iterations of the variable metric
#' # algorithm
#' fit32 <- fitGMVAR(data, p=3, M=2, maxit=10)
#' fit32
#'
#' # Iterate more:
#' fit32_2 <- iterate_more(fit32)
#' fit32_2
#' }
#' @export

iterate_more <- function(gmvar, maxit=100, calc_std_errors=TRUE) {
  check_gmvar(gmvar)
  stopifnot(maxit %% 1 == 0 & maxit >= 1)
  minval <- -(10^(ceiling(log10(nrow(gmvar$data))) + ncol(gmvar$data) + 1) - 1)

  fn <- function(params) {
    tryCatch(loglikelihood_int(data=gmvar$data, p=gmvar$model$p, M=gmvar$model$M, params=params,
                               conditional=gmvar$model$conditional, parametrization=gmvar$model$parametrization,
                               constraints=gmvar$model$constraints, check_params=TRUE, to_return="loglik",
                               minval=minval), error=function(e) minval)
  }
  gr <- function(params) {
    calc_gradient(x=params, fn=fn)
  }

  res <- optim(par=gmvar$params, fn=fn, gr=gr, method=c("BFGS"), control=list(fnscale=-1, maxit=maxit))
  if(res$convergence == 1) message("The maximum number of iterations was reached! Consired iterating more.")

  GMVAR(data=gmvar$data, p=gmvar$model$p, M=gmvar$model$M, params=res$par,
        conditional=gmvar$model$conditional, parametrization=gmvar$model$parametrization,
        constraints=gmvar$model$constraints, calc_std_errors=calc_std_errors)
}
