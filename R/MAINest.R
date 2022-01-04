#' @title Two-phase maximum likelihood estimation of a GMVAR, StMVAR, or G-StMVAR model
#'
#' @description \code{fitGSMVAR} estimates a GMVAR, StMVAR, or G-StMVAR model model in two phases:
#'   in the first phase it uses a genetic algorithm to find starting values for a gradient based
#'   variable metric algorithm, which it then uses to finalize the estimation in the second phase.
#'   Parallel computing is utilized to perform multiple rounds of estimations in parallel.
#'
#' @inheritParams GAfit
#' @param ncalls the number of estimation rounds that should be performed.
#' @param ncores the number CPU cores to be used in parallel computing.
#' @param maxit the maximum number of iterations in the variable metric algorithm.
#' @param seeds a length \code{ncalls} vector containing the random number generator seed for each call to the genetic algorithm,
#'  or \code{NULL} for not initializing the seed. Exists for creating reproducible results.
#' @param print_res should summaries of estimation results be printed?
#' @param ... additional settings passed to the function \code{GAfit} employing the genetic algorithm.
#' @details
#'  If you wish to estimate a structural model without overidentifying constraints that is identified statistically,
#'  specify your W matrix is \code{structural_pars} to be such that it contains the same sign constraints in a single row
#'  (e.g. a row of ones) and leave the other elements as \code{NA}. In this way, the genetic algorithm works the best.
#'  The ordering and signs of the columns of the W matrix can be changed afterwards with the functions
#'  \code{reorder_W_columns} and \code{swap_W_signs}.
#'
#'  Because of complexity and high multimodality of the log-likelihood function, it's \strong{not certain} that the estimation
#'  algorithms will end up in the global maximum point. It's expected that most of the estimation rounds will end up in
#'  some local maximum or saddle point instead. Therefore, a (sometimes large) number of estimation rounds is required
#'  for reliable results. Because of the nature of the model, the estimation may fail especially in the cases where the
#'  number of mixture components is chosen too large. \strong{With two regimes and couple hundred observations in a two-dimensional
#'  time series, 50 rounds is usually enough. Several hundred estimation rounds often suffices for reliably fitting two-regimes
#'  models to 3 or 4 dimensional time series. With more than two regimes and more than couple hundred
#'  observations, thousands of estimation rounds (or more) are often required to obtain reliable results.}
#'
#'  The estimation process is computationally heavy and it might take considerably long time for large models with
#'  large number of observations. If the iteration limit \code{maxit} in the variable metric algorithm is reached,
#'  one can continue the estimation by iterating more with the function \code{iterate_more}. Alternatively, one may
#'  use the found estimates as starting values for the genetic algorithm and and employ another round of estimation
#'  (see \code{?GAfit} how to set up an initial population with the dot parameters).
#'
#'  \strong{If the estimation algorithm fails to create an initial population for the genetic algorithm},
#'  it usually helps to scale the individual series so that the AR coefficients (of a VAR model) will be
#'  relative small, preferably less than one. Even if one is able to create an initial population, it should
#'  be preferred to scale the series so that most of the AR coefficients will not be very large, as the
#'  estimation algorithm works better with relatively small AR coefficients. If needed, another package can be used
#'  to fit linear VARs to the series to see which scaling of the series results in relatively small AR coefficients.
#'
#'  The code of the genetic algorithm is mostly based on the description by \emph{Dorsey and Mayer (1995)} but it
#'  includes some extra features that were found useful for this particular estimation problem. For instance,
#'  the genetic algorithm uses a slightly modified version of the individually adaptive crossover and mutation
#'  rates described by \emph{Patnaik and Srinivas (1994)} and employs (50\%) fitness inheritance discussed
#'  by \emph{Smith, Dike and Stegmann (1995)}.
#'
#'  The gradient based variable metric algorithm used in the second phase is implemented with function \code{optim}
#'  from the package \code{stats}.
#'
#'  Note that the structural models are even more difficult to estimate than the reduced form models due to
#'  the different parametrization of the covariance matrices, so larger number of estimation rounds should be considered.
#'  Also, be aware that if the lambda parameters are constrained in any other way than by restricting some of them to be
#'  identical, the parameter "lambda_scale" of the genetic algorithm (see \code{?GAfit}) needs to be carefully adjusted accordingly.
#'  \strong{When estimating a structural model that imposes overidentifiying constraints to a time series with \eqn{d>3},
#'  it is highly recommended to create an initial population based on the estimates of a statistically identified model
#'  (when \eqn{M=2}). This is because currently obtaining the ML estimate reliably to such a structural model seems
#'  difficult in many application.}
#'
#'  Finally, the function fails to calculate approximate standard errors and the parameter estimates are near the border
#'  of the parameter space, it might help to use smaller numerical tolerance for the stationarity and positive
#'  definiteness conditions. The numerical tolerance of an existing model can be changed with the function
#'  \code{update_numtols}.
#' @return Returns an object of class \code{'gsmvar'} defining the estimated (reduced form or structural) GMVAR, StMVAR, or G-StMVAR model.
#'   Multivariate quantile residuals (Kalliovirta and Saikkonen 2010) are also computed and included in the returned object.
#'   In addition, the returned object contains the estimates and log-likelihood values from all the estimation rounds performed.
#'   The estimated parameter vector can be obtained at \code{gsmvar$params} (and corresponding approximate standard errors
#'   at \code{gsmvar$std_errors}). See \code{?GSMVAR} for the form of the parameter vector, if needed.
#'
#'   Remark that the first autocovariance/correlation matrix in \code{$uncond_moments} is for the lag zero,
#'   the second one for the lag one, etc.
#' @section S3 methods:
#'   The following S3 methods are supported for class \code{'gsmvar'}: \code{logLik}, \code{residuals}, \code{print}, \code{summary},
#'    \code{predict}, \code{simulate}, and \code{plot}.
#' @seealso \code{\link{GSMVAR}}, \code{\link{iterate_more}}, \code{\link{stmvar_to_gstmvar}}, \code{\link{predict.gsmvar}},
#'   \code{\link{profile_logliks}}, \code{\link{simulate.gsmvar}}, \code{\link{quantile_residual_tests}}, \code{\link{print_std_errors}},
#'   \code{\link{swap_parametrization}}, \code{\link{get_gradient}}, \code{\link{GIRF}}, \code{\link{GFEVD}}, \code{\link{LR_test}},
#'   \code{\link{Wald_test}}, \code{\link{gsmvar_to_sgsmvar}}, \code{\link{stmvar_to_gstmvar}}, \code{\link{reorder_W_columns}},
#'    \code{\link{swap_W_signs}}, \code{\link{cond_moment_plot}}, \code{\link{update_numtols}}
#' @references
#'  \itemize{
#'    \item Dorsey R. E. and Mayer W. J. 1995. Genetic algorithms for estimation problems with multiple optima,
#'          nondifferentiability, and other irregular features. \emph{Journal of Business & Economic Statistics},
#'         \strong{13}, 53-66.
#'    \item Kalliovirta L., Meitz M. and Saikkonen P. 2016. Gaussian mixture vector autoregression.
#'          \emph{Journal of Econometrics}, \strong{192}, 485-498.
#'    \item Patnaik L.M. and Srinivas M. 1994. Adaptive Probabilities of Crossover and Mutation in Genetic Algorithms.
#'          \emph{Transactions on Systems, Man and Cybernetics} \strong{24}, 656-667.
#'    \item Smith R.E., Dike B.A., Stegmann S.A. 1995. Fitness inheritance in genetic algorithms.
#'          \emph{Proceedings of the 1995 ACM Symposium on Applied Computing}, 345-350.
#'    \item Virolainen S. 2021. Structural Gaussian mixture vector autoregressive model. Unpublished working
#'      paper, available as arXiv:2007.04713.
#'    \item Virolainen S. 2021. Gaussian and Student's t mixture vector autoregressive model. Unpublished working
#'      paper, available as arXiv:2109.13648.
#'  }
#' @examples
#' \donttest{
#' ## These are long running examples that use parallel computing!
#' # Running all the below examples will take approximately 3-4 minutes.
#'
#' # GMVAR(1,2) model: 10 estimation rounds with seeds set
#' # for reproducibility
#' fit12 <- fitGSMVAR(gdpdef, p=1, M=2, ncalls=10, seeds=1:10)
#' fit12
#' plot(fit12)
#' summary(fit12)
#' print_std_errors(fit12)
#' profile_logliks(fit12)
#'
#' # The rest of the examples only use a single estimation round with a given
#' # seed that produces the MLE to reduce running time of the examples. When
#' # estimating models for empirical applications, a large number of estimation
#' # rounds (ncalls = a large number) should be performed to ensure reliability
#' # of the estimates (see the section details).
#'
#' # StMVAR(2, 2) model
#' fit22t <- fitGSMVAR(gdpdef, p=2, M=2, model="StMVAR", ncalls=1, seeds=1)
#' fit22t # Overly large degrees of freedom estimate in the 2nd regime!
#' fit22gs <- stmvar_to_gstmvar(fit22t) # So switch it to GMVAR type!
#' fit22gs # This is the appropriate G-StMVAR model based on the above StMVAR model.
#' fit22gss <- gsmvar_to_sgsmvar(fit22gs) # Switch to structural model
#' fit22gss # This is the implied statistically identified structural model.
#'
#' # Structural GMVAR(1,2) model identified with sign
#' # constraints.
#' W_122 <- matrix(c(1, 1, -1, 1), nrow=2)
#' fit12s <- fitGSMVAR(gdpdef, p=1, M=2, structural_pars=list(W=W_122),
#'   ncalls=1, seeds=1)
#' fit12s
#' # A statistically identified structural model can also be obtained with
#' # gsmvar_to_sgsmvar(fit12)
#'
#'
#' # GMVAR(2,2) model with autoregressive parameters restricted
#' # to be the same for both regimes
#' C_mat <- rbind(diag(2*2^2), diag(2*2^2))
#' fit22c <- fitGSMVAR(gdpdef, p=2, M=2, constraints=C_mat, ncalls=1, seeds=1)
#' fit22c
#'
#' fit22gscm <- fitGSMVAR(gdpdef, p=2, M=c(1, 1), model="G-StMVAR", constraints=C_mat,
#'    parametrization="mean", same_means=list(1:2), ncalls=1, seeds=1)
#'
#' # GMVAR(2,2) model with autoregressive parameters restricted
#' # to be the same for both regimes and non-diagonal elements
#' # the coefficient matrices constrained to zero.
#' tmp <- matrix(c(1, rep(0, 10), 1, rep(0, 8), 1, rep(0, 10), 1),
#'  nrow=2*2^2, byrow=FALSE)
#' C_mat2 <- rbind(tmp, tmp)
#' fit22c2 <- fitGSMVAR(gdpdef, p=2, M=2, constraints=C_mat2, ncalls=1,
#'   seeds=1)
#' fit22c2
#' }
#' @export

fitGSMVAR <- function(data, p, M, model=c("GMVAR", "StMVAR", "G-StMVAR"), conditional=TRUE, parametrization=c("intercept", "mean"),
                      constraints=NULL, same_means=NULL, structural_pars=NULL, ncalls=(M + 1)^5, ncores=2, maxit=1000,
                      seeds=NULL, print_res=TRUE, ...) {

  model <- match.arg(model)
  parametrization <- match.arg(parametrization)
  check_pMd(p=p, M=M, model=model)
  if(!all_pos_ints(c(ncalls, ncores, maxit))) stop("Arguments ncalls, ncores, and maxit must be positive integers")
  stopifnot(length(ncalls) == 1)
  if(!is.null(seeds) && length(seeds) != ncalls) stop("The argument 'seeds' should be NULL or a vector of length 'ncalls'")
  data <- check_data(data=data, p=p)
  d <- ncol(data)
  n_obs <- nrow(data)
  check_same_means(parametrization=parametrization, same_means=same_means)
  check_constraints(p=p, M=M, d=d, constraints=constraints, same_means=same_means, structural_pars=structural_pars)
  npars <- n_params(p=p, M=M, d=d, model=model, constraints=constraints, same_means=same_means, structural_pars=structural_pars)
  if(npars >= d*nrow(data)) stop("There are at least as many parameters in the model as there are observations in the data")
  dot_params <- list(...)
  minval <- ifelse(is.null(dot_params$minval), get_minval(data), dot_params$minval)
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

  ### Optimization with the genetic algorithm ###
  cl <- parallel::makeCluster(ncores)
  on.exit(try(parallel::stopCluster(cl), silent=TRUE)) # Close the cluster on exit, if not already closed.
  parallel::clusterExport(cl, ls(environment(fitGSMVAR)), envir = environment(fitGSMVAR)) # assign all variables from package:gmvarkit
  parallel::clusterEvalQ(cl, c(library(Brobdingnag), library(mvnfast), library(pbapply)))

  cat("Optimizing with a genetic algorithm...\n")
  GAresults <- pbapply::pblapply(1:ncalls, function(i1) GAfit(data=data, p=p, M=M, model=model, conditional=conditional,
                                                              parametrization=parametrization, constraints=constraints,
                                                              same_means=same_means, structural_pars=structural_pars,
                                                              seed=seeds[i1], ...), cl=cl)

  loks <- vapply(1:ncalls, function(i1) loglikelihood_int(data=data, p=p, M=M, params=GAresults[[i1]], model=model,
                                                          conditional=conditional, parametrization=parametrization,
                                                          constraints=constraints, same_means=same_means,
                                                          structural_pars=structural_pars, check_params=TRUE,
                                                          to_return="loglik", minval=minval), numeric(1))

  if(print_res) {
    print_loks <- function() {
      printfun <- function(txt, FUN) cat(paste(txt, round(FUN(loks), 3)), "\n")
      printfun("The lowest loglik: ", min)
      printfun("The mean loglik:   ", mean)
      printfun("The largest loglik:", max)
    }
    cat("Results from the genetic algorithm:\n")
    print_loks()
  }

  ### Optimization with the variable metric algorithm ###

  # Logarithmize degrees of freedom parameters to get overly large degrees of freedom parameters
  # value to the same range as other parameters. This adjusts the difference 'h' to be larger
  # for larger df parameters in non-log scale to avoid the numerical problems associated with overly
  # large degrees of freedom parameters.
  manipulateDFS <- function(M, params, model, FUN) { # The function to log/exp the dfs
    FUN <- match.fun(FUN)
    M2 <- ifelse(model == "StMVAR", M, M[2])
    params[(npars - M2 + 1):npars] <- FUN(params[(npars - M2 + 1):npars])
    params
  }
  if(model == "StMVAR" | model == "G-StMVAR") { # Logarithmize the degrees of freedom parameters
    GAresults <- lapply(1:ncalls, function(i1) manipulateDFS(M=M, params=GAresults[[i1]], model=model, FUN=log))
  }

  loglik_fn <- function(params) {
    if(model == "StMVAR" | model == "G-StMVAR") {
      params <- manipulateDFS(M=M, params=params, model=model, FUN=exp) # Unlogarithmize dfs for calculating log-likelihood
    }
    tryCatch(loglikelihood_int(data=data, p=p, M=M, params=params, model=model,
                               conditional=conditional, parametrization=parametrization,
                               constraints=constraints, same_means=same_means,
                               structural_pars=structural_pars, check_params=TRUE,
                               to_return="loglik", minval=minval), error=function(e) minval)
  }

  h <- 6e-6
  I <- diag(rep(1, npars))
  loglik_grad <- function(params) {
    vapply(1:npars, function(i1) (loglik_fn(params + I[i1,]*h) - loglik_fn(params - I[i1,]*h))/(2*h), numeric(1))
  }

  cat("Optimizing with a variable metric algorithm...\n")
  NEWTONresults <- pbapply::pblapply(1:ncalls, function(i1) optim(par=GAresults[[i1]], fn=loglik_fn, gr=loglik_grad, method="BFGS",
                                                                  control=list(fnscale=-1, maxit=maxit)), cl=cl)
  parallel::stopCluster(cl=cl)

  loks <- vapply(1:ncalls, function(i1) NEWTONresults[[i1]]$value, numeric(1)) # Log-likelihoods
  converged <- vapply(1:ncalls, function(i1) NEWTONresults[[i1]]$convergence == 0, logical(1)) # Which coverged

  if(print_res) {
    cat("Results from the variable metric algorithm:\n")
    print_loks()
  }


  ### Obtain estimates and standard errors, calculate IC ###
  all_estimates <- lapply(NEWTONresults, function(x) x$par)
  if(model == "StMVAR" || model == "G-StMVAR") { # Unlogarithmize degrees of freedom parameter values
    all_estimates <- lapply(1:ncalls, function(i1) manipulateDFS(M=M, params=all_estimates[[i1]], model=model, FUN=exp))
  }
  which_best_fit <- which(loks == max(loks))[1]
  best_fit <- all_estimates[[which_best_fit]]
  params <- best_fit
  if(is.null(constraints) && is.null(structural_pars$C_lambda) && is.null(same_means)) {
    params <- sort_components(p=p, M=M, d=d, params=params, model=model, structural_pars=structural_pars)
    all_estimates <- lapply(all_estimates, function(pars) sort_components(p=p, M=M, d=d, params=pars, model=model, structural_pars=structural_pars))
  }
  if(NEWTONresults[[which_best_fit]]$convergence == 1) {
    message("Iteration limit was reached when estimating the best fitting individual! Consider further estimation with the function 'iterate_more'")
  }
  mixing_weights <- loglikelihood_int(data=data, p=p, M=M, params=params, model=model,
                                      conditional=conditional, parametrization=parametrization,
                                      constraints=constraints, same_means=same_means,
                                      structural_pars=structural_pars, to_return="mw",
                                      check_params=TRUE, minval=NULL)
  if(any(vapply(1:sum(M), function(i1) sum(mixing_weights[,i1] > red_criteria[1]) < red_criteria[2]*n_obs, logical(1)))) {
    message("At least one of the mixture components in the estimated model seems to be wasted!")
  }


  ### Wrap up ###
  cat("Calculating approximate standard errors...\n")
  ret <- GSMVAR(data=data, p=p, M=M, d=d, params=params, model=model,
                conditional=conditional, parametrization=parametrization,
                constraints=constraints, same_means=same_means,
                structural_pars=structural_pars, calc_std_errors=TRUE)
  ret$all_estimates <- all_estimates
  ret$all_logliks <- loks
  ret$which_converged <- converged
  ret$which_round <- which_best_fit # Which estimation round induced the largest log-likelihood?
  warn_eigens(ret)
  cat("Finished!\n")
  ret
}


#' @title Maximum likelihood estimation of a GMVAR, StMVAR, or G-StMVAR model with preliminary estimates
#'
#' @description \code{iterate_more} uses a variable metric algorithm to finalize maximum likelihood
#'  estimation of a GMVAR, StMVAR, or G-StMVAR model (object of class \code{'gsmvar'}) which already has preliminary estimates.
#'
#' @inheritParams quantile_residual_tests
#' @inheritParams fitGSMVAR
#' @inheritParams GSMVAR
#' @inheritParams standard_errors
#' @details The purpose of \code{iterate_more} is to provide a simple and convenient tool to finalize
#'   the estimation when the maximum number of iterations is reached when estimating a GMVAR, StMVAR, or G-StMVAR model
#'   with the main estimation function \code{fitGSMVAR}. \code{iterate_more} is essentially a wrapper
#'   around the function \code{optim} from the package \code{stats} and \code{GSMVAR} from the package
#'   \code{gmvarkit}.
#' @return Returns an object of class \code{'gsmvar'} defining the estimated GMVAR, StMVAR, or G-StMVAR model.
#' @seealso \code{\link{fitGSMVAR}}, \code{\link{GSMVAR}}, \code{\link[stats]{optim}},
#'  \code{\link{profile_logliks}}, \code{\link{update_numtols}}
#' @inherit GSMVAR references
#' @examples
#' \donttest{
#' ## These are long running examples that use parallel computing!
#' ## Running the below examples takes approximately 2 minutes
#'
#' # GMVAR(1,2) model, only 5 iterations of the variable metric
#' # algorithm
#' fit12 <- fitGSMVAR(gdpdef, p=1, M=2, ncalls=1, maxit=5, seeds=1)
#' fit12
#'
#' # Iterate more:
#' fit12_2 <- iterate_more(fit12)
#' fit12_2
#' }
#' @export

iterate_more <- function(gsmvar, maxit=100, calc_std_errors=TRUE, custom_h=NULL,
                         stat_tol=1e-3, posdef_tol=1e-8, df_tol=1e-8) {
  gsmvar <- gmvar_to_gsmvar(gsmvar) # Backward compatibility
  check_gsmvar(gsmvar)
  stopifnot(maxit %% 1 == 0 & maxit >= 1)
  if(is.null(custom_h)) { # Adjust h for overly large degrees of freedom parameters
    varying_h <- get_varying_h(M=gsmvar$model$M, params=gsmvar$params, model=gsmvar$model$model)
  } else { # Utilize user-specified h
    stopifnot(length(custom_h) == length(gsmvar$params))
    varying_h <- custom_h
  }
  minval <- get_minval(gsmvar$data)

  fn <- function(params) {
    tryCatch(loglikelihood_int(data=gsmvar$data, p=gsmvar$model$p, M=gsmvar$model$M, params=params, model=gsmvar$model$model,
                               conditional=gsmvar$model$conditional, parametrization=gsmvar$model$parametrization,
                               constraints=gsmvar$model$constraints, same_means=gsmvar$model$same_means,
                               structural_pars=gsmvar$model$structural_pars, check_params=TRUE,
                               to_return="loglik", minval=minval,
                               stat_tol=stat_tol, posdef_tol=posdef_tol, df_tol=df_tol),
             error=function(e) minval)
  }
  gr <- function(params) {
    calc_gradient(x=params, fn=fn, varying_h=varying_h)
  }

  res <- optim(par=gsmvar$params, fn=fn, gr=gr, method=c("BFGS"), control=list(fnscale=-1, maxit=maxit))
  if(res$convergence == 1) message("The maximum number of iterations was reached! Consired iterating more.")

  ret <- GSMVAR(data=gsmvar$data, p=gsmvar$model$p, M=gsmvar$model$M, params=res$par, model=gsmvar$model$model,
                conditional=gsmvar$model$conditional, parametrization=gsmvar$model$parametrization,
                constraints=gsmvar$model$constraints, same_means=gsmvar$model$same_means,
                structural_pars=gsmvar$model$structural_pars, calc_std_errors=calc_std_errors,
                stat_tol=stat_tol, posdef_tol=posdef_tol, df_tol=df_tol)

  ret$all_estimates <- gsmvar$all_estimates
  ret$all_logliks <- gsmvar$all_logliks
  ret$which_converged <- gsmvar$which_converged
  if(!is.null(gsmvar$which_round)) {
    ret$which_round <- gsmvar$which_round
    ret$all_estimates[[gsmvar$which_round]] <- ret$params
    ret$all_logliks[gsmvar$which_round] <- ret$loglik
    ret$which_converged[gsmvar$which_round] <- res$convergence == 0
  }
  warn_eigens(ret)
  ret
}


#' @title Returns the default smallest allowed log-likelihood for given data.
#'
#' @description \code{get_minval} returns the default smallest allowed log-likelihood for given data.
#'
#' @inheritParams GAfit
#' @details This function exists to avoid dublication inside the package.
#' @return Returns \code{-(10^(ceiling(log10(nrow(data)) + ncol(data))) - 1)}
#' @seealso \code{\link{fitGSMVAR}}, \code{\link{GAfit}}
#' @keywords internal

get_minval <- function(data) {
  -(10^(ceiling(log10(nrow(data)) + ncol(data))) - 1)
}



