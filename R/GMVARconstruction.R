#' @title Create a class 'gmvar' object defining a reduced form or structural GMVAR model
#'
#' @description \code{GMVAR} creates a class \code{'gmvar'} object that defines
#'  a reduced form or structural GMVAR model
#'
#' @inheritParams loglikelihood_int
#' @param data a matrix or class \code{'ts'} object with \code{d>1} columns. Each column is taken to represent
#'  a single times series. \code{NA} values are not supported. Ignore if defining a model without data is desired.
#' @param d number of times series in the system, i.e. \code{ncol(data)}. This can be
#'   used to define GMVAR models without data and can be ignored if \code{data} is provided.
#' @param calc_cond_moments should conditional means and covariance matrices should be calculated?
#'   Default is \code{TRUE} if the model contains data and \code{FALSE} otherwise.
#' @param calc_std_errors should approximate standard errors be calculated?
#' @details If data is provided, then also multivariate quantile residuals (\emph{Kalliovirta and Saikkonen 2010})
#'   are computed and included in the returned object.
#'
#'   If the function fails to calculate approximative standard errors and the parameter values are near the border
#'   of the parameter space, it might help to use smaller numerical tolerance for the stationarity and positive
#'   definiteness conditions.
#' @return Returns an object of class \code{'gmvar'} defining the specified reduced form or structural GMVAR model.
#'   Can be used to work with other functions provided in \code{gmvarkit}.
#'
#'   Remark that the first autocovariance/correlation matrix in \code{$uncond_moments} is for the lag zero,
#'   the second one for the lag one, etc.
#' @section About S3 methods:
#'   Only the \code{print} method is available if data is not provided.
#'   If data is provided, then in addition to the ones listed above, the \code{predict} method is also available.
#' @seealso \code{\link{fitGMVAR}}, \code{\link{add_data}}, \code{\link{swap_parametrization}}, \code{\link{GIRF}},
#'   \code{\link{gmvar_to_sgmvar}}, \code{\link{reorder_W_columns}}, \code{\link{swap_W_signs}}, \code{\link{update_numtols}}
#' @references
#'  \itemize{
#'    \item Kalliovirta L., Meitz M. and Saikkonen P. 2016. Gaussian mixture vector autoregression.
#'          \emph{Journal of Econometrics}, \strong{192}, 485-498.
#'    \item Kalliovirta L. and Saikkonen P. 2010. Reliable Residuals for Multivariate Nonlinear
#'          Time Series Models. \emph{Unpublished Revision of HECER Discussion Paper No. 247}.
#'    \item Virolainen S. 2020. Structural Gaussian mixture vector autoregressive model. Unpublished working
#'      paper, available as arXiv:2007.04713.
#'  }
#' @examples
#' # GMVAR(1, 2), d=2 model:
#' params12 <- c(0.55, 0.112, 0.344, 0.055, -0.009, 0.718, 0.319, 0.005,
#'   0.03, 0.619, 0.173, 0.255, 0.017, -0.136, 0.858, 1.185, -0.012,
#'   0.136, 0.674)
#' mod12 <- GMVAR(gdpdef, p=1, M=2, params=params12)
#' mod12
#'
#' # GMVAR(1, 2), d=2 model without data
#' mod12_2 <- GMVAR(p=1, M=2, d=2, params=params12)
#' mod12_2
#'
#' # GMVAR(2, 2), d=2 model with mean-parametrization:
#' params22 <- c(0.869, 0.549, 0.223, 0.059, -0.151, 0.395, 0.406,
#'  -0.005, 0.083, 0.299, 0.215, 0.002, 0.03, 0.576, 1.168, 0.218,
#'  0.02, -0.119, 0.722, 0.093, 0.032, 0.044, 0.191, 1.101, -0.004,
#'  0.105, 0.58)
#' mod22 <- GMVAR(gdpdef, p=2, M=2, params=params22, parametrization="mean")
#' mod22
#'
#' # Structural GMVAR(2, 2), d=2 model identified with sign-constraints:
#' params22s <- c(0.36, 0.121, 0.484, 0.072, 0.223, 0.059, -0.151, 0.395,
#'   0.406, -0.005, 0.083, 0.299, 0.218, 0.02, -0.119, 0.722, 0.093, 0.032,
#'   0.044, 0.191, 0.057, 0.172, -0.46, 0.016, 3.518, 5.154, 0.58)
#' W_22 <- matrix(c(1, 1, -1, 1), nrow=2, byrow=FALSE)
#' mod22s <- GMVAR(gdpdef, p=2, M=2, params=params22s,
#'  structural_pars=list(W=W_22))
#' mod22s
#' @export

GMVAR <- function(data, p, M, d, params, conditional=TRUE, parametrization=c("intercept", "mean"), constraints=NULL,
                  same_means=NULL, structural_pars=NULL, calc_cond_moments, calc_std_errors=FALSE,
                  stat_tol=1e-3, posdef_tol=1e-8) {

  parametrization <- match.arg(parametrization)
  if(missing(calc_cond_moments)) calc_cond_moments <- ifelse(missing(data) || is.null(data), FALSE, TRUE)
  if(!all_pos_ints(c(p, M))) stop("Arguments p and M must be positive integers")
  if(missing(data) & missing(d)) stop("data or d must be provided")
  if(missing(data) || is.null(data)) {
    data <- NULL
  } else {
    data <- check_data(data=data, p=p)
    if(missing(d)) {
      d <- ncol(data)
    } else if(ncol(data) != d) {
      warning("ncol(data) does not equal d. Using d = ncol(data)")
      d <- ncol(data)
    }
  }
  check_constraints(p=p, M=M, d=d, constraints=constraints, same_means=same_means, structural_pars=structural_pars)
  check_parameters(p=p, M=M, d=d, params=params, parametrization=parametrization, constraints=constraints,
                   same_means=same_means, structural_pars=structural_pars,
                   stat_tol=stat_tol, posdef_tol=posdef_tol)
  npars <- n_params(p=p, M=M, d=d, constraints=constraints, same_means=same_means, structural_pars=structural_pars)

  if(is.null(data)) {
    lok_and_mw <- list(loglik=NA, mw=NA)
    IC <- data.frame(AIC=NA, HQIC=NA, BIC=NA)
    qresiduals <- NA
  } else {
    if(npars >= d*nrow(data)) stop("There are at least as many parameters in the model than there are observations in the data")
    lok_and_mw <- loglikelihood_int(data=data, p=p, M=M, params=params,
                                    conditional=conditional, parametrization=parametrization,
                                    constraints=constraints, same_means=same_means,
                                    structural_pars=structural_pars,
                                    to_return="loglik_and_mw",
                                    check_params=FALSE, minval=NA,
                                    stat_tol=stat_tol, posdef_tol=posdef_tol)
    qresiduals <- quantile_residuals_int(data=data, p=p, M=M, params=params,
                                         conditional=conditional, parametrization=parametrization,
                                         constraints=constraints, same_means=same_means,
                                         structural_pars=structural_pars, stat_tol=stat_tol,
                                         posdef_tol=posdef_tol)
    obs <- ifelse(conditional, nrow(data) - p, nrow(data))
    IC <- get_IC(loglik=lok_and_mw$loglik, npars=npars, obs=obs)
  }
  if(calc_std_errors) {
    if(is.null(data)) {
      warning("Approximate standard errors can't be calculated")
      std_errors <- rep(NA, npars)
    } else {
      std_errors <- tryCatch(standard_errors(data=data, p=p, M=M, params=params,
                                             conditional=conditional, parametrization=parametrization,
                                             constraints=constraints, same_means=same_means,
                                             structural_pars=structural_pars,
                                             minval=-(10^(ceiling(log10(nrow(data))) + ncol(data) + 1) - 1),
                                             stat_tol=stat_tol, posdef_tol=posdef_tol),
                             error=function(e) {
                               warning("Approximate standard errors can't be calculated")
                               std_errors=rep(NA, npars)
                             })
    }
  } else {
    std_errors <- rep(NA, npars)
  }
  if(calc_cond_moments == FALSE || is.null(data)) {
    if(calc_cond_moments) warning("Conditional moments can't be calculated without data")
    regime_cmeans <- NA
    total_cmeans <- NA
    total_ccovs <- NA
  } else {
    get_cm <- function(to_return) loglikelihood_int(data=data, p=p, M=M, params=params,
                                                    conditional=conditional, parametrization=parametrization,
                                                    constraints=constraints, same_means=same_means,
                                                    structural_pars=structural_pars,
                                                    check_params=TRUE,
                                                    to_return=to_return, minval=NA,
                                                    stat_tol=stat_tol, posdef_tol=posdef_tol)
    regime_cmeans <- get_cm("regime_cmeans")
    total_cmeans <- get_cm("total_cmeans")
    total_ccovs <- get_cm("total_ccovs")
  }

  structure(list(data=data,
                 model=list(p=p,
                            M=M,
                            d=d,
                            conditional=conditional,
                            parametrization=parametrization,
                            constraints=constraints,
                            same_means=same_means,
                            structural_pars=structural_pars),
                 params=params,
                 std_errors=std_errors,
                 mixing_weights=lok_and_mw$mw,
                 regime_cmeans=regime_cmeans,
                 total_cmeans=total_cmeans,
                 total_ccovs=total_ccovs,
                 quantile_residuals=qresiduals,
                 loglik=structure(lok_and_mw$loglik,
                                  class="logLik",
                                  df=npars),
                 IC=IC,
                 uncond_moments=uncond_moments_int(p=p, M=M, d=d, params=params,
                                                   parametrization=parametrization,
                                                   constraints=constraints,
                                                   same_means=same_means,
                                                   structural_pars=structural_pars),
                 all_estimates=NULL,
                 all_logliks=NULL,
                 which_converged=NULL,
                 which_round=NULL,
                 num_tols=list(stat_tol=stat_tol,
                               posdef_tol=posdef_tol)),
            class="gmvar")
}


#' @title Add data to an object of class 'gmvar' defining a GMVAR model
#'
#' @description \code{add_data} adds or updates data to object of class '\code{gmvar}' that defines a GMVAR model.
#'  Also calculates mixing weights and quantile residuals accordingly.
#'
#' @inheritParams loglikelihood_int
#' @inheritParams simulateGMVAR
#' @inheritParams GMVAR
#' @return Returns an object of class 'gmvar' defining the specified GMVAR model with the data added to the model.
#'   If the object already contained data, the data will be updated.
#' @seealso \code{\link{fitGMVAR}}, \code{\link{GMVAR}}, \code{\link{iterate_more}}, \code{\link{update_numtols}}
#' @references
#'  \itemize{
#'    \item Kalliovirta L., Meitz M. and Saikkonen P. 2016. Gaussian mixture vector autoregression.
#'          \emph{Journal of Econometrics}, \strong{192}, 485-498.
#'    \item Virolainen S. 2020. Structural Gaussian mixture vector autoregressive model. Unpublished working
#'      paper, available as arXiv:2007.04713.
#'  }
#' @examples
#' # GMVAR(1,2), d=2 model:
#' params12 <- c(0.55, 0.112, 0.344, 0.055, -0.009, 0.718, 0.319, 0.005,
#'   0.03, 0.619, 0.173, 0.255, 0.017, -0.136, 0.858, 1.185, -0.012,
#'   0.136, 0.674)
#' mod12 <- GMVAR(p=1, M=2, d=2, params=params12)
#' mod12
#'
#' mod12_2 <- add_data(gdpdef, mod12)
#' mod12_2
#'
#'
#' # Structural GMVAR(2, 2), d=2 model identified with sign-constraints:
#' params22s <- c(0.36, 0.121, 0.484, 0.072, 0.223, 0.059, -0.151, 0.395,
#'   0.406, -0.005, 0.083, 0.299, 0.218, 0.02, -0.119, 0.722, 0.093, 0.032,
#'   0.044, 0.191, 0.057, 0.172, -0.46, 0.016, 3.518, 5.154, 0.58)
#' W_22 <- matrix(c(1, 1, -1, 1), nrow=2, byrow=FALSE)
#' mod22s <- GMVAR(p=2, M=2, d=2, params=params22s, structural_pars=list(W=W_22))
#' mod22s
#'
#' mod22s_2 <- add_data(gdpdef, mod22s)
#' mod22s_2
#' @export

add_data <- function(data, gmvar, calc_cond_moments=TRUE, calc_std_errors=FALSE) {
  check_gmvar(gmvar)
  GMVAR(data=data, p=gmvar$model$p, M=gmvar$model$M, params=gmvar$params, conditional=gmvar$model$conditional,
        parametrization=gmvar$model$parametrization, constraints=gmvar$model$constraints,
        same_means=gmvar$model$same_means, structural_pars=gmvar$model$structural_pars,
        calc_cond_moments=calc_cond_moments, calc_std_errors=calc_std_errors,
        stat_tol=gmvar$num_tols$stat_tol, posdef_tol=gmvar$num_tols$posdef_tol)
}


#' @title Swap the parametrization of a GMVAR model
#'
#' @description \code{swap_parametrization} swaps the parametrization of a GMVAR model
#'  to \code{"mean"} if the current parametrization is \code{"intercept"}, and vice versa.
#'
#' @inheritParams simulateGMVAR
#' @details \code{swap_parametrization} is a convenient tool if you have estimated the model in
#'  "intercept"-parametrization, but wish to work with "mean"-parametrization in the future, or vice versa.
#'  In \code{gmvarkit}, the approximate standard errors are only available for parametrized parameters.
#' @inherit GMVAR references return
#' @inherit add_data seealso
#' @examples
#' \donttest{
#' # GMVAR(2, 2), d=2 model with mean-parametrization:
#' params22 <- c(0.869, 0.549, 0.223, 0.059, -0.151, 0.395, 0.406,
#'  -0.005, 0.083, 0.299, 0.215, 0.002, 0.03, 0.576, 1.168, 0.218,
#'  0.02, -0.119, 0.722, 0.093, 0.032, 0.044, 0.191, 1.101, -0.004,
#'  0.105, 0.58)
#' mod22 <- GMVAR(gdpdef, p=2, M=2, params=params22, parametrization="mean")
#' mod22 # mean parametrization
#'
#' mod22_2 <- swap_parametrization(mod22)
#' mod22_2 # intercept parametrization
#'
#'
#' # Structural GMVAR(2, 2), d=2 model identified with sign-constraints:
#' params22s <- c(0.36, 0.121, 0.484, 0.072, 0.223, 0.059, -0.151, 0.395,
#'   0.406, -0.005, 0.083, 0.299, 0.218, 0.02, -0.119, 0.722, 0.093, 0.032,
#'   0.044, 0.191, 0.057, 0.172, -0.46, 0.016, 3.518, 5.154, 0.58)
#' W_22 <- matrix(c(1, 1, -1, 1), nrow=2, byrow=FALSE)
#' mod22s <- GMVAR(p=2, M=2, d=2, params=params22s, structural_pars=list(W=W_22))
#' mod22s # intercept parametrization
#'
#' mod22s_2 <- swap_parametrization(mod22s)
#' mod22s_2 # mean parametrization
#' }
#' @export

swap_parametrization <- function(gmvar) {
  check_gmvar(gmvar)
  if(!is.null(gmvar$model$same_means)) {
    stop("Cannot change parametrization to intercept is the mean parameters are constrained")
  }
  change_to <- ifelse(gmvar$model$parametrization == "intercept", "mean", "intercept")
  new_params <- change_parametrization(p=gmvar$model$p, M=gmvar$model$M, d=gmvar$model$d, params=gmvar$params,
                                       constraints=gmvar$model$constraints, structural_pars=gmvar$model$structural_pars,
                                       change_to=change_to)
  GMVAR(data=gmvar$data, p=gmvar$model$p, M=gmvar$model$M, d=gmvar$model$d, params=new_params, conditional=gmvar$model$conditional,
        parametrization=change_to, constraints=gmvar$model$constraints, structural_pars=gmvar$model$structural_pars,
        calc_std_errors=ifelse(is.null(gmvar$data), FALSE, TRUE), stat_tol=gmvar$num_tols$stat_tol,
        posdef_tol=gmvar$num_tols$posdef_tol)
}


#' @title Construct a GMVAR model based on results from an arbitrary estimation round of \code{fitGMVAR}
#'
#' @description \code{alt_gmvar} constructs a GMVAR model based on results from an arbitrary estimation round of \code{fitGMVAR}.
#'
#' @inheritParams simulateGMVAR
#' @inheritParams GMVAR
#' @param which_round based on which estimation round should the model be constructed? An integer value in 1,...,\code{ncalls}.
#' @param which_largest based on estimation round with which largest log-likelihood should the model be constructed?
#'   An integer value in 1,...,\code{ncalls}. For example, \code{which_largest=2} would take the second largest log-likelihood
#'   and construct the model based on the corresponding estimates. If used, then \code{which_round} is ignored.
#' @details It's sometimes useful to examine other estimates than the one with the highest log-likelihood. This function
#'   is wrapper around \code{GMVAR} that picks the correct estimates from an object returned by \code{fitGMVAR}.
#' @inherit GMVAR references return
#' @inherit add_data seealso
#' @examples
#' \donttest{
#' # GMVAR(1,2) model
#' fit12 <- fitGMVAR(gdpdef, p=1, M=2, ncalls=2, seeds=4:5)
#' fit12
#' fit12_2 <- alt_gmvar(fit12, which_largest=2)
#' fit12_2
#' }
#' @export

alt_gmvar <- function(gmvar, which_round=1, which_largest, calc_cond_moments=TRUE, calc_std_errors=TRUE) {
  stopifnot(!is.null(gmvar$all_estimates))
  stopifnot(which_round >= 1 && which_round <= length(gmvar$all_estimates))
  if(!missing(which_largest)) {
    stopifnot(which_largest >= 1 && which_largest <= length(gmvar$all_estimates))
    which_round <- order(gmvar$all_logliks, decreasing=TRUE)[which_largest]
  }
  ret <- GMVAR(data=gmvar$data, p=gmvar$model$p, M=gmvar$model$M, d=gmvar$model$d, params=gmvar$all_estimates[[which_round]],
               conditional=gmvar$model$conditional, parametrization=gmvar$model$parametrization,
               constraints=gmvar$model$constraints, same_means=gmvar$model$same_means, structural_pars=gmvar$model$structural_pars,
               calc_cond_moments=calc_cond_moments, calc_std_errors=calc_std_errors, stat_tol=gmvar$num_tols$stat_tol,
               posdef_tol=gmvar$num_tols$posdef_tol)

  # Pass the estimation results to the new object
  ret$all_estimates <- gmvar$all_estimates
  ret$all_logliks <- gmvar$all_logliks
  ret$which_converged <- gmvar$which_converged
  if(!is.null(gmvar$which_round)) {
   ret$which_round <- which_round
  }
  warn_eigens(ret)
  ret
}


#' @title Switch from two-regime reduced form GMVAR model to a structural GMVAR model.
#'
#' @description \code{gmvar_to_sgmvar} constructs SGMVAR model based on a reduced form GMVAR model.
#'
#' @inheritParams simulateGMVAR
#' @inheritParams GMVAR
#' @details The switch is made by simultaneously diagonalizing the two error term covariance matrices
#'   with a well known matrix decomposition (Muirhead, 1982, Theorem A9.9) and then normalizing the
#'   diagonal of the matrix W positive (which implies positive diagonal of the B-matrix). Models with
#'   more that two regimes are not supported because the matrix decomposition does not generally
#'   exists for more than two covariance matrices. If the model has only one regime (= regular SVAR model),
#'   a symmetric and pos. def. square root matrix of the error term covariance matrix is used.
#'
#'   The columns of \eqn{W} as well as the lambda parameters can be re-ordered (without changing the implied
#'   reduced form model) afterwards with the function \code{reorder_W_columns}. Also all signs in any column
#'   of \eqn{W} can be swapped (without changing the implied reduced form model) afterwards with the function
#'   \code{swap_W_signs}. These two functions work with models containing any number of regimes.
#' @return Returns an object of class \code{'gmvar'} defining a structural GMVAR model based on a
#'   two-regime reduced form GMVAR model with the main diagonal of the B-matrix normalized to be
#'   positive.
#' @seealso \code{\link{fitGMVAR}}, \code{\link{GMVAR}}, \code{\link{GIRF}}, \code{\link{reorder_W_columns}},
#'  \code{\link{swap_W_signs}}
#' @references
#'  \itemize{
#'    \item Muirhead R.J. 1982. Aspects of Multivariate Statistical Theory, \emph{Wiley}.
#'    \item Kalliovirta L., Meitz M. and Saikkonen P. 2016. Gaussian mixture vector autoregression.
#'          \emph{Journal of Econometrics}, \strong{192}, 485-498.
#'    \item Virolainen S. 2020. Structural Gaussian mixture vector autoregressive model. Unpublished working
#'      paper, available as arXiv:2007.04713.
#'  }
#' @examples
#' \donttest{
#' # Reduced form GMVAR(1,2) model
#' params12 <- c(0.55, 0.112, 0.344, 0.055, -0.009, 0.718, 0.319,
#'  0.005, 0.03, 0.619, 0.173, 0.255, 0.017, -0.136, 0.858, 1.185,
#'  -0.012, 0.136, 0.674)
#' mod12 <- GMVAR(gdpdef, p=1, M=2, params=params12)
#'
#' # Form a structural model based on the reduced form model:
#' mod12s <- gmvar_to_sgmvar(mod12)
#' mod12s
#' }
#' @export

gmvar_to_sgmvar <- function(gmvar, calc_std_errors=TRUE) {
  check_gmvar(gmvar)
  if(is.null(gmvar$data)) calc_std_errors <- FALSE
  if(!is.null(gmvar$model$structural_pars)) stop("Only reduced form models are supported!")
  p <- gmvar$model$p
  M <- gmvar$model$M
  if(M > 2) stop("Only models with at most two regimes are supported!")
  d <- gmvar$model$d
  constraints <- gmvar$model$constraints
  same_means <- gmvar$model$same_means
  params <- reform_constrained_pars(p=p, M=M, d=d, params=gmvar$params, constraints=constraints,
                                    same_means=same_means)
  all_Omega <- pick_Omegas(p=p, M=M, d=d, params=params, structural_pars=NULL)
  if(M == 1) {
    W <- matrix(diag_Omegas(Omega1=all_Omega[, , 1]), nrow=d, ncol=d, byrow=FALSE)
    lambdas <- numeric(0)
  } else { # M == 2
    tmp <- diag_Omegas(Omega1=all_Omega[, , 1], Omega2=all_Omega[, , 2])
    W <- matrix(tmp[1:(d^2)], nrow=d, ncol=d, byrow=FALSE)
    lambdas <- tmp[(d^2 + 1):(d^2 + d)]
  }

  # Normalize the main diagonal of W to be positive
  for(i1 in 1:d) {
    if(W[i1, i1] < 0) W[,i1] <- -W[,i1]
  }

  # Create SGMVAR parameter vector
  g <- ifelse(is.null(same_means), M, length(same_means)) # Number of groups of regimes with the same mean parameters
  if(is.null(same_means)) {
    all_phi0 <- pick_phi0(p=p, M=M, d=d, params=params, structural_pars=NULL)
  } else {
    all_phi0 <- gmvar$params[1:(d*g)]
  }
  alphas <- pick_alphas(p=p, M=M, d=d, params=params)
  if(is.null(constraints)) {
    all_A <- as.vector(pick_allA(p=p, M=M, d=d, params=params))
  } else {
    all_A <- gmvar$params[(d*g + 1):(d*g + ncol(constraints))] # \psi
  }
  new_params <- c(all_phi0, all_A, vec(W), lambdas, alphas[-M])
  new_W <- matrix(NA, nrow=d, ncol=d)
  diag(new_W) <- rep(1, times=d)

  # Construct the SGMVAR model based on the obtained structural parameters
  GMVAR(data=gmvar$data, p=p, M=M, d=d, params=new_params, conditional=gmvar$model$conditional,
        parametrization=gmvar$model$parametrization, constraints=constraints, same_means=same_means,
        structural_pars=list(W=new_W), calc_std_errors=calc_std_errors, stat_tol=gmvar$num_tols$stat_tol,
        posdef_tol=gmvar$num_tols$posdef_tol)
}


#' @title Reorder columns of the W-matrix and lambda parameters of a structural GMVAR model.
#'
#' @description \code{reorder_W_columns} reorder columns of the W-matrix and lambda parameters
#'   of a structural GMVAR model.
#'
#' @inheritParams simulateGMVAR
#' @param perm an integer vector of length \eqn{d} specifying the new order of the columns of \eqn{W}.
#'   Also lambda parameters of each regime will be reordered accordingly.
#' @details The order of the columns of \eqn{W} can be changed without changing the implied reduced
#'   form model as long as the order of lambda parameters is also changed accordingly. Note that the
#'   constraints imposed on \eqn{W} (or the B-matrix) will also be modified accordingly.
#'
#'   This function does not support models with constraints imposed on the lambda parameters!
#'
#'   Also all signs in any column of \eqn{W} can be swapped (without changing the implied reduced form model)
#'   with the function \code{swap_W_signs} but this obviously also swaps the sign constraints in the
#'   corresponding columns of \eqn{W}.
#' @return Returns an object of class \code{'gmvar'} defining a structural GMVAR model with the modified
#'   structural parameters and constraints.
#' @seealso \code{\link{fitGMVAR}}, \code{\link{GMVAR}}, \code{\link{GIRF}}, \code{\link{gmvar_to_sgmvar}},
#'  \code{\link{swap_W_signs}}
#' @inherit in_paramspace_int references
#' @examples
#' # Structural GMVAR(2, 2), d=2 model identified with sign-constraints:
#' params22s <- c(0.36, 0.121, 0.484, 0.072, 0.223, 0.059, -0.151, 0.395,
#'   0.406, -0.005, 0.083, 0.299, 0.218, 0.02, -0.119, 0.722, 0.093, 0.032,
#'   0.044, 0.191, 0.057, 0.172, -0.46, 0.016, 3.518, 5.154, 0.58)
#' W_22 <- matrix(c(1, 1, -1, 1), nrow=2, byrow=FALSE)
#' mod22s <- GMVAR(p=2, M=2, d=2, params=params22s, structural_pars=list(W=W_22))
#' mod22s
#'
#' # The same reduced form model, reordered W and lambda in the structual model:
#' mod22s_2 <- reorder_W_columns(mod22s, perm=2:1)
#' mod22s_2
#' @export

reorder_W_columns <- function(gmvar, perm) {
  check_gmvar(gmvar)
  if(is.null(gmvar$model$structural_pars)) stop("Only structural models are supported!")
  if(!is.null(gmvar$model$structural_pars$C_lambda)) stop("Models with constraints imposed on the lambda parameters are not supported!")
  p <- gmvar$model$p
  M <- gmvar$model$M
  d <- gmvar$model$d
  stopifnot(length(perm) == d && all(perm %in% 1:d) && length(unique(perm)) == d)
  constraints <- gmvar$model$constraints
  same_means <- gmvar$model$same_means
  structural_pars <- gmvar$model$structural_pars
  params <- reform_constrained_pars(p=p, M=M, d=d, params=gmvar$params, constraints=constraints,
                                    same_means=same_means, structural_pars=structural_pars)
  W <- pick_W(p=p, M=M, d=d, params=params, structural_pars=structural_pars)
  lambdas <- pick_lambdas(p=p, M=M, d=d, params=params, structural_pars=structural_pars)

  # Create the new parameter vector
  W <- W[, perm]
  W <- Wvec(W) # Zeros removed
  if(M > 1) {
    lambdas <- matrix(lambdas, nrow=d, ncol=M - 1, byrow=FALSE)
    lambdas <- vec(lambdas[perm,])
  }
  new_params <- gmvar$params
  new_params[(length(new_params) - (M - 1 + length(W) + length(lambdas)) + 1):(length(new_params) - (M - 1))] <- c(W, lambdas)
  new_W <- structural_pars$W[, perm]

  # Construct the SGMVAR model based on the obtained structural parameters
  calc_std_errors <- ifelse(all(is.na(gmvar$std_errors)) || is.null(gmvar$data), FALSE, TRUE)
  GMVAR(data=gmvar$data, p=p, M=M, d=d, params=new_params, conditional=gmvar$model$conditional,
        parametrization=gmvar$model$parametrization, constraints=constraints, same_means=same_means,
        structural_pars=list(W=new_W), calc_std_errors=calc_std_errors, stat_tol=gmvar$num_tols$stat_tol,
        posdef_tol=gmvar$num_tols$posdef_tol)
}


#' @title Swap all signs in pointed columns a the \eqn{W} matrix of a structural GMVAR model.
#'
#' @description \code{swap_W_signs} swaps all signs in pointed columns a the \eqn{W} matrix
#'  of a structural GMVAR model. Consequently, signs in the columns of the B-matrix are also swapped
#'  accordingly.
#'
#' @inheritParams simulateGMVAR
#' @param which_to_swap a numeric vector of length at most \eqn{d} and elemnts in \eqn{1,..,d}
#'   specifying the columns of \eqn{W} whose sign should be swapped.
#' @details All signs in any column of \eqn{W} can be swapped without changing the implied reduced form model.
#'   Consequently, also the signs in the columns of the B-matrix are swapped. Note that the sign constraints
#'   imposed on \eqn{W} (or the B-matrix) are also swapped in the corresponding columns accordingly.
#'
#'   Also the order of the columns of \eqn{W} can be changed (without changing the implied reduced
#'   form model) as long as the order of lambda parameters is also changed accordingly. This can be
#'   done with the function \code{reorder_W_columns}.
#' @return Returns an object of class \code{'gmvar'} defining a structural GMVAR model with the modified
#'   structural parameters and constraints.
#' @seealso \code{\link{fitGMVAR}}, \code{\link{GMVAR}}, \code{\link{GIRF}}, \code{\link{reorder_W_columns}},
#'  \code{\link{gmvar_to_sgmvar}}
#' @inherit reorder_W_columns references
#' @examples
#' # Structural GMVAR(2, 2), d=2 model identified with sign-constraints:
#' params22s <- c(0.36, 0.121, 0.484, 0.072, 0.223, 0.059, -0.151, 0.395,
#'   0.406, -0.005, 0.083, 0.299, 0.218, 0.02, -0.119, 0.722, 0.093, 0.032,
#'   0.044, 0.191, 0.057, 0.172, -0.46, 0.016, 3.518, 5.154, 0.58)
#' W_22 <- matrix(c(1, 1, -1, 1), nrow=2, byrow=FALSE)
#' mod22s <- GMVAR(p=2, M=2, d=2, params=params22s, structural_pars=list(W=W_22))
#' mod22s
#'
#' # The same reduced form model, with signs in the second column of W swapped:
#' swap_W_signs(mod22s, which_to_swap=2)
#'
#' # The same reduced form model, with signs in both column of W swapped:
#' swap_W_signs(mod22s, which_to_swap=1:2)
#' @export

swap_W_signs <- function(gmvar, which_to_swap) {
  check_gmvar(gmvar)
  if(is.null(gmvar$model$structural_pars)) stop("Only structural models are supported!")
  p <- gmvar$model$p
  M <- gmvar$model$M
  d <- gmvar$model$d
  stopifnot(length(which_to_swap) <= d && length(which_to_swap) >= 1 && all(which_to_swap %in% 1:d)
            && length(unique(which_to_swap)) == length(which_to_swap))
  constraints <- gmvar$model$constraints
  same_means <- gmvar$model$same_means
  structural_pars <- gmvar$model$structural_pars
  params <- reform_constrained_pars(p=p, M=M, d=d, params=gmvar$params, constraints=constraints,
                                    same_means=same_means, structural_pars=structural_pars)
  W <- pick_W(p=p, M=M, d=d, params=params, structural_pars=structural_pars)

  # Create the new parameter vector
  W[, which_to_swap] <- -W[, which_to_swap]
  W <- Wvec(W) # Zeros removed
  r <- ifelse(is.null(structural_pars$C_lambda), d*(M - 1), r <- ncol(structural_pars$C_lambda))

  new_params <- gmvar$params
  new_params[(length(new_params) - (M - 1 + length(W) + r) + 1):(length(new_params) - (M - 1 + r))] <- W
  new_W <- structural_pars$W
  new_W[, which_to_swap] <- -new_W[, which_to_swap]

  # Construct the SGMVAR model based on the obtained structural parameters
  calc_std_errors <- ifelse(all(is.na(gmvar$std_errors)) || is.null(gmvar$data), FALSE, TRUE)
  GMVAR(data=gmvar$data, p=p, M=M, d=d, params=new_params, conditional=gmvar$model$conditional,
        parametrization=gmvar$model$parametrization, constraints=constraints, same_means=same_means,
        structural_pars=list(W=new_W), calc_std_errors=calc_std_errors, stat_tol=gmvar$num_tols$stat_tol,
        posdef_tol=gmvar$num_tols$posdef_tol)
}



#' @title Update the stationarity and positive definiteness numerical tolerances of an
#'   existing class 'gmvar' model.
#'
#' @description \code{update_numtols} updates the stationarity and positive definiteness
#'   numerical tolerances of an existing class 'gmvar' model.
#'
#' @inheritParams simulateGMVAR
#' @inheritParams in_paramspace_int
#' @details All signs in any column of \eqn{W} can be swapped without changing the implied reduced form model.
#'   Consequently, also the signs in the columns of the B-matrix are swapped. Note that the sign constraints
#'   imposed on \eqn{W} (or the B-matrix) are also swapped in the corresponding columns accordingly.
#'
#'   Also the order of the columns of \eqn{W} can be changed (without changing the implied reduced
#'   form model) as long as the order of lambda parameters is also changed accordingly. This can be
#'   done with the function \code{reorder_W_columns}.
#' @return Returns an object of class \code{'gmvar'} defining a structural GMVAR model with the modified
#'   structural parameters and constraints.
#' @seealso \code{\link{fitGMVAR}}, \code{\link{GMVAR}}, \code{\link{GIRF}}, \code{\link{reorder_W_columns}},
#'  \code{\link{gmvar_to_sgmvar}}
#' @inherit reorder_W_columns references
#' @examples
#' # Structural GMVAR(2, 2), d=2 model identified with sign-constraints:
#' params22s <- c(0.36, 0.121, 0.484, 0.072, 0.223, 0.059, -0.151, 0.395,
#'   0.406, -0.005, 0.083, 0.299, 0.218, 0.02, -0.119, 0.722, 0.093, 0.032,
#'   0.044, 0.191, 0.057, 0.172, -0.46, 0.016, 3.518, 5.154, 0.58)
#' W_22 <- matrix(c(1, 1, -1, 1), nrow=2, byrow=FALSE)
#' mod22s <- GMVAR(p=2, M=2, d=2, params=params22s, structural_pars=list(W=W_22))
#' mod22s
#'
#' # Update numerical tolerances:
#' mod22s <- update_numtols(mod22s, stat_tol=1e-4, posdef_tol=1e-9)
#' mod22s # The same model
#' @export

update_numtols <- function(gmvar, stat_tol=1e-3, posdef_tol=1e-8) {
  GMVAR(data=gmvar$data, p=gmvar$model$p, M=gmvar$model$M, d=gmvar$model$d, params=gmvar$params,
        conditional=gmvar$model$conditional, parametrization=gmvar$model$parametrization,
        constraints=gmvar$model$constraints, same_means=gmvar$model$same_means,
        structural_pars=gmvar$model$structural_pars,
        calc_std_errors=ifelse(is.null(gmvar$data), FALSE, TRUE), stat_tol=stat_tol,
        posdef_tol=posdef_tol)
}
