#' @title Create a class 'gsmvar' object defining a reduced form or structural GMVAR, StMVAR, or G-StMVAR model
#'
#' @description \code{GSMVAR} creates a class \code{'gsmvar'} object that defines
#'  a reduced form or structural GMVAR, StMVAR, or G-StMVAR model
#'
#' @inheritParams loglikelihood_int
#' @param data a matrix or class \code{'ts'} object with \code{d>1} columns. Each column is taken to represent
#'  a single times series. \code{NA} values are not supported. Ignore if defining a model without data is desired.
#' @param d number of times series in the system, i.e. \code{ncol(data)}. This can be
#'   used to define GSMVAR models without data and can be ignored if \code{data} is provided.
#' @param calc_cond_moments should conditional means and covariance matrices should be calculated?
#'   Default is \code{TRUE} if the model contains data and \code{FALSE} otherwise.
#' @param calc_std_errors should approximate standard errors be calculated?
#' @details If data is provided, then also multivariate quantile residuals (\emph{Kalliovirta and Saikkonen 2010})
#'   are computed and included in the returned object.
#'
#'   If the function fails to calculate approximative standard errors and the parameter values are near the border
#'   of the parameter space, it might help to use smaller numerical tolerance for the stationarity and positive
#'   definiteness conditions.
#' @return Returns an object of class \code{'gsmvar'} defining the specified reduced form or structural GMVAR,
#'   StMVAR, or G-StMVAR model. Can be used to work with other functions provided in \code{gmvarkit}.
#'
#'   Note that the first autocovariance/correlation matrix in \code{$uncond_moments} is for the lag zero,
#'   the second one for the lag one, etc.
#' @section About S3 methods:
#'   If data is not provided, only the \code{print} and \code{simulate} methods are available.
#'   If data is provided, then in addition to the ones listed above, \code{predict} method is also available.
#'   See \code{?simulate.gsmvar} and \code{?predict.gsmvar} for details about the usage.
#' @seealso \code{\link{fitGSMVAR}}, \code{\link{add_data}}, \code{\link{swap_parametrization}}, \code{\link{GIRF}},
#'   \code{\link{gsmvar_to_sgsmvar}}, \code{\link{stmvar_to_gstmvar}}, \code{\link{reorder_W_columns}},
#'   \code{\link{swap_W_signs}}, \code{\link{update_numtols}}
#' @references
#'  \itemize{
#'    \item Kalliovirta L., Meitz M. and Saikkonen P. 2016. Gaussian mixture vector autoregression.
#'          \emph{Journal of Econometrics}, \strong{192}, 485-498.
#'    \item Kalliovirta L. and Saikkonen P. 2010. Reliable Residuals for Multivariate Nonlinear
#'          Time Series Models. \emph{Unpublished Revision of HECER Discussion Paper No. 247}.
#'    \item Virolainen S. 2022. Structural Gaussian mixture vector autoregressive model with application to the asymmetric
#'      effects of monetary policy shocks. Unpublished working paper, available as arXiv:2007.04713.
#'    \item Virolainen S. 2022. Gaussian and Student's t mixture vector autoregressive model with application to the
#'      asymmetric effects of monetary policy shocks in the Euro area. Unpublished working
#'      paper, available as arXiv:2109.13648.
#'  }
#' @examples
#' # GMVAR(1, 2), d=2 model:
#' params12 <- c(0.55, 0.112, 0.344, 0.055, -0.009, 0.718, 0.319, 0.005,
#'   0.03, 0.619, 0.173, 0.255, 0.017, -0.136, 0.858, 1.185, -0.012,
#'   0.136, 0.674)
#' mod12 <- GSMVAR(gdpdef, p=1, M=2, params=params12)
#' mod12
#'
#' # GMVAR(1, 2), d=2 model without data
#' mod12_2 <- GSMVAR(p=1, M=2, d=2, params=params12)
#' mod12_2
#'
#' # StMVAR(1, 2), d=2 model:
#' mod12t <- GSMVAR(gdpdef, p=1, M=2, params=c(params12, 10, 20),
#'                  model="StMVAR")
#' mod12t
#'
#' # G-StMVAR(1, 1, 1), d=2 model:
#' mod12gs <- GSMVAR(gdpdef, p=1, M=c(1, 1), params=c(params12, 20),
#'                   model="G-StMVAR")
#' mod12gs
#'
#' # GMVAR(2, 2), d=2 model with mean-parametrization:
#' params22 <- c(0.869, 0.549, 0.223, 0.059, -0.151, 0.395, 0.406,
#'  -0.005, 0.083, 0.299, 0.215, 0.002, 0.03, 0.576, 1.168, 0.218,
#'  0.02, -0.119, 0.722, 0.093, 0.032, 0.044, 0.191, 1.101, -0.004,
#'  0.105, 0.58)
#' mod22 <- GSMVAR(gdpdef, p=2, M=2, params=params22, parametrization="mean")
#' mod22
#'
#' # Structural GMVAR(2, 2), d=2 model identified with sign-constraints:
#' params22s <- c(0.36, 0.121, 0.484, 0.072, 0.223, 0.059, -0.151, 0.395,
#'   0.406, -0.005, 0.083, 0.299, 0.218, 0.02, -0.119, 0.722, 0.093, 0.032,
#'   0.044, 0.191, 0.057, 0.172, -0.46, 0.016, 3.518, 5.154, 0.58)
#' W_22 <- matrix(c(1, 1, -1, 1), nrow=2, byrow=FALSE)
#' mod22s <- GSMVAR(gdpdef, p=2, M=2, params=params22s,
#'  structural_pars=list(W=W_22))
#' mod22s
#' @export

GSMVAR <- function(data, p, M, d, params, conditional=TRUE, model=c("GMVAR", "StMVAR", "G-StMVAR"), parametrization=c("intercept", "mean"),
                  constraints=NULL, same_means=NULL, weight_constraints=NULL, structural_pars=NULL, calc_cond_moments, calc_std_errors=FALSE,
                  stat_tol=1e-3, posdef_tol=1e-8, df_tol=1e-8) {
  model <- match.arg(model)
  parametrization <- match.arg(parametrization)
  if(missing(calc_cond_moments)) calc_cond_moments <- ifelse(missing(data) || is.null(data), FALSE, TRUE)
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
  check_pMd(p=p, M=M, d=d, model=model)
  check_constraints(p=p, M=M, d=d, constraints=constraints, same_means=same_means, weight_constraints=weight_constraints,
                    structural_pars=structural_pars)
  check_parameters(p=p, M=M, d=d, params=params, model=model, parametrization=parametrization, constraints=constraints,
                   same_means=same_means, weight_constraints=weight_constraints, structural_pars=structural_pars,
                   stat_tol=stat_tol, posdef_tol=posdef_tol, df_tol=df_tol)
  npars <- n_params(p=p, M=M, d=d, model=model, constraints=constraints, same_means=same_means, weight_constraints=weight_constraints,
                    structural_pars=structural_pars)

  if(is.null(data)) {
    lok_and_mw <- list(loglik=NA, mw=NA)
    IC <- data.frame(AIC=NA, HQIC=NA, BIC=NA)
    qresiduals <- NA
  } else {
    if(npars >= d*nrow(data)) warning("There are at least as many parameters in the model as there are observations in the data")
    lok_and_mw <- loglikelihood_int(data=data, p=p, M=M, params=params, model=model,
                                    conditional=conditional, parametrization=parametrization,
                                    constraints=constraints, same_means=same_means,
                                    weight_constraints=weight_constraints,
                                    structural_pars=structural_pars,
                                    to_return="loglik_and_mw",
                                    check_params=FALSE, minval=NA,
                                    stat_tol=stat_tol, posdef_tol=posdef_tol, df_tol=df_tol)
    qresiduals <- quantile_residuals_int(data=data, p=p, M=M, params=params, model=model,
                                         conditional=conditional, parametrization=parametrization,
                                         constraints=constraints, same_means=same_means,
                                         weight_constraints=weight_constraints,
                                         structural_pars=structural_pars, stat_tol=stat_tol,
                                         posdef_tol=posdef_tol, df_tol=df_tol)
    IC <- get_IC(loglik=lok_and_mw$loglik, npars=npars, obs=ifelse(conditional, nrow(data) - p, nrow(data)))
  }
  warn_df(p=p, M=M, params=params, model=model)
  if(calc_std_errors) {
    if(is.null(data)) {
      warning("Approximate standard errors can't be calculated without data")
      std_errors <- rep(NA, npars)
    } else {
      std_errors <- tryCatch(standard_errors(data=data, p=p, M=M, params=params, model=model,
                                             conditional=conditional, parametrization=parametrization,
                                             constraints=constraints, same_means=same_means,
                                             weight_constraints=weight_constraints, structural_pars=structural_pars,
                                             minval=-(10^(ceiling(log10(nrow(data))) + ncol(data) + 1) - 1),
                                             stat_tol=stat_tol, posdef_tol=posdef_tol, df_tol=df_tol),
                             error=function(e) {
                               warning("Approximate standard errors can't be calculated:")
                               warning(e)
                               std_errors=rep(NA, npars)
                             })
    }
  } else {
    std_errors <- rep(NA, npars)
  }
  if(calc_cond_moments == FALSE || is.null(data)) {
    if(calc_cond_moments) warning("Conditional moments can't be calculated without data")
    regime_cmeans <- regime_ccovs <- total_cmeans <- total_ccovs <- arch_scalars <- NA
  } else {
    get_cm <- function(to_return) loglikelihood_int(data=data, p=p, M=M, params=params, model=model,
                                                    conditional=conditional, parametrization=parametrization,
                                                    constraints=constraints, same_means=same_means,
                                                    weight_constraints=weight_constraints, structural_pars=structural_pars,
                                                    check_params=TRUE,
                                                    to_return=to_return, minval=NA,
                                                    stat_tol=stat_tol, posdef_tol=posdef_tol, df_tol=df_tol)
    regime_cmeans <- get_cm("regime_cmeans")
    regime_ccovs <- get_cm("regime_ccovs")
    total_cmeans <- get_cm("total_cmeans")
    total_ccovs <- get_cm("total_ccovs")
    arch_scalars <- get_cm("arch_scalars")
  }

  structure(list(data=data,
                 model=list(p=p,
                            M=M,
                            d=d,
                            model=model,
                            conditional=conditional,
                            parametrization=parametrization,
                            constraints=constraints,
                            same_means=same_means,
                            weight_constraints=weight_constraints,
                            structural_pars=structural_pars),
                 params=params,
                 std_errors=std_errors,
                 mixing_weights=lok_and_mw$mw,
                 regime_cmeans=regime_cmeans,
                 regime_ccovs=regime_ccovs,
                 total_cmeans=total_cmeans,
                 total_ccovs=total_ccovs,
                 arch_scalars=arch_scalars,
                 quantile_residuals=qresiduals,
                 loglik=structure(lok_and_mw$loglik,
                                  class="logLik",
                                  df=npars),
                 IC=IC,
                 uncond_moments=uncond_moments_int(p=p, M=M, d=d, params=params, model=model,
                                                   parametrization=parametrization,
                                                   constraints=constraints,
                                                   same_means=same_means,
                                                   weight_constraints=weight_constraints,
                                                   structural_pars=structural_pars),
                 all_estimates=NULL,
                 all_logliks=NULL,
                 which_converged=NULL,
                 which_round=NULL,
                 num_tols=list(stat_tol=stat_tol,
                               posdef_tol=posdef_tol,
                               df_tol=df_tol)),
            class="gsmvar")
}


#' @title Add data to an object of class 'gsmvar' defining a GMVAR, StMVAR, or G-StMVAR model
#'
#' @description \code{add_data} adds or updates data to object of class '\code{gsmvar}' that defines
#'  a GMVAR, StMVAR, or G-StMVAR model. Also calculates mixing weights and quantile residuals accordingly.
#'
#' @inheritParams loglikelihood_int
#' @inheritParams quantile_residual_tests
#' @inheritParams GSMVAR
#' @return Returns an object of class 'gsmvar' defining the specified GSMVAR, StMVAR, or G-StMVAR model with the data added to the model.
#'   If the object already contained data, the data will be updated.
#' @seealso \code{\link{fitGSMVAR}}, \code{\link{GSMVAR}}, \code{\link{iterate_more}}, \code{\link{update_numtols}}
#' @references
#'  \itemize{
#'    \item Kalliovirta L., Meitz M. and Saikkonen P. 2016. Gaussian mixture vector autoregression.
#'          \emph{Journal of Econometrics}, \strong{192}, 485-498.
#'    \item Virolainen S. 2022. Structural Gaussian mixture vector autoregressive model with application to the asymmetric
#'      effects of monetary policy shocks. Unpublished working paper, available as arXiv:2007.04713.
#'    \item Virolainen S. 2022. Gaussian and Student's t mixture vector autoregressive model with application to the
#'      asymmetric effects of monetary policy shocks in the Euro area. Unpublished working
#'      paper, available as arXiv:2109.13648.
#'  }
#' @examples
#' # GMVAR(1, 2), d=2 model:
#' params12 <- c(0.55, 0.112, 0.344, 0.055, -0.009, 0.718, 0.319, 0.005,
#'   0.03, 0.619, 0.173, 0.255, 0.017, -0.136, 0.858, 1.185, -0.012,
#'   0.136, 0.674)
#' mod12 <- GSMVAR(p=1, M=2, d=2, params=params12)
#' mod12
#'
#' mod12_2 <- add_data(gdpdef, mod12)
#' mod12_2
#'
#' # StMVAR(1, 2), d=2 model:
#' mod12t <- GSMVAR(p=1, M=2, d=2, params=c(params12, 10, 12), model="StMVAR")
#' mod12t
#' mod12t_2 <- add_data(gdpdef, mod12t)
#' mod12t_2
#'
#' # Structural GMVAR(2, 2), d=2 model identified with sign-constraints:
#' params22s <- c(0.36, 0.121, 0.484, 0.072, 0.223, 0.059, -0.151, 0.395,
#'   0.406, -0.005, 0.083, 0.299, 0.218, 0.02, -0.119, 0.722, 0.093, 0.032,
#'   0.044, 0.191, 0.057, 0.172, -0.46, 0.016, 3.518, 5.154, 0.58)
#' W_22 <- matrix(c(1, 1, -1, 1), nrow=2, byrow=FALSE)
#' mod22s <- GSMVAR(p=2, M=2, d=2, params=params22s, structural_pars=list(W=W_22))
#' mod22s
#'
#' mod22s_2 <- add_data(gdpdef, mod22s)
#' mod22s_2
#' @export

add_data <- function(data, gsmvar, calc_cond_moments=TRUE, calc_std_errors=FALSE) {
  gsmvar <- gmvar_to_gsmvar(gsmvar) # Backward compatibility
  check_gsmvar(gsmvar)
  GSMVAR(data=data, p=gsmvar$model$p, M=gsmvar$model$M, params=gsmvar$params,
         model=gsmvar$model$model, conditional=gsmvar$model$conditional,
         parametrization=gsmvar$model$parametrization, constraints=gsmvar$model$constraints,
         same_means=gsmvar$model$same_means, weight_constraints=gsmvar$model$weight_constraints,
         structural_pars=gsmvar$model$structural_pars,
         calc_cond_moments=calc_cond_moments, calc_std_errors=calc_std_errors,
         stat_tol=gsmvar$num_tols$stat_tol, posdef_tol=gsmvar$num_tols$posdef_tol,
         df_tol=gsmvar$num_tols$df_tol)
}


#' @title Swap the parametrization of a GMVAR, StMVAR, or G-StMVAR model
#'
#' @description \code{swap_parametrization} swaps the parametrization of a GMVAR, StMVAR or G-StMVAR, model
#'  to \code{"mean"} if the current parametrization is \code{"intercept"}, and vice versa.
#'
#' @inheritParams quantile_residual_tests
#' @details \code{swap_parametrization} is a convenient tool if you have estimated the model in
#'  "intercept"-parametrization, but wish to work with "mean"-parametrization in the future, or vice versa.
#'  In \code{gmvarkit}, the approximate standard errors are only available for parametrized parameters.
#' @inherit GSMVAR references return
#' @inherit add_data seealso
#' @examples
#' \donttest{
#' # GMVAR(2, 2), d=2 model with mean-parametrization:
#' params22 <- c(0.869, 0.549, 0.223, 0.059, -0.151, 0.395, 0.406,
#'  -0.005, 0.083, 0.299, 0.215, 0.002, 0.03, 0.576, 1.168, 0.218,
#'  0.02, -0.119, 0.722, 0.093, 0.032, 0.044, 0.191, 1.101, -0.004,
#'  0.105, 0.58)
#' mod22 <- GSMVAR(gdpdef, p=2, M=2, params=params22, parametrization="mean")
#' mod22 # mean parametrization
#'
#' mod22_2 <- swap_parametrization(mod22)
#' mod22_2 # intercept parametrization
#'
#' # G-StMVAR(2, 1, 1), d=2 model with mean-parametrization:
#' mod22gs <- GSMVAR(gdpdef, p=2, M=c(1, 1), params=c(params22, 10), model="G-StMVAR",
#'                   parametrization="mean")
#' mod22gs # mean parametrization
#'
#' mod22gs_2 <- swap_parametrization(mod22gs)
#' mod22gs_2 # intercept parametrization
#'
#' # Structural GMVAR(2, 2), d=2 model identified with sign-constraints:
#' params22s <- c(0.36, 0.121, 0.484, 0.072, 0.223, 0.059, -0.151, 0.395,
#'   0.406, -0.005, 0.083, 0.299, 0.218, 0.02, -0.119, 0.722, 0.093, 0.032,
#'   0.044, 0.191, 0.057, 0.172, -0.46, 0.016, 3.518, 5.154, 0.58)
#' W_22 <- matrix(c(1, 1, -1, 1), nrow=2, byrow=FALSE)
#' mod22s <- GSMVAR(p=2, M=2, d=2, params=params22s, structural_pars=list(W=W_22))
#' mod22s # intercept parametrization
#'
#' mod22s_2 <- swap_parametrization(mod22s)
#' mod22s_2 # mean parametrization
#' }
#' @export

swap_parametrization <- function(gsmvar) {
  gsmvar <- gmvar_to_gsmvar(gsmvar) # Backward compatibility
  check_gsmvar(gsmvar)
  if(!is.null(gsmvar$model$same_means)) {
    stop("Cannot change parametrization to intercept if the mean parameters are constrained")
  }
  change_to <- ifelse(gsmvar$model$parametrization == "intercept", "mean", "intercept")
  new_params <- change_parametrization(p=gsmvar$model$p, M=gsmvar$model$M, d=gsmvar$model$d, params=gsmvar$params,
                                       model=gsmvar$model$model, constraints=gsmvar$model$constraints,
                                       weight_constraints=gsmvar$model$weight_constraints,
                                       structural_pars=gsmvar$model$structural_pars, change_to=change_to)
  GSMVAR(data=gsmvar$data, p=gsmvar$model$p, M=gsmvar$model$M, d=gsmvar$model$d, params=new_params,
         model=gsmvar$model$model, conditional=gsmvar$model$conditional, parametrization=change_to,
         constraints=gsmvar$model$constraints, weight_constraints=gsmvar$model$weight_constraints,
         structural_pars=gsmvar$model$structural_pars, calc_std_errors=ifelse(is.null(gsmvar$data), FALSE, TRUE),
         stat_tol=gsmvar$num_tols$stat_tol, posdef_tol=gsmvar$num_tols$posdef_tol, df_tol=gsmvar$num_tols$df_tol)
}


#' @title Construct a GMVAR, StMVAR, or G-StMVAR model based on results from an arbitrary estimation round of \code{fitGSMVAR}
#'
#' @description \code{alt_gsmvar} constructs a GMVAR, StMVAR, or G-StMVAR model based on results from
#'   an arbitrary estimation round of \code{fitGSMVAR}.
#'
#' @inheritParams quantile_residual_tests
#' @inheritParams GSMVAR
#' @param which_round based on which estimation round should the model be constructed? An integer value in 1,...,\code{ncalls}.
#' @param which_largest based on estimation round with which largest log-likelihood should the model be constructed?
#'   An integer value in 1,...,\code{ncalls}. For example, \code{which_largest=2} would take the second largest log-likelihood
#'   and construct the model based on the corresponding estimates. If used, then \code{which_round} is ignored.
#' @details It's sometimes useful to examine other estimates than the one with the highest log-likelihood. This function
#'   is wrapper around \code{GSMVAR} that picks the correct estimates from an object returned by \code{fitGSMVAR}.
#' @inherit GSMVAR references return
#' @inherit add_data seealso
#' @examples
#' \donttest{
#' # GMVAR(1,2) model
#' fit12 <- fitGSMVAR(gdpdef, p=1, M=2, ncalls=2, seeds=4:5)
#' fit12
#' fit12_2 <- alt_gsmvar(fit12, which_largest=2)
#' fit12_2
#' }
#' @export

alt_gsmvar <- function(gsmvar, which_round=1, which_largest, calc_cond_moments=TRUE, calc_std_errors=TRUE) {
  gsmvar <- gmvar_to_gsmvar(gsmvar) # Backward compatibility
  stopifnot(!is.null(gsmvar$all_estimates))
  stopifnot(which_round >= 1 && which_round <= length(gsmvar$all_estimates))
  if(!missing(which_largest)) {
    stopifnot(which_largest >= 1 && which_largest <= length(gsmvar$all_estimates))
    which_round <- order(gsmvar$all_logliks, decreasing=TRUE)[which_largest]
  }
  ret <- GSMVAR(data=gsmvar$data, p=gsmvar$model$p, M=gsmvar$model$M, d=gsmvar$model$d,
                params=gsmvar$all_estimates[[which_round]], model=gsmvar$model$model,
                conditional=gsmvar$model$conditional, parametrization=gsmvar$model$parametrization,
                constraints=gsmvar$model$constraints, same_means=gsmvar$model$same_means,
                weight_constraints=gsmvar$model$weight_constraints, structural_pars=gsmvar$model$structural_pars,
                calc_cond_moments=calc_cond_moments, calc_std_errors=calc_std_errors,
                stat_tol=gsmvar$num_tols$stat_tol, posdef_tol=gsmvar$num_tols$posdef_tol,
                df_tol=gsmvar$num_tols$df_tol)

  # Pass the estimation results to the new object
  ret$all_estimates <- gsmvar$all_estimates
  ret$all_logliks <- gsmvar$all_logliks
  ret$which_converged <- gsmvar$which_converged
  if(!is.null(gsmvar$which_round)) {
   ret$which_round <- which_round
  }
  warn_eigens(ret)
  ret
}


#' @title Switch from two-regime reduced form GMVAR, StMVAR, or G-StMVAR model to a structural model.
#'
#' @description \code{gsmvar_to_sgsmvar} constructs SGMVAR, SStMVAR, or SG-StMVAR model based on a reduced
#'   form GMVAR, StMVAR, or G-StMVAR model.
#'
#' @inheritParams quantile_residual_tests
#' @inheritParams GSMVAR
#' @param cholesky if \code{M == 1}, should the lower triangular Cholesky identification be employed?
#'    See details for using Cholesky identification with \code{M > 1}.
#' @details The switch is made by simultaneously diagonalizing the two error term covariance matrices
#'   with a well known matrix decomposition (Muirhead, 1982, Theorem A9.9) and then normalizing the
#'   diagonal of the matrix W positive (which implies positive diagonal of the B-matrix). Models with
#'   more that two regimes are not supported because the matrix decomposition does not generally
#'   exists for more than two covariance matrices. If the model has only one regime (= regular SVAR model),
#'   a symmetric and pos. def. square root matrix of the error term covariance matrix is used \strong{unless}
#'   \code{cholesky = TRUE} is set in the arguments, in which case Cholesky identification is employed.
#'
#'   In order to employ a structural model with Cholesky identification and multiple regimes (\code{M > 1}),
#'   use the function \code{GIRF} directly with a reduced form model (see \code{?GIRF}).
#'
#'   The columns of \eqn{W} as well as the lambda parameters can be re-ordered (without changing the implied
#'   reduced form model) afterwards with the function \code{reorder_W_columns}. Also all signs in any column
#'   of \eqn{W} can be swapped (without changing the implied reduced form model) afterwards with the function
#'   \code{swap_W_signs}. These two functions work with models containing any number of regimes.
#' @return Returns an object of class \code{'gsmvar'} defining a structural GMVAR, StMVAR, or G-StMVAR model based on a
#'   two-regime reduced form GMVAR, StMVAR, or G-StMVAR model, with the main diagonal of the B-matrix normalized to be
#'   positive.
#' @seealso \code{\link{fitGSMVAR}}, \code{\link{GSMVAR}}, \code{\link{GIRF}}, \code{\link{reorder_W_columns}},
#'  \code{\link{swap_W_signs}}, \code{\link{stmvar_to_gstmvar}}
#' @references
#'  \itemize{
#'    \item Muirhead R.J. 1982. Aspects of Multivariate Statistical Theory, \emph{Wiley}.
#'    \item Kalliovirta L., Meitz M. and Saikkonen P. 2016. Gaussian mixture vector autoregression.
#'          \emph{Journal of Econometrics}, \strong{192}, 485-498.
#'    \item Virolainen S. 2022. Structural Gaussian mixture vector autoregressive model with application to the asymmetric
#'      effects of monetary policy shocks. Unpublished working paper, available as arXiv:2007.04713.
#'    \item Virolainen S. 2022. Gaussian and Student's t mixture vector autoregressive model with application to the
#'      asymmetric effects of monetary policy shocks in the Euro area. Unpublished working
#'      paper, available as arXiv:2109.13648.
#'  }
#' @examples
#' \donttest{
#' # Reduced form GMVAR(1,2) model
#' params12 <- c(0.55, 0.112, 0.344, 0.055, -0.009, 0.718, 0.319,
#'  0.005, 0.03, 0.619, 0.173, 0.255, 0.017, -0.136, 0.858, 1.185,
#'  -0.012, 0.136, 0.674)
#' mod12 <- GSMVAR(gdpdef, p=1, M=2, params=params12)
#'
#' # Form a structural model based on the reduced form model:
#' mod12s <- gsmvar_to_sgsmvar(mod12)
#' mod12s
#'
#' #' # Reduced form StMVAR(1,2) model
#' mod12t <- GSMVAR(gdpdef, p=1, M=2, params=c(params12, 11, 12), model="StMVAR")
#'
#' # Form a structural model based on the reduced form model:
#' mod12ts <- gsmvar_to_sgsmvar(mod12t)
#' mod12ts
#' }
#' @export

gsmvar_to_sgsmvar <- function(gsmvar, calc_std_errors=TRUE, cholesky=FALSE) {
  gsmvar <- gmvar_to_gsmvar(gsmvar) # Backward compatibility
  check_gsmvar(gsmvar)
  if(is.null(gsmvar$data)) calc_std_errors <- FALSE
  if(!is.null(gsmvar$model$structural_pars)) stop("Only reduced form models are supported!")
  p <- gsmvar$model$p
  M <- gsmvar$model$M
  if(sum(M) > 2) stop("Only models with at most two regimes are supported!")
  d <- gsmvar$model$d
  model <- gsmvar$model$model
  constraints <- gsmvar$model$constraints
  same_means <- gsmvar$model$same_means
  weight_constraints <- gsmvar$model$weight_constraints
  params <- reform_constrained_pars(p=p, M=M, d=d, params=gsmvar$params, model=model,
                                    constraints=constraints, same_means=same_means,
                                    weight_constraints=weight_constraints, structural_pars=NULL)
  all_Omega <- pick_Omegas(p=p, M=M, d=d, params=params, structural_pars=NULL)
  if(sum(M) == 1) {
    lambdas <- numeric(0)
    if(cholesky) {
      W <- t(chol(all_Omega[, , 1]))
    } else {
      W <- matrix(diag_Omegas(Omega1=all_Omega[, , 1]), nrow=d, ncol=d, byrow=FALSE)
    }
  } else { # M == 2
    tmp <- diag_Omegas(Omega1=all_Omega[, , 1], Omega2=all_Omega[, , 2])
    W <- matrix(tmp[1:(d^2)], nrow=d, ncol=d, byrow=FALSE)
    lambdas <- tmp[(d^2 + 1):(d^2 + d)]
    if(cholesky) {
      warning("The argument 'cholesky=TRUE' is ignored for models with more than one mixture component")
    }
  }

  # Normalize the main diagonal of W to be positive
  for(i1 in 1:d) {
    if(W[i1, i1] < 0) W[,i1] <- -W[,i1]
  }

  # Create SGSMVAR model parameter vector
  g <- ifelse(is.null(same_means), sum(M), length(same_means)) # Number of groups of regimes with the same mean parameters
  if(is.null(same_means)) {
    all_phi0 <- pick_phi0(p=p, M=M, d=d, params=params, structural_pars=NULL)
  } else {
    all_phi0 <- gsmvar$params[1:(d*g)]
  }
  all_df <- pick_df(M=M, params=params, model=model)
  if(is.null(constraints)) {
    all_A <- as.vector(pick_allA(p=p, M=M, d=d, params=params))
  } else {
    all_A <- gsmvar$params[(d*g + 1):(d*g + ncol(constraints))] # \psi
  }
  if(is.null(weight_constraints)) {
    alphas <- pick_alphas(p=p, M=M, d=d, params=params, model=model)
    new_params <- c(all_phi0, all_A, Wvec(W), lambdas, alphas[-sum(M)], all_df)
  } else { # No alpha params in the parameter vector
    new_params <- c(all_phi0, all_A, Wvec(W), lambdas, all_df)
  }
  new_W <- matrix(NA, nrow=d, ncol=d)
  diag(new_W) <- rep(1, times=d)
  if(sum(M) == 1 && cholesky) {
    new_W[W == 0] <- 0 # Relevant for one regime Cholesky identification only
  }

  # Construct the SGSMVAR model based on the obtained structural parameters
  GSMVAR(data=gsmvar$data, p=p, M=M, d=d, params=new_params, model=model, conditional=gsmvar$model$conditional,
         parametrization=gsmvar$model$parametrization, constraints=constraints, same_means=same_means,
         weight_constraints=gsmvar$model$weight_constraints, structural_pars=list(W=new_W), calc_std_errors=calc_std_errors,
         stat_tol=gsmvar$num_tols$stat_tol, posdef_tol=gsmvar$num_tols$posdef_tol, df_tol=gsmvar$num_tols$df_tol)
}



#' @title Estimate a G-StMVAR model based on a StMVAR model that has large degrees of freedom parameters
#'
#' @description \code{stmvar_to_gstmvar} estimates a G-StMVAR model based on a StMVAR model that has
#'  large degrees of freedom parameters.
#'
#' @inheritParams quantile_residual_tests
#' @inheritParams GSMVAR
#' @inheritParams fitGSMVAR
#' @inheritParams stmvarpars_to_gstmvar
#' @param estimate set \code{TRUE} if the new model should be estimated with a variable metric algorithm
#'  using the StMAR model parameter value as the initial value. By default \code{TRUE} iff the model
#'  contains data.
#' @param calc_std_errors set \code{TRUE} if the approximate standard errors should be calculated.
#' @param maxit the maximum number of iterations for the variable metric algorithm. Ignored if \code{estimate==FALSE}.
#' @details If a StMVAR model contains large estimates for the degrees of freedom parameters,
#'   one should consider switching to the corresponding G-StMAR model that lets the corresponding
#'   regimes to be GMVAR type. \code{stmvar_to_gstmvar} does this switch conveniently. Also G-StMVAR models
#'   are supported if some of the StMVAR type regimes have large degrees of freedom paraters.
#'
#'   Note that if the model imposes constraints on the autoregressive parameters, or if a structural model imposes
#'   constraints on the lambda parameters, and the ordering the regimes changes, the constraints are removed from
#'   the model. This is because of the form of the constraints that does not generally allow to switch the ordering
#'   of the regimes. If you wish to keep the constraints, you may construct the resulting G-StMVAR model parameter
#'   vector by hand, redefine your constraints accordingly, build the model with the function \code{GSMVAR}, and then
#'   estimate it with the function \code{iterate_more}. Alternatively, you can always directly estimate the constrained
#'   G-StMVAR model with the function \code{fitGSMVAR}.
#' @return Returns an object of class \code{'gsmvar'} defining a G-StMVAR model based on the provided StMVAR (or G-StMVAR)
#'   model with the regimes that had large degrees of freedom parameters changed to GMVAR type.
#' @seealso \code{\link{fitGSMVAR}}, \code{\link{GSMVAR}}, \code{\link{GIRF}}, \code{\link{reorder_W_columns}},
#'  \code{\link{swap_W_signs}}, \code{\link{gsmvar_to_sgsmvar}}
#' @references
#'  \itemize{
#'    \item Muirhead R.J. 1982. Aspects of Multivariate Statistical Theory, \emph{Wiley}.
#'    \item Kalliovirta L., Meitz M. and Saikkonen P. 2016. Gaussian mixture vector autoregression.
#'          \emph{Journal of Econometrics}, \strong{192}, 485-498.
#'    \item Virolainen S. 2022. Structural Gaussian mixture vector autoregressive model with application to the asymmetric
#'      effects of monetary policy shocks. Unpublished working paper, available as arXiv:2007.04713.
#'    \item Virolainen S. 2022. Gaussian and Student's t mixture vector autoregressive model with application to the
#'      asymmetric effects of monetary policy shocks in the Euro area. Unpublished working
#'      paper, available as arXiv:2109.13648.
#'  }
#' @examples
#' \donttest{
#' # StMVAR(1, 2), d=2 model:
#' params12t <- c(0.5453, 0.1157, 0.331, 0.0537, -0.0422, 0.7089, 0.4181, 0.0018,
#'   0.0413, 1.6004, 0.4843, 0.1256, -0.0311, -0.6139, 0.7221, 1.2123, -0.0357,
#'   0.1381, 0.8337, 7.5564, 90000)
#' mod12t <- GSMVAR(gdpdef, p=1, M=2, params=params12t, model="StMVAR")
#' mod12t
#'
#' # Switch to the G-StMVAR model:
#' mod12gs <- stmvar_to_gstmvar(mod12t)
#' mod12gs
#' }
#' @export

stmvar_to_gstmvar <- function(gsmvar, estimate, calc_std_errors=estimate, maxdf=100, maxit=100) {
  check_gsmvar(gsmvar)
  if(missing(estimate)) estimate <- ifelse(is.null(gsmvar$data), FALSE, TRUE)
  stopifnot(all_pos_ints(c(maxdf, maxit)))
  M <- gsmvar$model$M

  # G-StMVAR model parameter vector with the large df regimes changed to GMVAR type
  new_params <- stmvarpars_to_gstmvar(p=gsmvar$model$p, M=M, d=gsmvar$model$d,
                                      params=gsmvar$params, model=gsmvar$model$model,
                                      constraints=gsmvar$model$constraints, same_means=gsmvar$model$same_means,
                                      weight_constraints=gsmvar$model$weight_constraints,
                                      structural_pars=gsmvar$model$structural_pars, maxdf=maxdf)

  # New constraints, if they we removed
  new_weight_constraints <- new_params$weight_constraints
  if(!all(new_params$reg_order == 1:sum(M))) {
    new_constraints <- NULL
    if(!is.null(gsmvar$model$structural_pars)) {
      new_structural_pars <- list(W=gsmvar$model$structural_pars$W, C_lambda=NULL, fixed_lambdas=new_params$fixed_lambdas)
    }
  } else {
    new_constraints <- gsmvar$constraints
    new_structural_pars <- gsmvar$model$structural_pars
  }
  if(is.null(gsmvar$model$structural_pars)) new_structural_pars <- NULL
  new_same_means <- new_params$same_means

 # print(new_structural_pars)

  #check_constraints(p=gsmvar$model$p, M=gsmvar$model$M, d=gsmvar$model$d, structural_pars = new_structural_pars)
  #print("lol")
  # The type and order M of the new model
  if(new_params$M[2] == 0) {
    new_model <- "GMVAR"
    new_M <- new_params$M[1]
  } else if(new_params$M[1] == 0) {
    new_model <- "StMVAR"
    new_M <- new_params$M[2]
  } else {
    new_model <- "G-StMVAR"
    new_M <- new_params$M
  }
  new_params <- new_params$params

  if(estimate) { # Build and estimate the new model
    tmp_mod <- suppressWarnings(GSMVAR(data=gsmvar$data, p=gsmvar$model$p, M=new_M, d=gsmvar$model$d, params=new_params,
                                       model=new_model, conditional=gsmvar$model$conditional,
                                       parametrization=gsmvar$model$parametrization, constraints=new_constraints,
                                       same_means=new_same_means, weight_constraints=new_weight_constraints,
                                       structural_pars=new_structural_pars, calc_std_errors=FALSE, calc_cond_moments=FALSE,
                                       stat_tol=gsmvar$num_tols$stat_tol, posdef_tol=gsmvar$num_tols$posdef_tol,
                                       df_tol=gsmvar$num_tols$df_tol))
    new_mod <- iterate_more(tmp_mod, maxit=maxit, calc_std_errors=calc_std_errors, stat_tol=gsmvar$num_tols$posdef_tol,
                            posdef_tol=gsmvar$num_tols$posdef_tol, df_tol=gsmvar$num_tols$df_tol)
  } else { # Build the new model without estimation
    new_mod <- GSMVAR(data=gsmvar$data, p=gsmvar$model$p, M=new_M, d=gsmvar$model$d, params=new_params, model=new_model,
                      conditional=gsmvar$model$conditional, parametrization=gsmvar$model$parametrization,
                      constraints=new_constraints, same_means=new_same_means,
                      weight_constraints=new_weight_constraints, structural_pars=new_structural_pars,
                      calc_std_errors=calc_std_errors, stat_tol=gsmvar$num_tols$stat_tol,
                      posdef_tol=gsmvar$num_tols$posdef_tol, df_tol=gsmvar$num_tols$df_tol)
  }
  new_mod
}


#' @title Reorder columns of the W-matrix and lambda parameters of a structural GMVAR, StMVAR, or G-StMVAR model.
#'
#' @description \code{reorder_W_columns} reorder columns of the W-matrix and lambda parameters
#'   of a structural GMVAR, StMVAR, or G-StMVAR model.
#'
#' @inheritParams quantile_residual_tests
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
#' @return Returns an object of class \code{'gsmvar'} defining a structural GMVAR, StMVAR, or G-StMVAR model with the modified
#'   structural parameters and constraints.
#' @seealso \code{\link{fitGSMVAR}}, \code{\link{GSMVAR}}, \code{\link{GIRF}}, \code{\link{gsmvar_to_sgsmvar}},
#'  \code{\link{stmvar_to_gstmvar}}, \code{\link{swap_W_signs}}
#' @inherit in_paramspace_int references
#' @examples
#' # Structural GMVAR(2, 2), d=2 model identified with sign-constraints:
#' params22s <- c(0.36, 0.121, 0.484, 0.072, 0.223, 0.059, -0.151, 0.395,
#'   0.406, -0.005, 0.083, 0.299, 0.218, 0.02, -0.119, 0.722, 0.093, 0.032,
#'   0.044, 0.191, 0.057, 0.172, -0.46, 0.016, 3.518, 5.154, 0.58)
#' W_22 <- matrix(c(1, 1, -1, 1), nrow=2, byrow=FALSE)
#' mod22s <- GSMVAR(p=2, M=2, d=2, params=params22s, structural_pars=list(W=W_22))
#' mod22s
#'
#' # The same reduced form model, reordered W and lambda in the structual model:
#' mod22s_2 <- reorder_W_columns(mod22s, perm=2:1)
#' mod22s_2
#'
#' # Structural StMVAR(2, 2), d=2 model identified with sign-constraints:
#' mod22ts <- GSMVAR(p=2, M=2, d=2, params=c(params22s, 10, 20), model="StMVAR",
#'                  structural_pars=list(W=W_22))
#' mod22ts
#'
#' # The same reduced form model, reordered W and lambda in the structual model:
#' mod22ts_2 <- reorder_W_columns(mod22ts, perm=2:1)
#' mod22ts_2
#' @export

reorder_W_columns <- function(gsmvar, perm) {
  gsmvar <- gmvar_to_gsmvar(gsmvar) # Backward compatibility
  check_gsmvar(gsmvar)
  if(is.null(gsmvar$model$structural_pars)) stop("Only structural models are supported!")
  if(!is.null(gsmvar$model$structural_pars$C_lambda)) stop("Models with C_lambda constraints are not supported!")
  p <- gsmvar$model$p
  M <- gsmvar$model$M
  d <- gsmvar$model$d
  model <- gsmvar$model$model
  stopifnot(length(perm) == d && all(perm %in% 1:d) && length(unique(perm)) == d)
  constraints <- gsmvar$model$constraints
  same_means <- gsmvar$model$same_means
  structural_pars <- gsmvar$model$structural_pars
  params <- reform_constrained_pars(p=p, M=M, d=d, params=gsmvar$params, model=model,
                                    constraints=constraints, same_means=same_means,
                                    weight_constraints=gsmvar$model$weight_constraints,
                                    structural_pars=structural_pars)
  unconstr_structural_pars <- get_unconstrained_structural_pars(structural_pars=structural_pars)
  W <- pick_W(p=p, M=M, d=d, params=params, structural_pars=unconstr_structural_pars)
  lambdas <- pick_lambdas(p=p, M=M, d=d, params=params, structural_pars=unconstr_structural_pars)
  all_df <- pick_df(M=M, params=params, model=model)

  n_alphas <- ifelse(is.null(gsmvar$model$weight_constraints), sum(M) - 1, 0)

  # Create the new parameter vector
  W <- W[, perm]
  W <- Wvec(W) # Zeros removed
  if(sum(M) > 1) {
    lambdas <- matrix(lambdas, nrow=d, ncol=sum(M) - 1, byrow=FALSE)
    lambdas <- vec(lambdas[perm,])
  }
  new_params <- gsmvar$params
  if(is.null(structural_pars$fixed_lambdas) || sum(M) == 1) { # Add W and lambdas (if M > 2)
    new_params[(length(new_params) - (n_alphas + length(W) + length(lambdas)
                                      + length(all_df)) + 1):(length(new_params) - (n_alphas + length(all_df)))] <- c(W, lambdas)
    new_fixed_lambdas <- NULL
  } else { # Fixed lambdas so only add W
    new_params[(length(new_params) - (n_alphas + length(W) + length(all_df)) + 1):(length(new_params) - (n_alphas + length(all_df)))] <- W
    new_fixed_lambdas <- lambdas
  }
  new_W <- structural_pars$W[, perm]

  # Construct the SGSMVAR model based on the obtained structural parameters
  calc_std_errors <- ifelse(all(is.na(gsmvar$std_errors)) || is.null(gsmvar$data), FALSE, TRUE)
  GSMVAR(data=gsmvar$data, p=p, M=M, d=d, params=new_params, model=model, conditional=gsmvar$model$conditional,
         parametrization=gsmvar$model$parametrization, constraints=constraints, weight_constraints=gsmvar$model$weight_constraints,
         same_means=same_means, structural_pars=list(W=new_W, fixed_lambdas=new_fixed_lambdas), calc_std_errors=calc_std_errors,
         stat_tol=gsmvar$num_tols$stat_tol, posdef_tol=gsmvar$num_tols$posdef_tol, df_tol=gsmvar$num_tols$df_tol)
}


#' @title Swap all signs in pointed columns a the \eqn{W} matrix of a structural GMVAR, StMVAR, or G-StMVAR model.
#'
#' @description \code{swap_W_signs} swaps all signs in pointed columns a the \eqn{W} matrix
#'  of a structural GMVAR, StMVAR, or G-StMVAR model. Consequently, signs in the columns of the B-matrix are also swapped
#'  accordingly.
#'
#' @inheritParams quantile_residual_tests
#' @param which_to_swap a numeric vector of length at most \eqn{d} and elemnts in \eqn{1,..,d}
#'   specifying the columns of \eqn{W} whose sign should be swapped.
#' @details All signs in any column of \eqn{W} can be swapped without changing the implied reduced form model.
#'   Consequently, also the signs in the columns of the B-matrix are swapped. Note that the sign constraints
#'   imposed on \eqn{W} (or the B-matrix) are also swapped in the corresponding columns accordingly.
#'
#'   Also the order of the columns of \eqn{W} can be changed (without changing the implied reduced
#'   form model) as long as the order of lambda parameters is also changed accordingly. This can be
#'   done with the function \code{reorder_W_columns}.
#' @return Returns an object of class \code{'gsmvar'} defining a structural GMVAR, StMVAR, or G-StMVAR model with the modified
#'   structural parameters and constraints.
#' @seealso \code{\link{fitGSMVAR}}, \code{\link{GSMVAR}}, \code{\link{GIRF}}, \code{\link{reorder_W_columns}},
#'  \code{\link{gsmvar_to_sgsmvar}}, \code{\link{stmvar_to_gstmvar}}
#' @inherit reorder_W_columns references
#' @examples
#' # Structural GMVAR(2, 2), d=2 model identified with sign-constraints:
#' params22s <- c(0.36, 0.121, 0.484, 0.072, 0.223, 0.059, -0.151, 0.395,
#'   0.406, -0.005, 0.083, 0.299, 0.218, 0.02, -0.119, 0.722, 0.093, 0.032,
#'   0.044, 0.191, 0.057, 0.172, -0.46, 0.016, 3.518, 5.154, 0.58)
#' W_22 <- matrix(c(1, 1, -1, 1), nrow=2, byrow=FALSE)
#' mod22s <- GSMVAR(p=2, M=2, d=2, params=params22s, structural_pars=list(W=W_22))
#' mod22s
#'
#' # The same reduced form model, with signs in the second column of W swapped:
#' swap_W_signs(mod22s, which_to_swap=2)
#'
#' # The same reduced form model, with signs in both column of W swapped:
#' swap_W_signs(mod22s, which_to_swap=1:2)
#'
#' #' # Structural G-StMVAR(2, 1, 1), d=2 model identified with sign-constraints:
#' mod22gss <- GSMVAR(p=2, M=c(1, 1), d=2, params=c(params22s, 10), model="G-StMVAR",
#'                    structural_pars=list(W=W_22))
#' mod22gss
#'
#' # The same reduced form model, with signs in the first column of W swapped:
#' swap_W_signs(mod22gss, which_to_swap=1)
#' @export

swap_W_signs <- function(gsmvar, which_to_swap) {
  gsmvar <- gmvar_to_gsmvar(gsmvar) # Backward compatibility
  check_gsmvar(gsmvar)
  if(is.null(gsmvar$model$structural_pars)) stop("Only structural models are supported!")
  p <- gsmvar$model$p
  M <- gsmvar$model$M
  d <- gsmvar$model$d
  model <- gsmvar$model$model
  stopifnot(length(which_to_swap) <= d && length(which_to_swap) >= 1 && all(which_to_swap %in% 1:d)
            && length(unique(which_to_swap)) == length(which_to_swap))
  constraints <- gsmvar$model$constraints
  same_means <- gsmvar$model$same_means
  structural_pars <- gsmvar$model$structural_pars
  params <- reform_constrained_pars(p=p, M=M, d=d, params=gsmvar$params, model=model,
                                    constraints=constraints, same_means=same_means,
                                    weight_constraints=gsmvar$model$weight_constraints,
                                    structural_pars=structural_pars)
  unconstr_structural_pars <- get_unconstrained_structural_pars(structural_pars=structural_pars)
  W <- pick_W(p=p, M=M, d=d, params=params, structural_pars=unconstr_structural_pars)
  all_df <- pick_df(M=M, params=params, model=model)
  n_alphas <- ifelse(is.null(gsmvar$model$weight_constraints), sum(M) - 1, 0)

  # Create the new parameter vector
  W[, which_to_swap] <- -W[, which_to_swap]
  W <- Wvec(W) # Zeros removed
  if(is.null(structural_pars$C_lambda) && is.null(structural_pars$fixed_lambdas)) {
    r <- d*(sum(M) - 1) # The number of lambda parameters
  } else if(!is.null(structural_pars$C_lambda)) {
    r <- ncol(structural_pars$C_lambda)
  } else { # !is.null(structural_pars$fixed_lambdas)
    r <- 0
  }

  new_params <- gsmvar$params
  new_params[(length(new_params) - (n_alphas + length(W) + r
                                    + length(all_df)) + 1):(length(new_params) - (n_alphas + r + length(all_df)))] <- W
  new_W <- structural_pars$W
  new_W[, which_to_swap] <- -new_W[, which_to_swap]

  # Construct the SGSMVAR model based on the obtained structural parameters
  calc_std_errors <- ifelse(all(is.na(gsmvar$std_errors)) || is.null(gsmvar$data), FALSE, TRUE)
  GSMVAR(data=gsmvar$data, p=p, M=M, d=d, params=new_params, model=model, conditional=gsmvar$model$conditional,
         parametrization=gsmvar$model$parametrization, constraints=constraints, same_means=same_means,
         weight_constraints=gsmvar$model$weight_constraints,
         structural_pars=list(W=new_W, C_lambda=structural_pars$C_lambda, fixed_lambdas=structural_pars$fixed_lambdas),
         calc_std_errors=calc_std_errors, stat_tol=gsmvar$num_tols$stat_tol,
         posdef_tol=gsmvar$num_tols$posdef_tol, df_tol=gsmvar$num_tols$df_tol)
}



#' @title Update the stationarity and positive definiteness numerical tolerances of an
#'   existing class 'gsmvar' model.
#'
#' @description \code{update_numtols} updates the stationarity and positive definiteness
#'   numerical tolerances of an existing class 'gsmvar' model.
#'
#' @inheritParams quantile_residual_tests
#' @inheritParams in_paramspace_int
#' @details All signs in any column of \eqn{W} can be swapped without changing the implied reduced form model.
#'   Consequently, also the signs in the columns of the B-matrix are swapped. Note that the sign constraints
#'   imposed on \eqn{W} (or the B-matrix) are also swapped in the corresponding columns accordingly.
#'
#'   Also the order of the columns of \eqn{W} can be changed (without changing the implied reduced
#'   form model) as long as the order of lambda parameters is also changed accordingly. This can be
#'   done with the function \code{reorder_W_columns}.
#' @return Returns an object of class \code{'gsmvar'} defining a structural GSMVAR model with the modified
#'   structural parameters and constraints.
#' @seealso \code{\link{fitGSMVAR}}, \code{\link{GSMVAR}}, \code{\link{GIRF}}, \code{\link{reorder_W_columns}},
#'  \code{\link{gsmvar_to_sgsmvar}}, \code{\link{stmvar_to_gstmvar}}
#' @inherit reorder_W_columns references
#' @examples
#' # Structural GMVAR(2, 2), d=2 model identified with sign-constraints:
#' params22s <- c(0.36, 0.121, 0.484, 0.072, 0.223, 0.059, -0.151, 0.395,
#'   0.406, -0.005, 0.083, 0.299, 0.218, 0.02, -0.119, 0.722, 0.093, 0.032,
#'   0.044, 0.191, 0.057, 0.172, -0.46, 0.016, 3.518, 5.154, 0.58)
#' W_22 <- matrix(c(1, 1, -1, 1), nrow=2, byrow=FALSE)
#' mod22s <- GSMVAR(p=2, M=2, d=2, params=params22s, structural_pars=list(W=W_22))
#' mod22s
#'
#' # Update numerical tolerances:
#' mod22s <- update_numtols(mod22s, stat_tol=1e-4, posdef_tol=1e-9, df_tol=1e-10)
#' mod22s # The same model
#' @export

update_numtols <- function(gsmvar, stat_tol=1e-3, posdef_tol=1e-8, df_tol=1e-8) {
  gsmvar <- gmvar_to_gsmvar(gsmvar) # Backward compatibility
  GSMVAR(data=gsmvar$data, p=gsmvar$model$p, M=gsmvar$model$M, d=gsmvar$model$d,
         params=gsmvar$params, model=gsmvar$model$model, conditional=gsmvar$model$conditional,
         parametrization=gsmvar$model$parametrization, constraints=gsmvar$model$constraints,
         same_means=gsmvar$model$same_means, weight_constraints=gsmvar$model$weight_constraints,
         structural_pars=gsmvar$model$structural_pars,
         calc_std_errors=ifelse(is.null(gsmvar$data), FALSE, TRUE),
         stat_tol=stat_tol, posdef_tol=posdef_tol, df_tol=df_tol)
}
