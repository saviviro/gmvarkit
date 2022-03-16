
#' @title Perform Wald test for a GMVAR, StMVAR, or G-StMVAR model
#'
#' @description \code{Wald_test} performs a Wald test for a GMVAR, StMVAR, or G-StMVAR model
#'
#' @inheritParams quantile_residual_tests
#' @inheritParams calc_gradient
#' @param A a size \eqn{(k x n_params)} matrix with full row rank specifying part of the null hypothesis
#'   where \eqn{n_params} is the number of parameters in the (unconstrained) model.
#'   See details for more information.
#' @param c a length \eqn{k} vector specifying part of the null hypothesis. See details for more information.
#' @param custom_h a numeric vector with the same length as \code{x} specifying the difference \code{h}
#'  for each dimension separately. If \code{NULL} (default), then the difference \code{1e-6} used for
#'   all but overly large degrees of freedom parameters. For them, the difference is adjusted to avoid
#'   numerical problems.
#' @details Denoting the true parameter value by \eqn{\theta_{0}}, we test the null hypothesis \eqn{A\theta_{0}=c}.
#'   Under the null, the test statistic is asymptotically \eqn{\chi^2}-distributed with \eqn{k}
#'   (\code{=nrow(A)}) degrees of freedom. The parameter \eqn{\theta_{0}} is assumed to have the same form as in
#'   the model supplied in the argument \code{gsmvar} and it is presented in the documentation of the argument
#'   \code{params} in the function \code{GSMVAR} (see \code{?GSMVAR}).
#'
#'   Finally, note that this function does \strong{not} check whether the specified constraints are feasible (e.g. whether
#'   the implied constrained model would be stationary or have positive definite error term covariance matrices).
#' @return A list with class "hypotest" containing the test results and arguments used to calculate the test.
#' @seealso \code{\link{LR_test}}, \code{\link{fitGSMVAR}}, \code{\link{GSMVAR}}, \code{\link{diagnostic_plot}},
#'  \code{\link{profile_logliks}}, \code{\link{quantile_residual_tests}}, \code{\link{cond_moment_plot}}
#' @inherit in_paramspace_int references
#' @examples
#' \donttest{
#'  # Structural GMVAR(2, 2), d=2 model with recursive identification
#'  W22 <- matrix(c(1, NA, 0, 1), nrow=2, byrow=FALSE)
#'  fit22s <- fitGSMVAR(gdpdef, p=2, M=2, structural_pars=list(W=W22),
#'                     ncalls=1, seeds=2)
#'  fit22s
#'
#'  # Test whether the lambda parameters (of the second regime) are identical
#'  # (due to the zero constraint, the model is identified under the null):
#'  # fit22s has parameter vector of length 26 with the lambda parameters
#'  # in elements 24 and 25.
#'  A <- matrix(c(rep(0, times=23), 1, -1, 0), nrow=1, ncol=26)
#'  c <- 0
#'  Wald_test(fit22s, A=A, c=c)
#'
#'  # Test whether the off-diagonal elements of the first regime's first
#'  # AR coefficient matrix (A_11) are both zero:
#'  # fit22s has parameter vector of length 26 and the off-diagonal elements
#'  # of the 1st regime's 1st AR coefficient matrix are in the elements 6 and 7.
#'  A <- rbind(c(rep(0, times=5), 1, rep(0, times=20)),
#'             c(rep(0, times=6), 1, rep(0, times=19)))
#'  c <- c(0, 0)
#'  Wald_test(fit22s, A=A, c=c)
#' }
#' @export

Wald_test <- function(gsmvar, A, c, custom_h=NULL) {
  gsmvar <- gmvar_to_gsmvar(gsmvar) # Backward compatibility
  params <- gsmvar$params
  stopifnot(is.matrix(A) && ncol(A) == length(params) && nrow(A) <= ncol(A))
  stopifnot(length(c) == nrow(A))
  stopifnot(!is.null(gsmvar$data))
  if(qr(A)$rank != nrow(A)) stop("The constraint matrix 'A' should have full row rank")

  # Calculate Hessian matrix at the estimate
  Hess <- tryCatch(get_hessian(gsmvar, custom_h=custom_h),
                   error=function(e) {
                     print(paste("Failed to calculate Hessian matrix:", e))
                     return(NA)})

  # Invert the Hessian matrix
  inv_Hess <- tryCatch(solve(Hess), error=function(e) {
    print(paste("Failed to invert Hessian matrix:", e))
    return(NA)
  })
  if(anyNA(inv_Hess)) stop("Couldn't invert Hessian matrix of the log-likelihood function. This might happen when the mixing weights are very close to zero for some regime (if so, reduce the redundant regime from the model).")

  # Calculate the test statistic
  test_stat <- as.numeric(crossprod(A%*%params - c, solve(-tcrossprod(A%*%inv_Hess, A), A%*%params - c))) # t(A%*%params - c)%*%solve(-A%*%inv_Hess%*%t(A))%*%(A%*%params - c)

  # Calculate the p-value
  df <- nrow(A)
  p_value <- pchisq(test_stat, df=df, lower.tail=FALSE)

  # Return
  #dname <- paste0(deparse(substitute(gsmvar)),", ", deparse(substitute(A)), ", ", deparse(substitute(c)))
  #structure(list(statistic=c("W"=test_stat),
  #               parameter=c("df"=df),
  #               p.value=p_value,
  #               alternative="the true parameter theta does not satisfy A%*%theta = c",
  #               data.name=dname,
  #               method="Wald test",
  #               gsmvar=gsmvar,
  #               A=A,
  #               c=c,
  #               h=h),
  #          class="htest")
  structure(list(gsmvar=gsmvar,
                 A=A,
                 c=c,
                 df=df,
                 test_stat=test_stat,
                 df=df,
                 p_value=p_value,
                 type="Wald test"),
            class="hypotest")
}


#' @title Perform likelihood ratio test for a GMVAR, StMVAR, or G-StMVAR model
#'
#' @description \code{LR_test} performs a likelihood ratio test for a GMVAR, StMVAR, or G-StMVAR model
#'
#' @param gsmvar1 an object of class \code{'gsmvar'} generated by \code{fitGSMVAR} or \code{GSMVAR}, containing
#'   the \strong{freely estimated} model.
#' @param gsmvar2 an object of class \code{'gsmvar'} generated by \code{fitGSMVAR} or \code{GSMVAR}, containing
#'   the \strong{constrained} model.
#' @details Performs a likelihood ratio test, testing the null hypothesis that the true parameter value lies
#'   in the constrained parameter space. Under the null, the test statistic is asymptotically
#'   \eqn{\chi^2}-distributed with \eqn{k} degrees of freedom, \eqn{k} being the difference in the dimensions
#'   of the unconstrained and constrained parameter spaces.
#'
#'   Note that this function does \strong{not} verify that the two models are actually nested.
#' @seealso \code{\link{Wald_test}}, \code{\link{fitGSMVAR}}, \code{\link{GSMVAR}}, \code{\link{diagnostic_plot}},
#'  \code{\link{profile_logliks}}, \code{\link{quantile_residual_tests}}, \code{\link{cond_moment_plot}}
#' @inherit in_paramspace_int references
#' @inherit Wald_test return
#' @examples
#' \donttest{
#'  ## These are long running examples that use parallel computing!
#'  ## The below examples take around 1 minute to run.
#'
#'  # Structural GMVAR(2, 2), d=2 model with recursive identification
#'  W22 <- matrix(c(1, NA, 0, 1), nrow=2, byrow=FALSE)
#'  fit22s <- fitGSMVAR(gdpdef, p=2, M=2, structural_pars=list(W=W22),
#'                      ncalls=1, seeds=2)
#'
#'  # The same model but the AR coefficients restricted to be the same
#'  # in both regimes:
#'  C_mat <- rbind(diag(2*2^2), diag(2*2^2))
#'  fit22sc <- fitGSMVAR(gdpdef, p=2, M=2, constraints=C_mat,
#'                       structural_pars=list(W=W22), ncalls=1, seeds=1)
#'
#'  # Test the AR constraints with likelihood ratio test:
#'  LR_test(fit22s, fit22sc)
#'  }
#' @export

LR_test <- function(gsmvar1, gsmvar2) {
  gsmvar1 <- gmvar_to_gsmvar(gsmvar1) # Backward compatibility
  gsmvar2 <- gmvar_to_gsmvar(gsmvar2) # Backward compatibility

  check_gsmvar(gsmvar1, object_name="gsmvar1")
  check_gsmvar(gsmvar2, object_name="gsmvar2")
  stopifnot(length(gsmvar1$params) > length(gsmvar2$params))
  stopifnot(gsmvar1$loglik >= gsmvar2$loglik)

  test_stat <- as.numeric(2*(gsmvar1$loglik - gsmvar2$loglik))
  df <- length(gsmvar1$params) - length(gsmvar2$params)
  p_value <- pchisq(test_stat, df=df, lower.tail=FALSE)

  # Return
  #dname <- paste(deparse(substitute(gsmvar1)), "and", deparse(substitute(gsmvar2)))
  #structure(list(statistic=c("LR"=test_stat),
  #               parameter=c("df"=df),
  #               p.value=p_value,
  #               alternative=paste("the true parameter does not satisfy the constraints imposed in", deparse(substitute(gsmvar2))),
  #               data.name=dname,
  #               method="Likelihood ratio test",
  #               gsmvar1=gsmvar1,
  #               gsmvar2=gsmvar2),
  #          class="htest")
  structure(list(gsmvar1=gsmvar1,
                 gsmvar2=gsmvar2,
                 test_stat=test_stat,
                 df=df,
                 p_value=p_value,
                 type="Likelihood ratio test"),
            class="hypotest")
}


