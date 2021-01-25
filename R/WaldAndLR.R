
#' @title Perform Wald test for a GMVAR or SGMVAR model
#'
#' @description \code{Wald_test} performs a Wald test for a GMVAR or SGMVAR model
#'
#' @inheritParams simulateGMVAR
#' @inheritParams calc_gradient
#' @param A a size \eqn{(k x n_params)} matrix with full row rank specifying part of the null hypothesis
#'   where \eqn{n_params} is the number of parameters in the (unconstrained) model.
#'   See details for more information.
#' @param c a length \eqn{k} vector specifying part of the null hypothesis. See details for more information.
#' @details Denoting the true parameter value by \eqn{\theta_{0}}, we test the null hypothesis \eqn{A\theta_{0}=c}.
#'   Under the null, the test statistic is asymptotically \eqn{\chi^2}-distributed with \eqn{k}
#'   (\code{=nrow(A)}) degrees of freedom. The parameter \eqn{\theta_{0}} is assumed to have the same form as in
#'   the model supplied in the argument \code{gmvar} and it is presented in the documentation of the argument
#'   \code{params} in the function \code{GMVAR} (see \code{?GMVAR}).
#'
#'   Finally, note that this function does \strong{not} check whether the specified constraints are feasible (e.g. whether
#'   the implied constrained model would be stationary or have positive definite error term covariance matrices).
#' @return Returns an object of class \eqn{'wald'} containing the test statistic and the related p-value.
#' @seealso \code{\link{LR_test}}, \code{\link{fitGMVAR}}, \code{\link{GMVAR}}, \code{\link{diagnostic_plot}},
#'  \code{\link{profile_logliks}}, \code{\link{quantile_residual_tests}}, \code{\link{cond_moment_plot}}
#' @inherit in_paramspace_int references
#' @examples
#' \donttest{
#'  ## These are long running examples that use parallel computing!
#'  ## The below examples take around 40 seconds to run.
#'  # Load the data
#'  data(eurusd, package="gmvarkit")
#'  data <- cbind(10*eurusd[,1], 100*eurusd[,2])
#'  colnames(data) <- colnames(eurusd)
#'
#'  # Structural GMVAR(2, 2), d=2 model identified similarly to Cholesky:
#'  W22 <- matrix(c(1, NA, 0, 1), nrow=2, byrow=FALSE)
#'  fit22s <- fitGMVAR(data, p=2, M=2, structural_pars=list(W=W22),
#'                     ncalls=1, seeds=4)
#'  fit22s
#'
#'  # Test whether the lambda parameters (of the second regime) are identical:
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

Wald_test <- function(gmvar, A, c, h=6e-6) {
  params <- gmvar$params
  stopifnot(is.matrix(A) && ncol(A) == length(params) && nrow(A) <= ncol(A))
  stopifnot(length(c) == nrow(A))
  stopifnot(!is.null(gmvar$data))
  if(qr(A)$rank != nrow(A)) stop("The constraint matrix 'A' should have full row rank")

  # Calculate Hessian matrix at the estimate
  minval <- get_minval(gmvar$data)
  loglik_fn <- function(pars) {
    tryCatch(loglikelihood_int(data=gmvar$data, p=gmvar$model$p, M=gmvar$model$M, params=pars,
                               conditional=gmvar$model$conditional,
                               parametrization=gmvar$model$parametrization,
                               constraints=gmvar$model$constraints,
                               same_means=gmvar$model$same_means,
                               structural_pars=gmvar$model$structural_pars,
                               check_params=TRUE,
                               to_return="loglik", minval=minval,
                               stat_tol=gmvar$num_tols$stat_tol,
                               posdef_tol=gmvar$num_tols$posdef_tol),
             error=function(e) {
               print(paste("Failed to evualuate log-likelihood function in the approximation of Hessian matrix:", e))
               return(NA)
               })
  }
  Hess <- calc_hessian(x=params, fn=loglik_fn, h=h)
  if(anyNA(Hess)) stop("Unable to fully calculate Hessian matrix of the log-likelihood function using central difference numerical approximation. Check whether there is something funny in the estimates or maybe try another difference 'h'?")

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
  structure(list(gmvar=gmvar,
                 A=A,
                 c=c,
                 df=df,
                 test_stat=test_stat,
                 df=df,
                 p_value=p_value),
            class="wald")
}


#' @title Perform likelihood ratio test for a GMVAR or SGMVAR model
#'
#' @description \code{LR_test} performs a likelihood ratio test for a GMVAR or SGMVAR model
#'
#' @param gmvar1 an object of class \code{'gmvar'} generated by \code{fitGMVAR} or \code{GMVAR}, containing
#'   the \strong{freely estimated} model.
#' @param gmvar2 an object of class \code{'gmvar'} generated by \code{fitGMVAR} or \code{GMVAR}, containing
#'   the \strong{constrained} model.
#' @details Performs a likelihood ratio test, testing the null hypothesis that the true parameter value lies
#'   in the constrained parameter space. Under the null, the test statistic is asymptotically
#'   \eqn{\chi^2}-distributed with \eqn{k} degrees of freedom, \eqn{k} being the difference in the dimensions
#'   of the unconstrained and constrained parameter spaces.
#'
#'   Note that this function does \strong{not} verify that the two models are actually nested.
#' @return Returns an object of class \eqn{'lr'} containing the test statistic and the related p-value.
#' @seealso \code{\link{Wald_test}}, \code{\link{fitGMVAR}}, \code{\link{GMVAR}}, \code{\link{diagnostic_plot}},
#'  \code{\link{profile_logliks}}, \code{\link{quantile_residual_tests}}, \code{\link{cond_moment_plot}}
#' @inherit in_paramspace_int references
#' @examples
#' \donttest{
#'  ## These are long running examples that use parallel computing!
#'  ## The below examples take around 1 minute to run.
#'
#'  # Load the data
#'  data(eurusd, package="gmvarkit")
#'  data <- cbind(10*eurusd[,1], 100*eurusd[,2])
#'  colnames(data) <- colnames(eurusd)
#'
#'  # Structural GMVAR(2, 2), d=2 model identified similarly to Cholesky:
#'  W22 <- matrix(c(1, NA, 0, 1), nrow=2, byrow=FALSE)
#'  fit22s <- fitGMVAR(data, p=2, M=2, structural_pars=list(W=W22),
#'                     ncalls=1, seeds=4)
#'
#'  # The same model but the AR coefficients restricted to be the same
#'  # in both regimes:
#'  C_mat <- rbind(diag(2*2^2), diag(2*2^2))
#'  fit22sc <- fitGMVAR(data, p=2, M=2, constraints=C_mat,
#'                      structural_pars=list(W=W22), ncalls=1, seeds=1)
#'
#'  # Test the AR constraints with likelihood ratio test:
#'  LR_test(fit22s, fit22sc)
#'  }
#' @export

LR_test <- function(gmvar1, gmvar2) {
  check_gmvar(gmvar1, object_name="gmvar1")
  check_gmvar(gmvar2, object_name="gmvar2")
  stopifnot(length(gmvar1$params) > length(gmvar2$params))
  stopifnot(gmvar1$loglik >= gmvar2$loglik)

  test_stat <- as.numeric(2*(gmvar1$loglik - gmvar2$loglik))
  df <- length(gmvar1$params) - length(gmvar2$params)
  p_value <- pchisq(test_stat, df=df, lower.tail=FALSE)

  structure(list(gmvar1=gmvar1,
                 gmvar2=gmvar2,
                 test_stat=test_stat,
                 df=df,
                 p_value=p_value),
            class="lr")
}


