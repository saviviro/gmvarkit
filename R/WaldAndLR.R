
#' @title Perform Wald test for a GMvAR or SGMVAR model
#'
#' @description \code{Wald_test} performs a Wald test for a GMvAR or SGMVAR model
#'
#' @inheritParams simulateGMVAR
#' @param gmvar
#' @param A a size \eqn(k x n_params) matrix with full row rank specifying part of the null hypothesis
#'   where \eqn{n_params} is the number of parameters in the (unconstrained) model.
#'   See details for more information.
#' @param c a length \eqn{k} vector specifying part of the null hypothesis. See details for more information.
#' @details Denoting the true parameter value by \eqn{\theta_{0}}, we test the null hypothesis \eqn{A\theta_{0}=c}.
#'   Under the null hypothesis, the test statistic is asymptotically \eqn{\chi^2} distributed with \eqn{k}
#'    (\code{=nrow(A)}) degrees of freedom. The parameter \eqn{\theta_{0}} is assumed to have the same form as in
#'    the model supplied in the argument \code{gmvar} and it is presented in the documentation of the argument
#'    \code{params} in the function \code{GMVAR} (see \eqn{?GMVAR}).
#'
#'    Finally, note that this function does not check whether the specified constraints are feasible (e.g. whether
#'    the implied constrained model would be stationary or have positive definite error term covariance matrices).
#' @return Returns an object of class \eqn{'wald'} containing the test statistic and the related p-value.
#' @seealso \code{\link{fitGMVAR}}, \code{\link{GMVAR}}, \code{\link{diagnostic_plot}},
#'  \code{\link{profile_logliks}}, \code{\link{quantile_residual_tests}}
#' @inherit in_paramspace_int references
#' @examples
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
#'  # FILL IN
#' @export

Wald_test <- function(gmvar) {
  NULL
}
