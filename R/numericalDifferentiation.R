#' @title Calculate gradient or Hessian matrix
#'
#' @description \code{calc_gradient} or \code{calc_hessian} calculates the gradient or Hessian matrix
#'   of the given function at the given point using central difference numerical approximation.
#'   \code{get_gradient} or \code{get_hessian} calculates the gradient or Hessian matrix of the
#'   log-likelihood function at the parameter estimates of class \code{'gmvar'} object. \code{get_soc}
#'   returns eigenvalues of the Hessian matrix.
#'
#' @inheritParams quantile_residuals
#' @param x a numeric vector specifying the point where the gradient or Hessian should be calculated.
#' @param fn a function that takes in argument \code{x} as the \strong{first} argument.
#' @param h difference used to approximate the derivatives.
#' @param ... other arguments passed to \code{fn}
#' @details Especially the functions \code{get_gradient} or \code{get_hessian} can be used to check whether
#'   the found estimates denote a (local) maximum point, a saddle point or something else.
#' @return Gradient functions return numerical approximation of the gradient, and Hessian functions return
#'   numerical approximation of the Hessian. \code{get_soc} returns eigenvalues of the Hessian matrix.
#' @section Warning:
#'   No argument checks!
#' @examples
#'   # Simple function
#'   foo <- function(x) x^2 + x
#'   calc_gradient(x=1, fn=foo)
#'   calc_gradient(x=-0.5, fn=foo)
#'
#'   # More complicated function
#'   foo <- function(x, a, b) a*x[1]^2 - b*x[2]^2
#'   calc_gradient(x=c(1, 2), fn=foo, a=0.3, b=0.1)
#'
#'   # These examples below use the data 'eurusd' which comes
#'   # with the package, but in a scaled form.
#'   data <- cbind(10*eurusd[,1], 100*eurusd[,2])
#'   colnames(data) <- colnames(eurusd)
#'
#'   \donttest{
#'   # GMVAR(1,2), d=2 model:
#'   params122 <- c(0.623, -0.129, 0.959, 0.089, -0.006, 1.006, 1.746,
#'     0.804, 5.804, 3.245, 7.913, 0.952, -0.037, -0.019, 0.943, 6.926,
#'     3.982, 12.135, 0.789)
#'   mod122 <- GMVAR(data, p=1, M=2, params=params122)
#'   get_gradient(mod122)
#'   get_hessian(mod122)
#'   get_soc(mod122)
#'   }
#' @export

calc_gradient <- function(x, fn, h=6e-06, ...) {
  fn <- match.fun(fn)
  n <- length(x)
  I <- diag(1, nrow=n, ncol=n)
  vapply(1:n, function(i1) (fn(x + h*I[i1,], ...) - fn(x - h*I[i1,], ...))/(2*h), numeric(1))
}


#' @rdname calc_gradient
#' @export
calc_hessian <- function(x, fn, h=6e-06, ...) {
  fn <- match.fun(fn)
  n <- length(x)
  I <- diag(1, nrow=n, ncol=n)
  Hess <- matrix(ncol=n, nrow=n)
  for(i1 in 1:n) {
    for(i2 in i1:n) {
      dr1 <- (fn(x + h*I[i1,] + h*I[i2,], ...) - fn(x - h*I[i1,] + h*I[i2,], ...))/(2*h)
      dr2 <- (fn(x + h*I[i1,] - h*I[i2,], ...) - fn(x - h*I[i1,] - h*I[i2,], ...))/(2*h)
      Hess[i1, i2] <- (dr1 - dr2)/(2*h)
      Hess[i2, i1] <- Hess[i1, i2] # Take use of symmetry
    }
  }
  Hess
}


#' @rdname calc_gradient
#' @export
get_gradient <- function(gmvar, h=6e-06) {
  check_gmvar(gmvar)
  foo <- function(x) {
    loglikelihood(data=gmvar$data, p=gmvar$model$p, M=gmvar$model$M, params=x,
                  conditional=gmvar$model$conditional, parametrization=gmvar$model$parametrization,
                  constraints=gmvar$model$constraints, minval = NA)
  }
  calc_gradient(x=gmvar$params, fn=foo, h=h)
}

#' @rdname calc_gradient
#' @export
get_hessian <- function(gmvar, h=6e-06) {
  check_gmvar(gmvar)
  foo <- function(x) {
    loglikelihood(data=gmvar$data, p=gmvar$model$p, M=gmvar$model$M, params=x,
                  conditional=gmvar$model$conditional, parametrization=gmvar$model$parametrization,
                  constraints=gmvar$model$constraints, minval = NA)
  }
  calc_hessian(x=gmvar$params, fn=foo, h=h)
}


#' @rdname calc_gradient
#' @export
get_soc <- function(gmvar, h=6e-6) {
  eigen(get_hessian(gmvar, h))$value
}
