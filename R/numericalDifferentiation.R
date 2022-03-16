#' @title Calculate gradient or Hessian matrix
#'
#' @description \code{calc_gradient} or \code{calc_hessian} calculates the gradient or Hessian matrix
#'   of the given function at the given point using central difference numerical approximation.
#'   \code{get_gradient} or \code{get_hessian} calculates the gradient or Hessian matrix of the
#'   log-likelihood function at the parameter estimates of a class \code{'gsmvar'} object. \code{get_soc}
#'   returns eigenvalues of the Hessian matrix, and \code{get_foc} is the same as \code{get_gradient}
#'   but named conveniently.
#'
#' @inheritParams quantile_residuals
#' @param x a numeric vector specifying the point where the gradient or Hessian should be calculated.
#' @param fn a function that takes in argument \code{x} as the \strong{first} argument.
#' @param h difference used to approximate the derivatives.
#' @param varying_h a numeric vector with the same length as \code{x} specifying the difference \code{h}
#'  for each dimension separately. If \code{NULL} (default), then the difference given as parameter \code{h}
#'  will be used for all dimensions.
#' @param custom_h same as \code{varying_h} except thaqt if \code{NULL} (default), then the difference used
#'   for differentiating overly large degrees of freedom parameters is adjusted to avoid numerical problems,
#'   and the difference \code{6e-6} is used for the other parameters.
#' @param ... other arguments passed to \code{fn}
#' @details In particular, the functions \code{get_foc} and \code{get_soc} can be used to check whether
#'   the found estimates denote a (local) maximum point, a saddle point, or something else. Note that
#'   profile log-likelihood functions can be conveniently plotted with the function \code{profile_logliks}.
#' @return Gradient functions return numerical approximation of the gradient and Hessian functions return
#'   numerical approximation of the Hessian. \code{get_soc} returns eigenvalues of the Hessian matrix.
#' @seealso \code{\link{profile_logliks}}
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
#'   \donttest{
#'   # GMVAR(1,2), d=2 model:
#'   params12 <- c(0.55, 0.112, 0.344, 0.055, -0.009, 0.718, 0.319, 0.005,
#'    0.03, 0.619, 0.173, 0.255, 0.017, -0.136, 0.858, 1.185, -0.012,
#'    0.136, 0.674)
#'   mod12 <- GSMVAR(gdpdef, p=1, M=2, params=params12)
#'   get_gradient(mod12)
#'   get_hessian(mod12)
#'   get_soc(mod12)
#'   }
#' @export

calc_gradient <- function(x, fn, h=6e-06, varying_h=NULL, ...) {
  fn <- match.fun(fn)
  n <- length(x)
  I <- diag(1, nrow=n, ncol=n)
  if(is.null(varying_h)) { # The same difference h for all parameters
    h <- rep(h, times=n)
  } else { # Varying h
    stopifnot(length(varying_h) == length(x))
    h <- varying_h
  }
  vapply(1:n, function(i1) (fn(x + h[i1]*I[i1,], ...) - fn(x - h[i1]*I[i1,], ...))/(2*h[i1]), numeric(1))
}


#' @rdname calc_gradient
#' @export

calc_hessian <- function(x, fn, h=6e-06, varying_h=NULL, ...) {
  fn <- match.fun(fn)
  n <- length(x)
  I <- diag(1, nrow=n, ncol=n)
  if(is.null(varying_h)) { # The same difference h for all parameters
    h <- rep(h, times=n)
  } else { # Varying h
    stopifnot(length(varying_h) == length(x))
    h <- varying_h
  }
  Hess <- matrix(ncol=n, nrow=n)
  for(i1 in 1:n) {
    for(i2 in i1:n) {
      dr1 <- (fn(x + h[i1]*I[i1,] + h[i2]*I[i2,], ...) - fn(x - h[i1]*I[i1,] + h[i2]*I[i2,], ...))/(2*h[i1])
      dr2 <- (fn(x + h[i1]*I[i1,] - h[i2]*I[i2,], ...) - fn(x - h[i1]*I[i1,] - h[i2]*I[i2,], ...))/(2*h[i1])
      Hess[i1, i2] <- (dr1 - dr2)/(2*h[i2])
      Hess[i2, i1] <- Hess[i1, i2] # Take use of symmetry
    }
  }
  Hess
}


#' @rdname calc_gradient
#' @export

get_gradient <- function(gsmvar, custom_h=NULL) {
  check_gsmvar(gsmvar)
  if(is.null(custom_h)) { # Adjust h for overly large degrees of freedom parameters
    varying_h <- get_varying_h(M=gsmvar$model$M, params=gsmvar$params, model=gsmvar$model$model)
  } else { # Utilize user-specified h
    stopifnot(length(custom_h) == length(gsmvar$params))
    varying_h <- custom_h
  }
  foo <- function(x) {
    loglikelihood_int(data=gsmvar$data, p=gsmvar$model$p, M=gsmvar$model$M, params=x, model=gsmvar$model$model,
                      conditional=gsmvar$model$conditional, parametrization=gsmvar$model$parametrization,
                      structural_pars=gsmvar$model$structural_pars, constraints=gsmvar$model$constraints,
                      to_return="loglik", check_params=TRUE, minval=NA, stat_tol=gsmvar$num_tols$stat_tol,
                      posdef_tol=gsmvar$num_tols$posdef_tol, df_tol=gsmvar$num_tols$df_tol)
  }
  calc_gradient(x=gsmvar$params, fn=foo, varying_h=varying_h)
}



#' @rdname calc_gradient
#' @export

get_hessian <- function(gsmvar, custom_h=NULL) {
  check_gsmvar(gsmvar)
  if(is.null(custom_h)) { # Adjust h for overly large degrees of freedom parameters
    varying_h <- get_varying_h(M=gsmvar$model$M, params=gsmvar$params, model=gsmvar$model$model)
  } else { # Utilize user-specified h
    stopifnot(length(custom_h) == length(gsmvar$params))
    varying_h <- custom_h
  }
  foo <- function(x) {
    loglikelihood_int(data=gsmvar$data, p=gsmvar$model$p, M=gsmvar$model$M, params=x, model=gsmvar$model$model,
                      conditional=gsmvar$model$conditional, parametrization=gsmvar$model$parametrization,
                      constraints=gsmvar$model$constraints, structural_pars=gsmvar$model$structural_pars,
                      to_return="loglik", check_params=TRUE, minval=NA, stat_tol=gsmvar$num_tols$stat_tol,
                      posdef_tol=gsmvar$num_tols$posdef_tol, df_tol=gsmvar$num_tols$df_tol)
  }
  calc_hessian(x=gsmvar$params, fn=foo, varying_h=varying_h)
}

#' @rdname calc_gradient
#' @export
get_foc <- function(gsmvar, custom_h=NULL) {
  get_gradient(gsmvar, custom_h=custom_h)
}

#' @rdname calc_gradient
#' @export
get_soc <- function(gsmvar, custom_h=NULL) {
  eigen(get_hessian(gsmvar, custom_h=custom_h))$value
}



#' @title Get differences 'h' which are adjusted for overly large degrees of freedom parameters
#'
#' @description \code{get_varying_h} adjusts differences for overly large degrees of freedom parameters
#'   for finite difference approximation of the derivatives of the log-likelihood function.
#'
#' @inheritParams loglikelihood_int
#' @details This function is used for approximating gradient and Hessian of a StMVAR or G-StMVAR model.
#'   Very large degrees of freedom parameters cause significant numerical error if too small differences
#'   are used.
#' @return Returns a vector with the same length as \code{params}. For other parameters than degrees
#'   of freedom parameters larger than 100, the differences will be \code{6e-6}. For the large degrees of
#'   freedom parameters, the difference will be \code{signif(df/1000, digits=2)}.
#' @keywords internal

get_varying_h <- function(M, params, model=c("GMVAR", "StMVAR", "G-StMVAR")) {
  model <- match.arg(model)
  if(model != "GMVAR") {
    all_df <- pick_df(M=M, params=params, model=model)
    adj_diffs <- numeric(length(all_df))
    adj_diffs[all_df <= 100] <- 6e-6 # h is not adjusted for dfs not larger than hundred
    adj_diffs[all_df > 100] <- signif(all_df[all_df > 100]/1000, digits=2) # The adjusted differences h
    varying_h <- c(rep(6e-6, times=length(params) - length(all_df)), adj_diffs) # Difference h for all parameters
  } else { # No degrees of freedom parameters in GMVAR model, so the default difference is used
    varying_h <- rep(6e-6, times=length(params))
  }
  varying_h
}
