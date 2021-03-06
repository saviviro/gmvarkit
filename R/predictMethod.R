#' @title Predict method for class 'gmvar' objects
#'
#' @description Forecast GMVAR process defined as a class \code{'gmvar'} object. The forecasts are
#'   computed by performing independent simulations and using the sample medians or means as point
#'   forecasts and empirical quantiles as prediction intervals. For one-step-ahead predictions
#'   using the exact conditional mean is also supported.
#'
#' @param object an object of class \code{'gmvar'}, generated by function \code{fitGMVAR} or \code{GMVAR}.
#' @param n_ahead how many steps ahead should be predicted?
#' @param n_simu to how many independent simulations should the forecast be based on?
#' @param pi a numeric vector specifying the confidence levels of the prediction intervals.
#' @param pi_type should the prediction intervals be "two-sided", "upper", or "lower"?
#' @param pred_type should the prediction be based on sample "median" or "mean"? Or should it
#'   be one-step-ahead forecast based on the exact conditional mean (\code{"cond_mean"})?
#'   Prediction intervals won't be calculated if the exact conditional mean is used.
#' @param plot_res should the results be plotted?
#' @param mix_weights \code{TRUE} if forecasts for mixing weights should be plotted,
#'   \code{FALSE} in not.
#' @param nt a positive integer specifying the number of observations to be plotted
#'   along with the prediction (ignored if \code{plot_res==FALSE}). Default is \code{round(nrow(data)*0.15)}.
#' @param ... additional arguments passed to \code{grid} (ignored if \code{plot_res==FALSE}) which plots
#'   grid to the figure.
#' @return Returns a class '\code{gmvarpred}' object containing, among the specifications,...
#'  \describe{
#'    \item{$pred}{Point forecasts}
#'    \item{$pred_int}{Prediction intervals, as \code{[, , d]}.}
#'    \item{$mix_pred}{Point forecasts for the mixing weights}
#'    \item{mix_pred_int}{Individual prediction intervals for mixing weights, as \code{[, , m]}, m=1,..,M.}
#'  }
#' @seealso \code{\link{GIRF}}, \code{\link{GFEVD}}, \code{\link{simulateGMVAR}}
#' @inherit in_paramspace_int references
#' @examples
#' # These examples use the data 'eurusd' which comes with the
#' # package, but in a scaled form.
#' data <- cbind(10*eurusd[,1], 100*eurusd[,2])
#' colnames(data) <- colnames(eurusd)
#'
#' # GMVAR(2,2) model
#' params22 <- c(1.386, -0.767, 1.314, 0.145, 0.094, 1.292, -0.389, -0.07,
#'  -0.109, -0.281, 0.92, -0.025, 4.839, 0.998, 5.916, 1.248, 0.077, -0.04,
#'  1.266, -0.272, -0.074, 0.034, -0.313, 5.855, 3.569, 9.837, 0.741)
#' fit22 <- GMVAR(data, p=2, M=2, params=params22)
#' p1 <- predict(fit22, n_ahead=10, pred_type="median", n_simu=500)
#' p1
#' p2 <- predict(fit22, n_ahead=10, nt=20, lty=1, n_simu=500)
#' p2
#' p3 <- predict(fit22, n_ahead=10, pi=c(0.99, 0.90, 0.80, 0.70),
#'               nt=30, lty=0, n_simu=500)
#' p3
#'
#' # Structural GMVAR(2, 2), d=2 model identified with sign-constraints:
#' params222s <- c(-11.964, 155.024, 11.636, 124.988, 1.314, 0.145, 0.094, 1.292,
#'  -0.389, -0.07, -0.109, -0.281, 1.248, 0.077, -0.04, 1.266, -0.272, -0.074,
#'   0.034, -0.313, 0.903, 0.718, -0.324, 2.079, 7.00, 1.44, 0.742)
#' W_222 <- matrix(c(1, 1, -1, 1), nrow=2, byrow=FALSE)
#' mod222s <- GMVAR(data, p=2, M=2, params=params222s, parametrization="mean",
#'  structural_pars=list(W=W_222))
#' p1 <- predict(mod222s, n_ahead=10, n_simu=500)
#' @export

predict.gmvar <- function(object, ..., n_ahead, n_simu=2000, pi=c(0.95, 0.80), pi_type=c("two-sided", "upper", "lower", "none"),
                          pred_type=c("median", "mean", "cond_mean"), plot_res=TRUE, mix_weights=TRUE, nt) {
  gmvar <- object
  check_gmvar(gmvar)
  check_null_data(gmvar)
  if(is.null(gmvar$data)) {
    stop("The model needs to contain data as forecasting requires initial values. Data can be added to the model with the function 'add_data'.")
  }
  data <- gmvar$data
  if(missing(n_ahead)) {
    warning("Argument n_ahead is missing. Using n_ahead = 1.")
    n_ahead <- 1
  }
  if(!all_pos_ints(c(n_ahead, n_simu))) stop("Arguments n_ahead and n_simu must be positive integers")
  stopifnot(all(pi > 0 & pi < 1))
  pi_type <- match.arg(pi_type)
  pred_type <- match.arg(pred_type)
  stopifnot(pi_type %in% c("two-sided", "upper", "lower", "none"))
  stopifnot(pred_type %in% c("mean", "median", "cond_mean"))
  if(missing(nt)) {
    nt <- round(nrow(data)*0.15)
  } else {
    stopifnot(nt > 0 & nt %% 1 == 0)
    if(nt > nrow(data)) {
      warning("nt > nrow(data), using nt = nrow(data)")
      nt <- nrow(data)
    }
  }

  if(pred_type == "cond_mean") {
    if(n_ahead > 1) warning("Only one-step-ahead forecasts are supported for prediction type 'cond_mean'")
    p <- gmvar$model$p
    M <- gmvar$model$M
    d <- gmvar$model$d
    constraints <- gmvar$model$constraints
    same_means <- gmvar$model$same_means
    structural_pars <- gmvar$model$structural_pars

    params <- gmvar$params
    n_obs <- nrow(data)
    mw <- loglikelihood_int(data, p, M,
                            params=params,
                            conditional=gmvar$model$conditional,
                            parametrization=gmvar$model$parametrization,
                            constraints=constraints,
                            same_means=same_means,
                            structural_pars=structural_pars,
                            to_return="mw_tplus1",
                            stat_tol=gmvar$num_tols$stat_tol,
                            posdef_tol=gmvar$num_tols$posdef_tol)
    mw <- mw[nrow(mw),]

    # Collect parameter values
    params <- reform_constrained_pars(p=p, M=M, d=d,
                                      params=params,
                                      constraints=constraints,
                                      same_means=same_means,
                                      structural_pars=structural_pars)
    structural_pars <- get_unconstrained_structural_pars(structural_pars)
    if(gmvar$model$parametrization == "mean") {
      params <- change_parametrization(p=p, M=M, d=d, params=params, constraints=NULL,
                                       structural_pars=structural_pars, change_to="intercept")
    }
    all_phi0 <- pick_phi0(p=p, M=M, d=d, params=params, structural_pars=structural_pars)
    all_A <- pick_allA(p=p, M=M, d=d, params=params, structural_pars=structural_pars)

    # Calculate the conditional mean
    pred <- rowSums(vapply(1:M, function(m) mw[m]*(all_phi0[, m] + rowSums(vapply(1:p, function(i1) all_A[, , i1, m]%*%data[n_obs + 1 - i1,],
                                                                                  numeric(d)))), numeric(d)))
    pi <- NULL
    pi_type <- "none"
    pred_ints <- NULL
    q_tocalc <- numeric(0)
    mix_weights <- FALSE
    mix_pred <- NULL
    mix_pred_ints <- NULL
  } else { # pred_type != cond_mean

    # Simulations
    simulations <- simulateGMVAR(gmvar, nsimu=n_ahead, init_values=gmvar$data, ntimes=n_simu)
    sample <- simulations$sample
    alpha_mt <- simulations$mixing_weights
    colnames(sample) <- colnames(data)
    colnames(alpha_mt) <- vapply(1:gmvar$model$M, function(m) paste("regime", m), character(1))

    # Calculate quantiles from the third dimension of 3D simulation array
    dim3_quantiles <- function(x, q) {
      apply(x, c(1, 2), quantile, probs=q, names=FALSE)
    }

    # Predictions
    if(pred_type == "mean") {
      pred <- rowMeans(sample, dims=2)
      mix_pred <- rowMeans(alpha_mt, dims=2)
    } else {
      pred <- dim3_quantiles(sample, q=0.5)
      mix_pred <- dim3_quantiles(alpha_mt, q=0.5)
    }
    if(is.null(colnames(pred))) colnames(pred) <- vapply(1:gmvar$model$d, function(m) paste("Comp.", m), character(1))

    # Prediction intervals
    if(pi_type == "upper") {
      q_tocalc <- pi
    } else if(pi_type == "lower") {
      q_tocalc <- 1 - pi
    } else if(pi_type == "two-sided") {
      lower <- (1 - pi)/2
      upper <- rev(1 - lower)
      q_tocalc <- c(lower, upper)
    } else {
      q_tocalc <- numeric(0)
      pi <- NULL
    }

    q_tocalc <- sort(q_tocalc, decreasing=FALSE)
    pred_ints <- dim3_quantiles(sample, q_tocalc)
    mix_pred_ints <- dim3_quantiles(alpha_mt, q_tocalc)

    if(pi_type != "none") {
      if(length(q_tocalc) == 1) {
        pred_ints <- array(pred_ints, dim=c(n_ahead, ncol(data), length(q_tocalc)), dimnames=list(NULL, colnames(sample), q_tocalc)) # Make it an array with length(q_tocalc) slices
        mix_pred_ints <- array(mix_pred_ints, dim=c(n_ahead, gmvar$model$M, length(q_tocalc)), dimnames=list(NULL, colnames(alpha_mt), q_tocalc))
        pred_ints <- aperm(pred_ints, perm=c(1, 3, 2))
        mix_pred_ints <- aperm(mix_pred_ints, perm=c(1, 3, 2))
      } else {
        pred_ints <- aperm(pred_ints, perm=c(2, 1, 3))
        mix_pred_ints <- aperm(mix_pred_ints, perm=c(2, 1, 3))
      }
      colnames(pred_ints) <- colnames(mix_pred_ints) <- q_tocalc
    }
  }

  ret <- structure(list(gmvar=gmvar,
                        pred=pred,
                        pred_ints=pred_ints,
                        mix_pred=mix_pred,
                        mix_pred_ints=mix_pred_ints,
                        n_ahead=n_ahead,
                        n_simu=n_simu,
                        pi=pi,
                        pi_type=pi_type,
                        pred_type=pred_type,
                        q=q_tocalc),
                   class="gmvarpred")
  if(plot_res) plot.gmvarpred(x=ret, nt=nt, mix_weights=mix_weights, ...)
  ret
}

