#' @title Predict method for class 'gsmvar' objects
#'
#' @description \code{predict.gsmvar} is a predict method for class \code{'gsmvar'} objects. The forecasts of
#'   the GMVAR, StMVAR, and G-StMVAR models are computed by performing independent simulations and using the
#'   sample medians or means as point forecasts and empirical quantiles as prediction intervals.
#'   For one-step-ahead predictions using the exact conditional mean is also supported.
#'
#' @inheritParams simulate.gsmvar
#' @param n_ahead how many steps ahead should be predicted?
#' @param nsim to how many independent simulations should the forecast be based on?
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
#' @return Returns a class '\code{gsmvarpred}' object containing, among the specifications,...
#'  \describe{
#'    \item{$pred}{Point forecasts}
#'    \item{$pred_int}{Prediction intervals, as \code{[, , d]}.}
#'    \item{$mix_pred}{Point forecasts for the mixing weights}
#'    \item{mix_pred_int}{Individual prediction intervals for mixing weights, as \code{[, , m]}, m=1,..,M.}
#'  }
#' @seealso \code{\link{GIRF}}, \code{\link{GFEVD}}, \code{\link{simulate.gsmvar}}
#' @inherit in_paramspace_int references
#' @examples
#' # GMVAR(2, 2), d=2 model
#' params22 <- c(0.36, 0.121, 0.223, 0.059, -0.151, 0.395, 0.406, -0.005,
#'  0.083, 0.299, 0.215, 0.002, 0.03, 0.484, 0.072, 0.218, 0.02, -0.119,
#'   0.722, 0.093, 0.032, 0.044, 0.191, 1.101, -0.004, 0.105, 0.58)
#' mod22 <- GSMVAR(gdpdef, p=2, M=2, d=2, params=params22)
#' p1 <- predict(mod22, n_ahead=10, pred_type="median", nsim=500)
#' p1
#' p2 <- predict(mod22, n_ahead=10, nt=20, lty=1, nsim=500)
#' p2
#' p3 <- predict(mod22, n_ahead=10, pi=c(0.99, 0.90, 0.80, 0.70),
#'               nt=30, lty=0, nsim=500)
#' p3
#'
#' # StMVAR(2, 2), d=2 model
#' params22t <- c(0.36, 0.121, 0.223, 0.059, -0.151, 0.395, 0.406, -0.005,
#'  0.083, 0.299, 0.215, 0.002, 0.03, 0.484, 0.072, 0.218, 0.02, -0.119,
#'   0.722, 0.093, 0.032, 0.044, 0.191, 1.101, -0.004, 0.105, 0.58, 3, 4)
#' mod22t <- GSMVAR(gdpdef, p=2, M=2, d=2, params=params22t, model="StMVAR")
#' p1 <- predict(mod22t, n_ahead=12, pred_type="median", nsim=500, pi=0.9)
#' p1
#' @export

predict.gsmvar <- function(object, ..., n_ahead, nsim=2000, pi=c(0.95, 0.80), pi_type=c("two-sided", "upper", "lower", "none"),
                          pred_type=c("median", "mean", "cond_mean"), plot_res=TRUE, mix_weights=TRUE, nt) {
  gsmvar <- object
  check_gsmvar(gsmvar)
  check_null_data(gsmvar)
  if(is.null(gsmvar$data)) {
    stop("The model needs to contain data as forecasting requires initial values. Data can be added to the model with the function 'add_data'.")
  }
  data <- gsmvar$data
  if(missing(n_ahead)) {
    warning("Argument n_ahead is missing. Using n_ahead = 1.")
    n_ahead <- 1
  }
  if(!all_pos_ints(c(n_ahead, nsim))) stop("Arguments n_ahead and nsim must be positive integers")
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
    p <- gsmvar$model$p
    M <- gsmvar$model$M
    d <- gsmvar$model$d
    model <- gsmvar$model$model
    constraints <- gsmvar$model$constraints
    same_means <- gsmvar$model$same_means
    structural_pars <- gsmvar$model$structural_pars

    params <- gsmvar$params
    n_obs <- nrow(data)
    mw <- loglikelihood_int(data=data, p=p, M=M,
                            params=params, model=model,
                            conditional=gsmvar$model$conditional,
                            parametrization=gsmvar$model$parametrization,
                            constraints=constraints,
                            same_means=same_means,
                            structural_pars=structural_pars,
                            to_return="mw_tplus1",
                            stat_tol=gsmvar$num_tols$stat_tol,
                            posdef_tol=gsmvar$num_tols$posdef_tol)
    mw <- mw[nrow(mw),]

    # Collect parameter values
    params <- reform_constrained_pars(p=p, M=M, d=d,
                                      params=params, model=model,
                                      constraints=constraints,
                                      same_means=same_means,
                                      structural_pars=structural_pars)
    structural_pars <- get_unconstrained_structural_pars(structural_pars)
    if(gsmvar$model$parametrization == "mean") {
      params <- change_parametrization(p=p, M=M, d=d, params=params, model=model, constraints=NULL,
                                       structural_pars=structural_pars, change_to="intercept")
    }
    all_phi0 <- pick_phi0(p=p, M=M, d=d, params=params, structural_pars=structural_pars)
    all_A <- pick_allA(p=p, M=M, d=d, params=params, structural_pars=structural_pars)

    # Calculate the conditional mean
    pred <- rowSums(vapply(1:sum(M), function(m) mw[m]*(all_phi0[, m] + rowSums(vapply(1:p, function(i1) all_A[, , i1, m]%*%data[n_obs + 1 - i1,],
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
    simulations <- simulate.gsmvar(gsmvar, nsim=n_ahead, init_values=gsmvar$data, ntimes=nsim)
    sample <- simulations$sample
    alpha_mt <- simulations$mixing_weights
    colnames(sample) <- colnames(data)
    colnames(alpha_mt) <- vapply(1:sum(gsmvar$model$M), function(m) paste("regime", m), character(1))

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
    if(is.null(colnames(pred))) colnames(pred) <- vapply(1:gsmvar$model$d, function(m) paste("Comp.", m), character(1))

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
        mix_pred_ints <- array(mix_pred_ints, dim=c(n_ahead, gsmvar$model$M, length(q_tocalc)), dimnames=list(NULL, colnames(alpha_mt), q_tocalc))
        pred_ints <- aperm(pred_ints, perm=c(1, 3, 2))
        mix_pred_ints <- aperm(mix_pred_ints, perm=c(1, 3, 2))
      } else {
        pred_ints <- aperm(pred_ints, perm=c(2, 1, 3))
        mix_pred_ints <- aperm(mix_pred_ints, perm=c(2, 1, 3))
      }
      colnames(pred_ints) <- colnames(mix_pred_ints) <- q_tocalc
    }
  }

  ret <- structure(list(gsmvar=gsmvar,
                        pred=pred,
                        pred_ints=pred_ints,
                        mix_pred=mix_pred,
                        mix_pred_ints=mix_pred_ints,
                        n_ahead=n_ahead,
                        nsim=nsim,
                        pi=pi,
                        pi_type=pi_type,
                        pred_type=pred_type,
                        q=q_tocalc),
                   class="gsmvarpred")
  if(plot_res) plot.gsmvarpred(x=ret, nt=nt, mix_weights=mix_weights, ...)
  ret
}

