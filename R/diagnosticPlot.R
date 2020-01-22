#' @import graphics
#'
#' @title Quantile residual diagnostic plot for a GMVAR model
#'
#' @description \code{diagnostic_plot} plots a multivariate quantile residual diagnostic plot
#'   for either autocorrelation, conditional heteroskedasticity, or normality, or simply draws
#'   the quantile residual time series.
#'
#' @inheritParams simulateGMVAR
#' @param type which type of diagnostic plot should be plotted?
#'   \itemize{
#'     \item{\code{"series"} the quantile residual time series.}
#'     \item{\code{"ac"} the quantile residual autocorrelation and cross-correlation functions.}
#'     \item{\code{"ch"} the squared quantile residual autocorrelation and cross-correlation functions.}
#'     \item{\code{"norm"} the quantile residual histogram with theoretical standard normal
#'       density (dashed line) and standard normal QQ-plots.}
#'   }
#' @param maxlag the maximum lag considered in types \code{"ac"} and \code{"ch"}.
#' @details Auto- and cross-correlations (types \code{"ac"} and \code{"ch"}) are calculated with the function
#'  \code{acf} from the package \code{stats} and the plot method for class \code{'acf'} objects is employed.
#' @inherit quantile_residual_tests references
#' @seealso \code{\link{profile_logliks}}, \code{\link{fitGMVAR}}, \code{\link{GMVAR}}, \code{\link{quantile_residual_tests}},
#'  \code{\link[stats]{acf}}, \code{\link[stats]{density}}, \code{\link{predict.gmvar}}
#' @examples
#' # These examples use the data 'eurusd' which comes with the
#' # package, but in a scaled form.
#' data <- cbind(10*eurusd[,1], 100*eurusd[,2])
#' colnames(data) <- colnames(eurusd)
#'
#' # GMVAR(1,2), d=2 model:
#' params122 <- c(0.623, -0.129, 0.959, 0.089, -0.006, 1.006, 1.746,
#'  0.804, 5.804, 3.245, 7.913, 0.952, -0.037, -0.019, 0.943, 6.926,
#'  3.982, 12.135, 0.789)
#' mod122 <- GMVAR(data, p=1, M=2, params=params122)
#' diagnostic_plot(mod122, type="series")
#' diagnostic_plot(mod122, type="ac")
#'
#' # GMVAR(2,2), d=2 model:
#' params222 <-  c(1.386, -0.765, 1.314, 0.145, 0.094, 1.292, -0.389,
#'  -0.070, -0.109, -0.281, 0.920, -0.025, 4.839, 1.005, 5.928, 1.248,
#'   0.077, -0.040, 1.266, -0.272, -0.074, 0.034, -0.313, 5.855, 3.570,
#'   9.838, 0.740)
#' mod222 <- GMVAR(data, p=2, M=2, params=params222)
#' diagnostic_plot(mod222, type="ch")
#' diagnostic_plot(mod222, type="norm")
#'
#' # GMVAR(2,2), d=2 model with AR-parameters restricted to be
#' # the same for both regimes:
#' C_mat <- rbind(diag(2*2^2), diag(2*2^2))
#' params222c <- c(1.031, 2.356, 1.786, 3.000, 1.250, 0.060, 0.036,
#'  1.335, -0.290, -0.083, -0.047, -0.356, 0.934, -0.152, 5.201, 5.883,
#'  3.560, 9.799, 0.368)
#' mod222c <- GMVAR(data, p=2, M=2, params=params222c, constraints=C_mat)
#' diagnostic_plot(mod222c)
#' diagnostic_plot(mod222c, type="ac", maxlag=12)
#' @export

diagnostic_plot <- function(gmvar, type=c("series", "ac", "ch", "norm"), maxlag=10) {
  check_gmvar(gmvar)
  check_null_data(gmvar)
  type <- match.arg(type)
  qres <- gmvar$quantile_residuals
  colnames(qres) <- colnames(as.ts(gmvar$data))
  if(type == "series") {
    plot.ts(qres, plot.type="multiple", main="Quantile residual time series", xlab=NULL)
  } else if(type == "ac") {
    acf(qres, lag.max=maxlag, plot=TRUE)
  } else if(type == "ch") {
    acf(qres^2, lag.max=maxlag, plot=TRUE)
  } else if(type == "norm") {
    old_par <- par(no.readonly=TRUE)
    on.exit(par(old_par))
    d <- gmvar$model$d
    par(mfrow=c(2, d), mar=c(2.5, 2.5, 2.1, 1.1))
    for(i1 in 1:d) {
      hs <- hist(qres[,i1], breaks="Scott", probability=TRUE, col="skyblue", plot=TRUE,
                 main=colnames(qres)[i1], ylim=c(0, 0.5))
      x <- seq(from=min(hs$breaks), to=max(hs$breaks), length.out=1000)
      lines(x=x, y=dnorm(x), lty=2, col="darkred")
    }
    for(i1 in 1:d) {
      qqnorm(qres[,i1], main="")
      qqline(qres[,i1], col="darkred")
    }
  }
}


#' @title Plot profile log-likehoods around the estimates
#'
#' @description \code{profile_logliks} plots profile log-likelihoods around the estimates.
#'
#' @inheritParams simulateGMVAR
#' @param which_pars the profile log-likelihood function of which parameters should be plotted? An integer
#'  vector specifying the positions of the parameters in the parameter vector. The parameter vector has the
#'  form...
#'   \describe{
#'     \item{\strong{For unconstrained models:}}{
#'       Should be size \eqn{((M(pd^2+d+d(d+1)/2+1)-1)x1)} and have form
#'       \strong{\eqn{\theta}}\eqn{ = }(\strong{\eqn{\upsilon}}\eqn{_{1}},
#'       ...,\strong{\eqn{\upsilon}}\eqn{_{M}}, \eqn{\alpha_{1},...,\alpha_{M-1}}), where:
#'       \itemize{
#'         \item \strong{\eqn{\upsilon}}\eqn{_{m}} \eqn{ = (\phi_{m,0},}\strong{\eqn{\phi}}\eqn{_{m}}\eqn{,\sigma_{m})}
#'         \item \strong{\eqn{\phi}}\eqn{_{m}}\eqn{ = (vec(A_{m,1}),...,vec(A_{m,p})}
#'         \item and \eqn{\sigma_{m} = vech(\Omega_{m})}, m=1,...,M.
#'       }
#'     }
#'     \item{\strong{For constrained models:}}{
#'       Should be size \eqn{((M(d+d(d+1)/2+1)+q-1)x1)} and have form
#'       \strong{\eqn{\theta}}\eqn{ = (\phi_{1,0},...,\phi_{M,0},}\strong{\eqn{\psi}}
#'       \eqn{,\sigma_{1},...,\sigma_{M},\alpha_{1},...,\alpha_{M-1})}, where:
#'       \itemize{
#'         \item \strong{\eqn{\psi}} \eqn{(qx1)} satisfies (\strong{\eqn{\phi}}\eqn{_{1}}\eqn{,...,}
#'         \strong{\eqn{\phi}}\eqn{_{M}) =} \strong{\eqn{C \psi}}. Here \strong{\eqn{C}} is \eqn{(Mpd^2xq)}
#'         constraint matrix.
#'       }
#'     }
#'   }
#'   Above, \eqn{\phi_{m,0}} is the intercept parameter, \eqn{A_{m,i}} denotes the \eqn{i}:th coefficient matrix of the \eqn{m}:th
#'   mixture component, \eqn{\Omega_{m}} denotes the error term covariance matrix of the \eqn{m}:th mixture component, and
#'   \eqn{\alpha_{m}} is the mixing weight parameter.
#'
#'   The default is that profile log-likelihood functions for all parameters are plotted.
#' @param scale a numeric scalar specifying the interval plotted for each estimate:
#'  the estimate plus-minus \code{abs(scale*estimate)}.
#' @param nrows how many rows should be in the plot-matrix? The default is \code{max(ceiling(log2(length(which_pars)) - 1), 1)}.
#' @param ncols how many columns should be in the plot-matrix? The default is \code{ceiling(length(which_pars)/nrows)}.
#'   Note that \code{nrows*ncols} should not be smaller than the length of \code{which_pars}.
#' @param precission at how many points should each profile log-likelihood be evaluated at?
#' @details When the number of parameters is large, it might be better to plot a smaller number of profile
#'  log-likelihood functions at a time using the argument \code{which_pars}.
#'
#' The red vertical line points the estimate.
#' @return  Only plots to a graphical device and doesn't return anything.
#' @inherit loglikelihood references
#' @seealso  \code{\link{get_soc}}, \code{\link{diagnostic_plot}}, \code{\link{fitGMVAR}}, \code{\link{GMVAR}}
#' @examples
#' \donttest{
#' # These examples use the data 'eurusd' which comes with the
#' # package, but in a scaled form (similar to Kalliovirta et al. 2016).
#' data(eurusd, package="gmvarkit")
#' data <- cbind(10*eurusd[,1], 100*eurusd[,2])
#' colnames(data) <- colnames(eurusd)
#'
#' # GMVAR(1,2) model: 10 estimation rounds with seeds set
#' # for reproducibility
#' fit12 <- fitGMVAR(data, p=1, M=2, ncalls=10, seeds=1:10)
#' fit12
#' profile_logliks(fit12)

#' # GMVAR(2,2) model with mean parametrization
#' fit22 <- fitGMVAR(data, p=2, M=2, parametrization="mean",
#'                   ncalls=16, seeds=11:26)
#' profile_logliks(fit22)
#'
#' # GMVAR(2,2) model with autoregressive parameters restricted
#' # to be the same for both regimes
#' C_mat <- rbind(diag(2*2^2), diag(2*2^2))
#' fit22c <- fitGMVAR(data, p=2, M=2, constraints=C_mat)
#' profile_logliks(fit22c)
#'
#' # GMVAR(2,2) model with autoregressive parameters restricted
#' # to be the same for both regimes and non-diagonl elements
#' # the coefficient matrices constrained to zero.
#' tmp <- matrix(c(1, rep(0, 10), 1, rep(0, 8), 1, rep(0, 10), 1),
#'  nrow=2*2^2, byrow=FALSE)
#' C_mat2 <- rbind(tmp, tmp)
#' fit22c2 <- fitGMVAR(data, p=2, M=2, constraints=C_mat2)
#' profile_logliks(fit22c2)
#' }
#' @export

profile_logliks <- function(gmvar, which_pars, scale=0.02, nrows, ncols, precission=200) {
  check_gmvar(gmvar)
  check_null_data(gmvar)
  p <- gmvar$model$p
  M <- gmvar$model$M
  d <- gmvar$model$d
  params <- gmvar$params
  parametrization <- gmvar$model$parametrization
  if(missing(which_pars)) which_pars <- 1:length(params)
  if(!all_pos_ints(which_pars) || any(which_pars > length(params))) {
    stop("The argument 'which_pars' should contain strictly positive integers not larger than length of the parameter vector.")
  } else if(anyDuplicated(which_pars) != 0) {
    stop("There are dublicates in which_pars")
  }
  constraints <- gmvar$model$constraints
  npars <- length(which_pars)

  if(missing(nrows)) nrows <- max(ceiling(log2(npars) - 1), 1)
  if(missing(ncols)) ncols <- ceiling(npars/nrows)
  stopifnot(all_pos_ints(c(nrows, ncols)) && nrows*ncols >= npars)

  # Graphical settings: restore on exit.
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  par(mar=c(2.1, 2.1, 1.6, 1.1), mfrow=c(nrows, ncols))

  # In order to get the labels right, we first determine which indeces in params
  # correspond to which parameters: different procedure for constrained models.
  if(is.null(constraints)) {
    all_q <- rep(d^2*p + d + d*(d + 1)/2, M) # Length of (phi_0m, \bold{phi_m}, sigma_m) for each m.
    cum_q <- c(0, cumsum(all_q)) # After this index, a new regime starts; after the last one the mixing weights start
  } else {
    q <- ncol(constraints)
  }

  for(i1 in which_pars) { # Go though the parameters

    pars <- params
    range <- abs(scale*pars[i1])
    vals <- seq(from=pars[i1] - range, to=pars[i1] + range, length.out=precission) # Loglik to be evaluated at these values of the parameter considered
    logliks <- vapply(vals, function(val) {
      new_pars <- pars
      new_pars[i1] <- val # Change the single parameter value
      loglikelihood_int(data=gmvar$data, p=p, M=M, params=new_pars, constraints=constraints,
                        conditional=gmvar$model$conditional, parametrization=parametrization,
                        check_params=TRUE, minval=NA)
    }, numeric(1))

    # In order to get the labels right, we first determine which parameter is in question.
    # We consider constrained models separately.
    if(is.null(constraints)) {
      if(i1 <= max(cum_q)) { # phi and sigma parameters first
        m <- sum(i1 > cum_q) # Which regime are we considering
        if(i1 > cum_q[m + 1] - d*(d + 1)/2 && i1 <= cum_q[m + 1]) { # Omega params
          pos <- i1 - (cum_q[m] + d + p*d^2) # Position in vech(Omega_m)
          cum_d <- c(0, cumsum(d - 0:(d - 1))) # Index after which column changes in the current vech(Omega_m)
          col_ind <- sum(pos > cum_d) # Which column in the current Omega
          row_inds <- unlist(lapply(1:d, function(i2) i2:d)) # At which row are we in the current Omega_m for each pos?
          row_ind <- row_inds[pos] # At which row of Omega_m we are
          main <- substitute(Omega[foo](foo2), list(foo=m, foo2=paste0(row_ind, ",", col_ind)))
        } else if(i1 <= max(cum_q)) { # The phi parameters (or mean + AR parameters)
          if(i1 > cum_q[m] && i1 <= cum_q[m] + d) { # phi_{m,0} or mean parameters
            pos <- i1 - cum_q[m] # Which time series? 1,..,d
            mylist <- list(foo=paste0(m, ",0"), foo2=pos)
            if(parametrization == "intercept") {
              main <- substitute(phi[foo](foo2), mylist)
            } else {
              main <- substitute(mu[foo](foo2), mylist)
            }
          } else {  # The elements of A_m1,...,A_mp
            pos1 <- i1 - (cum_q[m] + d) # Position in vec(A_m1),...,vec(A_mp)
            cum_a <- c(0, cumsum(rep(d^2, times=p))) # the index after which new matrix A_m,p starts
            which_mat <- sum(pos1 > cum_a) # in which matrix A_m,j, j=1,..,p we are?
            pos2 <- pos1 - (which_mat - 1)*d^2 # Position in the current matrix A_m,j
            cum_a_cols <- c(0, cumsum(rep(d, times=d))) # The index in currect matrix after which a new column in A_m,j starts
            col_ind <- sum(pos2 > cum_a_cols) # The column where we are in the current matrix
            row_ind <- rep(1:d, times=d)[pos2] # The row where we are in the current matrix
            main <- substitute(A[foo](foo2), list(foo=paste0(m, ",", which_mat), foo2=paste0(row_ind, ",", col_ind)))
          }
        }
      } else { # alphas; we know that M > 1 by the fact that we are here
            m <- i1 - max(cum_q)
            main <- substitute(alpha[foo], list(foo=m))
      }
    } else { ## Linear constraints are employed

      if(i1 <= M*d + q) { # phi_{m,0} and AR parameters
        if(i1 <= M*d) { # phi_{m,0}
          cum_d <- c(0, cumsum(rep(d, M))) # The index after which regime changes
          m <- sum(i1 > cum_d)
          pos <- i1 - cum_d[m] # Which Time series?
          mylist <- list(foo=paste0(m, ",0"), foo2=pos)
          if(parametrization == "intercept") {
            main <- substitute(phi[foo](foo2), mylist)
          } else {
            main <- substitute(mu[foo](foo2), mylist)
          }
        } else { # The AR parameters
          pos <- i1 - M*d
          main <- substitute(AR(foo), list(foo=pos))
        }
      } else if(i1 <= M*(d + d*(d + 1)/2) + q) { # Omega parameters
        cum_s <- M*d + q + c(0, cumsum(rep(d*(d + 1)/2, times=M))) # Index after which regime changes
        m <- sum(i1 > cum_s)
        i1 - cum_s[1] # Position in vech(Omega_1),...,vech(Omega_M)
        pos <- i1 - cum_s[m] # position in vech(Omega_m)
        cum_d <- c(0, cumsum(d - 0:(d - 1))) # Index after which column changes in the current vech(Omega_m)
        col_ind <- sum(pos > cum_d) # Which column in the current Omega
        row_inds <- unlist(lapply(1:d, function(i2) i2:d)) # At which row are we in the current Omega_m for each pos?
        row_ind <- row_inds[pos] # At which row of Omega_m we are
        main <- substitute(Omega[foo](foo2), list(foo=m, foo2=paste0(row_ind, ", ", col_ind)))
      } else { # alphas; we know M > 1 since we ended up here
        m <- i1 - (M*(d + d*(d + 1)/2) + q)
        main <- substitute(alpha[foo], list(foo=m))
      }
    }
    plot(x=vals, y=logliks, type="l", main=main)
    abline(v=pars[i1], col="red") # Points the estimate
  }
}
