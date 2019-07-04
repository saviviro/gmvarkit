#' @import graphics
#'
#' @title Quantile residual diagnostic plot for GMVAR model
#'
#' @description \code{diagnostic_plot} plots a multivariate quantile residual diagnostic plot
#'   for either autocorrelation, conditional heteroskedasticity, normality or simply draws the
#'   quantile residual time series.
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
#'  \code{acf} from the package \code{stats} and the plot method for class \code{'acf'} objects is used.
#' @inherit quantile_residual_tests references
#' @seealso \code{\link{fitGMVAR}}, \code{\link{GMVAR}}, \code{\link{quantile_residual_tests}},
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
