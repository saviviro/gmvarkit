#' @import graphics
#'
#' @title  Conditional mean or variance plot for a GMVAR model
#'
#' @description \code{cond_moment_plot} plots the one-step in-sample conditional means/variances of the model along with
#'  the individual time series contained in the model (e.g. the time series the model was fitted to). Also plots
#'  the regimewise conditional means/variances multiplied with mixing weights.
#'
#' @inheritParams simulateGMVAR
#' @param which_moment should conditional means or variances be plotted?
#' @param grid add grid to the plots?
#' @param ... additional paramters passed to \code{grid(...)} plotting the grid if \code{grid == TRUE}.
#' @details The conditional mean plot works best if the data contains positive values only.
#'  \code{acf} from the package \code{stats} and the plot method for class \code{'acf'} objects is employed.
#' @inherit simulateGMVAR references
#' @seealso \code{\link{profile_logliks}}, \code{\link{fitGMVAR}}, \code{\link{GMVAR}}, \code{\link{quantile_residual_tests}},
#'  \code{\link{LR_test}}, \code{\link{Wald_test}}, \code{\link{diagnostic_plot}}
#' @examples
#' # GMVAR(2, 2), d=2 model;
#' params22 <- c(0.36, 0.121, 0.223, 0.059, -0.151, 0.395, 0.406, -0.005,
#'  0.083, 0.299, 0.215, 0.002, 0.03, 0.484, 0.072, 0.218, 0.02, -0.119,
#'  0.722, 0.093, 0.032, 0.044, 0.191, 1.101, -0.004, 0.105, 0.58)
#' mod22 <- GMVAR(gdpdef, p=2, M=2, params=params22)
#'
#' cond_moment_plot(mod22, which_moment="mean")
#' cond_moment_plot(mod22, which_moment="variance")
#' cond_moment_plot(mod22, which_moment="mean", grid=TRUE, lty=3)
#' @export

cond_moment_plot <- function(gmvar, which_moment=c("mean", "variance"), grid=FALSE, ...) {
  check_gmvar(gmvar)
  stopifnot(!is.null(gmvar$data))
  which_moment <- match.arg(which_moment)
  p <- gmvar$model$p
  M <- gmvar$model$M
  d <- gmvar$model$d
  data <- gmvar$data
  if(is.null(gmvar$regime_cmeans)) stop("Conditional moments were not calculated when building this model")

  if(which_moment == "mean") {
    total_moments <- gmvar$total_cmeans # [t, d]
    mw_x_reg <- lapply(1:d, function(d1) gmvar$mixing_weights*gmvar$regime_cmeans[, d1, ]) # [[d]][t, m]
    vals <- lapply(1:d, function(d1) c(total_moments[,d1], vec(mw_x_reg[[d1]]), data[,d1]))
  } else {
    total_moments <- t(vapply(1:dim(gmvar$total_ccovs)[3], function(i1) diag(gmvar$total_ccovs[, ,i1, drop=TRUE]), numeric(d))) # [t, d]
    params <- reform_constrained_pars(p=p, M=M, d=d, params=gmvar$params, constraints=gmvar$model$constraints,
                                      same_means=gmvar$model$same_means, structural_pars=gmvar$model$structural_pars)
    omegas <- pick_Omegas(p=p, M=M, d=d, params=params, structural_pars=get_unconstrained_structural_pars(gmvar$model$structural_pars))
    vars <- vapply(1:M, function(m) diag(omegas[, , m]), numeric(d)) # Regs in cols, d in rows
    mw_x_reg <- lapply(1:d, function(d1) t(t(gmvar$mixing_weights)*vars[d1,])) # [[d]][t, m]
    vals <- lapply(1:d, function(d1) c(total_moments[,d1], vec(mw_x_reg[[d1]])))
  }

  old_par <- par(no.readonly=TRUE)
  on.exit(par(old_par))
  right_marg <- ifelse(which_moment == "mean", 1, 2.5)
  graphics::par(mfrow=c(d, 1), mar=c(0.5, 2.5, 2.1, right_marg))
  colpal_reg <- grDevices::colorRampPalette(c("blue", "turquoise1", "green", "red"))(M)
  names_ts <- colnames(as.ts(data))
  names_reg <- paste0("mix.comp.", 1:M)
  make_ts <- function(dat) ts(c(rep(NA, p), dat), start=start(data), frequency=frequency(data))

  for(d1 in 1:d) {
    xaxt <- "n"
    if(d1 == d) {
      xaxt <- "s"
      par(mar=c(2.5, 2.5, 0.5, right_marg))
    } else if(d1 > 1) {
      par(mar=c(0.5, 2.5, 0.5, right_marg))
    }
    ymin <- floor(min(vals[[d1]]))
    ymax <- ceiling(max(vals[[d1]]))
    if(which_moment == "mean") {
      main <- ifelse(d1 == 1, "Conditional means", "")
      plot(data[,d1], ylim=c(ymin, ymax), xlab="", ylab="", xaxt=xaxt, main=main)
      lines(make_ts(total_moments[,d1]), col="grey", lty=2, lwd=2)
    } else { # Plot conditional variances
      main <- ifelse(d1 == 1, "Conditional variances", "")
      plot(data[,d1], xlab="", ylab="", xaxt=xaxt, main=main)
      par(new=TRUE)
      plot(make_ts(total_moments[,d1]), ylim=c(ymin, ymax), col="grey", lty=2, lwd=2, xlab="", ylab="", yaxt="n", xaxt="n")
      axis(side=4, col="grey", lwd=2)
    }
    if(grid) grid(...)
    for(m1 in 1:M) {
      lines(make_ts(mw_x_reg[[d1]][,m1]), col=colpal_reg[m1], lty=3)
    }
    legend("topleft", legend=names_ts[d1], bty="n", col="black", text.font=2, cex=0.65, x.intersp=0.5, y.intersp=1)
    if(d1 == 1) {
      legend("topright", legend=c("total", paste0("regime ", 1:M)), bty="n", col=c("grey", colpal_reg),
             lty=c(2, rep(3, M)), lwd=2, text.font=2, cex=0.65, x.intersp=0.5, y.intersp=1)
    }
  }
}


