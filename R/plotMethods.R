#' @import graphics
#'
#' @title plot method for class 'gsmvarpred' objects
#'
#' @description \code{plot.gsmvarpred} is plot method for gsmvarpred objects.
#'
#' @inheritParams predict.gsmvar
#' @param x object of class \code{'gsmvarpred'} generated by \code{predict.gsmvar}.
#' @param add_grid should grid be added to the plots?
#' @param ... arguments passed to \code{grid} which plots grid to the figure.
#' @details This method is used plot forecasts of GSMVAR processes
#' @inherit in_paramspace_int references
#' @export

plot.gsmvarpred <- function(x, ..., nt, mix_weights=TRUE, add_grid=TRUE) {
  old_par <- par(no.readonly=TRUE)
  on.exit(par(old_par))
  gsmvarpred <- x
  data <- as.ts(gsmvarpred$gsmvar$data)
  n_obs <- nrow(data)
  d <- ncol(data)
  q <- gsmvarpred$q
  M <- sum(gsmvarpred$gsmvar$model$M)
  mixing_weights <- gsmvarpred$gsmvar$mixing_weights
  n_mix <- nrow(mixing_weights)
  mix_weights <- mix_weights && M > 1 # Don't plot mixing weights if M == 1
  if(missing(nt)) {
    nt <- round(nrow(data)*0.15)
  } else {
    stopifnot(nt > 0 & nt %% 1 == 0)
    if(nt > nrow(data)) {
      warning("nt > nrow(data), using nt = nrow(data)")
      nt <- nrow(data)
    }
  }
  if(mix_weights) {
    par(mfrow=c(d + 1, 1), mar=c(2.5, 2.5, 2.1, 1))
  } else {
    par(mfrow=c(d, 1), mar=c(2.5, 2.5, 2.1, 1))
  }
  make_ts <- function(x, mix=FALSE) { # Make ts that has the first value the same as the last value of the observed series/estim. m.weights.
    if(mix) {
      last_obs <- mixing_weights[n_mix,]
    } else {
      last_obs <- data[n_obs,]
    }
    ts(rbind(last_obs, x), start=time(data)[n_obs], frequency=frequency(data))
  }
  ts_pred <- make_ts(gsmvarpred$pred)
  ts_mix_pred <- make_ts(gsmvarpred$mix_pred, mix=TRUE)
  ts_dat <- ts(data[(n_obs - nt + 1):n_obs,], start=time(data)[n_obs - nt + 1], frequency=frequency(data))
  ts_mix <- ts(mixing_weights[(n_mix - nt + 1):n_mix,], start=time(data)[n_obs - nt + 1], frequency=frequency(data))
  t0 <- time(ts_pred)
  ts_names <- attributes(data)$dimnames[[2]]
  reg_names <- attributes(gsmvarpred$mix_pred)$dimnames[[2]]
  if(gsmvarpred$pi_type != "none") {
    pred_ints <- aperm(gsmvarpred$pred_ints, perm=c(1, 3, 2)) # [step, series, quantiles]
    mix_pred_ints <- aperm(gsmvarpred$mix_pred_ints, perm=c(1, 3, 2)) # [step, series, quantiles]
  }

  # All values to indicate ylims
  if(gsmvarpred$pi_type == "none") {
    all_val <- lapply(1:d, function(j) c(ts_dat[,j], ts_pred[,j]))
  } else {
    all_val <- lapply(1:d, function(j) c(ts_dat[,j], ts_pred[,j], simplify2array(pred_ints, higher=TRUE)[, j, ]))
  }

  # Prediction intervals, we lapply through quantiles [, , q]
  ts_fun_fact <- function(inds) function(pred_ints, mix=FALSE) lapply(inds, function(i1) make_ts(pred_ints[, , i1], mix))
  if(gsmvarpred$pi_type == "two-sided") {
    ts1_lapply <- ts_fun_fact(1:(length(q)/2)) # Lower bounds
    ts2_lapply <- ts_fun_fact((length(q)/2 + 1):length(q)) # Upper bounds

    ints1 <- pred_ints
    ints1_mix <- mix_pred_ints

  } else { # Lower or upper
    ts1_lapply <- function(pred_ints, mix=FALSE) lapply(1:length(q), function(i1) make_ts(pred_ints, mix)[-1,]) # Lower or upper bound (dummy)
    ts2_lapply <- ts_fun_fact(seq_along(q))

    myfun <- ifelse(gsmvarpred$pi_type == "upper", min, max)
    myfun2 <- ifelse(gsmvarpred$pi_type == "upper", function(x) x - 3, function(x) x + 3)
    ints1 <- vapply(1:d, function(j) rep(myfun2(round(myfun(all_val[[j]]))), times=nrow(ts_pred)), numeric(nrow(ts_pred)))
    ints1_mix <- matrix(ifelse(gsmvarpred$pi_type == "upper", 0, 1), nrow=nrow(ts_pred), ncol=M)
  }

  ts1 <- ts1_lapply(ints1)
  ts2 <- ts2_lapply(pred_ints)

  if(mix_weights) {
    ts1_mix <- ts1_lapply(ints1_mix, mix=TRUE)
    ts2_mix <- ts2_lapply(mix_pred_ints, mix=TRUE)
  }

  # Plot forecasts for the series
  draw_poly <- function(ts1_or_ts2, pred_ts, col) polygon(x=c(t0, rev(t0)), y=c(ts1_or_ts2, rev(pred_ts)), col=col, border=NA)
  col_pred <- grDevices::rgb(0, 0, 1, 0.2)
    for(i1 in 1:d) {
    ts.plot(ts_dat[,i1], ts_pred[,i1], gpars=list(col=c("black", "blue"),
                                                  lty=1:2,
                                                  ylim=c(floor(min(all_val[[i1]])),
                                                         ceiling(max(all_val[[i1]]))),
                                                  main=ts_names[i1]))
    if(add_grid) grid(...)
    if(gsmvarpred$pi_type %in% c("two-sided", "upper", "lower")) {
      for(i2 in 1:length(gsmvarpred$pi)) {
        draw_poly(ts1[[i2]][,i1], ts_pred[,i1], col=col_pred)
        draw_poly(ts2[[i2]][,i1], ts_pred[,i1], col=col_pred)
      }
    }
  }

  # Plot forecasts for the mixing weights
  if(mix_weights) {

    # Point forecasts
    colpal_mw <- grDevices::colorRampPalette(c("blue", "turquoise1", "green", "red"))(M)
    colpal_mw2 <- grDevices::adjustcolor(colpal_mw, alpha.f=0.5)
    ts.plot(ts_mix, ts_mix_pred, gpars=list(col=c(colpal_mw2, colpal_mw),
                                            ylim=c(0, 1), lty=c(rep(1, M), rep(2, M)),
                                            main="Mixing weights"))
    legend("topleft", legend=paste0("regime ", 1:M), bty="n", col=colpal_mw, lty=1, lwd=2,
           text.font=2, cex=0.9, x.intersp=0.5, y.intersp=1)
    if(add_grid) grid(...)

    # Individual prediction intervals as for the mixing weights
    colpal_mw3 <- grDevices::adjustcolor(colpal_mw, alpha.f=0.2)
    if(gsmvarpred$pi_type %in% c("two-sided", "upper", "lower")) {
      for(m in 1:M) { # Go through regimes
        for(i2 in 1:length(gsmvarpred$pi)) { # Go through the prediction intervals
          draw_poly(ts1_mix[[i2]][,m], ts_mix_pred[,m], col=colpal_mw3[m])
          draw_poly(ts2_mix[[i2]][,m], ts_mix_pred[,m], col=colpal_mw3[m])
        }
      }
    }

  }

  invisible(gsmvarpred)
}


#' @describeIn GIRF plot method
#' @inheritParams print.girf
#' @param add_grid should grid be added to the plots?
#' @param margs numeric vector of length four that adjusts the
#'  \code{[bottom_marginal, left_marginal, top_marginal, right_marginal]}
#'  as the relative sizes of the marginals to the figures of the responses.
#' @param ... arguments passed to \code{grid} which plots grid to the figure.
#' @export

plot.girf <- function(x, add_grid=FALSE, margs, ...) {

  # Relevant statistics etc
  girf <- x
  girf_res <- girf$girf_res
  nresps <- ncol(girf_res[[1]]$point_est)
  resp_names <- colnames(girf_res[[1]]$point_est)
  ngirfs <- length(girf_res)

  # Graphical settings
  if(missing(margs)) {
    margs <- c(max(0.4, 0.4 + log(0.31 + log(nresps))/6),
               max(0.4, 0.4 + log(0.31 + log(ngirfs))),
               max(0.35, 0.35 + log(0.31 + log(nresps))/6),
               max(0.1, 0.1 + log(0.31 + log(ngirfs))/10))
    if(ngirfs == 1) margs[2] <- 0.3
    margs <- vapply(1:length(margs), function(i1) min(margs[i1], 1), numeric(1))
  } else {
    stopifnot(all(margs > 0))
  }
  old_par <- par(no.readonly=TRUE)
  on.exit(par(old_par))
  par(las=1, mar=c(0, 0, 0, 0))
  nrows <- nresps + 2 # + 2 for bottom and top marginals
  ncols <- 3*ngirfs # 3x for the left and right margin in each column of figures
  nfigs <- nrows*ncols
  layoutmat <- matrix(seq_len(nfigs), nrow=nrows, ncol=ncols, byrow=FALSE)
  layout(layoutmat, # Below -0.2 for not including the ylab and also adding to the right marginals
         widths=c(margs[2], 1, margs[4], rep(c(margs[2] - 0.2, 1, margs[4]), times=ngirfs - 1)),
         heights=c(margs[3], rep(1, times=nrows - 2), margs[1]))

  # Function to plot empty plots (for the marginals)
  empty_plot <- function() plot(0, xaxt='n', yaxt='n', bty='n', pch='', ylab='', xlab='')

  # Function to plot the GIRF for each response separately
  plot_girf <- function(resp_ind, main="", xaxt="n", ylab="", first_col=FALSE) {

    # Plot point estimate
    point_est <- girf_i1$point_est[, resp_ind]
    conf_ints <- girf_i1$conf_ints[, , resp_ind]
    plot(x=0:(length(point_est) - 1), y=point_est, type="l", ylim=c(min(0, min(conf_ints)), max(0, max(conf_ints))),
         main="", ylab="", xlab="", xaxt=xaxt, lwd=2, col="blue")
    if(first_col) mtext(resp_names[resp_ind], side=2, cex=0.8, font=2, las=0, padj=-4) # Add yaxis label to the first column of responses
    if(resp_ind == 1) mtext(main, padj=-0.5, cex=1, font=2)
    if(add_grid) grid(...)

    # Plot confidence intervals
    inds <- 0:girf$N
    draw_poly <- function(up_or_low) polygon(x=c(inds, rev(inds)), y=c(up_or_low, rev(point_est)), col=grDevices::rgb(0, 0, 1, 0.2), border=NA)

    for(i1 in 1:length(girf$ci)) {
      draw_poly(conf_ints[, i1]) # lower
      draw_poly(conf_ints[, ncol(conf_ints) + 1 - i1]) # upper
    }

    abline(h=0, lty=3, col="red")
  }

  # Loop through the shocks
  for(i1 in 1:ngirfs) {
    # Plot a column of empty plots as the left margins
    for(i2 in 1:(nrows + 1)) { # + 1 for the top margin of the first row of responses
      empty_plot()
    }

    # Plot the responses of each variable to shock i1
    girf_i1 <- girf_res[[i1]]

    # Plot the GIRF of shocks i1
    first_col <- i1 == 1
    plot_girf(resp_ind=1, main=paste("Shock", girf$shocks[i1]),
              ylab=resp_names[1], first_col=first_col)
    if(nresps > 2) {
      for(i2 in 2:(nresps - 1)) {
        plot_girf(resp_ind=i2, ylab=resp_names[i2], first_col=first_col)
      }
    }
    plot_girf(resp_ind=nresps, xaxt="s", ylab=resp_names[nresps], first_col=first_col)
    empty_plot() # To bottom margin of the last row of responses

    # Plot a column of empty plots as the right margins
    for(i2 in 1:nrows) {
      empty_plot()
    }
  }
}



#' @describeIn GFEVD plot method
#' @inheritParams print.gfevd
#' @param ... currently not used.
#' @export

plot.gfevd <- function(x, ...) {
  old_par <- par(no.readonly=TRUE)
  on.exit(par(old_par))

  gfevd <- x
  gfevd_res <- gfevd$gfevd_res
  n_gfevds <- dim(gfevd_res)[3]
  varnames <- dimnames(gfevd_res)[[3]]
  n_shocks <- dim(gfevd_res)[2]
  graphics::par(las=1, mfrow=c(1, 1), mar=c(2.6, 2.6, 2.6, 3.1))

  # Function to plot the GFEVD for each variable separately
  plot_gfevd <- function(var_ind, main) {

    one_gfevd <- gfevd_res[, , var_ind] # [horizon, shock]
    mycums <- as.matrix(1 - apply(one_gfevd[, 1:(ncol(one_gfevd) - 1), drop=FALSE], MARGIN=1, FUN=cumsum))
    if(ncol(mycums) > 1) mycums <- t(mycums)
    upper_ints <- cbind(rep(1, times=nrow(one_gfevd)), mycums)
    lower_ints <- cbind(mycums,rep(0, times=nrow(one_gfevd)))
    x_points <- seq(from=-0.5, to=gfevd$N + 0.5, by=1)

    # Plot the template
    colpal <- grDevices::adjustcolor(grDevices::topo.colors(ncol(one_gfevd)), alpha.f=0.3)
    xlim_adj <- ifelse(gfevd$N < 13, 0, 0.04*(gfevd$N - 12))
    plot(NA, ylim=c(0, 1), xlim=c(0 + xlim_adj, gfevd$N - xlim_adj),
         main=main, ylab="", xlab="")
    x_bars <- (0:(gfevd$N + 1) - 0.5)
    segments(x0=x_bars, y0=rep(0, times=length(x_bars)), x1=x_bars, y1=rep(1, times=length(x_bars)))

    # Go through the shocks
    for(i1 in 1:ncol(one_gfevd)) {
      for(i2 in 1:(length(x_points) - 1)) {
        polygon(x=c(x_points[c(i2, i2 + 1)], x_points[c(i2 + 1, i2)]),
                y=c(upper_ints[c(i2, i2), i1], lower_ints[c(i2, i2), i1]),
                col=colpal[i1], border=NA)
      }
    }

    # Add shock legengs
    ylim_ajd <- ifelse(n_shocks < 500, 0.02*n_shocks, 0.01*n_shocks)
    colpal2 <- grDevices::adjustcolor(colpal, alpha.f=3)
    xtext_adj <- ifelse(gfevd$N < 10, 0.75 - gfevd$N^(-1/3), 0.4 - gfevd$N/150)
    text(x=rep(gfevd$N + xtext_adj, n_shocks), y=seq(from=1, to=1 - 0.02*n_shocks, length.out=n_shocks),
         labels=paste("Shock", 1:n_shocks),
         col=colpal2, pos=4, font=2, cex=0.8, xpd=TRUE)

  }

  # Loop through the GFEVDs
  for(i1 in 1:n_gfevds) {
    if(i1 > 1) grDevices::devAskNewPage(TRUE)
    plot_gfevd(var_ind=i1, main=ifelse(i1 <= n_shocks,
                                       paste("GFEVD for ", varnames[i1]),
                                       paste("GFEVD for regime", i1 - n_shocks, "mix. weight")))

  }
}



#' @import graphics
#' @describeIn GSMVAR plot method for class 'gsmvar'
#' @param x object of class \code{'gsmvar'} generated by \code{fitGSMVAR} or \code{GSMVAR}.
#' @param ... currently not used.
#' @param type which type figure should be produced? Or both?
#' @details The first plot displays the time series together with estimated mixing weights.
#'   The second plot displays a (Gaussian) kernel density estimates of the individual series
#'   together with the marginal stationary density implied by the model. The colored regimewise
#'   stationary densities are multiplied with the mixing weight parameter estimates.
#' @export

plot.gsmvar <- function(x, ..., type=c("both", "series", "density")) {
  type <- match.arg(type)
  gsmvar <- x
  check_null_data(gsmvar)
  data <- as.ts(gsmvar$data)
  n_obs <- nrow(data)
  p <- gsmvar$model$p
  M <- gsmvar$model$M
  d <- ncol(data)
  model <- gsmvar$model$model
  params <- reform_constrained_pars(p=p, M=M, d=d, params=gsmvar$params, model=model,
                                    constraints=gsmvar$model$constraints,
                                    same_means=gsmvar$model$same_means,
                                    weight_constraints=gsmvar$model$weight_constraints,
                                    structural_pars=gsmvar$model$structural_pars)
  all_df <- pick_df(M=M, params=params, model=model)
  alphas <- pick_alphas(p=p, M=M, d=d, params=params, model=model)
  if(model == "GMVAR") {
    M1 <- M
  } else if(model == "StMVAR") {
    M1 <- 0
  } else { # model == "G-StMVAR"
    M1 <- M[1]
  }
  M <- sum(M)
  ts_mw <- ts(rbind(matrix(NA, nrow=p, ncol=M), gsmvar$mixing_weights),
              start=start(data), frequency=frequency(data)) # First p observations are starting values

  # The 1-dimensional marginal Student's density function for StMVAR type regimes (in log first because non-log gamma function may produce Inf)
  my_dt <- function(y, mean, var, df) {
    exp(lgamma(0.5*(1 + df)) - lgamma(0.5*df) - 0.5*log(pi*(df - 2)) - 0.5*log(var) -
      0.5*(1 + df)*log(1 + (y - mean)^2/(var*(df - 2))))
  }

  old_par <- par(no.readonly=TRUE)
  on.exit(par(old_par))
  graphics::par(mfrow=c(2, 1), mar=c(2.5, 3.0, 2.1, 1), las=1)
  colpal_ts <- grDevices::colorRampPalette(c("darkgreen", "darkblue", "darkmagenta", "red3"))(d)
  colpal_mw <- grDevices::colorRampPalette(c("blue", "turquoise1", "green", "red"))(M)
  names_ts <- colnames(data)
  names_mw <- paste0("mix.comp.", 1:M)
  draw_legend <- function(nams, cols) {
    legend("topleft", legend=nams, bty="n", col=cols, lty=1, lwd=2, text.font=2, cex=0.6, x.intersp=0.5, y.intersp=1)
  }

  if(type %in% c("both", "series")) {
    # Time series and mixing weights
    ts.plot(data, gpars=list(main="Time series", col=colpal_ts, lty=1:d))
    draw_legend(names_ts, cols=colpal_ts)
    ts.plot(ts_mw, gpars=list(main="Mixing weights", ylim=c(0, 1), col=colpal_mw, lty=2))
    draw_legend(names_mw, cols=colpal_mw)
  }

  if(type == "both") {
    grDevices::devAskNewPage(TRUE)
  }


  if(type %in% c("both", "density")) {
    # Marginal stationary distributions
    nrows <- max(ceiling(log2(d) - 1), 1)
    ncols <- ceiling(d/nrows)
    graphics::par(mfrow=c(nrows, ncols), mar=c(2.5, 2.5, 2.1, 1))

    means <- get_regime_means(gsmvar) # [d, m]
    vars <- vapply(1:M, function(m) diag(get_regime_autocovs(gsmvar)[, , 1, m]), numeric(d)) # [d, m]
    reg_dens <- function(d1, m, xx) { # Marginal stat dens of d1:th series and m:th regime multiplied with mix weight parameter
      if(m <= M1) { # GMVAR type regime
        alphas[m]*dnorm(xx, mean=means[d1, m], sd=sqrt(vars[d1, m]))
      } else { # StMVAR type regime
        alphas[m]*my_dt(y=xx, mean=means[d1, m], var=vars[d1, m], df=all_df[m - M1])
      }
    }
    d1_dens_f <- function(d1, xx) rowSums(vapply(1:M, function(m) reg_dens(d1, m, xx), numeric(length(xx)))) # Marginal stat dens of d1:th series

    for(d1 in 1:d) { # Go through marginal series
      data_dens <- density(data[,d1], kernel="gaussian") # Kernel estimate
      mod_mean <- gsmvar$uncond_moments$uncond_mean[d1]
      mod_sd <- sqrt(gsmvar$uncond_moments$autocovs[d1, d1, 1])
      x0 <- min(mod_mean - 3*mod_sd, min(data_dens$x)) # Plot limits on horizontal axis
      x1 <- max(mod_mean + 3*mod_sd, max(data_dens$x))
      xpp <- seq(from=x0, to=x1, length.out=500)
      mod_dens <- d1_dens_f(d1=d1, xx=xpp) # Model implied marginal density
      y0 <- 0
      y1 <- max(c(data_dens$y, mod_dens))

      # Plot the densities
      plot(x=data_dens$x, y=data_dens$y, xlim=c(x0, x1), ylim=c(y0, y1), main=paste("Density:", names_ts[d1]),
           ylab="", xlab="", cex.axis=0.8, font.axis=2, type="l")
      lines(x=xpp, y=mod_dens, type="l", lty=2, lwd=2, col="darkgrey")
      for(m in 1:M) {
        lines(x=xpp, y=reg_dens(d1=d1, m=m, xx=xpp), type="l", lty=3, col=colpal_mw[m])
      }
      if(d1 == 1) {
        legend("topleft", legend=c("series", "model", names_mw), col=c("black", "darkgrey", colpal_mw),
               bty="n", lty=c(1, 2, 2, 2), lwd=2, text.font=2, cex=0.6, x.intersp=0.5, y.intersp=1)#,
               #seg.len=1, inset=c(-0.025, -0.01))
      }
    }
  }
}



#' @import graphics
#' @describeIn quantile_residual_tests Plot p-values of the autocorrelation and conditional
#'  heteroskedasticity tests.
#' @inheritParams print.qrtest
#' @export

plot.qrtest <- function(x, ...) {
  old_par <- par(no.readonly=TRUE)
  on.exit(par(old_par))
  qrtest <- x
  par(mfrow=c(1, 2), mar=c(5.1, 3.1, 3.1, 1.1))

  plot_pvalues <- function(which_ones) { # ac_res, ch_res
    res <- qrtest[[which(names(qrtest) == which_ones)]]
    pvals <- res$test_results$p_val
    seq_pvals <- seq_along(pvals)
    plot(pvals, ylim=c(0, 1), xlim=c(min(seq_pvals) - 0.2, max(seq_pvals) + 0.2), ylab="", xlab="lags",
         main=ifelse(which_ones == "ac_res", "Autocorrelation", "Cond. h.skedasticity"),
         xaxt="n", yaxt="n", pch=16, col="blue")
    axis(side=1, at=seq_pvals, labels=res$test_results$lags)
    levels <- c(0.01, 0.05, 0.10, seq(from=0.20, to=1.00, by=0.20))
    axis(side=2, at=levels, las=1, cex.axis=0.8)
    abline(h=0, lwd=2)
    abline(h=c(0.01, 0.05, 0.10, 1.00), lty=2, col=c(rep("red", 3), "green4"))
    segments(x0=seq_pvals, y0=0, y1=pvals, x1=seq_pvals, ...)
    points(pvals)
  }

  plot_pvalues("ac_res")
  plot_pvalues("ch_res")
}



#' @describeIn linear_IRF plot method
#' @inheritParams print.irf
#' @param shocks_to_plot IRFs of which shocks should be plotted? A numeric vector
#'   with elements in \code{1,...,d}.
#' @export

plot.irf <- function(x, shocks_to_plot, ...) {

  # Relevant statistics etc
  irf <- x
  point_est <- irf$point_est
  conf_ints <- irf$conf_ints
  d <- irf$gsmvar$model$d
  if(missing(shocks_to_plot)) {
    shocks_to_plot <- 1:d
  } else {
    stopifnot(all(shocks_to_plot %in% 1:irf$gsmvar$model$d))
  }
  var_names <- dimnames(point_est)[[1]]
  shock_names <- dimnames(point_est)[[2]]


  # Graphical settings
  old_par <- par(no.readonly=TRUE)
  on.exit(par(old_par))
  margs <- c(max(0.4, 0.4 + log(0.31 + log(d))/6), 0.3,
             max(0.35, 0.35 + log(0.31 + log(d))/6), 0.1)
  margs <- vapply(1:length(margs), function(i1) min(margs[i1], 1), numeric(1))
  par(las=1, mar=c(0, 0, 0, 0))
  nrows <- d + 2 # +2 for bottom and top marginals
  ncols <- 3 # +2 for the left and right margin; IRF of each shock in separate figure
  nfigs <- nrows*ncols
  layoutmat <- matrix(seq_len(nfigs), nrow=nrows, ncol=ncols, byrow=FALSE)
  layout(layoutmat,
         widths=c(margs[2], 1, margs[4]),
         heights=c(margs[3], rep(1, times=nrows - 2), margs[1]))

  # Function to plot empty plots (for the marginals)
  empty_plot <- function() plot(0, xaxt='n', yaxt='n', bty='n', pch='', ylab='', xlab='')

  # Function to plot the IRF for each variable separately (for a given shock)
  plot_irf <- function(var_ind, shock_ind, main="", xaxt="n", ylab="") {

    # Plot point estimate
    pe_var_shock <- point_est[var_ind, shock_ind, ]
    if(is.null(conf_ints)) {
      ylim <- c(min(0, min(pe_var_shock)), max(0, max(pe_var_shock)))
    } else {
      ci_var_shock <- conf_ints[var_ind, shock_ind, , ] # [horizon, bound]
      ylim <- c(min(0, min(cbind(ci_var_shock, pe_var_shock))),
                max(0, max(cbind(ci_var_shock, pe_var_shock))))
    }
    plot(x=0:(length(pe_var_shock) - 1), y=pe_var_shock, type="l", ylim=ylim,
         main="", ylab="", xlab="", xaxt=xaxt, lwd=2, col="black")
    mtext(var_names[var_ind], side=2, cex=0.8, font=2, las=0, padj=-4) # Add yaxis label to the first column of responses
    if(var_ind == 1) mtext(main, padj=-0.5, cex=1, font=2)

    # Plot confidence intervals
    if(!is.null(conf_ints)) {
      lines(x=0:(length(pe_var_shock) - 1), y=ci_var_shock[,1], lty=2)
      lines(x=0:(length(pe_var_shock) - 1), y=ci_var_shock[,2], lty=2)
    }

    abline(h=0, lty=3, col="lightgrey")
  }

  # Loop through the shocks
  for(i1 in shocks_to_plot) {
    if(i1 != shocks_to_plot[1]) {
      grDevices::devAskNewPage(TRUE)
    }

    # Plot a column of empty plots as the left margins
    for(i2 in 1:(nrows + 1)) { # + 1 for the top margin of the first row of responses
      empty_plot()
    }

    # Plot the responses of each variable to shock i1
    plot_irf(var_ind=1, shock_ind=i1, main=shock_names[i1], ylab=var_names[1])
    if(d > 2) {
      for(i2 in 2:(d - 1)) {
        plot_irf(var_ind=i2, shock_ind=i1, ylab=var_names[i2])
      }
    }
    plot_irf(var_ind=d, shock_ind=i1, xaxt="s", ylab=var_names[var])
    empty_plot() # To bottom margin of the last row of responses

    # Plot a column of empty plots as the right margins
    for(i2 in 1:nrows) {
      empty_plot()
    }
  }
}
