#' @title Function factory for value formatting
#'
#' @description \code{format_valuef} is a function factory for
#'   formatting values with certain number of digits.
#'
#' @param digits the number of decimals to print
#' @return Returns a function that takes an atomic vector as argument
#'   and returns it formatted to character with \code{digits} decimals.
#' @keywords internal

format_valuef <- function(digits) {
  function(x) tryCatch(format(round(x, digits), nsmall=digits), error=function(e) x)
}


#' @describeIn GSMVAR print method
#' @inheritParams plot.gsmvar
#' @param digits number of digits to be printed.
#' @param summary_print if set to \code{TRUE} then the print
#'   will include log-likelihood and information criteria values.
#' @export

print.gsmvar <- function(x, ..., digits=2, summary_print=FALSE) {
  gsmvar <- x
  stopifnot(digits >= 0 & digits%%1 == 0)
  format_value <- format_valuef(digits)
  p <- gsmvar$model$p
  M <- gsmvar$model$M
  d <- gsmvar$model$d
  model <- gsmvar$model$model
  IC <- gsmvar$IC
  constraints <- gsmvar$model$constraints
  same_means <- gsmvar$model$same_means
  structural_pars <- gsmvar$model$structural_pars
  all_mu <- round(get_regime_means(gsmvar), digits)
  params <- gsmvar$params
  npars <- length(params)
  T_obs <- ifelse(is.null(gsmvar$data), NA, nrow(gsmvar$data))
  params <- reform_constrained_pars(p=p, M=M, d=d, params=params, model=model,
                                    constraints=constraints, same_means=same_means,
                                    structural_pars=structural_pars)
  if(gsmvar$model$parametrization == "mean") {
    params <- change_parametrization(p=p, M=M, d=d, params=params, model=model,
                                     constraints=NULL, same_means=NULL,
                                     structural_pars=structural_pars, change_to="intercept")
  }
  structural_pars <- get_unconstrained_structural_pars(structural_pars=structural_pars)
  all_phi0 <- pick_phi0(p=p, M=M, d=d, params=params, structural_pars=structural_pars)
  all_A <- pick_allA(p=p, M=M, d=d, params=params, structural_pars=structural_pars)
  all_Omega <- pick_Omegas(p=p, M=M, d=d, params=params, structural_pars=structural_pars)
  alphas <- pick_alphas(p=p, M=M, d=d, params=params, model=model)
  all_df <- pick_df(M=M, params=params, model=model)
  cat(ifelse(is.null(structural_pars), "Reduced form", "Structural"), model, "model:\n")
  cat(paste0(" p = ", p, ", "))
  if(model == "G-StMVAR") {
    cat(paste0("M1 = ", M[1], ", M2 = ", M[2], ", "))
  } else { # model == "GMVAR" or "StMVAR"
    cat(paste0("M = ", M, ", "))
  }
  cat(paste0("d = ", d, ", #parameters = " , npars, ","),
      ifelse(is.na(T_obs), "\n", paste0("#observations = ", T_obs, " x ", d, ",\n")),
      ifelse(gsmvar$model$conditional, "conditional", "exact"),
      "log-likelihood,",
      ifelse(gsmvar$model$parametrization == "mean", "mean parametrization,", "intercept parametrization,"),
      ifelse(is.null(same_means),
             ifelse(is.null(constraints), "no AR parameter constraints", "AR parameters constrained"),
             ifelse(is.null(constraints), "mean paremeters constrained, no AR parameter constraints",
                                          "mean parameters constrained, AR parameters constrained")), "\n")
  cat("\n")

  if(summary_print) {
    all_boldA_eigens <- get_boldA_eigens(gsmvar)
    all_omega_eigens <- get_omega_eigens(gsmvar)
    form_val2 <- function(txt, val) paste(txt, format_value(val))
    cat(paste(form_val2(" log-likelihood:", gsmvar$loglik),
                    form_val2("AIC:", IC$AIC),
                    form_val2("HQIC:", IC$HQIC),
                    form_val2("BIC:", IC$BIC),
                    sep=", "), "\n\n")
  }

  plus <- c("+", rep(" ", times=d-1))
  arch_scalar <- c(rep(" ", times=d-1), "ARCH_mt")
  round_lbrackets <- rep("(", times=d)
  round_rbrackets <- rep(")", times=d)
  Y <- paste0("y", 1:d)
  tmp_names <- paste0("tmp", 1:(p*(d + 2) + d + 2))

  for(m in seq_len(sum(M))) {
    count <- 1
    if(model == "GMVAR") {
      regime_type <- "GMVAR"
    } else if(model == "StMVAR") {
      regime_type <- "StMVAR"
      M1 <- 0
    } else {
      M1 <- M[1]
      regime_type <- ifelse(m <= M1, "GMVAR", "StMVAR")
    }
    cat(paste("Regime", m))
    if(model == "G-StMVAR") cat(paste0(" (", regime_type, " type)"))
    cat("\n")
    if(summary_print) {
      cat(paste("Moduli of 'bold A' eigenvalues: ", paste0(format_value(all_boldA_eigens[,m]), collapse=", ")),"\n")
      cat(paste("Cov. matrix 'Omega' eigenvalues:", paste0(format_value(all_omega_eigens[,m]), collapse=", ")),"\n")
    }
    cat(paste("Mixing weight:", format_value(alphas[m])), "\n")
    cat("Regime means:", paste0(format_value(all_mu[,m]), collapse=", "))
    if(regime_type == "StMVAR") { # Print degrees of freedom parameter for StMVAR type regimes
      cat("\nDf parameter: ", format_value(all_df[m - M1]))
    }
    cat("\n\n")

    left_brackets <- rep("[", times=d)
    right_brackets <- rep("]", times=d)
    df <- data.frame(Y=Y,
                     eq=c("=", rep(" ", d - 1)),
                     eq=left_brackets,
                     phi0=format_value(all_phi0[, m, drop=FALSE]),
                     eq=rep("]", times=d),
                     plus)
    for(i1 in seq_len(p)) {
      Amp_colnames <- c(paste0("A", i1), tmp_names[count:(count + d - 1 - 1)]); count <- count + d - 1
      df[, tmp_names[count]] <- left_brackets; count <- count + 1
      df[, Amp_colnames] <- format_value(all_A[, ,i1 , m])
      df[, tmp_names[count]] <- rep("]", times=d); count <- count + 1
      df[, tmp_names[count]] <- paste0(Y, ".", i1); count <- count + 1
      df <- cbind(df, plus)
    }
    if(regime_type == "StMVAR") { # Time varying ARCH scalar multiplying the constant part of error term covariance matrix
      df <- cbind(df, round_lbrackets, arch_scalar)
    }
    df[, tmp_names[p*(d + 2) + 1]] <- left_brackets
    df[, c("Omega", tmp_names[(p*(d + 2) + 2):(p*(d + 2) + d)])] <- format_value(all_Omega[, , m])
    df[, tmp_names[p*(d + 2) + d + 1]] <- rep("]", times=d)
    if(regime_type == "StMVAR") {
      df <- cbind(df, round_rbrackets)
    }
    df[, "1/2"] <- rep(" ", d)
    df[, tmp_names[p*(d + 2) + d + 2]] <- paste0("eps", 1:d)
    names_to_omit <- unlist(lapply(c("plus", "eq", "arch_scalar", "round_lbrackets", "round_rbrackets", tmp_names),
                                   function(nam) grep(nam, colnames(df))))
    colnames(df)[names_to_omit] <- " "
    print(df)
    cat("\n")
    if(summary_print) {
      cat("Error term correlation matrix:\n")
      print(cov2cor(all_Omega[, , m]), digits=digits)
      cat("\n")
    }
  }
  if(!is.null(structural_pars)) {
    cat("Structural parameters:\n")
    W <- format_value(pick_W(p=p, M=M, d=d, params=params, structural_pars=structural_pars))

    tmp <- c(rep(" ", times=d - 1), ",")
    df2 <- data.frame(left_brackets, W=W[,1])
    for(i1 in 2:d) {
      df2 <- cbind(df2, W[, i1])
      colnames(df2)[1 + i1] <- "tmp"
    }
    df2 <- cbind(df2, right_brackets)
    if(sum(M) > 1) {
      lambdas <- format_value(pick_lambdas(p=p, M=M, d=d, params=params, structural_pars=structural_pars))
      tmp <- c(rep(" ", times=d - 1), ",")
      lambdas <- matrix(lambdas, nrow=d, ncol=sum(M) - 1, byrow=FALSE) # Column for each regime
      for(i1 in 1:(sum(M) - 1)) {
        lmb <- lambdas[,i1]
        df2 <- cbind(df2, tmp, left_brackets, lmb, right_brackets)
        colnames(df2)[grep("lmb", colnames(df2))] <- paste0("lamb", i1 + 1)
      }
    }
    names_to_omit <- unlist(lapply(c("left_brackets", "right_brackets", "tmp"), function(nam) grep(nam, colnames(df2))))
    colnames(df2)[names_to_omit] <- " "
    print(df2)
    cat("\n")
    W_orig <- gsmvar$model$structural_pars$W
    n_zero <- sum(W_orig == 0, na.rm=TRUE)
    n_free <- sum(is.na(W_orig))
    n_sign <- d^2 - n_zero - n_free
    cat("The B-matrix (or equally W) is subject to", n_zero, "zero constraints and", n_sign, "sign constraints.\n")
    cat("The eigenvalues lambda_{mi} are", ifelse(is.null(gsmvar$model$structural_pars$C_lambda), "not subject to linear constraints.",
                                                  "subject to linear constraints."))
    cat("\n")
  }

  if(summary_print) {
    cat("Print approximate standard errors with the function 'print_std_errors'.\n")
  }
  invisible(gsmvar)
}


#' @title Summary print method from objects of class 'gsmvarsum'
#'
#' @description \code{print.gsmvarsum} is a print method for object \code{'gsmvarsum'} generated
#'   by \code{summary.gsmvar}.
#'
#' @param x object of class 'gsmvarsum' generated by \code{summary.gsmvar}.
#' @param ... currently not used.
#' @param digits the number of digits to be printed.
#' @examples
#' # GMVAR(2, 2), d=2 model;
#' params22 <- c(0.36, 0.121, 0.223, 0.059, -0.151, 0.395, 0.406, -0.005,
#'  0.083, 0.299, 0.215, 0.002, 0.03, 0.484, 0.072, 0.218, 0.02, -0.119,
#'  0.722, 0.093, 0.032, 0.044, 0.191, 1.101, -0.004, 0.105, 0.58)
#' mod22 <- GSMVAR(gdpdef, p=2, M=2, params=params22)
#' sumry22 <- summary(mod22)
#' print(sumry22)
#' @export

print.gsmvarsum <- function(x, ..., digits) {
  gsmvarsum <- x
  if(missing(digits)) digits <- gsmvarsum$digits
  print.gsmvar(gsmvarsum$gsmvar, ..., digits=digits, summary_print=TRUE)
  invisible(gsmvarsum)
}


#' @title Print method for class 'gsmvarpred' objects
#'
#' @description \code{print.gsmvarpred} is a print method for object generated
#'  by \code{predict.gsmvar}.
#'
#' @inheritParams plot.gsmvarpred
#' @param digits the number of decimals to print
#' @param ... currently not used.
#' @examples
#' # GMVAR(2, 2), d=2 model;
#' params22 <- c(0.36, 0.121, 0.223, 0.059, -0.151, 0.395, 0.406, -0.005,
#'  0.083, 0.299, 0.215, 0.002, 0.03, 0.484, 0.072, 0.218, 0.02, -0.119,
#'  0.722, 0.093, 0.032, 0.044, 0.191, 1.101, -0.004, 0.105, 0.58)
#' mod22 <- GSMVAR(gdpdef, p=2, M=2, params=params22)
#' pred22 <- predict(mod22, n_ahead=3, plot_res=FALSE)
#' print(pred22)
#' print(pred22, digits=3)
#' @export

print.gsmvarpred <- function(x, ..., digits=2) {
  gsmvarpred <- x
  stopifnot(digits >= 0 & digits%%1 == 0)
  format_value <- format_valuef(digits)

  if(gsmvarpred$pred_type == "cond_mean") {
    cat("One-step-ahead forecast by exact conditional mean, no prediction intervals.\n")
    cat("Forecast:", paste0(format_value(gsmvarpred$pred), collapse=", "), "\n")

  } else if(gsmvarpred$pi_type == "none") {
    cat(paste0("Point forecast by ", gsmvarpred$pred_type, ", no prediction intervals."), "\n")
    cat(paste0("Forecast ", gsmvarpred$n_ahead, " steps ahead, based on ", gsmvarpred$nsim, " Monte Carlo repetitions.\n"))
    print(gsmvarpred$pred)

  } else {
    cat(paste0("Point forecast by ", gsmvarpred$pred_type, ", ", gsmvarpred$pi_type,
               " prediction intervals with levels ", paste(gsmvarpred$pi, collapse=", "), "."), "\n")
    cat(paste0("Forecast ", gsmvarpred$n_ahead, " steps ahead, based on ", gsmvarpred$nsim, " Monte Carlo repetitions.\n"))

    cat("\n")
    q <- gsmvarpred$q
    pred_ints <- gsmvarpred$pred_ints
    pred <- gsmvarpred$pred
    pred_type <- gsmvarpred$pred_type
    series_names <- colnames(gsmvarpred$pred)
    for(i1 in seq_len(gsmvarpred$gsmvar$model$d)) {
      cat(paste0(series_names[i1], ":"), "\n")
      df <- as.data.frame(lapply(1:length(gsmvarpred$q), function(i2) format_value(pred_ints[, i2, i1])))
      names(df) <- q
      df[, pred_type] <- format_value(pred[,i1])
      if(gsmvarpred$pi_type == "two-sided") {
        new_order <- as.character(c(q[1:(length(q)/2)], pred_type, q[(length(q)/2 + 1):length(q)]))
      } else if(gsmvarpred$pi_type == "upper") {
        new_order <- as.character(c(pred_type, q))
      } else {
        new_order <- names(df)
      }
      print(df[, new_order])
      cat("\n")
    }
    if(gsmvarpred$pred_type != "cond_mean") {
      cat("Point forecasts and prediction intervals for mixing weights can be obtained with $mix_pred and $mix_pred_ints, respectively.\n")
    }
  }
  invisible(gsmvarpred)
}


#' @describeIn quantile_residual_tests Print method for class 'qrtest'
#' @inheritParams print.gsmvarpred
#' @param x object of class \code{'qrtest'} generated by the function \code{quantile_residual_tests)}.
#' @param ... currently not used.
#' @export

print.qrtest <- function(x, ..., digits=3) {
  qrtest <- x
  format_value <- format_valuef(digits)
  format_lag <- format_valuef(0)
  cat(paste("Normality test p_value:", format_value(qrtest$norm_res$p_val)), "\n\n")

  cat("Autocorrelation tests:\nlags | p_value\n")
  for(i1 in seq_along(qrtest$ac_res$test_results$lags)) {
    if(qrtest$ac_res$test_results$lags[i1] < 10) {
      cat(" ", format_lag(qrtest$ac_res$test_results$lags[i1]), " | ", format_value(qrtest$ac_res$test_results$p_val[i1]), "\n")
    } else {
      cat(" ", format_lag(qrtest$ac_res$test_results$lags[i1]), "| ", format_value(qrtest$ac_res$test_results$p_val[i1]), "\n")
    }
  }
  cat("\nConditional hetetoskedasticity tests:\nlags | p_value\n")
  for(i1 in seq_along(qrtest$ch_res$test_results$lags)) {
    if(qrtest$ch_res$test_results$lags[i1] < 10) {
      cat(" ", format_lag(qrtest$ch_res$test_results$lags[i1]), " | ", format_value(qrtest$ch_res$test_results$p_val[i1]), "\n")
    } else {
      cat(" ", format_lag(qrtest$ch_res$test_results$lags[i1]), "| ", format_value(qrtest$ch_res$test_results$p_val[i1]), "\n")
    }
  }
  invisible(qrtest)
}


#' @describeIn GIRF print method
#' @inheritParams print.gsmvarpred
#' @param x object of class \code{'girf'} generated by the function \code{GIRF}.
#' @param N_to_print an integer specifying the horizon how far to print the estimates and
#'   confidence intervals. The default is that all the values are printed.
#' @export

print.girf <- function(x, ..., digits=2, N_to_print) {
  girf <- x
  girf_res <- girf$girf_res
  stopifnot(digits >= 0 & digits%%1 == 0)
  format_value <- format_valuef(digits)
  if(missing(N_to_print)) {
    N_to_print <- nrow(girf_res[[1]]$point_est)
  } else {
    stopifnot(N_to_print %in% 1:nrow(girf_res[[1]]$point_est))
  }
  if(length(girf$which_cumulative) > 0) {
    cat(paste0("The responses of the variables ",
               paste0(dimnames(girf_res[[1]]$point_est)[[2]][girf$which_cumulative], collapse=", "),
               " were cumulated."), "\n\n")
  }

  for(i1 in 1:length(girf_res)) {
    if(i1 > 1) cat("------------------------\n")
    cat(paste0("The GIRF of shock ", girf$shocks[i1], ":"), "\n")
    girf_i1 <- girf_res[[i1]]
    for(i2 in 1:dim(girf_i1$conf_ints)[3]) {
      cat(paste0("The response of ", dimnames(girf_i1$conf_ints)[[3]][i2], ":"), "\n")
      df <- as.data.frame(lapply(1:ncol(girf_i1$conf_ints[, , i2]), function(i3) format_value(girf_i1$conf_ints[, i3, i2])))
      q <- dimnames(girf_i1$conf_ints)[[2]]
      names(df) <- q

      df[, "mean"] <- format_value(girf_i1$point_est[, i2])
      new_order <- as.character(c(q[1:(length(q)/2)], "mean", q[(length(q)/2 + 1):length(q)]))
      print(utils::head(df[, new_order], n=N_to_print + 1))
      cat("\n")
    }
  }
  invisible(girf)
}



#' @describeIn GFEVD print method
#' @inheritParams print.gsmvarpred
#' @param x object of class \code{'gfevd'} generated by the function \code{GFEVD}.
#' @param N_to_print an integer specifying the horizon how far to print the estimates.
#'   The default is that all the values are printed.
#' @export

print.gfevd <- function(x, ..., digits=2, N_to_print) {
  gfevd <- x
  gfevd_res <- gfevd$gfevd_res
  stopifnot(digits >= 0 & digits%%1 == 0)
  format_value <- format_valuef(digits)
  if(missing(N_to_print)) {
    N_to_print <- nrow(gfevd_res[, , 1]) - 1
  } else {
    stopifnot(N_to_print %in% 1:nrow(gfevd_res[, , 1]))
  }
  if(length(gfevd$which_cumulative) > 0) {
    cat(paste0("The responses of the variables ",
               paste0(dimnames(gfevd_res)[[3]][gfevd$which_cumulative], collapse=", "),
               " were cumulated."), "\n\n")
  }

  for(i1 in 1:dim(gfevd_res)[3]) { # Go through GFEVDs of each variable and possibly mixing weights
    if(i1 > 1) cat("------------------------\n")
    cat(paste0("The GFEVD for ", dimnames(gfevd_res)[[3]][i1], ":"), "\n")
    print(round(gfevd_res[1:(N_to_print  + 1), , i1], digits=digits))
    cat("\n")
  }
  invisible(gfevd)
}


#' @title Print method for the class hypotest
#'
#' @description \code{print.hypotest} is the print method for the class hypotest
#'  objects.
#' @param digits how many significant digits to print?
#' @param x object of class \code{'hypotest'} generated by the function \code{Wald_test} or \code{LR_test}.
#' @param ... currently not in use.
#' @export

print.hypotest <- function(x, ..., digits=4) {
  stopifnot(digits >= 0 & digits%%1 == 0)
  format_value <- function(a) format(a, digits=digits)

  cat(paste0(x$type, ":"), "\n",
      paste0("test stat = ", format_value(x$test_stat),
             ", df = ", x$df,
             ", p-value = ", format_value(x$p_value)))
  invisible(x)
}
