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


#' @describeIn GMVAR print method
#' @inheritParams plot.gmvar
#' @param digits number of digits to be printed.
#' @param summary_print if set to \code{TRUE} then the print
#'   will include log-likelihood and information criteria values.
#' @export

print.gmvar <- function(x, ..., digits=2, summary_print=FALSE) {
  gmvar <- x
  stopifnot(digits >= 0 & digits%%1 == 0)
  format_value <- format_valuef(digits)
  p <- gmvar$model$p
  M <- gmvar$model$M
  d <- gmvar$model$d
  IC <- gmvar$IC
  constraints <- gmvar$model$constraints
  same_means <- gmvar$model$same_means
  structural_pars <- gmvar$model$structural_pars
  all_mu <- round(get_regime_means(gmvar), digits)
  params <- gmvar$params
  npars <- length(params)
  T_obs <- ifelse(is.null(gmvar$data), NA, nrow(gmvar$data))
  params <- reform_constrained_pars(p=p, M=M, d=d, params=params, constraints=constraints,
                                   same_means=same_means, structural_pars=structural_pars)
  if(gmvar$model$parametrization == "mean") {
    params <- change_parametrization(p=p, M=M, d=d, params=params, constraints=NULL,
                                     same_means=NULL, structural_pars=structural_pars,
                                     change_to="intercept")
  }
  structural_pars <- get_unconstrained_structural_pars(structural_pars=structural_pars)
  all_phi0 <- pick_phi0(p=p, M=M, d=d, params=params, structural_pars=structural_pars)
  all_A <- pick_allA(p=p, M=M, d=d, params=params, structural_pars=structural_pars)
  all_Omega <- pick_Omegas(p=p, M=M, d=d, params=params, structural_pars=structural_pars)
  alphas <- pick_alphas(p=p, M=M, d=d, params=params)
  cat(ifelse(is.null(structural_pars), "Reduced form", "Structural"), "model:\n")
  cat(paste0(" p = ", p, ", M = ", M, ", d = ", d, ","),
      paste0("#parameters = " , npars, ","),
      ifelse(is.na(T_obs), "\n", paste0("#observations = ", T_obs, " x ", d, ",\n")),
      ifelse(gmvar$model$conditional, "conditional", "exact"),
      "log-likelihood,",
      ifelse(gmvar$model$parametrization == "mean", "mean parametrization,", "intercept parametrization,"),
      ifelse(is.null(same_means),
             ifelse(is.null(constraints), "no AR parameter constraints", "AR parameters constrained"),
             ifelse(is.null(constraints), "mean paremeters constrained, no AR parameter constraints",
                                          "mean parameters constrained, AR parameters constrained")), "\n")
  cat("\n")

  if(summary_print) {
    all_boldA_eigens <- get_boldA_eigens(gmvar)
    all_omega_eigens <- get_omega_eigens(gmvar)
    form_val2 <- function(txt, val) paste(txt, format_value(val))
    cat(paste(form_val2(" log-likelihood:", gmvar$loglik),
                    form_val2("AIC:", IC$AIC),
                    form_val2("HQIC:", IC$HQIC),
                    form_val2("BIC:", IC$BIC),
                    sep=", "), "\n\n")
  }

  plus <- c("+", rep(" ", d-1))
  Y <- paste0("y", 1:d)
  tmp_names <- paste0("tmp", 1:(p*(d + 2) + d + 2))

  for(m in seq_len(M)) {
    count <- 1
    cat(paste("Regime", m), "\n")
    if(summary_print) {
      cat(paste("Moduli of 'bold A' eigenvalues: ", paste0(format_value(all_boldA_eigens[,m]), collapse=", ")),"\n")
      cat(paste("Cov. matrix 'Omega' eigenvalues:", paste0(format_value(all_omega_eigens[,m]), collapse=", ")),"\n")
    }
    cat(paste("Mixing weight:", format_value(alphas[m])), "\n")
    cat("Regime means:", paste0(format_value(all_mu[,m]), collapse=", "), "\n\n")
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
    df[, tmp_names[p*(d + 2) + 1]] <- left_brackets
    df[, c("Omega", tmp_names[(p*(d + 2) + 2):(p*(d + 2) + d)])] <- format_value(all_Omega[, , m])
    df[, tmp_names[p*(d + 2) + d + 1]] <- rep("]", times=d)
    df[, "1/2"] <- rep(" ", d)
    df[, tmp_names[p*(d + 2) + d + 2]] <- paste0("eps", 1:d)
    names_to_omit <- unlist(lapply(c("plus", "eq", tmp_names), function(nam) grep(nam, colnames(df))))
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
    if(M > 1) {
      lambdas <- format_value(pick_lambdas(p=p, M=M, d=d, params=params, structural_pars=structural_pars))
      tmp <- c(rep(" ", times=d - 1), ",")
      lambdas <- matrix(lambdas, nrow=d, ncol=M - 1, byrow=FALSE) # Column for each regime
      for(i1 in 1:(M - 1)) {
        lmb <- lambdas[,i1]
        df2 <- cbind(df2, tmp, left_brackets, lmb, right_brackets)
        colnames(df2)[grep("lmb", colnames(df2))] <- paste0("lamb", i1 + 1)
      }
    }
    names_to_omit <- unlist(lapply(c("left_brackets", "right_brackets", "tmp"), function(nam) grep(nam, colnames(df2))))
    colnames(df2)[names_to_omit] <- " "
    print(df2)
    cat("\n")
    W_orig <- gmvar$model$structural_pars$W
    n_zero <- sum(W_orig == 0, na.rm=TRUE)
    n_free <- sum(is.na(W_orig))
    n_sign <- d^2 - n_zero - n_free
    cat("The B-matrix (or equally W) is subject to", n_zero, "zero constraints and", n_sign, "sign constraints.\n")
    cat("The eigenvalues lambda_{mi} are", ifelse(is.null(gmvar$model$structural_pars$C_lambda), "not subject to linear constraints.",
                                                  "subject to linear constraints."))
    cat("\n")
  }

  if(summary_print) {
    cat("Print approximate standard errors with the function 'print_std_errors'.\n")
  }
  invisible(gmvar)
}


#' @title Summary print method from objects of class 'gmvarsum'
#'
#' @description \code{print.gmvarsum} is a print method for object \code{'gmvarsum'} generated
#'   by \code{summary.gmvar}.
#'
#' @param x object of class 'gmvarsum' generated by \code{summary.gmvar}.
#' @param ... currently not used.
#' @param digits the number of digits to be printed.
#' @examples
#' # GMVAR(2, 2), d=2 model;
#' params22 <- c(0.36, 0.121, 0.223, 0.059, -0.151, 0.395, 0.406, -0.005,
#'  0.083, 0.299, 0.215, 0.002, 0.03, 0.484, 0.072, 0.218, 0.02, -0.119,
#'  0.722, 0.093, 0.032, 0.044, 0.191, 1.101, -0.004, 0.105, 0.58)
#' mod22 <- GMVAR(gdpdef, p=2, M=2, params=params22)
#' sumry22 <- summary(mod22)
#' print(sumry22)
#' @export

print.gmvarsum <- function(x, ..., digits) {
  gmvarsum <- x
  if(missing(digits)) digits <- gmvarsum$digits
  print.gmvar(gmvarsum$gmvar, ..., digits=digits, summary_print=TRUE)
  if(!is.null(gmvarsum$qrtest)) {
    cat("_____________________________________\n")
    cat("Quantile residual tests based on data\n\n")
    print.qrtest(gmvarsum$qrtest)
  }
  invisible(gmvarsum)
}


#' @title Print method for class 'gmvarpred' objects
#'
#' @description \code{print.gmvarpred} is a print method for object generated
#'  by \code{predict.gmvar}.
#'
#' @inheritParams plot.gmvarpred
#' @param digits the number of decimals to print
#' @param ... currently not used.
#' @examples
#' # GMVAR(2, 2), d=2 model;
#' params22 <- c(0.36, 0.121, 0.223, 0.059, -0.151, 0.395, 0.406, -0.005,
#'  0.083, 0.299, 0.215, 0.002, 0.03, 0.484, 0.072, 0.218, 0.02, -0.119,
#'  0.722, 0.093, 0.032, 0.044, 0.191, 1.101, -0.004, 0.105, 0.58)
#' mod22 <- GMVAR(gdpdef, p=2, M=2, params=params22)
#' pred22 <- predict(mod22, n_ahead=3, plot_res=FALSE)
#' print(pred22)
#' print(pred22, digits=3)
#' @export

print.gmvarpred <- function(x, ..., digits=2) {
  gmvarpred <- x
  stopifnot(digits >= 0 & digits%%1 == 0)
  format_value <- format_valuef(digits)

  if(gmvarpred$pred_type == "cond_mean") {
    cat("One-step-ahead forecast by exact conditional mean, no prediction intervals.\n")
    cat("Forecast:", paste0(format_value(gmvarpred$pred), collapse=", "), "\n")

  } else if(gmvarpred$pi_type == "none") {
    cat(paste0("Point forecast by ", gmvarpred$pred_type, ", no prediction intervals."), "\n")
    cat(paste0("Forecast ", gmvarpred$n_ahead, " steps ahead, based on ", gmvarpred$n_simu, " simulations.\n"))
    print(gmvarpred$pred)

  } else {
    cat(paste0("Point forecast by ", gmvarpred$pred_type, ", ", gmvarpred$pi_type,
               " prediction intervals with levels ", paste(gmvarpred$pi, collapse=", "), "."), "\n")
    cat(paste0("Forecast ", gmvarpred$n_ahead, " steps ahead, based on ", gmvarpred$n_simu, " simulations.\n"))

    cat("\n")
    q <- gmvarpred$q
    pred_ints <- gmvarpred$pred_ints
    pred <- gmvarpred$pred
    pred_type <- gmvarpred$pred_type
    for(i1 in seq_len(gmvarpred$gmvar$model$d)) {
      cat(paste0("Component ", i1, ":"), "\n")
      df <- as.data.frame(lapply(1:length(gmvarpred$q), function(i2) format_value(pred_ints[, i2, i1])))
      names(df) <- q
      df[, pred_type] <- format_value(pred[,i1])
      if(gmvarpred$pi_type == "two-sided") {
        new_order <- as.character(c(q[1:(length(q)/2)], pred_type, q[(length(q)/2 + 1):length(q)]))
      } else if(gmvarpred$pi_type == "upper") {
        new_order <- as.character(c(pred_type, q))
      } else {
        new_order <- names(df)
      }
      print(df[, new_order])
      cat("\n")
    }
    if(gmvarpred$pred_type != "cond_mean") {
      cat("Point forecasts and prediction intervals for mixing weights can be obtained with $mix_pred and $mix_pred_ints, respectively.\n")
    }
  }
  invisible(gmvarpred)
}


#' @describeIn quantile_residual_tests Print method for class 'qrtest'
#' @inheritParams print.gmvarpred
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
#' @inheritParams print.gmvarpred
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
#' @inheritParams print.gmvarpred
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

