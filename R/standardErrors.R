
#' @title Calculate standard errors for estimates of a GMVAR, StMVAR, or G-StMVAR model
#'
#' @description \code{standard_errors} calculates approximate standard errors for the GMVAR,
#'   StMVAR, or G-StMVAR model using square roots of the diagonal of inverse of observed information matrix
#'   and central-difference approximation for the differentiation.
#'
#' @inheritParams loglikelihood_int
#' @param custom_h A numeric vector with same the length as the parameter vector: i:th element of custom_h is the difference
#'  used in central difference approximation for partial differentials of the log-likelihood function for the i:th parameter.
#'  If \code{NULL} (default), then the difference used for differentiating overly large degrees of freedom parameters
#'  is adjusted to avoid numerical problems, and the difference is \code{6e-6} for the other parameters.
#' @return A vector containing the approximate standard errors of the estimates.
#' @inherit in_paramspace_int references
#' @keywords internal

standard_errors <- function(data, p, M, params, model=c("GMVAR", "StMVAR", "G-StMVAR"), conditional=TRUE,
                            parametrization=c("intercept", "mean"), constraints=NULL, same_means=NULL,
                            weight_constraints=NULL, structural_pars=NULL, minval, custom_h=NULL,
                            stat_tol=1e-3, posdef_tol=1e-8, df_tol=1e-8) {
  model <- match.arg(model)
  parametrization <- match.arg(parametrization)

  # The log-likelihood function to differentiate
  loglik_fn <- function(params) {
    tryCatch(loglikelihood_int(data=data, p=p, M=M, params=params, model=model, conditional=conditional,
                               parametrization=parametrization, constraints=constraints, same_means=same_means,
                               weight_constraints=weight_constraints, structural_pars=structural_pars,
                               check_params=TRUE, to_return="loglik", minval=minval,
                               stat_tol=stat_tol, posdef_tol=posdef_tol, df_tol=df_tol),
             error=function(e) NA)
  }

  # Calculate Hessian
  if(is.null(custom_h)) { # Adjust h for overly large degrees of freedom parameters
    varying_h <- get_varying_h(M=M, params=params, model=model)
  } else { # Utilize user-specified h
    stopifnot(length(custom_h) == length(params))
    varying_h <- custom_h
  }
  Hess <- calc_hessian(x=params, fn=loglik_fn, varying_h=varying_h)

  # Inverse of the observed information matrix
  inv_obs_inf <- tryCatch(solve(-Hess), error=function(e) matrix(NA, nrow=length(params), ncol=length(params)))

  # Calculate the standard errors
  unlist(lapply(diag(inv_obs_inf), function(x) ifelse(is.na(x) | x < 0, NA, sqrt(x))))
}


#' @title Print standard errors of a GMVAR, StMVAR, or G-StMVAR model in the same form as the model estimates are printed
#'
#' @description \code{print_std_errors} prints the approximate standard errors of a GMVAR, StMVAR, or G-StMVAR model in the
#'   same form as the parameters of objects of class \code{'gsmvar'} are printed.
#'
#' @inheritParams quantile_residual_tests
#' @param digits how many digits should be printed?
#' @details The main purpose of \code{print_std_errors} is to provide a convenient tool to match the standard
#'   errors to certain parameter estimates. Note that if the model is intercept parametrized, there won't
#'   be standard errors for the unconditional means, and vice versa. Also, there is no standard error for the
#'   last mixing weight alpha_M because it is not parametrized.
#'
#'   Note that if linear constraints are imposed and they involve summations or multiplications, then the AR
#'   parameter standard errors are printed separately as they don't correspond one-to-one to the model parameter
#'   standard errors.
#' @seealso \code{\link{profile_logliks}}, \code{\link{fitGSMVAR}}, \code{\link{GSMVAR}}, \code{\link{print.gsmvar}},
#'  \code{\link{swap_parametrization}}
#' @inherit GSMVAR references
#' @examples
#' \donttest{
#' # GMVAR(1,2) model
#' fit12 <- fitGSMVAR(gdpdef, p=1, M=2, ncalls=1, seeds=1)
#' fit12
#' print_std_errors(fit12)
#' }
#' @export

print_std_errors <- function(gsmvar, digits=3) {
  gsmvar <- gmvar_to_gsmvar(gsmvar) # Backward compatibility
  if(!all_pos_ints(digits)) stop("Argument digits must be positive integer")
  format_value <- format_valuef(digits)
  p <- gsmvar$model$p
  M <- gsmvar$model$M
  d <- gsmvar$model$d
  model <- gsmvar$model$model
  npars <- length(gsmvar$params)
  T_obs <- ifelse(is.null(gsmvar$data), NA, nrow(gsmvar$data))
  constraints <- gsmvar$model$constraints
  parametrization <- gsmvar$model$parametrization
  same_means <- gsmvar$model$same_means
  weight_constraints <- gsmvar$model$weight_constraints
  pars <- reform_constrained_pars(p=p, M=M, d=d, params=gsmvar$std_errors, model=model, constraints=constraints,
                                  same_means=gsmvar$model$same_means, weight_constraints=weight_constraints,
                                  structural_pars=gsmvar$model$structural_pars, change_na=TRUE)
  structural_pars <- get_unconstrained_structural_pars(structural_pars=gsmvar$model$structural_pars)
  all_phi0_or_mu <- pick_phi0(p=p, M=M, d=d, params=pars, structural_pars=structural_pars)
  all_A <- pick_allA(p=p, M=M, d=d, params=pars, structural_pars=structural_pars)
  if(is.null(structural_pars)) {
    all_Omega <- pick_Omegas(p=p, M=M, d=d, params=pars, structural_pars=structural_pars)
  } else {
    # No standard errors for cov. mats. as the model is parametrized with W and lambdas
    all_Omega <- array(" ", dim=c(d, d, sum(M)))
  }
  alphas <- pick_alphas(p=p, M=M, d=d, params=pars, model=model)
  all_df <- pick_df(M=M, params=pars, model=model)
  M_orig <- M
  M <- sum(M)
  alphas[M] <- NA # No standard error for the last alpha (it is not displayed anyway, though)
  if(!is.null(weight_constraints)) {
    alphas <- rep(NA, times=M)
  }

  if(parametrization == "mean") {
    all_mu <- all_phi0_or_mu
    all_phi0 <- matrix(" ", nrow=d, ncol=M)
  } else {
    all_mu <- matrix(NA, nrow=d, ncol=M)
    all_phi0 <- all_phi0_or_mu
  }
  if(!is.null(constraints)) {
    # The constrained AR parameter standard errors multiplied open in 'pars' are valid iff
    # the constraint matrix contains zeros and ones only, and there is at most one one in
    # each row (no multiplications or summations).
    if(any(constraints != 1 & constraints != 0) | any(rowSums(constraints) > 1)) {
      sep_AR <- TRUE # The AR parameter std errors must be printed separately
      all_A <- array(" ", dim=c(d, d, p, M))
      AR_stds <- gsmvar$std_errors[(M*d + 1):(M*d + ncol(constraints))] # Constrained AR param std errors
    } else {
      sep_AR <- FALSE
    }
  } else {
    sep_AR <- FALSE # No constraints imposed
  }

  cat(ifelse(is.null(structural_pars), "Reduced form", "Structural"), model, "model:\n")
  cat(paste0(" p = ", p, ", "))
  if(model == "G-StMVAR") {
    cat(paste0("M1 = ", M[1], ", M2 = ", M[2], ", "))
  } else { # model == "GMVAR" or "StMVAR"
    cat(paste0("M = ", M, ", "))
  }
  cat(paste0("d = ", d, ", #parameters = " , npars, ","),
      ifelse(is.na(T_obs), "\n", paste0("#observations = ", T_obs, " x ", d, ",\n")),
      ifelse(gsmvar$model$conditional, "conditional", "exact"), "log-likelihood,",
      paste0(ifelse(gsmvar$model$parametrization == "mean", "mean parametrization", "intercept parametrization"),
             ifelse(is.null(same_means), "", ", mean parameters constrained"),
             ifelse(is.null(constraints), "", ", AR matrices constrained"),
             ifelse(is.null(weight_constraints), "", ", alphas constrained")), "\n")
  cat("\n")
  cat("APPROXIMATE STANDARD ERRORS\n\n")

  left_brackets <- rep("[", times=d)
  right_brackets <- rep("]", times=d)
  plus <- c("+", rep(" ", d - 1))
  arch_scalar <- c(rep(" ", times=d-1), "ARCH_mt")
  round_lbrackets <- rep("(", times=d)
  round_rbrackets <- rep(")", times=d)
  Y <- paste0("Y", 1:d)
  tmp_names <- paste0("tmp", 1:(p*(d + 2) + d + 2))

  for(m in seq_len(M)) {
    count <- 1
    if(model == "GMVAR") {
      regime_type <- "GMVAR"
    } else if(model == "StMVAR") {
      regime_type <- "StMVAR"
      M1 <- 0
    } else {
      M1 <- M_orig[1]
      regime_type <- ifelse(m <= M1, "GMVAR", "StMVAR")
    }

    cat(paste("Regime", m))
    if(model == "G-StMVAR") cat(paste0(" (", regime_type, " type)"))
    cat("\n")
    if(m < M) cat(paste("Mixing weight:", format_value(alphas[m])), "\n")
    if(parametrization == "mean") cat("Regime means:", paste0(format_value(all_mu[,m]), collapse=", "), "\n")
    if(regime_type == "StMVAR") { # Print degrees of freedom parameter for StMVAR type regimes
      cat("Df parameter: ", format_value(all_df[m - M1]), "\n")
    }
    cat("\n")
    df <- data.frame(Y=Y,
                     eq=c("=", rep(" ", d - 1)),
                     eq=left_brackets,
                     phi0=format_value(all_phi0[, m, drop=FALSE]),
                     eq=right_brackets,
                     plus)
    for(i1 in seq_len(p)) {
      Amp_colnames <- c(paste0("A", i1), tmp_names[count:(count + d - 1 - 1)]); count <- count + d - 1
      df[, tmp_names[count]] <- left_brackets; count <- count + 1
      df[, Amp_colnames] <- format_value(all_A[, ,i1 , m])
      df[, tmp_names[count]] <- right_brackets; count <- count + 1
      df[, tmp_names[count]] <- paste0(Y, ".", i1); count <- count + 1
      df <- cbind(df, plus)
    }
    if(regime_type == "StMVAR") { # Time varying ARCH scalar multiplying the constant part of error term covariance matrix
      df <- cbind(df, round_lbrackets, arch_scalar)
    }
    df[, tmp_names[p*(d + 2) + 1]] <- left_brackets
    df[, c("Omega", tmp_names[(p*(d + 2) + 2):(p*(d + 2) + d)])] <- format_value(all_Omega[, , m])
    df[, tmp_names[p*(d + 2) + d + 1]] <- right_brackets
    df[, "1/2"] <- rep(" ", d)
    df[, tmp_names[p*(d + 2) + d + 2]] <- paste0("eps", 1:d)
    names_to_omit <- unlist(lapply(c("plus", "eq", "arch_scalar", "round_lbrackets", "round_rbrackets", tmp_names),
                                   function(nam) grep(nam, colnames(df))))
    colnames(df)[names_to_omit] <- " "
    print(df)
    cat("\n")
  }
  if(sep_AR) cat(paste0("AR parameters: ", paste0(format_value(AR_stds), collapse=", ")), "\n\n")

  if(!is.null(structural_pars)) {
    cat("Structural parameters:\n")
    W <- format_value(pick_W(p=p, M=M, d=d, params=pars, structural_pars=structural_pars))

    if(M > 1) {
      lambdas <- format_value(pick_lambdas(p=p, M=M, d=d, params=pars, structural_pars=structural_pars))
      lambdas <- matrix(lambdas, nrow=d, ncol=M - 1, byrow=FALSE) # Column for each regime

      # Similarly to the AR parameters, constrained lambda parameter standard errors multiplied open
      # in 'pars' are valid iff the constraint matrix "C_lambda" contains zeros and ones only, and
      # there is at most one one in each row (no multiplications or summations).
      C_lambda <- gsmvar$model$structural_pars$C_lambda
      if(!is.null(C_lambda)) {
        if(any(C_lambda != 1 & C_lambda != 0) | any(rowSums(C_lambda) > 1)) {
          sep_lambda <- TRUE # The lambda parameter std errors must be printed separately
          lambdas <- matrix(NA, nrow=d, ncol=M - 1)
          n_zeros <- sum(W == 0, na.rm=TRUE)
          lambda_stds <- gsmvar$std_errors[(M*d + M*d^2*p + d^2 - n_zeros + 1):(M*d + M*d^2*p + d^2 - n_zeros + ncol(C_lambda))]
        } else {
          sep_lambda <- FALSE
        }
      } else {
        sep_lambda <- FALSE
      }
    }
    if(!is.null(gsmvar$model$structural_pars$fixed_lambdas)) {
      sep_lambda <- TRUE
    }

    tmp <- c(rep(" ", times=d - 1), ",")
    df2 <- data.frame(left_brackets, W=W[,1])
    for(i1 in 2:d) {
      df2 <- cbind(df2, W[, i1])
      colnames(df2)[1 + i1] <- "tmp"
    }
    df2 <- cbind(df2, right_brackets)
    if(M > 1) {
      tmp <- c(rep(" ", times=d - 1), ",")
      for(i1 in 1:(M - 1)) {
        if(sep_lambda) {
          lmb <- rep(NA, times=d)
        } else {
          lmb <- lambdas[,i1]
        }
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
    if(sep_lambda && is.null(gsmvar$model$structural_pars$fixed_lambdas)) {
      cat(paste0("lambda parameters: ", paste0(format_value(lambda_stds), collapse=", ")), "\n\n")
    }
    cat("The B-matrix (or equally W) is subject to", n_zero, "zero constraints and", n_sign, "sign constraints.\n")
    cat("The eigenvalues lambda_{mi} are", ifelse(is.null(gsmvar$model$structural_pars$C_lambda), "not subject to linear constraints.",
                                                  "subject to linear constraints."))
    cat("\n")
  }
  invisible(gsmvar)
}

