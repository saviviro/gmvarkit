
#' @title Calculate standard errors for estimates of GMVAR model
#'
#' @description \code{standard_errors} numerically calculates approximate standard errors for the GMVAR model using square
#'   roots of the diagonal of inverse of observed information matrix.
#'
#' @inheritParams loglikelihood_int
#' @return A vector containing the approximate standard errors of the estimates.
#' @inherit in_paramspace_int references

standard_errors <- function(data, p, M, params, conditional=TRUE, parametrization=c("intercept", "mean"),
                            constraints=NULL, structural_pars=NULL, minval, stat_tol=1e-3, posdef_tol=1e-8) {

  parametrization <- match.arg(parametrization)

  loglik_fn <- function(params) {
    tryCatch(loglikelihood_int(data=data, p=p, M=M, params=params, conditional=conditional, parametrization=parametrization,
                               constraints=constraints, structural_pars=structural_pars, check_params=TRUE,
                               to_return="loglik", minval=minval, stat_tol=stat_tol, posdef_tol=posdef_tol),
             error=function(e) NA)
  }

  # Calculate Hessian
  Hess <- calc_hessian(x=params, fn=loglik_fn, h=6e-6)

  # Inverse of the observed information matrix
  inv_obs_inf <- tryCatch(solve(-Hess), error=function(e) matrix(NA, nrow=length(params), ncol=length(params)))

  # Calculate the standard errors
  diag_inv_obs_inf <- diag(inv_obs_inf)
  unlist(lapply(diag_inv_obs_inf, function(x) ifelse(is.na(x) | x < 0, NA, sqrt(x))))
}


#' @title Print standard errors of GMVAR model in the same form as the model estimates are printed
#'
#' @description \code{print_std_errors} prints the approximate standard errors of GMVAR model in the
#'   same form as the parameters of objects of class \code{'gmvar'} are printed.
#'
#' @inheritParams simulateGMVAR
#' @param digits how many digits should be printed?
#' @details The main purpose of \code{print_std_errors} is to provide a convenient tool to match the standard
#'   errors to certain parameter estimates.
#'
#'   Note that if linear constraints are imposed and they involve summations or multiplications, then the AR
#'   parameter standard errors are printed separately as they don't correspond one-to-one to the model parameter
#'   standard errors.
#' @seealso \code{\link{profile_logliks}}, \code{\link{fitGMVAR}}, \code{\link{GMVAR}}, \code{\link{print.gmvar}},
#'  \code{\link{swap_parametrization}}
#' @inherit GMVAR references
#' @examples
#' \donttest{
#' ## These are long running examples that use parallel computing!
#'
#' # These examples use the data 'eurusd' which comes with the
#' # package, but in a scaled form.
#' data <- cbind(10*eurusd[,1], 100*eurusd[,2])
#' colnames(data) <- colnames(eurusd)
#'
#' # GMVAR(1,2) model with default settings
#' fit12 <- fitGMVAR(data, p=1, M=2)
#' fit12
#' print_std_errors(fit12)
#'
#' # GMVAR(2,2) model with mean parametrization
#' fit22 <- fitGMVAR(data, p=2, M=2, parametrization="mean")
#' fit22
#' print_std_errors(fit22)
#'
#' # GMVAR(2,2) model with autoregressive parameters restricted
#' # to be the same for all regimes
#' C_mat <- rbind(diag(2*2^2), diag(2*2^2))
#' fit22c <- fitGMVAR(data, p=2, M=2, constraints=C_mat)
#' fit22c
#' print_std_errors(fit22c)
#'
#' # GMVAR(2,2) model with autoregressive parameters restricted
#' # to be the same for all regimes and non-diagonl elements
#' # the coefficient matrices constrained to zero.
#' tmp <- matrix(c(1, rep(0, 10), 1, rep(0, 8), 1, rep(0, 10), 1),
#'  nrow=2*2^2, byrow=FALSE)
#' C_mat2 <- rbind(tmp, tmp)
#' fit22c2 <- fitGMVAR(data, p=2, M=2, constraints=C_mat2, ncalls=10)
#' fit22c2
#' print_std_errors(fit22c2)
#' }
#' @export

print_std_errors <- function(gmvar, digits=3) {
  if(!all_pos_ints(digits)) stop("Argument digits must be positive integer")
  format_value <- format_valuef(digits)
  p <- gmvar$model$p
  M <- gmvar$model$M
  d <- gmvar$model$d
  constraints <- gmvar$model$constraints
  pars <- reform_constrained_pars(p=p, M=M, d=d, params=gmvar$std_errors, constraints=constraints,
                                  structural_pars=gmvar$model$structural_pars, change_na=TRUE)
  structural_pars <- get_unconstrained_structural_pars(structural_pars=gmvar$model$structural_pars)
  all_phi0_or_mu <- pick_phi0(p=p, M=M, d=d, params=pars, structural_pars=structural_pars)
  all_A <- pick_allA(p=p, M=M, d=d, params=pars, structural_pars=structural_pars)
  if(is.null(structural_pars)) {
    all_Omega <- pick_Omegas(p=p, M=M, d=d, params=pars, structural_pars=structural_pars)
  } else {
    # No standard errors for cov. mats. as the model is parametrized with W and lambdas
    all_Omega <- array(NA, dim=c(d, d, M))
  }
  alphas <- pick_alphas(p=p, M=M, d=d, params=pars)
  alphas[M] <- NA
  if(gmvar$model$parametrization == "mean") {
    all_mu <- all_phi0_or_mu
    all_phi0 <- matrix(NA, nrow=d, ncol=M)
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
      all_A <- array(NA, dim=c(d, d, p, M))
      AR_stds <- gmvar$std_errors[(M*d + 1):(M*d + ncol(constraints))] # Constrained AR param std errors
    } else {
      sep_AR <- FALSE
    }
  } else {
    sep_AR <- FALSE # No constraints imposed
  }

  cat(ifelse(is.null(structural_pars), "Reduced form", "Structural"), "model:\n")
  cat(paste0("p = ", p, ", M = ", M, ","),
      ifelse(gmvar$model$conditional, "conditional", "exact"),
      "log-likelihood,",
      ifelse(gmvar$model$parametrization == "mean", "mean parametrization,", "intercept parametrization,"),
      ifelse(is.null(constraints), "no AR parameter constraints", "linear constraints imposed on AR parameters"), "\n")
  cat("\n")
  cat("APPROXIMATE STANDARD ERRORS\n\n")

  left_brackets <- rep("[", times=d)
  right_brackets <- rep("]", times=d)
  plus <- c("+", rep(" ", d - 1))
  Y <- paste0("Y", 1:d)
  tmp_names <- paste0("tmp", 1:(p*(d + 2) + d + 2))

  for(m in seq_len(M)) {
    count <- 1
    cat(paste("Regime", m), "\n")
    cat(paste("Mixing weight:", format_value(alphas[m])), "\n")
    cat("Regime means:", paste0(format_value(all_mu[,m]), collapse=", "), "\n\n")
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
    df[, tmp_names[p*(d + 2) + 1]] <- left_brackets
    df[, c("Omega", tmp_names[(p*(d + 2) + 2):(p*(d + 2) + d)])] <- format_value(all_Omega[, , m])
    df[, tmp_names[p*(d + 2) + d + 1]] <- right_brackets
    df[, "1/2"] <- rep(" ", d)
    df[, tmp_names[p*(d + 2) + d + 2]] <- paste0("eps", 1:d)
    names_to_omit <- unlist(lapply(c("plus", "eq", tmp_names), function(nam) grep(nam, colnames(df))))
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
      C_lambda <- gmvar$model$structural_pars$C_lambda
      if(!is.null(C_lambda)) {
        if(any(C_lambda != 1 & C_lambda != 0) | any(rowSums(C_lambda) > 1)) {
          sep_lambda <- TRUE # The lambda parameter std errors must be printed separately
          lambdas <- matrix(NA, nrow=d, ncol=M - 1)
          n_zeros <- sum(W == 0, na.rm=TRUE)
          lambda_stds <- gmvar$std_errors[(M*d + M*d^2*p + d^2 - n_zeros + 1):(M*d + M*d^2*p + d^2 - n_zeros + ncol(C_lambda))]
        } else {
          sep_lambda <- FALSE
        }
      } else {
        sep_lambda <- FALSE
      }
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
    W_orig <- gmvar$model$structural_pars$W
    n_zero <- sum(W_orig == 0, na.rm=TRUE)
    n_free <- sum(is.na(W_orig))
    n_sign <- d^2 - n_zero - n_free
    if(sep_lambda) cat(paste0("lambda parameters: ", paste0(format_value(lambda_stds), collapse=", ")), "\n\n")
    cat("The B-matrix (or equally W) is subject to", n_zero, "zero constraints and", n_sign, "sign constraints.\n")
    cat("The eigenvalues lambda_{mi} are", ifelse(is.null(gmvar$model$structural_pars$C_lambda), "not subject to linear constraints.",
                                                  "subject to linear constraints."))
    cat("\n")
  }
  invisible(gmvar)
}

