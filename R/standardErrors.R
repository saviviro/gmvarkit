

#' @title Calculate standard errors for estimates of GMVAR model
#'
#' @description \code{standard_errors} numerically approximates standard errors for the given estimates of GMVAR model using square
#'   roots of the diagonal of inverse of observed information matrix.
#'
#' @inheritParams loglikelihood_int
#' @return a vector containing the approximate standard errors of the estimates
#' @inherit in_paramspace_int references

standard_errors <- function(data, p, M, params, conditional=TRUE, parametrization=c("intercept", "mean"), constraints=NULL, minval) {

  parametrization <- match.arg(parametrization)

  loglik_fn <- function(params) {
    tryCatch(loglikelihood_int(data, p, M, params=params, conditional=conditional, parametrization=parametrization,
                      constraints=constraints, check_params=TRUE, to_return="loglik", minval=minval),
             error=function(e) NA)
  }

  npars <- length(params)
  I <- diag(1, ncol=npars, nrow=npars) # Indicates which parameter is derivated

  # Calculate Hessian
  Hess <- calc_hessian(x=params, fn=loglik_fn, h=6e-6)

  # Inverse of the observed information matrix
  inv_obs_inf <- tryCatch(solve(-Hess), error=function(e) matrix(NA, nrow=npars, ncol=npars))

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
#' @seealso \code{\link{fitGMVAR}}, \code{\link{GMVAR}}, \code{\link{print.gmvar}}, \code{\link{swap_parametrization}}
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
  pars <- reform_constrained_pars(p=p, M=M, d=d, params=gmvar$std_errors, constraints=constraints, change_na=TRUE)
  all_phi0_or_mu <- pick_phi0(p=p, M=M, d=d, params=pars)
  all_A <- pick_allA(p=p, M=M, d=d, params=pars)
  all_Omega <- pick_Omegas(p=p, M=M, d=d, params=pars)
  alphas <- pick_alphas(p=p, M=M, d=d, params=pars)
  alphas[M] <- NA
  if(gmvar$model$parametrization == "mean") {
    all_mu <- all_phi0_or_mu
    all_phi0 <- matrix(NA, nrow=d, ncol=M)
  } else {
    all_mu <- matrix(NA, nrow=d, ncol=M)
    all_phi0 <- all_phi0_or_mu
  }

  cat("Model:\n")
  cat(paste0("p = ", p, ", M = ", M, ","),
      ifelse(gmvar$model$conditional, "conditional,", "exact,"),
      ifelse(gmvar$model$parametrization=="mean", "mean parametrization,", "intercept parametrization,"),
      ifelse(is.null(constraints), "no constraints", "linear constraints employed"), "\n")
  cat("\n")
  cat("APPROXIMATE STANDARD ERRORS\n\n")

  plus <- c("+", rep(" ", d-1))
  Y <- paste0("Y", 1:d)
  tmp_names <- paste0("tmp", 1:(p*(d + 2) + d + 2))

  for(m in seq_len(M)) {
    count <- 1
    cat(paste("Regime", m), "\n")
    cat(paste("Mixing weight:", format_value(alphas[m])), "\n")
    cat("Regime means:", paste0(format_value(all_mu[,m]), collapse=", "), "\n\n")
    df <- data.frame(Y=Y,
                     eq=c("=", rep(" ", d-1)),
                     eq=rep("[", d),
                     phi0=format_value(all_phi0[, m, drop=FALSE]),
                     eq=rep("]", d),
                     plus)
    for(i1 in seq_len(p)) {
      Amp_colnames <- c(paste0("A", i1), tmp_names[count:(count + d - 1 - 1)]); count <- count + d - 1
      df[, tmp_names[count]] <- rep("[", d); count <- count + 1
      df[, Amp_colnames] <- format_value(all_A[, ,i1 , m])
      df[, tmp_names[count]] <- rep("]", d); count <- count + 1
      df[, tmp_names[count]] <- paste0(Y, ".l", i1); count <- count + 1
      df <- cbind(df, plus)
    }
    df[, tmp_names[p*(d + 2) + 1]] <- rep("[", d)
    df[, c("Omega", tmp_names[(p*(d + 2) + 2):(p*(d + 2) + d)])] <- format_value(all_Omega[, , m])
    df[, tmp_names[p*(d + 2) + d + 1]] <- rep("]", d)
    df[, "1/2"] <- rep(" ", d)
    df[, tmp_names[p*(d + 2) + d + 2]] <- paste0("eps", 1:d)
    names_to_omit <- unlist(lapply(c("plus", "eq", tmp_names), function(nam) grep(nam, colnames(df))))
    colnames(df)[names_to_omit] <- " "
    print(df)
    cat("\n")
  }
  invisible(gmvar)
}

