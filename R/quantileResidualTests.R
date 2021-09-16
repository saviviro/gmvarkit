#' @title Quantile residual tests
#'
#' @description \code{quantile_residual_tests} performs quantile residual tests described
#'  by \emph{Kalliovirta and Saikkonen 2010}, testing autocorrelation, conditional heteroskedasticity,
#'  and normality.
#'
#' @inheritParams quantile_residuals
#' @inheritParams loglikelihood_int
#' @param lags_ac a positive integer vector specifying the lags used to test autocorrelation.
#' @param lags_ch a positive integer vector specifying the lags used to test conditional heteroskedasticity.
#' @param nsimu to how many simulations should the covariance matrix Omega used in the qr-tests be based on?
#'   If smaller than sample size, then the covariance matrix will be evaluated from the sample. Larger number
#'   of simulations might improve the tests size properties but it increases the computation time.
#' @param print_res should the test results be printed while computing the tests?
#' @details If the function fails to calculate the tests because of numerical problems and the parameter values
#'   are near the border of the parameter space, it might help to use smaller numerical tolerance for the
#'   stationarity and positeve definiteness conditions. The numerical tolerance of an existing model
#'   can be changed with the function \code{update_numtols} or you can set it directly with the arguments
#'   \code{stat_tol} and \code{posdef_tol}.
#' @return Returns an object of class \code{'qrtest'} which has its own print method. The returned object
#'   is a list containing the quantile residual test results for normality, autocorrelation, and conditional
#'   heteroskedasticity. The autocorrelation and conditional heteroskedasticity results also contain the
#'   associated (vectorized) individual statistics divided by their standard errors
#'   (see \emph{Kalliovirta and Saikkonen 2010}, s.17-20) under the label \code{$ind_stats}.
#' @inherit quantile_residuals references
#' @seealso \code{\link{fitGMVAR}}, \code{\link{GMVAR}}, \code{\link{quantile_residuals}}, \code{\link{GIRF}},
#'   \code{\link{diagnostic_plot}}, \code{\link{predict.gmvar}}, \code{\link{profile_logliks}},
#'   \code{\link{LR_test}}, \code{\link{Wald_test}}, \code{\link{cond_moment_plot}}, \code{\link{update_numtols}}
#' @examples
#' \donttest{
#' # GMVAR(3,2) model
#' fit32 <- fitGMVAR(gdpdef, p=3, M=2, ncalls=1, seeds=2)
#' qrtests32 <- quantile_residual_tests(fit32)
#' qrtests32
#' plot(qrtests32)
#'
#' # Structural GMVAR(1,2) model identified with sign
#' # constraints and build with hand-specified parameter values.
#' # Tests based on simulation procedure with nsimu=1000:
#' params12s <- c(0.55, 0.112, 0.619, 0.173, 0.344, 0.055, -0.009, 0.718,
#'  0.255, 0.017, -0.136, 0.858, 0.541, 0.057, -0.162, 0.162, 3.623,
#'  4.726, 0.674)
#' W_12 <- matrix(c(1, 1, -1, 1), nrow=2)
#' mod12s <- GMVAR(gdpdef, p=1, M=2, params=params12s,
#'                 structural_pars=list(W=W_12))
#' qrtests12s <- quantile_residual_tests(mod12s, nsimu=1000)
#' qrtests12s
#' }
#' @export

quantile_residual_tests <- function(gmvar, lags_ac=c(1, 3, 6, 12), lags_ch=lags_ac, nsimu=1, print_res=TRUE,
                                    stat_tol, posdef_tol) {
  check_gmvar(gmvar)
  check_null_data(gmvar)
  if(!all_pos_ints(c(lags_ac, lags_ch))) stop("arguments 'lags_ac' and 'lags_ch' must be strictly positive integer vectors")
  p <- gmvar$model$p
  M <- gmvar$model$M
  d <- gmvar$model$d
  data <- gmvar$data
  n_obs <- nrow(data)
  T_obs <- n_obs - p
  if(missing(stat_tol)) stat_tol <- gmvar$num_tols$stat_tol
  if(missing(posdef_tol)) posdef_tol <- gmvar$num_tols$posdef_tol
  if(max(c(lags_ac, lags_ch)) >= T_obs) stop("The lags are too large compared to the data size")

  qresiduals <- gmvar$quantile_residuals
  if(nsimu > n_obs) {
    omega_data <- simulateGMVAR(gmvar, nsimu=nsimu, init_values=NULL, ntimes=1)$sample
  } else {
    omega_data <- data
  }

  # Function to calculate covariance matrix omega
  try_to_get_omega <- function(g, dim_g, which_test, which_lag=NA) {
    print_message <- function(which_test, which_lag, because_of) {
      if(which_test == "norm") {
        message(paste("Can't perform normality test", because_of))
      } else if(which_test == "ac") {
        message(paste("Can't perform autocorrelation test for lag", which_lag, because_of))
      } else if(which_test == "ch") {
        message(paste("Can't perform conditional heteroskedasticity test for lag", which_lag, because_of))
      }
    }
    omg <- tryCatch(get_test_Omega(data=omega_data, p=p, M=M,
                                   params=gmvar$params,
                                   conditional=gmvar$model$conditional,
                                   parametrization=gmvar$model$parametrization,
                                   constraints=gmvar$model$constraints,
                                   same_means=gmvar$model$same_means,
                                   structural_pars=gmvar$model$structural_pars,
                                   g=g, dim_g=dim_g, stat_tol=stat_tol,
                                   posdef_tol=posdef_tol),
                    error=function(e) {
                      print_message(which_test, which_lag, because_of="because of numerical problems")
                      return(NA)
                    })
    if(is.matrix(omg) & anyNA(omg)) {
      print_message(which_test, which_lag, because_of="- probably because the model fits too poorly")
    } else if(length(omg) == 1) {
      if(is.na(omg)) return(matrix(NA, nrow=dim_g, ncol=dim_g))
    }
    omg
  }

  # Function to calculate general test statistic
  calc_test_stat <- function(g, m_dim, Omega) {
    g_qres <- g(qresiduals)
    sumg <- colSums(g_qres[m_dim:nrow(g_qres),])
    t(sumg)%*%solve(Omega, sumg)/(T_obs - m_dim + 1)
  }

  # Function to calculate individual statistics divided by their standard errors
  calc_ind_stats <- function(g, Omega) {
    g_qres <- g(qresiduals)
    g_qres2 <- g_qres[,(ncol(g_qres) - d^2 + 1):ncol(g_qres)] # Take the last d^2 columns: r_t*r_{t-K}' or v_t*v_{t-K}'
    c_lag <- rowMeans(array(t(g_qres2), dim=c(d, d, nrow(g_qres2))), dims=2)
    c_stderr <- sqrt(diag(Omega)[(ncol(Omega) - d^2 + 1):ncol(Omega)]/T_obs)
    vec(c_lag)/c_stderr # individual statistic divided by it's standard error
  }

  format_value0 <- format_valuef(0)
  format_value3 <- format_valuef(3)
  print_resf <- function(lag, p_val) {
    if(lag < 10) {
      cat(" ", format_value0(lag), " | ", format_value3(p_val), "\n")
    } else {
      cat(" ", format_value0(lag), "| ", format_value3(p_val), "\n")
    }
  }

  ######################
  # Test for normality # (Kalliovirta and Saikkonen 2010, sec. 3.3)
  ######################

  dim_g <- 3*d
  g <- function(r) { # "r" should be (T x d) quantile residual matrix
    T0 <- nrow(r)
    matrix((vapply(1:d, function(j) c(r[,j]^2 - 1, r[,j]^3, r[,j]^4 - 3), numeric(T0*3))), nrow=T0, ncol=dim_g, byrow=FALSE)
  } # Returns (T x dim_g) matrix with values of g_t at each row

  # Get estimated Omega based on large sample and calculate the test statistic and p-value
  Omega <- try_to_get_omega(g=g, dim_g=dim_g, which_test="norm", which_lag=NA)
  N <- calc_test_stat(g=g, m_dim=1, Omega=Omega)
  p_val <- 1 - pchisq(N, df=dim_g)

  if(print_res) cat(paste0("Normality test p-value: ", format_value3(p_val)), "\n\n")
  norm_res <- data.frame(test_stat=N, df=dim_g, p_val=p_val)


  ############################
  # Test for autocorrelation # (Kalliovirta and Saikkonen 2010, sec. 3.1)
  ############################

  tmp <- rep(NA, length(lags_ac))
  ac_res <- list(test_results=data.frame(lags=lags_ac, test_stat=tmp, df=tmp, p_val=tmp),
                 ind_stats=data.frame(row.names=1:d^2))

  # Function factory to produce function g for different lags
  get_g <- function(lag) {
    function(r) {
      t(vapply((lag + 1):nrow(r), function(t) vapply(1:lag, function(i1) tcrossprod(r[t,], r[t - i1,]), numeric(d^2)), numeric(lag*d^2)))
    }
  } # Returns (T - lag x dim_g) matrix with values of g_t at each row, starting from t=lag+1 at the first row

  if(print_res) cat("Autocorrelation tests:\nlags | p-value\n")

  for(i1 in seq_along(lags_ac)) {
    lag <- lags_ac[i1]
    dim_g <- lag*d^2
    g <- get_g(lag)
    m_dim <- lag + 1
    Omega <- try_to_get_omega(g=g, dim_g=dim_g, which_test="ac", which_lag=lag)
    A <- calc_test_stat(g=g, m_dim=m_dim, Omega=Omega)
    p_val <- 1 - pchisq(A, df=dim_g)

    # Calculate the individual statistics c_s (Kalliovirta and Saikkonen 2010, s.17)
    # in vectorised form and obtain their standard errors from relevant diagonal of Omega.
    ac_res$ind_stats[, paste0("lag", lag)] <- calc_ind_stats(g=g, Omega=Omega) # individual statistic divided by it's standard error

    if(print_res) print_resf(lag=lag, p_val=p_val)
    ac_res$test_results[i1, 2:4] <- c(A, dim_g, p_val)
  }


  ###########################################
  # Test for conditional heteroskedasticity # (Kalliovirta and Saikkonen 2010, sec. 3.2)
  ###########################################

  tmp <- rep(NA, length(lags_ch))
  ch_res <- list(test_results=data.frame(lags=lags_ch, test_stat=tmp, df=tmp, p_val=tmp),
                 ind_stats=data.frame(row.names=1:d^2))

  # Function factory to produce function g for different lags
  get_g <- function(lag) {
    function(r) {
      v <- r^2 - 1
      t(vapply((lag + 1):nrow(v), function(t) vapply(1:lag, function(i1) tcrossprod(v[t,], v[t - i1,]), numeric(d^2)), numeric(lag*d^2)))
    }
  } # Returns (T - lag x dim_g) matrix with values of g_t at each row, starting from t=lag+1 at the first row

  if(print_res) cat("\nConditional heteroskedasticity tests:\nlags | p-value\n")

  for(i1 in seq_along(lags_ch)) {
    lag <- lags_ch[i1]
    dim_g <- lag*d^2
    g <- get_g(lag)
    m_dim <- lag + 1
    Omega <- try_to_get_omega(g=g, dim_g=dim_g, which_test="ch", which_lag=lag)
    H <- calc_test_stat(g=g, m_dim=m_dim, Omega=Omega)
    p_val <- 1 - pchisq(H, df=dim_g)

    # Calculate the individual statistics d_s (Kalliovirta and Saikkonen 2010, s.19)
    # in vectorised form and obtain their standard errors from relevant diagonal of Omega.
    ch_res$ind_stats[, paste0("lag", lag)] <- calc_ind_stats(g=g, Omega=Omega) # individual statistic divided by it's standard error

    if(print_res) print_resf(lag=lag, p_val=p_val)
    ch_res$test_results[i1, 2:4] <- c(H, dim_g, p_val)
  }

  structure(list(norm_res=norm_res,
                 ac_res=ac_res,
                 ch_res=ch_res),
            class="qrtest")
}



#' @title Compute covariance matrix Omega used in quantile residual tests
#'
#' @description \code{get_test_Omega} computes the covariance matrix Omega used in the
#'  quantile residuals tests described by \emph{Kalliovirta and Saikkonen 2010}.
#'
#' @inheritParams loglikelihood
#' @param g function g specifying the transformation.
#' @param dim_g output dimension of the transformation \code{g}.
#' @return Returns the covariance matrix Omega described by \emph{Kalliovirta and Saikkonen 2010}.
#' @inherit quantile_residuals references
#' @keywords internal

get_test_Omega <- function(data, p, M, params, conditional, parametrization, constraints, same_means, structural_pars=NULL, g, dim_g,
                           stat_tol=1e-3, posdef_tol=1e-8) {

  n_obs <- nrow(data)
  T_obs <- n_obs - p
  d <- ncol(data)
  minval <- get_minval(data)

  # Function used to to calculate gradient for function g
  g_fn <- function(pars) {
    qresiduals <- quantile_residuals_int(data=data, p=p, M=M, params=pars, conditional=conditional,
                                         parametrization=parametrization, constraints=constraints,
                                         same_means=same_means, structural_pars=structural_pars,
                                         stat_tol=stat_tol, posdef_tol=posdef_tol)
    g(qresiduals) # a row for each t=1,...,T and column for each output of g
  }

  # Function used to calculate gradient for log-likelihood
  loglik_fn <- function(pars) {
    loglikelihood_int(data, p, M, params=pars, conditional=conditional, parametrization=parametrization,
                      constraints=constraints, same_means=same_means, structural_pars=structural_pars,
                      check_params=TRUE, to_return="terms", minval=minval, stat_tol=stat_tol, posdef_tol=posdef_tol)
  }

  npars <- length(params)
  I <- diag(1, nrow=npars, ncol=npars)
  h <- 6e-06
  central_diff <- function(params, fn, i1) (fn(params + h*I[i1,]) - fn(params - h*I[i1,]))/(2*h)

  g_qres <- g_fn(params) # Function g applied to model quantile residuals, row for each t and column for each output of g().
  T0 <- nrow(g_qres)

  # Calculate matrix G (Kalliovirta ja Saikkonen 2010, s.13)
  dg <- array(dim=c(T0, dim_g, npars))
  for(i1 in 1:npars) {
    dg[, , i1] <- central_diff(params, g_fn, i1)
  }
  G <- colMeans(dg)

  # Calculate gradients of the terms l_t
  dl <- vapply(1:npars, function(i1) central_diff(params, loglik_fn, i1), numeric(T_obs)) # (T x npars)

  # Approximate Fisher information matrix, calculate Psi matrix and H matrix
  diff0 <- nrow(dl) - T0
  Fish_inf <- crossprod(dl, dl)/nrow(dl)
  Psi <- crossprod(g_qres, dl[(1 + diff0):nrow(dl),])/nrow(g_qres)
  H <- crossprod(g_qres, g_qres)/nrow(g_qres)

  inv_Fish <- solve(Fish_inf) # Can cause error sometimes since this is not always (numerically) invertible

  # Calculate covariance matrix Omega
  FG <- tcrossprod(inv_Fish, G)
  G%*%FG + Psi%*%FG + G%*%tcrossprod(inv_Fish, Psi) + H
}

