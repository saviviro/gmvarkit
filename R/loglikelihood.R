#' @import stats
#'
#' @title Compute log-likelihood of a GMVAR, StMVAR, and G-StMVAR models
#'
#' @description \code{loglikelihood_int} computes log-likelihoodof a GMVAR, StMVAR, and G-StMVAR models.
#'
#' @inheritParams in_paramspace_int
#' @param data a matrix or class \code{'ts'} object with \code{d>1} columns. Each column is taken to represent
#'  a univariate time series. \code{NA} values are not supported.
#' @param p a positive integer specifying the autoregressive order of the model.
#' @param M \describe{
#'   \item{For \strong{GMVAR} and \strong{StMVAR} models:}{a positive integer specifying the number of mixture components.}
#'   \item{For \strong{G-StMVAR} models:}{a size (2x1) integer vector specifying the number of \emph{GMVAR type} components \code{M1}
#'    in the first element and \emph{StMVAR type} components \code{M2} in the second element. The total number of mixture components
#'    is \code{M=M1+M2}.}
#' }
#' @param params a real valued vector specifying the parameter values.
#'   \describe{
#'     \item{\strong{For unconstrained models:}}{
#'       Should be size \eqn{((M(pd^2+d+d(d+1)/2+2)-M1-1)x1)} and have the form
#'       \strong{\eqn{\theta}}\eqn{ = }(\strong{\eqn{\upsilon}}\eqn{_{1}},
#'       ...,\strong{\eqn{\upsilon}}\eqn{_{M}}, \eqn{\alpha_{1},...,\alpha_{M-1},}\strong{\eqn{\nu}}\eqn{)}, where
#'       \itemize{
#'         \item \strong{\eqn{\upsilon}}\eqn{_{m}} \eqn{ = (\phi_{m,0},}\strong{\eqn{\phi}}\eqn{_{m}}\eqn{,\sigma_{m})}
#'         \item \strong{\eqn{\phi}}\eqn{_{m}}\eqn{ = (vec(A_{m,1}),...,vec(A_{m,p})}
#'         \item and \eqn{\sigma_{m} = vech(\Omega_{m})}, m=1,...,M,
#'         \item \strong{\eqn{\nu}}\eqn{=(\nu_{M1+1},...,\nu_{M})}
#'         \item \eqn{M1} is the number of GMVAR type regimes.
#'       }
#'     }
#'     \item{\strong{For constrained models:}}{
#'       Should be size \eqn{((M(d+d(d+1)/2+2)+q-M1-1)x1)} and have the form
#'       \strong{\eqn{\theta}}\eqn{ = (\phi_{1,0},...,\phi_{M,0},}\strong{\eqn{\psi},}
#'       \eqn{\sigma_{1},...,\sigma_{M},\alpha_{1},...,\alpha_{M-1},}\strong{\eqn{\nu}}), where
#'       \itemize{
#'         \item \strong{\eqn{\psi}} \eqn{(qx1)} satisfies (\strong{\eqn{\phi}}\eqn{_{1}}\eqn{,...,}
#'         \strong{\eqn{\phi}}\eqn{_{M}) =} \strong{\eqn{C \psi}} where \strong{\eqn{C}} is a \eqn{(Mpd^2xq)}
#'         constraint matrix.
#'       }
#'     }
#'     \item{\strong{For same_means models:}}{
#'       Should have the form
#'       \strong{\eqn{\theta}}\eqn{ = (}\strong{\eqn{\mu},}\strong{\eqn{\psi},}
#'       \eqn{\sigma_{1},...,\sigma_{M},\alpha_{1},...,\alpha_{M-1},}\strong{\eqn{\nu}}\eqn{)}, where
#'       \itemize{
#'         \item \strong{\eqn{\mu}}\eqn{= (\mu_{1},...,\mu_{g})} where
#'           \eqn{\mu_{i}} is the mean parameter for group \eqn{i} and
#'           \eqn{g} is the number of groups.
#'         \item If AR constraints are employed, \strong{\eqn{\psi}} is as for constrained
#'           models, and if AR constraints are not employed, \strong{\eqn{\psi}}\eqn{ = }
#'           (\strong{\eqn{\phi}}\eqn{_{1}}\eqn{,...,}\strong{\eqn{\phi}}\eqn{_{M})}.
#'       }
#'     }
#'      \item{\strong{For models with weight_constraints:}}{Drop \eqn{\alpha_1,...,\alpha_{M-1}} from
#'       the parameter vector.}
#'     \item{\strong{For structural models:}}{
#'       Reduced form models can be directly used as recursively identified structural models. If the structural model is
#'       identified by conditional heteroskedasticity, the parameter vector should have the form
#'       \strong{\eqn{\theta}}\eqn{ = (\phi_{1,0},...,\phi_{M,0},}\strong{\eqn{\phi}}\eqn{_{1},...,}\strong{\eqn{\phi}}\eqn{_{M},
#'       vec(W),}\strong{\eqn{\lambda}}\eqn{_{2},...,}\strong{\eqn{\lambda}}\eqn{_{M},\alpha_{1},...,\alpha_{M-1},}\strong{\eqn{\nu}}\eqn{)},
#'        where
#'       \itemize{
#'         \item\strong{\eqn{\lambda}}\eqn{_{m}=(\lambda_{m1},...,\lambda_{md})} contains the eigenvalues of the \eqn{m}th mixture component.
#'       }
#'       \describe{
#'         \item{\strong{If AR parameters are constrained: }}{Replace \strong{\eqn{\phi}}\eqn{_{1}}\eqn{,...,}
#'         \strong{\eqn{\phi}}\eqn{_{M}} with \strong{\eqn{\psi}} \eqn{(qx1)} that satisfies (\strong{\eqn{\phi}}\eqn{_{1}}\eqn{,...,}
#'         \strong{\eqn{\phi}}\eqn{_{M}) =} \strong{\eqn{C \psi}}, as above.}
#'         \item{\strong{If same_means: }}{Replace \eqn{(\phi_{1,0},...,\phi_{M,0})} with \eqn{(\mu_{1},...,\mu_{g})},
#'           as above.}
#'         \item{\strong{If \eqn{W} is constrained:}}{Remove the zeros from \eqn{vec(W)} and make sure the other entries satisfy
#'          the sign constraints.}
#'         \item{\strong{If \eqn{\lambda_{mi}} are constrained via \code{C_lambda}:}}{Replace \strong{\eqn{\lambda}}\eqn{_{2},...,}
#'         \strong{\eqn{\lambda}}\eqn{_{M}} with \strong{\eqn{\gamma}} \eqn{(rx1)} that satisfies (\strong{\eqn{\lambda}}\eqn{_{2}}
#'         \eqn{,...,} \strong{\eqn{\lambda}}\eqn{_{M}) =} \strong{\eqn{C_{\lambda} \gamma}} where \eqn{C_{\lambda}} is
#'          a \eqn{(d(M-1) x r)} constraint matrix.}
#'         \item{\strong{If \eqn{\lambda_{mi}} are constrained via \code{fixed_lambdas}:}}{Drop \strong{\eqn{\lambda}}\eqn{_{2},...,}
#'         \strong{\eqn{\lambda}}\eqn{_{M}} from the parameter vector.}
#'       }
#'     }
#'   }
#'   Above, \eqn{\phi_{m,0}} is the intercept parameter, \eqn{A_{m,i}} denotes the \eqn{i}th coefficient matrix of the \eqn{m}th
#'   mixture component, \eqn{\Omega_{m}} denotes the error term covariance matrix of the \eqn{m}:th mixture component, and
#'   \eqn{\alpha_{m}} is the mixing weight parameter. The \eqn{W} and \eqn{\lambda_{mi}} are structural parameters replacing the
#'   error term covariance matrices (see Virolainen, 2022). If \eqn{M=1}, \eqn{\alpha_{m}} and \eqn{\lambda_{mi}} are dropped.
#'   If \code{parametrization=="mean"}, just replace each \eqn{\phi_{m,0}} with regimewise mean \eqn{\mu_{m}}.
#'   \eqn{vec()} is vectorization operator that stacks columns of a given matrix into a vector. \eqn{vech()} stacks columns
#'   of a given matrix from the principal diagonal downwards (including elements on the diagonal) into a vector.
#'
#'   In the \strong{GMVAR model}, \eqn{M1=M} and \strong{\eqn{\nu}} is dropped from the parameter vector. In the \strong{StMVAR} model,
#'   \eqn{M1=0}. In the \strong{G-StMVAR} model, the first \code{M1} regimes are \emph{GMVAR type} and the rest \code{M2} regimes are
#'   \emph{StMVAR type}. In \strong{StMVAR} and \strong{G-StMVAR} models, the degrees of freedom parameters in \strong{\eqn{\nu}} should
#'   be strictly larger than two.
#'
#'   The notation is similar to the cited literature.
#' @param model is "GMVAR", "StMVAR", or "G-StMVAR" model considered? In the G-StMVAR model, the first \code{M1} components
#'  are GMVAR type and the rest \code{M2} components are StMVAR type.
#' @param conditional a logical argument specifying whether the conditional or exact log-likelihood function
#'  should be used.
#' @param parametrization \code{"intercept"} or \code{"mean"} determining whether the model is parametrized with intercept
#'   parameters \eqn{\phi_{m,0}} or regime means \eqn{\mu_{m}}, m=1,...,M.
#' @param constraints a size \eqn{(Mpd^2 x q)} constraint matrix \strong{\eqn{C}} specifying general linear constraints
#'   to the autoregressive parameters. We consider constraints of form
#'   (\strong{\eqn{\phi}}\eqn{_{1}}\eqn{,...,}\strong{\eqn{\phi}}\eqn{_{M}) = }\strong{\eqn{C \psi}},
#'   where \strong{\eqn{\phi}}\eqn{_{m}}\eqn{ = (vec(A_{m,1}),...,vec(A_{m,p}) (pd^2 x 1), m=1,...,M},
#'   contains the coefficient matrices and \strong{\eqn{\psi}} \eqn{(q x 1)} contains the related parameters.
#'   For example, to restrict the AR-parameters to be the same for all regimes, set \strong{\eqn{C}}=
#'   [\code{I:...:I}]\strong{'} \eqn{(Mpd^2 x pd^2)} where \code{I = diag(p*d^2)}.
#'   Ignore (or set to \code{NULL}) if linear constraints should \strong{not} be employed.
#' @param same_means Restrict the mean parameters of some regimes to be the same? Provide a list of numeric vectors
#'   such that each numeric vector contains the regimes that should share the common mean parameters. For instance, if
#'   \code{M=3}, the argument \code{list(1, 2:3)} restricts the mean parameters of the second and third regime to be
#'   the same but the first regime has freely estimated (unconditional) mean. Ignore or set to \code{NULL} if mean parameters
#'   should not be restricted to be the same among any regimes. \strong{This constraint is available only for mean parametrized models;
#'   that is, when \code{parametrization="mean"}.}
#' @param weight_constraints a numeric vector of length \eqn{M-1} specifying fixed parameter values for the mixing weight parameters
#'   \eqn{\alpha_m, \ m=1,...,M-1}. Each element should be strictly between zero and one, and the sum of all the elements should
#'   be strictly less than one.
#' @param structural_pars If \code{NULL} a reduced form model is considered. Reduced models can be used directly as recursively
#'   identified structural models. For a structural model identified by conditional heteroskedasticity, should be a list containing
#'   at least the first one of the following elements:
#'   \itemize{
#'     \item \code{W} - a \eqn{(dxd)} matrix with its entries imposing constraints on \eqn{W}: \code{NA} indicating that the element is
#'       unconstrained, a positive value indicating strict positive sign constraint, a negative value indicating strict
#'       negative sign constraint, and zero indicating that the element is constrained to zero.
#'     \item \code{C_lambda} - a \eqn{(d(M-1) x r)} constraint matrix that satisfies (\strong{\eqn{\lambda}}\eqn{_{2}}\eqn{,...,}
#'       \strong{\eqn{\lambda}}\eqn{_{M}) =} \strong{\eqn{C_{\lambda} \gamma}} where \strong{\eqn{\gamma}} is the new \eqn{(r x 1)}
#'       parameter subject to which the model is estimated (similarly to AR parameter constraints). The entries of \code{C_lambda}
#'       must be either \strong{positive} or \strong{zero}. Ignore (or set to \code{NULL}) if the eigenvalues \eqn{\lambda_{mi}}
#'       should not be constrained.
#'     \item \code{fixed_lambdas} - a length \eqn{d(M-1)} numeric vector (\strong{\eqn{\lambda}}\eqn{_{2}}\eqn{,...,}
#'       \strong{\eqn{\lambda}}\eqn{_{M})} with elements strictly larger than zero specifying the fixed parameter values for the
#'       parameters \eqn{\lambda_{mi}} should be constrained to. This constraint is alternative \code{C_lambda}.
#'       Ignore (or set to \code{NULL}) if the eigenvalues \eqn{\lambda_{mi}} should not be constrained.
#'   }
#'   See Virolainen (2022) for the conditions required to identify the shocks and for the B-matrix as well (it is \eqn{W} times
#'   a time-varying diagonal matrix with positive diagonal entries).
#' @param check_params should it be checked that the parameter vector satisfies the model assumptions? Can be skipped to save
#'   computation time if it does for sure.
#' @param minval the value that will be returned if the parameter vector does not lie in the parameter space
#'   (excluding the identification condition).
#' @param to_return should the returned object be the log-likelihood value, which is default, or something else?
#'  See the section "Return" for all the options.
#' @details \code{loglikelihood_int} takes use of the function \code{dmvn} from the package \code{mvnfast}.
#' @return
#'  \describe{
#'   \item{By default:}{log-likelihood value of the specified GMVAR, StMVAR, or G-StMVAR model,}
#'   \item{If \code{to_return=="mw"}:}{a size ((n_obs-p)xM) matrix containing the mixing weights: for m:th component in m:th column.}
#'   \item{If \code{to_return=="mw_tplus1"}:}{a size ((n_obs-p+1)xM) matrix containing the mixing weights: for m:th component in m:th column.
#'     The last row is for \eqn{\alpha_{m,T+1}}.}
#'   \item{If \code{to_return=="terms"}:}{a size ((n_obs-p)x1) numeric vector containing the terms \eqn{l_{t}}.}
#'   \item{if \code{to_return=="loglik_and_mw"}:}{a list of two elements. The first element contains the log-likelihood value and the
#'     second element contains the mixing weights.}
#'   \item{If \code{to_return=="regime_cmeans"}:}{an \code{[T-p, d, M]} array containing the regimewise conditional means
#'    (the first p values are used as the initial values).}
#'    \item{If \code{to_return=="regime_ccovs"}:}{an \code{[d, d, T-p, M]} array containing the regimewise conditional
#'    covariance matrices (the first p values are used as the initial values). The index \code{[ , , t, m]} gives the time
#'    \code{t} conditional covariance matrix for the regime \code{m}.}
#'   \item{If \code{to_return=="total_cmeans"}:}{a \code{[T-p, d]} matrix containing the conditional means of the process
#'    (the first p values are used as the initial values).}
#'   \item{If \code{to_return=="total_ccov"}:}{an \code{[d, d, T-p]} array containing the conditional covariance matrices of the process
#'    (the first p values are used as the initial values).}
#'   \item{If \code{to_return=="arch_scalars"}:}{a \code{[T-p, M]} matrix containing the regimewise arch scalars
#'    multiplying error term covariance matrix in the conditional covariance matrix of the regime. For GMVAR type regimes, these
#'    are all ones (the first p values are used as the initial values).}
#'   \item{if \code{to_return=="loglik_mw_archscalars"}:}{a list of three elements. The first element contains the log-likelihood value, the
#'     second element contains the mixing weights, the third element contains the arch scalars
#'     (this is used in \code{quantile_residuals_int}).}
#'  }
#' @references
#'  \itemize{
#'    \item Kalliovirta L., Meitz M. and Saikkonen P. 2016. Gaussian mixture vector autoregression.
#'          \emph{Journal of Econometrics}, \strong{192}, 485-498.
#'    \item Lütkepohl H. 2005. New Introduction to Multiple Time Series Analysis,
#'          \emph{Springer}.
#'    \item McElroy T. 2017. Computation of vector ARMA autocovariances.
#'          \emph{Statistics and Probability Letters}, \strong{124}, 92-96.
#'    \item Virolainen S. 2022. Structural Gaussian mixture vector autoregressive model with application to the asymmetric
#'      effects of monetary policy shocks. Unpublished working paper, available as arXiv:2007.04713.
#'    \item Virolainen S. 2022. Gaussian and Student's t mixture vector autoregressive model with application to the
#'      asymmetric effects of monetary policy shocks in the Euro area. Unpublished working
#'      paper, available as arXiv:2109.13648.
#'  }
#' @keywords internal

loglikelihood_int <- function(data, p, M, params, model=c("GMVAR", "StMVAR", "G-StMVAR"), conditional=TRUE,
                              parametrization=c("intercept", "mean"), constraints=NULL, same_means=NULL,
                              weight_constraints=NULL, structural_pars=NULL,
                              to_return=c("loglik", "mw", "mw_tplus1", "loglik_and_mw", "terms",
                                          "regime_cmeans", "regime_ccovs", "total_cmeans", "total_ccovs",
                                          "arch_scalars", "loglik_mw_archscalars"),
                              check_params=TRUE, minval=NULL, stat_tol=1e-3, posdef_tol=1e-8, df_tol=1e-8) {

  # Compute required values
  epsilon <- round(log(.Machine$double.xmin) + 10) # Logarithm of the smallest value that can be handled normally
  d <- ncol(data)
  n_obs <- nrow(data)
  T_obs <- n_obs - p
  to_return <- match.arg(to_return)

  # Collect parameter values
  model <- match.arg(model)
  parametrization <- match.arg(parametrization)
  check_same_means(parametrization=parametrization, same_means=same_means)
  params <- reform_constrained_pars(p=p, M=M, d=d, params=params, model=model, constraints=constraints,
                                    same_means=same_means, weight_constraints=weight_constraints,
                                    structural_pars=structural_pars) # All constraints are expanded and removed from the parameter vector
  W_constraints <- structural_pars$W
  structural_pars <- get_unconstrained_structural_pars(structural_pars=structural_pars)
  if(parametrization == "intercept") {
    all_phi0 <- pick_phi0(p=p, M=M, d=d, params=params, structural_pars=structural_pars)
  } else {
    mu <- pick_phi0(p=p, M=M, d=d, params=params, structural_pars=structural_pars) # mean parameters instead of phi0
  }
  all_A <- pick_allA(p=p, M=M, d=d, params=params, structural_pars=structural_pars) # A_{m,i}, m=1,...,M, i=1,..,p
  all_Omega <- pick_Omegas(p=p, M=M, d=d, params=params, structural_pars=structural_pars) # Omega_m
  all_boldA <- form_boldA(p=p, M=M, d=d, all_A=all_A) # The 'bold A' for each m=1,..,M, Lütkepohl 2005, eq.(2.1.8)
  alphas <- pick_alphas(p=p, M=M, d=d, params=params, model=model) # Mixing weight parameters
  all_df <- pick_df(M=M, params=params, model=model) # Degrees of freedom parameters (numeric(0) for GMVAR models)

  # Check that the parameter vector lies in the parameter space (excluding indentifiability)
  if(check_params) {
    if(!in_paramspace_int(p=p, M=M, d=d, params=params, model=model, all_boldA=all_boldA, alphas=alphas, all_Omega=all_Omega,
                          W_constraints=W_constraints, stat_tol=stat_tol, posdef_tol=posdef_tol, df_tol=df_tol)) {
      return(minval)
    }
  }

  if(model == "GMVAR") {
    M1 <- M
    M2 <- 0
  } else if(model == "StMVAR") {
    M1 <- 0
    M2 <- M
  } else { # model == "G-StMVAR"
    M1 <- M[1]
    M2 <- M[2]
  }
  M <- sum(M) # The total number of mixture components


  # i:th row denotes the vector \bold{y_{i-1}} = (y_{i-1},...,y_{i-p}) (dpx1),
  # assuming the observed data is y_{-p+1},...,y_0,y_1,...,y_{T}
  Y <- reform_data(data, p)

  # Calculate expected values (column per component) or phi0-parameters if using mean-parametrization
  Id <- diag(nrow=d)
  if(parametrization == "intercept") {
    mu <- vapply(1:M, function(m) solve(Id - rowSums(all_A[, , , m, drop=FALSE], dims=2), all_phi0[,m]), numeric(d)) # rowSums: sum over dims+1=3
  } else {
    all_phi0 <- vapply(1:M, function(m) (Id - rowSums(all_A[, , , m, drop=FALSE], dims=2))%*%mu[,m], numeric(d))
  }

  # Calculate the covariance matrices Sigma_{m,p} (Lutkepohl 2005, eq. (2.1.39) or the algorithm proposed by McElroy 2017)
  Sigmas <- get_Sigmas(p=p, M=M, d=d, all_A=all_A, all_boldA=all_boldA, all_Omega=all_Omega) # Store the (dpxdp) covariance matrices
  inv_Sigmas <- array(dim=c(d*p, d*p, M)) # Only used for StMVAR type regimes
  chol_Sigmas <- array(dim=c(d*p, d*p, M))
  for(m in 1:M) {
    # Take Cholesky here to avoid unnecessary warnings from mvnfast::dmvn, also used in inverting for m > M1:
    chol_Sigmas[, , m] <- chol(Sigmas[, , m])
    if(m > M1) {
      inv_Sigmas[, , m] <- chol2inv(chol_Sigmas[, , m])
    }
  }

  # Calculate the dp-dimensional multinormal densities (KMS 2016, eq.(6)) or log Students t densities (Virolainen 2022, eq. (3.4)):
  # i:th row for index i-1 etc, m:th column for m:th component.
  # We calculate in logarithm because the non-log values may be too close to zero for machine accuracy (if they are too close to zero
  # for all regimes and computer handles them as zero, we would divide by zero when calculating the mixing weights)
  log_mvdvalues <- matprods <- matrix(nrow=n_obs - p + 1, ncol=M) # The quadratic forms in Student's t density
  if(M1 > 0) { # Multinormals
    log_mvdvalues[,1:M1] <- vapply(1:M1, function(m) mvnfast::dmvn(X=Y, mu=rep(mu[,m], p), sigma=chol_Sigmas[, , m],
                                                                   log=TRUE, ncores=1, isChol=TRUE), numeric(T_obs + 1))
  }
  if(M2 > 0) { # Multistudents
    for(m in (M1 + 1):M) { # Go through the StMVAR type regimes
      tmp_mat <- Y - matrix(rep(mu[,m], p), nrow=nrow(Y), ncol=ncol(Y), byrow=TRUE)
      matprods[,m] <- rowSums(tmp_mat%*%inv_Sigmas[, , m]*tmp_mat)
      logC <- lgamma(0.5*(d*p + all_df[m - M1])) - 0.5*d*p*log(base::pi) - 0.5*d*p*log(all_df[m - M1] - 2) - lgamma(0.5*all_df[m - M1])
      log_det_Sigma <- 2*log(prod(diag(chol_Sigmas[, , m])))  #log(det(Sigmas[, , m]))
      log_mvdvalues[,m] <- logC - 0.5*log_det_Sigma - 0.5*(d*p + all_df[m - M1])*log(1 + matprods[,m]/(all_df[m - M1] - 2))
    }
  }

  ## Calculate the mixing weights alpha_{m,t} (KMS 2016, eq.(7))
  if(to_return != "mw_tplus1") {
    log_mvdvalues <- log_mvdvalues[1:T_obs, , drop=FALSE] # alpha_mt uses y_{t-1} so the last row is not needed
  }
  alpha_mt_and_l_0 <- get_alpha_mt(M=M, log_mvdvalues=log_mvdvalues, alphas=alphas,
                                   epsilon=epsilon, conditional=conditional, also_l_0=TRUE)
  alpha_mt <- alpha_mt_and_l_0$alpha_mt
  l_0 <- alpha_mt_and_l_0$l_0 # The first term in the exact log-likelihood function (=0 for conditional)

  if(to_return == "mw" | to_return == "mw_tplus1") {
    return(alpha_mt)
  }

  # Calculate the conditional means mu_{m,t} (KMS 2016, Condition 1 (a))
  # The dimensions of mu_mt will be: [t, p, m]
  all_A2 <- array(all_A, dim=c(d, d*p, M)) # cbind coefficient matrices of each component: m:th component is obtained at [, , m]
  Y2 <- Y[1:T_obs,] # Last row is not needed because mu_mt uses lagged values
  mu_mt <- array(vapply(1:M, function(m) t(all_phi0[, m] + tcrossprod(all_A2[, , m], Y2)), numeric(d*T_obs)), dim=c(T_obs, d, M)) # [, , m]

  # For the StMVAR type regimes, calculate the time-varying ARC-type scalar that multiplies the conditional covariance matrix:
  # i:th row calculated from y_{i-1} observation when the indexing is y_{-p+1},..,y_0,y_1,...,y_T - note that the last row of
  # arch_scalars is for time T+1
  arch_scalars <- matrix(1, nrow=nrow(matprods) - 1, ncol=M) # The arch-scalars are all 1 for the GMVAR type regimes
  if(M2 > 0) {
    # The last row is not needed because for the time t, we use the previous value y_{t-1}:
    matprods0 <- as.matrix(matprods[1:(nrow(matprods) - 1), (M1 + 1):M])
    tmp_mat <- matrix(all_df, nrow=nrow(matprods0), ncol=ncol(matprods0), byrow=TRUE)
    arch_scalars[,(M1 + 1):M] <- (tmp_mat - 2 + matprods0)/(tmp_mat - 2 + d*p)
  }


  if(to_return == "regime_cmeans") {
    return(mu_mt)
  } else if(to_return == "arch_scalars") {
    return(arch_scalars)
  } else if(to_return == "total_cmeans") { # KMS 2016, eq.(3)
    return(matrix(rowSums(vapply(1:M, function(m) alpha_mt[,m]*mu_mt[, , m], numeric(d*T_obs))), nrow=T_obs, ncol=d, byrow=FALSE))
  } else if(to_return == "regime_ccovs") {
    regime_ccovs <- array(vapply(1:M, function(m) t(matrix(all_Omega[, , m], nrow=nrow(arch_scalars), ncol=d^2, byrow=TRUE)*arch_scalars[, m]),
                                 numeric(d*d*nrow(arch_scalars))), dim=c(d, d, nrow(arch_scalars), M))
    return(regime_ccovs)
  } else if(to_return == "total_ccovs") { # KMS 2016, eq.(4)
    first_term <- array(rowSums(vapply(1:M, function(m) rep(alpha_mt[, m]*arch_scalars[, m], each=d*d)*as.vector(all_Omega[, , m]),
                                       numeric(d*d*T_obs))), dim=c(d, d, T_obs))
    sum_alpha_mu <- matrix(rowSums(vapply(1:M, function(m) alpha_mt[, m]*mu_mt[, , m],
                                          numeric(d*T_obs))), nrow=T_obs, ncol=d, byrow=FALSE)
    second_term <- array(rowSums(vapply(1:M, function(m) rep(alpha_mt[, m],
                                                             each=d*d)*as.vector(vapply(1:nrow(alpha_mt),
                                                                                        function(i1) tcrossprod((mu_mt[, , m] - sum_alpha_mu)[i1,]),
                                                                                        numeric(d*d))), numeric(d*d*T_obs))), dim=c(d, d, T_obs))
    return(first_term + second_term)
  }

  # Calculate the second term of the log-likelihood (KMS 2016 eq.(10)), also see Virolainen (2022) for the StMVAR type regimes
  dat <- data[(p + 1):n_obs,] # Initial values are not used here (conditional means and variances are already calculated)
  mvd_vals <- matrix(nrow=T_obs, ncol=M)
  if(M1 > 0) { # GMVAR type regimes
    mvd_vals[,1:M1] <- vapply(1:M1, function(m) mvnfast::dmvn(X=dat - mu_mt[, , m], mu=rep(0, times=d), sigma=all_Omega[, , m],
                                                              log=FALSE, ncores=1, isChol=FALSE), numeric(T_obs))
  }
  if(M2 > 0) { # StMVAR type regimes
    for(m in (M1 + 1):M) {
      df_m <- all_df[m - M1] + d*p # Degrees of freedom in regime m
      chol_Omega_m <- chol(all_Omega[, , m]) # Faster determinant and matrix inversion useing Cholesky decomposition
      det_Omega_m <- prod(diag(chol_Omega_m))^2
      inv_Omega_m <- chol2inv(chol_Omega_m)
      # Below, we calculate the d-dimensional conditional t-densities for the regime m, for t=1,...,T.
      tmp_mat <- dat - mu_mt[, , m]
      mvd_vals[, m] <- exp(lgamma(0.5*(d + df_m)) - 0.5*d*log(base::pi) - 0.5*d*log(df_m - 2) - lgamma(0.5*df_m))*
        arch_scalars[, m]^(-d/2)*1/sqrt(det_Omega_m)*(1 + rowSums(tmp_mat%*%inv_Omega_m*tmp_mat)/(arch_scalars[, m]*(df_m - 2)))^(-0.5*(d + df_m))
    }
  }

  weighted_mvd <- rowSums(alpha_mt*mvd_vals)
  weighted_mvd[weighted_mvd == 0] <- exp(epsilon)
  l_t <- log(weighted_mvd)
  loglik <- l_0 + sum(l_t)

  if(to_return == "terms") {
     return(l_t)
  } else if(to_return == "loglik_and_mw") {
    return(list(loglik=loglik,
                mw=alpha_mt))
  } else if(to_return == "loglik_mw_archscalars") {
    return(list(loglik=loglik,
                mw=alpha_mt,
                arch_scalars=arch_scalars))
  } else { # to_return == "loglik"
    return(loglik)
  }
}


#' @title Get mixing weights alpha_mt (this function is for internal use)
#'
#' @description \code{get_alpha_mt} computes the mixing weights based on
#'   the logarithm of the multivariate normal densities in the definition of
#'   the mixing weights.
#'
#' @inheritParams loglikelihood_int
#' @param log_mvdvalues \eqn{T x M} matrix containing the log multivariate normal densities.
#' @param alphas \eqn{M x 1} vector containing the mixing weight pa
#' @param epsilon the smallest number such that its exponent is wont classified as numerically zero
#'   (around \code{-698} is used).
#' @param also_l_0 return also l_0 (the first term in the exact log-likelihood function)?
#' @details Note that we index the time series as \eqn{-p+1,...,0,1,...,T} as in Kalliovirta et al. (2016).
#' @return Returns the mixing weights a matrix of the same dimension as \code{log_mvdvalues} so
#'   that the t:th row is for the time point t and m:th column is for the regime m.
#' @inherit in_paramspace_int references
#' @seealso \code{\link{loglikelihood_int}}
#' @keywords internal

get_alpha_mt <- function(M, log_mvdvalues, alphas, epsilon, conditional, also_l_0=FALSE) {
  if(M == 1) {
    if(!is.matrix(log_mvdvalues)) log_mvdvalues <- as.matrix(log_mvdvalues) # Possibly many time points but only one regime
    alpha_mt <- as.matrix(rep(1, nrow(log_mvdvalues)))
  } else {
    if(!is.matrix(log_mvdvalues)) log_mvdvalues <- t(as.matrix(log_mvdvalues)) # Only one time point but multiple regimes

    log_mvdvalues_orig <- log_mvdvalues
    small_logmvns <- log_mvdvalues < epsilon
    if(any(small_logmvns)) {
      # If too small or large non-log-density values are present (i.e., that would yield -Inf or Inf),
      # we replace them with ones that are not too small or large but imply the same mixing weights
      # up to negligible numerical tolerance.
      which_change <- rowSums(small_logmvns) > 0 # Which rows contain too small  values
      to_change <- log_mvdvalues[which_change, , drop=FALSE]
      largest_vals <- do.call(pmax, split(to_change, f=rep(1:ncol(to_change), each=nrow(to_change)))) # The largest values of those rows
      diff_to_largest <- to_change - largest_vals # Differences to the largest value of the row

      # For each element in each row, check the (negative) distance from the largest value of the row. If the difference
      # is smaller than epsilon, replace the with epsilon. The results are then the new log_mvn values.
      diff_to_largest[diff_to_largest < epsilon] <- epsilon

      # Replace the old log_mvdvalues with the new ones
      log_mvdvalues[which_change,] <- diff_to_largest
    }

    mvnvalues <- exp(log_mvdvalues)
    denominator <- as.vector(mvnvalues%*%alphas)
    alpha_mt <- (mvnvalues/denominator)%*%diag(alphas)
  }
  if(!also_l_0) {
    return(alpha_mt)
  } else {
    # First term of the exact log-likelihood (Kalliovirta et al. 2016, eq.(9))
    l_0 <- 0
    if(M == 1 && conditional == FALSE) {
      l_0 <- log_mvdvalues[1,]
    } else if(M > 1 && conditional == FALSE) {
      if(any(log_mvdvalues_orig[1,] < epsilon)) { # Need to use Brobdingnag
        l_0 <- log(Reduce("+", lapply(1:M, function(i1) alphas[i1]*exp(Brobdingnag::as.brob(log_mvdvalues_orig[1, i1])))))
      } else {
        l_0 <- log(sum(alphas*mvnvalues[1,]))
      }
    }
    return(list(alpha_mt=alpha_mt,
                l_0=l_0))
  }
}



#' @title Compute log-likelihood of a GMVAR, StMVAR, or G-StMVAR model using parameter vector
#'
#' @description \code{loglikelihood} computes log-likelihood of a GMVAR, StMVAR, or G-StMVAR model using parameter vector
#'   instead of an object of class 'gsmvar'. Exists for convenience if one wants to for example
#'   employ other estimation algorithms than the ones used in \code{fitGSMVAR}. Use \code{minval} to
#'   control what happens when the parameter vector is outside the parameter space.
#'
#' @inheritParams loglikelihood_int
#' @return Returns log-likelihood if \code{params} is in the parameters space and \code{minval} if not.
#' @inherit loglikelihood_int details references
#' @seealso \code{\link{fitGSMVAR}}, \code{\link{GSMVAR}}, \code{\link{calc_gradient}}
#' @examples
#' # GMVAR(2, 2), d=2 model;
#' params22 <- c(0.36, 0.121, 0.223, 0.059, -0.151, 0.395, 0.406, -0.005,
#'  0.083, 0.299, 0.215, 0.002, 0.03, 0.484, 0.072, 0.218, 0.02, -0.119,
#'  0.722, 0.093, 0.032, 0.044, 0.191, 1.101, -0.004, 0.105, 0.58)
#' loglikelihood(data=gdpdef, p=2, M=2, params=params22)
#'
#' # Structural GMVAR(2, 2), d=2 model identified with sign-constraints:
#' params22s <- c(0.36, 0.121, 0.484, 0.072, 0.223, 0.059, -0.151, 0.395,
#'  0.406, -0.005, 0.083, 0.299, 0.218, 0.02, -0.119, 0.722, 0.093, 0.032,
#'  0.044, 0.191, 0.057, 0.172, -0.46, 0.016, 3.518, 5.154, 0.58)
#' W_22 <- matrix(c(1, 1, -1, 1), nrow=2, byrow=FALSE)
#' loglikelihood(data=gdpdef, p=2, M=2, params=params22s, structural_pars=list(W=W_22))
#' @export

loglikelihood <- function(data, p, M, params, model=c("GMVAR", "StMVAR", "G-StMVAR"), conditional=TRUE, parametrization=c("intercept", "mean"),
                          constraints=NULL, same_means=NULL, weight_constraints=NULL, structural_pars=NULL, minval=NA,
                          stat_tol=1e-3, posdef_tol=1e-8, df_tol=1e-8) {
  parametrization <- match.arg(parametrization)
  model <- match.arg(model)
  check_same_means(parametrization=parametrization, same_means=same_means)
  data <- check_data(data, p)
  d <- ncol(data)
  check_pMd(p=p, M=M, d=d, model=model)
  check_constraints(p=p, M=M, d=d, constraints=constraints, same_means=same_means, weight_constraints=weight_constraints,
                    structural_pars=structural_pars)
  if(length(params) != n_params(p=p, M=M, d=d, model=model, constraints=constraints, same_means=same_means,
                                weight_constraints=weight_constraints, structural_pars=structural_pars)) {
    stop("Parameter vector has wrong dimension")
  }
  loglikelihood_int(data=data, p=p, M=M, params=params, model=model, conditional=conditional, parametrization=parametrization,
                    constraints=constraints, same_means=same_means, weight_constraints=weight_constraints,
                    structural_pars=structural_pars, to_return="loglik", check_params=TRUE, minval=minval,
                    stat_tol=stat_tol, posdef_tol=posdef_tol, df_tol=df_tol)
}


#' @title Compute conditional moments of a GMVAR, StMVAR, or G-StMVAR model
#'
#' @description \code{cond_moments} compute conditional regimewise means, conditional means, and conditional covariance matrices
#'  of a GMVAR, StMVAR, or G-StMVAR model.
#'
#' @inheritParams loglikelihood_int
#' @param to_return should the regimewise conditional means, total conditional means, or total conditional covariance matrices
#'   be returned?
#' @details The first p values are used as the initial values, and by conditional we mean conditioning on the past. Formulas
#'   for the conditional means and covariance matrices are given in equations (3) and (4) of KMS (2016).
#' @return
#'  \describe{
#'   \item{If \code{to_return=="regime_cmeans"}:}{an \code{[T-p, d, M]} array containing the regimewise conditional means
#'     (the first p values are used as the initial values).}
#'   \item{If \code{to_return=="regime_ccovs"}:}{an \code{[d, d, T-p, M]} array containing the regimewise conditional
#'     covariance matrices (the first p values are used as the initial values). The index \code{[ , , t, m]} gives the time
#'     \code{t} conditional covariance matrix for the regime \code{m}.}
#'   \item{If \code{to_return=="total_cmeans"}:}{a \code{[T-p, d]} matrix containing the conditional means of the process
#'     (the first p values are used as the initial values).}
#'   \item{If \code{to_return=="total_ccov"}:}{an \code{[d, d, T-p]} array containing the conditional covariance matrices of the process
#'     (the first p values are used as the initial values).}
#'   \item{If \code{to_return=="arch_scalars"}:}{a \code{[T-p, M]} matrix containing the regimewise arch scalars
#'    multiplying error term covariance matrix in the conditional covariance matrix of the regime. For GMVAR type regimes, these
#'    are all ones (the first p values are used as the initial values).}
#'  }
#' @inherit loglikelihood_int references
#' @family moment functions
#' @examples
#' # GMVAR(2, 2), d=2 model;
#' params22 <- c(0.36, 0.121, 0.223, 0.059, -0.151, 0.395, 0.406, -0.005,
#'  0.083, 0.299, 0.215, 0.002, 0.03, 0.484, 0.072, 0.218, 0.02, -0.119,
#'  0.722, 0.093, 0.032, 0.044, 0.191, 1.101, -0.004, 0.105, 0.58)
#' cond_moments(data=gdpdef, p=2, M=2, params=params22, to_return="regime_cmeans")
#' cond_moments(data=gdpdef, p=2, M=2, params=params22, to_return="total_cmeans")
#' cond_moments(data=gdpdef, p=2, M=2, params=params22, to_return="total_ccovs")
#' @export

cond_moments <- function(data, p, M, params, model=c("GMVAR", "StMVAR", "G-StMVAR"), parametrization=c("intercept", "mean"),
                         constraints=NULL, same_means=NULL, weight_constraints=NULL, structural_pars=NULL,
                         to_return=c("regime_cmeans", "regime_ccovs", "total_cmeans", "total_ccovs", "arch_scalars"),
                         minval=NA, stat_tol=1e-3, posdef_tol=1e-8, df_tol=1e-8) {
  parametrization <- match.arg(parametrization)
  model <- match.arg(model)
  to_return <- match.arg(to_return)
  check_same_means(parametrization=parametrization, same_means=same_means)
  data <- check_data(data, p)
  d <- ncol(data)
  check_pMd(p=p, M=M, d=d, model=model)
  check_constraints(p=p, M=M, d=d, constraints=constraints, same_means=same_means,
                    weight_constraints=weight_constraints, structural_pars=structural_pars)
  if(length(params) != n_params(p=p, M=M, d=d, model=model, constraints=constraints, same_means=same_means,
                                weight_constraints=weight_constraints, structural_pars=structural_pars)) {
    stop("Parameter vector has wrong dimension")
  }
  loglikelihood_int(data=data, p=p, M=M, params=params, model=model, parametrization=parametrization,
                    constraints=constraints, same_means=same_means, weight_constraints=weight_constraints,
                    structural_pars=structural_pars, to_return=to_return, check_params=TRUE, minval=minval,
                    stat_tol=stat_tol, posdef_tol=posdef_tol, df_tol=df_tol)
}


#' @title Calculate AIC, HQIC, and BIC
#'
#' @description \code{get_IC} calculates the information criteria values AIC, HQIC, and BIC.
#'
#' @param loglik log-likelihood value
#' @param npars number of (freely estimated) parameters in the model
#' @param obs numbers of observations with starting values excluded for conditional models.
#' @details Note that for conditional models with different autoregressive order p the
#'  information criteria values are \strong{NOT} comparable.
#' @return Returns a data frame containing the information criteria values.
#' @keywords internal

get_IC <- function(loglik, npars, obs) {
  AIC <- -2*loglik + 2*npars
  HQIC <- -2*loglik + 2*npars*log(log(obs))
  BIC <- -2*loglik + npars*log(obs)
  data.frame(AIC=AIC, HQIC=HQIC, BIC=BIC)
}




#' @title Calculate multivariate Pearson residuals of a GMVAR, StMVAR, or G-StMVAR model
#'
#' @description \code{Pearson_residuals} calculates multivariate Pearson residuals for a GMVAR, StMVAR, or G-StMVAR model.
#'
#' @inheritParams quantile_residual_tests
#' @param standardize Should the residuals be standardized? Use \code{FALSE} to obtain raw residuals.
#' @return Returns \eqn{((n_obs-p) x d)} matrix containing the residuals,
#'   \eqn{j}:th column corresponds to the time series in the \eqn{j}:th column of the data.
#' @inherit GSMVAR references
#' @seealso \code{\link{fitGSMVAR}}, \code{\link{GSMVAR}}, \code{\link{quantile_residuals}},
#'   \code{\link{diagnostic_plot}}
#' @examples
#' # GMVAR(1,2), d=2 model:
#' params12 <- c(0.55, 0.112, 0.344, 0.055, -0.009, 0.718, 0.319, 0.005, 0.03,
#'  0.619, 0.173, 0.255, 0.017, -0.136, 0.858, 1.185, -0.012, 0.136, 0.674)
#' mod12 <- GSMVAR(gdpdef, p=1, M=2, params=params12)
#' Pearson_residuals(mod12, standardize=FALSE) # Raw residuals
#' Pearson_residuals(mod12, standardize=TRUE) # Standardized to identity cov.matrix.
#'
#' # Structural GMVAR(2, 2), d=2 model identified with sign-constraints:
#' params22s <- c(0.36, 0.121, 0.484, 0.072, 0.223, 0.059, -0.151, 0.395,
#'  0.406, -0.005, 0.083, 0.299, 0.218, 0.02, -0.119, 0.722, 0.093, 0.032,
#'  0.044, 0.191, 0.057, 0.172, -0.46, 0.016, 3.518, 5.154, 0.58)
#' W_22 <- matrix(c(1, 1, -1, 1), nrow=2, byrow=FALSE)
#' mod22s <- GSMVAR(gdpdef, p=2, M=2, params=params22s, structural_pars=list(W=W_22))
#' Pearson_residuals(mod22s, standardize=FALSE) # Raw residuals
#' Pearson_residuals(mod22s, standardize=TRUE) # Standardized to identity cov.matrix.
#' @export

Pearson_residuals <- function(gsmvar, standardize=TRUE) {
  # Checks, etc
  gsmvar <- gmvar_to_gsmvar(gsmvar) # Backward compatibility
  check_gsmvar(gsmvar)
  p <- gsmvar$model$p
  M <- gsmvar$model$M
  d <- gsmvar$model$d
  data <- gsmvar$data
  n_obs <- nrow(data)
  T_obs <- n_obs - p

  # Conditional means
  mu_t <- loglikelihood_int(data=data, p=gsmvar$model$p, M=gsmvar$model$M, params=gsmvar$params,
                            model=gsmvar$model$model, conditional=gsmvar$model$conditional,
                            parametrization=gsmvar$model$parametrization,
                            constraints=gsmvar$model$constraints,
                            same_means=gsmvar$model$same_means,
                            weight_constraints=gsmvar$model$weight_constraints,
                            structural_pars=gsmvar$model$structural_pars,
                            to_return="total_cmeans", check_params=TRUE,
                            stat_tol=gsmvar$num_tols$stat_tol,
                            posdef_tol=gsmvar$num_tols$posdef_tol,
                            df_tol=gsmvar$num_tols$df_tol)

  # Nonstandardized residuals
  y_minus_mu <- data[(p + 1):nrow(data),] - mu_t # [T_obs, d]
  if(!standardize) {
    return(y_minus_mu)
  }

  # Conditional covariance matrices
  Omega_t <- loglikelihood_int(data=data, p=gsmvar$model$p, M=gsmvar$model$M, params=gsmvar$params,
                               model=gsmvar$model$model, conditional=gsmvar$model$conditional,
                               parametrization=gsmvar$model$parametrization,
                               constraints=gsmvar$model$constraints,
                               same_means=gsmvar$model$same_means,
                               weight_constraints=gsmvar$model$weight_constraints,
                               structural_pars=gsmvar$model$structural_pars,
                               to_return="total_ccovs", check_params=TRUE,
                               stat_tol=gsmvar$num_tols$stat_tol,
                               posdef_tol=gsmvar$num_tols$posdef_tol,
                               df_tol=gsmvar$num_tols$df_tol)

  all_residuals <- matrix(nrow=T_obs, ncol=d)

  # Calculate the Pearson residuals
  for(i1 in 1:T_obs) {
    all_residuals[i1,] <- solve(unvec(d=d, a=get_symmetric_sqrt(Omega_t[, , i1])), y_minus_mu[i1,])
  }
  all_residuals
}
