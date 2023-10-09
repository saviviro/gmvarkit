#' @title Estimate linear impulse response function based on a single regime of a structural GMVAR,
#'   StMVAR, or G-StMVAR model.
#'
#' @description \code{linear_IRF} estimates linear impulse response function based on a single regime
#'   of a structural GMVAR, StMVAR, or G-StMVAR model.
#'
#' @param gsmvar an object of class \code{'gsmvar'} defining a structural or reduced form
#'   GSMVAR model. For a reduced form model, the shocks are automatically identified by
#'   the lower triangular Cholesky decomposition.
#' @param which_shocks a numeric vector of length at most \eqn{d}
#'   (\code{=ncol(data)}) and elements in \eqn{1,...,d} specifying the
#'   structural shocks for which the IRF should be estimated.
#' @param N a positive integer specifying the horizon how far ahead should the
#'   linear impulse responses be calculated.
#' @param regime Based on which regime the linear IRF should be calculated?
#'   An integer in \eqn{1,...,M}.
#' @param which_cumulative a numeric vector with values in \eqn{1,...,d}
#'   (\code{d=ncol(data)}) specifying which the variables for which the linear impulse
#'   responses should be cumulative. Default is none.
#' @param scale should the linear IRFs to some of the shocks be scaled so that they
#'   correspond to a specific magnitude of instantaneous or peak response
#'   of some specific variable (see the argument \code{scale_type})?
#'   Provide a length three vector where the shock of interest
#'   is given in the first element (an integer in \eqn{1,...,d}), the variable of
#'   interest is given in the second element (an integer in \eqn{1,...,d}), and
#'   the magnitude of its instantaneous or peak response in the third element
#'   (a non-zero real number). If the linear IRFs of multiple shocks should be scaled,
#'   provide a matrix which has one column for each of the shocks with the columns being
#'   the length three vectors described above.
#' @param scale_type If argument \code{scale} is specified, should the linear IRFs be
#'   scaled to match an instantaneous response (\code{"instant"}) or peak response
#'   (\code{"peak"}). If \code{"peak"}, the scale is based on the largest magnitude
#'   of peak response in absolute value. Ignored if \code{scale} is not specified.
#' @param scale_horizon If \code{scale_type == "peak"} what the maximum horizon up
#'   to which peak response is expected? Scaling won't based on values after this.
#' @param ci a numeric vector with elements in \eqn{(0, 1)} specifying the
#'   confidence levels of the confidence intervals calculated via a bootstrap
#'   method, see the details section.
#' @param seed a length one numeric vector initializing the seed for the random generator.
#' @param ... parameters passed to the plot method \code{plot.irf} that plots
#'   the results.
#' @details The model DOES NOT need to be structural in order for this function to be
#'   applicable. When an identified structural GMVAR, StMVAR, or G-StMVAR model is
#'   provided in the argument \code{gsmvar}, the identification imposed by the model
#'   is used. When a reduced form model is provided in the argument \code{gsmvar},
#'   lower triangular Cholesky identification is used to identify the shocks.
#'
#'   FILL IN DETAILS ABOUT CONFIDENCE INTERVALS
#' @return Returns a class \code{'irf'} list with the linear IRFs in ... FILL IN!
#' @seealso \code{\link{GIRF}}, \code{\link{GFEVD}}, \code{\link{fitGSMVAR}}, \code{\link{GSMVAR}},
#'   \code{\link{gsmvar_to_sgsmvar}}, \code{\link{reorder_W_columns}}, \code{\link{swap_W_signs}}
#' @references
#'  \itemize{
#'    \item Kilian L. and Lütkepohl H. 2017. Structural Vectors Autoregressive Analysis.
#'          \emph{Cambridge University Press}, Cambridge.
#'  }
#' @examples
#'  # FILL IN
#' @export

linear_IRF <- function(gsmvar, which_shocks, N=30, regime=1, which_cumulative=numeric(0),
                       scale=NULL, scale_type=c("instant", "peak"), scale_horizon=N,
                       ci=c(0.95, 0.80), seed=NULL, ...) {
  # Get the parameter values etc
  stopifnot(all_pos_ints(c(N, regime)))
  p <- gsmvar$model$p
  M_orig <- gsmvar$model$M
  M <- sum(M_orig)
  stopifnot(regime <= M)
  d <- gsmvar$model$d
  model <- gsmvar$model$model
  constraints <- gsmvar$model$constraints
  same_means <- gsmvar$model$same_means
  structural_pars <- gsmvar$model$structural_pars
  params <- gsmvar$params
  params <- reform_constrained_pars(p=p, M=M_orig, d=d, params=params, model=model,
                                    constraints=gsmvar$model$constraints,
                                    same_means=gsmvar$model$same_means,
                                    structural_pars=structural_pars)
  structural_pars <- get_unconstrained_structural_pars(structural_pars=structural_pars)
  if(gsmvar$model$parametrization == "mean") {
    params <- change_parametrization(p=p, M=M, d=d, params=params, constraints=NULL,
                                     structural_pars=structural_pars, change_to="intercept")
  }
  all_mu <- get_regime_means(gsmvar)
  all_phi0 <- pick_phi0(p=p, M=M_orig, d=d, params=params, structural_pars=structural_pars)
  all_A <- pick_allA(p=p, M=M_orig, d=d, params=params, structural_pars=structural_pars)
  all_Omega <- pick_Omegas(p=p, M=M_orig, d=d, params=params, structural_pars=structural_pars)
  all_boldA <- form_boldA(p=p, M=M_orig, d=d, all_A=all_A)
  alphas <- pick_alphas(p=p, M=M_orig, d=d, params=params, model=model)
  all_df <- pick_df(M=M_orig, params=params, model=model)

  # Obtain the impact matrix of the regime the IRF is to calculated for
  if(!is.null(structural_pars)) { # Shocks identified by heteroskedasticity
    W <- pick_W(p=p, M=M_orig, d=d, params=params, structural_pars=structural_pars)
    if(regime > 1) { # Include lambdas
      lambdas <- matrix(pick_lambdas(p=p, M=M_orig, d=d, params=params, structural_pars=structural_pars), nrow=d, byrow=FALSE)
      Lambda_m <- diag(lambdas[, regime - 1])
      B_matrix <- W%*%sqrt(Lambda_m)
    } else { # regime == 1
      B_matrix <- W
    }
  } else { # Shocks identified by lower-triangular Cholesky decomposition
    B_matrix <- t(chol(all_Omega[, , regime]))
  }

  # Calculate the impulse response functions
  boldA <- all_boldA[, , regime] # Companion form AR matrix for the selected regime
  J_matrix <- create_J_matrix(d=d, p=p)
  all_boldA_powers <- array(NA, dim=c(d*p, d*p, N+1)) # The first [, , i1] is for the impact period, i1+1 for period i1
  all_phi_i <- array(NA, dim=c(d*p, d*p, N+1)) # JA^iJ' matrices; [, , 1] is for the impact period, i1+1 for period i1
  all_Phi_i <- array(NA, dim=c(d*p, d*p, N+1)) # IR-matrices [, , 1] is for the impact period, i1+1 for period i1
  for(i1 in 1:(N + 1)) { # Go through the periods, i1=1 for the impact period, i1+1 for the period i1 after the impact
    if(i1 == 1) {
      all_boldA_powers[, , i1] <- diag(d*p) # Diagonal matrix for the power 0
    } else {
      all_boldA_powers[, , i1] <- all_boldA_powers[, , i1 - 1]%*%boldA # boldA^{i1-1} because i1=1 is for the zero period
    }
    all_phi_i[, , i1] <- J_matrix%*%all_boldA_powers[, , i1]%*%t(J_matrix)
    all_Phi_i[, , i1] <- all_phi_i[, , i1]%*%B_matrix
  }
  # all_Phi_i[variable, shock, horizon] -> all_Phi_i[variable, shock, ] subsets the IRF!

  # NOTE which cumulative does not do anything yet! Maybe original + cumulative for all?
  # Depends on how the confidence bounds are created!

  # Maybe also remove "which_shocks"? IRF calculated for all anyway - are there implications to the ci?
  # All the IRFs are calculate for UNIT SHOCKS i.e., one standard error struct shocks
  # -> responses can be scaled without loss of generality as the IRFs are symmetric w.r.t the size of the shock

  # Return the results
  structure(list(all_irfs=all_Phi_i,
                 N=N,
                 ci=ci,
                 which_cumulative=which_cumulative,
                 seed=seed,
                 gsmvar=gsmvar),
            class="irf")
}
