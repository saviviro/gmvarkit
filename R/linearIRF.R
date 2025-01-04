#' @title Estimate linear impulse response function based on a single regime of a structural GMVAR,
#'   StMVAR, or G-StMVAR model.
#'
#' @description \code{linear_IRF} estimates linear impulse response function based on a single regime
#'   of a structural GMVAR, StMVAR, or G-StMVAR model.
#'
#' @param gsmvar an object of class \code{'gsmvar'} defining a structural or reduced form
#'   GSMVAR model. For a reduced form model, the shocks are automatically identified by
#'   the lower triangular Cholesky decomposition.
#' @param N a positive integer specifying the horizon how far ahead should the
#'   linear impulse responses be calculated.
#' @param regime Based on which regime the linear IRF should be calculated?
#'   An integer in \eqn{1,...,M}.
#' @param which_cumulative a numeric vector with values in \eqn{1,...,d}
#'   (\code{d=ncol(data)}) specifying which the variables for which the linear impulse
#'   responses should be cumulative. Default is none.
#' @param scale should the linear IRFs to some of the shocks be scaled so that they
#'   correspond to a specific instantaneous response of some specific
#'   variable? Provide a length three vector where the shock of interest
#'   is given in the first element (an integer in \eqn{1,...,d}), the variable of
#'   interest is given in the second element (an integer in \eqn{1,...,d}), and
#'   its instantaneous response in the third element (a non-zero real number).
#'   If the linear IRFs of multiple shocks should be scaled, provide a matrix which has one
#'   column for each of the shocks with the columns being the length three vectors described above.
#' @param ci a real number in \eqn{(0, 1)} specifying the confidence level of the
#'   confidence intervals calculated via a fixed-design wild residual bootstrap method.
#'   Available only for models that impose linear autoregressive dynamics
#'   (excluding changes in the volatility regime).
#' @param bootstrap_reps the number of bootstrap repetitions for estimating confidence bounds.
#' @param ncores the number of CPU cores to be used in parallel computing when bootstrapping confidence bounds.
#' @param ncalls on how many estimation rounds should each bootstrap estimation be based on?
#'   Does not have to be very large since initial estimates used are based on the initially fitted model.
#'   Larger number of rounds gives more reliable results but is computationally more demanding.
#' @param seeds a numeric vector of length \code{bootstrap_reps} initializing the seed for the random
#'   generator for each bootstrap replication.
#' @param ... parameters passed to the plot method \code{plot.irf} that plots
#'   the results.
#' @details The model DOES NOT need to be structural in order for this function to be
#'   applicable. When an identified structural GMVAR, StMVAR, or G-StMVAR model is
#'   provided in the argument \code{gsmvar}, the identification imposed by the model
#'   is used. When a reduced form model is provided in the argument \code{gsmvar},
#'   lower triangular Cholesky identification is used to identify the shocks.
#'
#'   If the autoregressive dynamics of the model are linear (i.e., either M == 1 or mean and AR parameters
#'   are constrained identical across the regimes), confidence bounds can be calculated based on a type of fixed-design
#'   wild residual bootstrap method. See Virolainen (forthcoming) for a related discussion. We employ the method described
#'   in Herwartz and Lütkepohl (2014); see also the relevant chapters in Kilian and Lütkepohl (2017).
#' @return Returns a class \code{'irf'} list with the following elements:
#'   \describe{
#'     \item{\code{$point_est}:}{a 3D array \code{[variables, shock, horizon]} containing the point estimates of the IRFs.
#'        Note that the first slice is for the impact responses and the slice i+1 for the period i. The response of the
#'        variable 'i1' to the shock 'i2' is subsetted as \code{$point_est[i1, i2, ]}.}
#'     \item{\code{$conf_ints}:}{bootstrapped confidence intervals for the IRFs in a \code{[variables, shock, horizon, bound]}
#'        4D array. The lower bound is obtained as \code{$conf_ints[, , , 1]}, and similarly the upper bound as
#'         \code{$conf_ints[, , , 2]}. The subsetted 3D array is then the bound in a form similar to \code{$point_est}.}
#'     \item{\code{$all_bootstrap_reps}:}{IRFs from all of the bootstrap replications in a \code{[variables, shock, horizon, rep]}.
#'        4D array. The IRF from replication i1 is obtained as \code{$all_bootstrap_reps[, , , i1]}, and the subsetted 3D array
#'        is then the in a form similar to \code{$point_est}.}
#'     \item{Other elements:}{contains some of the arguments the \code{linear_IRF} was called with.}
#'   }
#' @seealso \code{\link{GIRF}}, \code{\link{GFEVD}}, \code{\link{fitGSMVAR}}, \code{\link{GSMVAR}},
#'   \code{\link{gsmvar_to_sgsmvar}}, \code{\link{reorder_W_columns}}, \code{\link{swap_W_signs}}
#' @references
#'  \itemize{
#'    \item Herwartz H. and Lütkepohl H. 2014. Structural vector autoregressions with Markov switching:
#'      Combining conventional with statistical identification of shocks. \emph{Journal of Econometrics},
#'      183, pp. 104-116.
#'    \item Kilian L. and Lütkepohl H. 2017. Structural Vectors Autoregressive Analysis.
#'          \emph{Cambridge University Press}, Cambridge.
#'    \item Virolainen S. 2025. A statistically identified structural vector autoregression with endogenously
#'          switching volatility regime. \emph{Journal of Business & Economic Statistics}, \strong{43}, 1, 44-54.
#'  }
#' @examples
#'   \donttest{
#'  # These are long running examples that take a few minutes to run
#'
#'  ## GMVAR, p=5, M=2, d=2 model with linear AR dynamics.
#'  # recursive identification, IRF based on the first regime:
#'  params52cm <- c(0.788, 0.559, 0.277, 0.038, -0.061, 0.463, 0.286, 0,
#'                0.035, 0.161, -0.112, 0.031, -0.313, 0.183, 0.103, 0.014,
#'                0.002, 0.143, -0.089, -0.013, 0.182, -0.04, 1.3, 0.008,
#'                0.139, 0.277, -0.005, 0.032, 0.118)
#'  mod52cm <- GSMVAR(data=gdpdef, p=5, M=2, params=params52cm,
#'                    constraints=rbind(diag(5*2^2), diag(5*2^2)),
#'                    same_means=list(1:2), parametrization="mean")
#'  irf1 <- linear_IRF(mod52cm, regime=1, N=20, scale=cbind(c(1, 1, 1), c(2, 2, 1)))
#'  print(irf1, digits=3)
#'  plot(irf1)
#'
#'  # Identification by heteroskedasticity, bootstrapped confidence intervals and
#'  # and scaled instantaneous effects of the shocks. Note that in actual
#'  # empirical application, a larger number of bootstrap reps should be used.
#'  mod52cms <- gsmvar_to_sgsmvar(mod52cm)
#'  irf2 <- linear_IRF(mod52cms, regime=1, N=20, ci=0.68, bootstrap_reps=10,
#'                     ncalls=1, seeds=1:10, ncores=1)
#'  plot(irf2)
#'  }
#' @export

linear_IRF <- function(gsmvar, N=30, regime=1, which_cumulative=numeric(0),
                       scale=NULL, ci=NULL, bootstrap_reps=100, ncores=2, ncalls=1, seeds=NULL, ...) {
  # Get the parameter values etc
  stopifnot(all_pos_ints(c(N, regime, ncores, bootstrap_reps, ncalls)))
  if(!is.null(seeds)) stopifnot(length(seeds) == bootstrap_reps)
  p <- gsmvar$model$p
  M_orig <- gsmvar$model$M
  M <- sum(M_orig)
  stopifnot(regime <= M)
  d <- gsmvar$model$d
  model <- gsmvar$model$model
  constraints <- gsmvar$model$constraints
  same_means <- gsmvar$model$same_means
  weight_constraints <- gsmvar$model$weight_constraints
  structural_pars <- gsmvar$model$structural_pars
  data <- gsmvar$data
  n_obs <- nrow(data)
  T_obs <- n_obs - p

  # Check the argument scale and which_cumulative
  if(is.null(gsmvar$model$structural_pars)) { # Reduced form model
    B_constrs <- matrix(NA, nrow=d, ncol=d)
    B_constrs[upper.tri(B_constrs)] <- 0
    diag(B_constrs) <- 1 # Lower triangular Cholesky constraints
  } else {
    B_constrs <- gsmvar$model$structural_pars$W
  }
  if(!is.null(scale)) {
    scale <- as.matrix(scale)
    stopifnot(all(scale[1,] %in% 1:d)) # All shocks in 1,...,d
    stopifnot(length(unique(scale[1,])) == length(scale[1,])) # No duplicate scales for the same shock
    stopifnot(all(scale[2,] %in% 1:d)) # All variables in 1,...,d
    stopifnot(all(scale[3,] != 0)) # No zero initial magnitudes

    # For the considered shocks, check that there are not zero-constraints for the variable
    # whose initial response is scaled.
    for(i1 in 1:ncol(scale)) {
      if(!is.na(B_constrs[scale[2, i1], scale[1, i1]]) && B_constrs[scale[2, i1], scale[1, i1]] == 0) {
        if(is.null(is.null(gsmvar$model$structural_pars))) {
          stop(paste("Instantaneous response of the variable that has a zero constraint for",
                     "the considered shock cannot be scaled"))
        } else { # Reduced form
          stop(paste("Instantaneous response of the variable that has a zero constraint for the considered",
                     "shock cannot be scaled",
                     "(lower triangular recursive identification is assumed for reduced form models)"))
        }
      }
    }
  }
  if(length(which_cumulative) > 0) {
    which_cumulative <- unique(which_cumulative)
    stopifnot(all(which_cumulative %in% 1:d))
  }

  # Pick params etc
  params <- reform_constrained_pars(p=p, M=M_orig, d=d, params=gsmvar$params, model=model,
                                    constraints=constraints, same_means=same_means,
                                    weight_constraints=weight_constraints,
                                    structural_pars=structural_pars)
  structural_pars <- get_unconstrained_structural_pars(structural_pars=structural_pars)
  if(gsmvar$model$parametrization == "mean") {
    params <- change_parametrization(p=p, M=M, d=d, params=params, constraints=NULL, same_means=NULL,
                                     weight_constraints=NULL, structural_pars=structural_pars, change_to="intercept")
  }
  all_mu <- get_regime_means(gsmvar)
  all_phi0 <- pick_phi0(p=p, M=M_orig, d=d, params=params, structural_pars=structural_pars)
  all_A <- pick_allA(p=p, M=M_orig, d=d, params=params, structural_pars=structural_pars)
  all_Omega <- pick_Omegas(p=p, M=M_orig, d=d, params=params, structural_pars=structural_pars)
  all_boldA <- form_boldA(p=p, M=M_orig, d=d, all_A=all_A)
  alphas <- pick_alphas(p=p, M=M_orig, d=d, params=params, model=model)
  all_df <- pick_df(M=M_orig, params=params, model=model)
  all_lambdas <- pick_lambdas(p=p, M=M_orig, d=d, params=params,
                              structural_pars=structural_pars) # numeric(0) if reduced form model

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

  # Function to calculate IRF
  get_IRF <- function(p, d, N, boldA, B_matrix) {
    J_matrix <- create_J_matrix(d=d, p=p)
    all_boldA_powers <- array(NA, dim=c(d*p, d*p, N+1)) # The first [, , i1] is for the impact period, i1+1 for period i1
    all_Phi_i <- array(NA, dim=c(d, d, N+1)) # JA^iJ' matrices; [, , 1] is for the impact period, i1+1 for period i1
    all_Theta_i <- array(NA, dim=c(d, d, N+1)) # IR-matrices [, , 1] is for the impact period, i1+1 for period i1
    for(i1 in 1:(N + 1)) { # Go through the periods, i1=1 for the impact period, i1+1 for the period i1 after the impact
      if(i1 == 1) {
        all_boldA_powers[, , i1] <- diag(d*p) # Diagonal matrix for the power 0
      } else {
        all_boldA_powers[, , i1] <- all_boldA_powers[, , i1 - 1]%*%boldA # boldA^{i1-1} because i1=1 is for the zero period
      }
      all_Phi_i[, , i1] <- J_matrix%*%all_boldA_powers[, , i1]%*%t(J_matrix)
      all_Theta_i[, , i1] <- all_Phi_i[, , i1]%*%B_matrix
    }
    all_Theta_i # all_Theta_i[variable, shock, horizon] -> all_Theta_i[variable, shock, ] subsets the IRF!
  }

  ## Calculate the impulse response functions:
  point_est <- get_IRF(p=p, d=d, N=N,
                       boldA=all_boldA[, , regime], # boldA= Companion form AR matrix of the selected regime
                       B_matrix=B_matrix)
  dimnames(point_est)[[1]] <- colnames(gsmvar$data)
  dimnames(point_est)[[2]] <- paste("Shock", 1:gsmvar$model$d)

  ## Confidence bounds by fixed design wild residual bootstrap
  AR_mats_identical <- all(apply(all_boldA, MARGIN=3, FUN=function(x) identical(x, all_boldA[,,1])))
  means_identical <- !is.null(same_means) && length(same_means) == 1 && all(same_means[[1]] == 1:M)
  ci_possible <- (means_identical && AR_mats_identical) || M == 1

  if(!is.null(ci) && !ci_possible) {
    warning("Confidence bounds are not available as the autoregressive dynamics are not linear")
    all_bootstrap_IRF <- NULL
  } else if(!is.null(ci) && ci_possible) { # Bootstrap confidence bounds

    ## Create initial values for the two-phase estimation algorithm: does not vary across the bootstrap reps
    new_params <- gsmvar$params

    # For all models, bootstrapping conditions on the estimated mw parameters, so they need to be removed:
    if(is.null(weight_constraints) && M > 1) {
      new_params <- c(new_params[1:(length(new_params) - (M - 1) - length(all_df))], all_df) # Removes alphas
      new_weight_constraints <- alphas[-M]
    } else { # If weight constraints used or M == 1, no modifications are required
      new_weight_constraints <- weight_constraints
    }

    # For structural models, bootstrapping also conditions on the estimated lambda parameters to
    # keep the shocks in a fixed ordering (which is given for recursively identified models).
    # Also make sure that each column of W has a strict sign constraint: if not normalize diagonal elements to positive.
    if(!is.null(gsmvar$model$structural_pars)) {
      new_fixed_lambdas <- all_lambdas # If fixed_lambdas already used, they don't change
      new_W <- gsmvar$model$structural_pars$W
      for(i1 in 1:ncol(new_W)) { # Iterate through each column of W
        col_vec <- gsmvar$model$structural_pars$W[, i1]
        if(all(is.na(col_vec) | col_vec == 0, na.rm=TRUE)) { # Are all elements in col_vec NA or zero?
          if(is.na(new_W[i1, i1])) { # Check if the diagonal element is NA
            new_W[i1, i1] <- 1 # Impose a positive sign constraints to the diagonal
          } else { # Zero constraint in the diagonal elements
            new_W[which(is.na(col_vec))[1], i1] <- 1 # Impose positive sign constraint on the first non-zero element
          }
        }
      }
      new_structural_pars <- list(W=new_W, fixed_lambdas=new_fixed_lambdas)

      # Finally, we need to make sure that the W params in new_params are in line with the constraints new_W.
      # This amounts checking the strict sign constraints and swapping the signs of the columns that don't
      # match the sign constraints.
      if(sum(W == 0, na.rm=TRUE) != sum(gsmvar$model$structural_pars$W == 0, na.rm=TRUE)) {
        # Throws an error since Wvec wont work properly if W contains exact zeros that are not constrained to zeros.
        stop(paste("A parameter value in W exactly zero but not constrained to zero.",
                   "Please adjust gsmvar$model$structural_pars$W so that the exact zeros match the constraints"))
      }
      # Determine which columns to swap: compare the first non-NA and non-zero element of the column of new_W to the
      # corresponding element of the corresponding column of W, and swap the signs of the column if the signs don't match.
      for(i1 in 1:ncol(new_W)) { # Loop through the columns
        col_new_W <- new_W[,i1]
        col_old_W <- W[,i1]
        which_to_compare <- which(!is.na(col_new_W) & col_new_W != 0)[1] # The first element that imposes a sign constraint
        if(sign(col_new_W[which_to_compare]) != sign(col_old_W[which_to_compare])) { # Different sign than the constrained one
          W[,i1] <- -W[,i1] # Swap the signs of the column
        }
      }
      # New params with the new W that corresponds to new_W constraints. Note that AR parameters are assumed
      # identical across the regimes here.
      new_params <- c(gsmvar$params[1:(d + p*d^2)], Wvec(W), all_df) # No lambdas or alphas
    } else { # is.null(structural_pars), i.e., recursive identification
      new_structural_pars <- NULL
    }

    ## Obtain residuals
    # Each y_t fixed, so the initial values y_{-p+1},...,y_0 are fixed in any case.
    # For y_1,...,y_T, new residuals are drawn at each bootstrap rep.
    # First, obtain the original residuals:
    mu_mt <- loglikelihood_int(data=data, p=gsmvar$model$p, M=gsmvar$model$M,
                               params=gsmvar$params, model=gsmvar$model$model,
                               conditional=gsmvar$model$conditional,
                               parametrization=gsmvar$model$parametrization,
                               constraints=gsmvar$model$constraints,
                               same_means=gsmvar$model$same_means,
                               weight_constraints=gsmvar$model$weight_constraints,
                               structural_pars=gsmvar$model$structural_pars,
                               to_return="total_cmeans")
    u_t <- data[(p+1):nrow(data),] - mu_mt # Residuals [T_obs, d]

    ## Functions that performs one bootstrap replication and return the IRF
    get_one_bootstrap_IRF <- function(seed) { # Take rest of the arguments from parent environment
      set.seed(seed) # Set seed for data generation
      estim_seeds <- sample.int(n=1e+6, size=ncalls) # Seeds for estimation

      ## Create new data
     eta_t <- sample(c(-1, 1), size=nrow(u_t), replace=TRUE, prob=c(0.5, 0.5))
     new_resid <- eta_t*u_t # each row of u_t multiplied by -1 or 1 based on eta_t
     new_data <- rbind(data[1:p,], # Fixed initial values
                       mu_mt + new_resid) # Bootstrapped data

      ## Estimate the model to the new data
      new_mod <- suppressMessages(fitGSMVAR(data=new_data, p=gsmvar$model$p, M=gsmvar$model$M,
                                                             model=gsmvar$model$model,
                                                             conditional=gsmvar$model$conditional,
                                                             parametrization=gsmvar$model$parametrization,
                                                             constraints=gsmvar$model$constraints,
                                                             same_means=gsmvar$model$same_means,
                                                             weight_constraints=new_weight_constraints, # Condition on alphas
                                                             structural_pars=new_structural_pars, # Condition on lambdas
                                                             ncalls=ncalls, seeds=estim_seeds,
                                                             print_res=FALSE, use_parallel=FALSE,
                                                             filter_estimates=TRUE, calc_std_errors=FALSE,
                                                             initpop=list(new_params))) # Initial values

      ## Get the IRF from the bootstrap replication
      tmp_params <- reform_constrained_pars(p=p, M=M_orig, d=d, params=new_mod$params, model=model,
                                            constraints=new_mod$model$constraints,
                                            same_means=new_mod$model$same_means,
                                            weight_constraints=new_mod$model$weight_constraints,
                                            structural_pars=new_mod$model$structural_pars)
      tmp_all_A <- pick_allA(p=p, M=M_orig, d=d, params=tmp_params,
                             structural_pars=structural_pars) # unconstrained struct pars
      tmp_all_boldA <- form_boldA(p=p, M=M_orig, d=d, all_A=tmp_all_A)
      tmp_all_Omega <- pick_Omegas(p=p, M=M_orig, d=d, params=tmp_params, structural_pars=structural_pars)

      # Obtain the impact matrix of the regime the IRF is to calculated for
      if(!is.null(structural_pars)) { # Shocks identified by heteroskedasticity
        tmp_W <- pick_W(p=p, M=M_orig, d=d, params=tmp_params, structural_pars=structural_pars)
        if(regime > 1) { # Include lambdas
          tmp_lambdas <- matrix(pick_lambdas(p=p, M=M_orig, d=d, params=tmp_params, structural_pars=structural_pars),
                                nrow=d, byrow=FALSE)
          tmp_Lambda_m <- diag(tmp_lambdas[, regime - 1])
          tmp_B_matrix <- tmp_W%*%sqrt(tmp_Lambda_m)
        } else { # regime == 1
          tmp_B_matrix <- tmp_W
        }
      } else { # Shocks identified by lower-triangular Cholesky decomposition
        tmp_B_matrix <- t(chol(tmp_all_Omega[, , regime]))
      }
      # Calculate and return the IRF
      get_IRF(p=p, d=d, N=N, boldA=tmp_all_boldA[, , regime], B_matrix=tmp_B_matrix)
    }

    ## Calculate the bootstrap replications using parallel computing
    if(ncores > parallel::detectCores()) {
      ncores <- parallel::detectCores()
      message("ncores was set to be larger than the number of cores detected")
    }
    cat(paste("Using", ncores, "cores for", bootstrap_reps, "bootstrap replications..."), "\n")
    cl <- parallel::makeCluster(ncores)
    on.exit(try(parallel::stopCluster(cl), silent=TRUE)) # Close the cluster on exit, if not already closed.
    parallel::clusterExport(cl, ls(environment(fitGSMVAR)),
                            envir = environment(fitGSMVAR)) # assign all variables from package:gmvarkit
    parallel::clusterEvalQ(cl, c(library(Brobdingnag), library(mvnfast), library(pbapply)))
    all_bootstrap_IRF <- pbapply::pblapply(1:bootstrap_reps, function(i1) get_one_bootstrap_IRF(seed=seeds[i1]), cl=cl)
    parallel::stopCluster(cl=cl)

  } else { # is.null(ci), no bounds
    all_bootstrap_IRF <- NULL
  }

  ## Accumulate IRF based on which_cumulative
  if(length(which_cumulative) > 0) {
    for(which_var in which_cumulative) {
      # Accumulate the impulse responses of the variables in which_cumulative
      point_est[which_var, , ] <- t(apply(point_est[which_var, , , drop=FALSE], MARGIN=2, FUN=cumsum))
      if(!is.null(all_bootstrap_IRF)) { # Do the same accumulation for each bootstrap replication:
        for(i2 in 1:length(all_bootstrap_IRF)) {
          all_bootstrap_IRF[[i2]][which_var, , ] <- t(apply(all_bootstrap_IRF[[i2]][which_var, , , drop=FALSE],
                                                            MARGIN=2, FUN=cumsum))
        }
      }
    }
  }

  ## Scale the IRFs
  if(!is.null(scale)) {
    for(i1 in 1:ncol(scale)) {
      which_shock <- scale[1, i1]
      which_var <- scale[2, i1]
      scale_size <- scale[3, i1]
      # Scale the IRFs of which_shock to correspond scale_size impact response of the variable which_var:
      multiplier <- scale_size/point_est[which_var, which_shock, 1]
      point_est[, which_shock, ] <- multiplier*point_est[, which_shock, ] # Impact response to scale_size
      if(!is.null(all_bootstrap_IRF)) { # Do the same scaling for each bootstrap replication:
        for(i2 in 1:length(all_bootstrap_IRF)) {
          multiplier <- scale_size/all_bootstrap_IRF[[i2]][which_var, which_shock, 1]
          all_bootstrap_IRF[[i2]][, which_shock, ] <- multiplier*all_bootstrap_IRF[[i2]][, which_shock, ]
        }
      }
    }
  }

  ## Calculate the confidence bounds
  if(!is.null(all_bootstrap_IRF)) {
    # First we convert the list into a 4D array in order to use apply to calculate empirical qunatiles
    all_bootstrap_IRF_4Darray <- array(NA, dim = c(dim(all_bootstrap_IRF[[1]]), length(all_bootstrap_IRF)))
    for (i1 in 1:length(all_bootstrap_IRF)) { # Fill the arrays
      all_bootstrap_IRF_4Darray[, , , i1] <- all_bootstrap_IRF[[i1]]
    }

    # Calculate empirical quantiles to obtain ci
    quantile_fun <- function(x) quantile(x, probs=c((1 - ci)/2, 1 - (1 - ci)/2), na.rm=TRUE)
    conf_ints <- apply(all_bootstrap_IRF_4Darray, MARGIN=1:3, FUN=quantile_fun)
    conf_ints <- aperm(conf_ints, perm=c(2, 3, 4, 1))
    dimnames(conf_ints)[[1]] <- colnames(gsmvar$data)
    dimnames(conf_ints)[[2]] <- paste("Shock", 1:gsmvar$model$d)
  } else {
    conf_ints <- NULL
    all_bootstrap_IRF_4Darray <- NULL
  }

  # Return the results
  structure(list(point_est=point_est,
                 conf_ints=conf_ints,
                 all_bootstrap_reps=all_bootstrap_IRF_4Darray,
                 N=N,
                 ci=ci,
                 scale=scale,
                 which_cumulative=which_cumulative,
                 seeds=seeds,
                 gsmvar=gsmvar),
            class="irf")
}
