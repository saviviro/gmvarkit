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
#'   correspond to a specific magnitude of instantaneous response of some specific
#'   variable? Provide a length three vector where the shock of interest
#'   is given in the first element (an integer in \eqn{1,...,d}), the variable of
#'   interest is given in the second element (an integer in \eqn{1,...,d}), and
#'   the magnitude of its instantaneous response in the third element (a non-zero real number).
#'   If the linear IRFs of multiple shocks should be scaled, provide a matrix which has one
#'   column for each of the shocks with the columns being the length three vectors described above.
#' @param ci a real number in \eqn{(0, 1)} specifying the confidence level of the
#'   confidence intervals calculated via a fixed-design wild residual bootstrap method.
#'   Available only for models that impose linear autoregressive dynamics
#'   (excluding changes in the volatility regime).
#' @param bootstrap_reps the number of bootstrap repetitions for estimating confidence bounds.
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

linear_IRF <- function(gsmvar, N=30, regime=1, which_cumulative=numeric(0),
                       scale=NULL, ci=NULL, bootstrap_reps=NULL, seeds=NULL, ...) {
  # Get the parameter values etc
  stopifnot(all_pos_ints(c(N, regime, bootstrap_reps)))
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
    all_phi_i <- array(NA, dim=c(d, d, N+1)) # JA^iJ' matrices; [, , 1] is for the impact period, i1+1 for period i1
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
    all_Phi_i # all_Phi_i[variable, shock, horizon] -> all_Phi_i[variable, shock, ] subsets the IRF!
  }

  ## Calculate the impulse response functions:
  point_est_IRF <- get_IRF(p=p, d=d, N=N,
                           boldA=all_boldA[, , regime], # boldA= Companion form AR matrix of the selected regime
                           B_matrix=B_matrix)

  ## Confidence bounds by fixed design wild residual bootstrap
  AR_mats_identical <- all(apply(all_boldA, MARGIN=3, FUN=function(x) identical(x, all_boldA[,,1])))
  means_identical <- !is.null(same_means) && length(same_means) == 1 && all(same_means[[1]] == 1:M)
  ci_possible <- (means_identical && AR_mats_identical) || M == 1

  if(!is.null(ci) && !ci_possible) {
    warning("Confidence bounds are not available as the autoregressive dynamics are not linear")
  } else if(!is.null(ci) && ci_possible) { # Bootstrap confidence bounds

    ## Create initial values for the two-phase estimation algorithm: does not vary across the bootstrap reps
    new_params <- gsmvar$params

    # For all models, bootstrapping conditions on the estimated mw parameters, so they need to be removed:
    if(is.null(weight_constraints) && M > 1) {
      new_params <- c(new_params[1:(length(new_params) - (M - 1) - length(all_df))], all_df) # Removes alphas
      weight_constraints <- alphas[-M]
    } # If weight constraints used or M == 1, no modifications are required

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
      if(is.null(gsmvar$model$structural_pars$fixed_lambdas)) {
        # Remove lambda parameters from the parameter vector (note: alphas already removed)
        new_params <- c(new_params[1:(length(new_params) - d*(M - 1) - length(all_df))], all_df)
      }
    } else {
      new_structural_pars <- NULL
    }

    # TARVITAAN SWAP_W_SIGNSIA SITEN ETTÄ UUSISSA PARAMETREISSA W:N MERKIT VASTAA W:N RAJOITTEITA!
    # TEE SE SEURAAVAKSI!!

    ## Obtain residuals
    # Each y_t fixed, so the initial values y_{-p+1},...,y_0 are fixed in any case.
    # For y_1,...,y_T, new residuals are drawn at each bootstrap rep
    # First, obtain the residuals:
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

    # Tämä funktioksi ja parallel computingin sisään kun toimivuus on testattu!
    all_new_IRF <- list()
    for(i1 in 1:bootstrap_reps) {
      set.seed(seeds[i1]) # Set seed for data generation
      ncalls <- 1 # The number of estimation rounds in each bootstrap replication
      estim_seeds <- sample.int(n=1e+6, size=ncalls) # Seeds for estimation

      ## Create new data
      eta_t <- sample(c(-1, 1), size=nrow(u_t), replace=TRUE, prob=c(0.5, 0.5))
      new_resid <- eta_t*u_t # each row of u_t multiplied by -1 or 1 based on eta_t
      new_data <- rbind(data[1:p,], # Fixed initial values
                        mu_mt + new_resid) # Bootstrapped data

      ## Estimate the model to the new data ADD SUPPRESS WARNIGNS HERE AFTER TESTING!!!! warn_dfs pointless here
      new_mod <- suppressMessages(fitGSMVAR(data=new_data, p=gsmvar$model$p, M=gsmvar$model$M,
                                            model=gsmvar$model$model,
                                            conditional=gsmvar$model$conditional,
                                            parametrization=gsmvar$model$parametrization,
                                            constraints=gsmvar$model$constraints,
                                            same_means=gsmvar$model$same_means,
                                            weight_constraints=new_weight_constraints, # Condition on alphas
                                            structural_pars=new_structural_pars, # Condition on lambdas
                                            ncalls=ncalls, seeds=estim_seems,
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

      new_IRF <- get_IRF(p=p, d=d, N=N, boldA=tmp_all_boldA[, , regime], B_matrix=tmp_B_matrix)
      all_new_IRF[i1] <- new_IRF
    }

    # HUOM! Parallel computing seed to GA + seed before data gen.
  } else { # is.null(ci), no bounds

  }


  # Maybe also remove "which_shocks"? IRF calculated for all anyway - are there implications to the ci?
  # All the IRFs are calculate for UNIT SHOCKS i.e., one standard error struct shocks
  # -> responses can be scaled without loss of generality as the IRFs are symmetric w.r.t the size of the shock

  # Tuleeko CI samaan funktioon? Sitä varten olisi varmaan syytä ensin tarkistaa, että onko AR-parametrit rajoitettu lineaarisiksi!

  # HUOM wild bootstrapissa käytetään forecast erroreita eikä virhetermien empiirisiä vastineita.
  # IDEA (omaan paperiinsa??): arvo eta_i:n sijasta s_{m,t} uudestaan ja sen perusteella virhe käyttäen sekoitussuhteita
  # (ehkä tyhmä idea mutta tuli mieleen)

  # HUOM! Entä sokkien järjestys ja merkki bootstrapatessa? Herwards and Lutkepohl (2014) käyttää samoja Lambda-parametreja kuin
  # alkuperäisessä estimoinnissa! Eli pitää ehkä implentoida rajoite, jossa lambda-parametrit on rajoitettu tietyiksi luvuiksi?
  # Senhän voi määritellä kokonaan erikseen C_lambda rajoitteesta, jolloin vanhoja yksikkötestejä jne ei tarvitse muuttaa.
  # Samoin ne käyttää alkuperäisiä transition-probabiliteja! alpha-parametrit asettamalla tietyiksi luvuiksi päästään vähän
  # saman tyyppiseen.
  # Lutkepohl and Netsunajev (2014) käytti samaa kuin Herwards and Lutkepohl (2014). Meillähän tuo nyt eroaa niin, että koska
  # transition probabilityt riippuu AR-matriiseista ja kov.matriiseista, niin ihan samoihin transition probiliteihin ei voida
  # ehdollistaa by construction.
  # Netsunajev (2013) ei ehdollista samoihin Lambdoihin, mutta taitaa näyttää yli-identifoivia rajoitteita siten ettei
  # sokkien järjestys muutu.
  # Simukokeessahan sokin vain uudelleenjärjesteltiin niin että ne oli mahdollisimman lähellä todellisia parametriarvoja
  # HUOM! Merkkien kääntyminen pitää huomioida bootstrappauksessa jotenkin; joko rajoittamalla merkit alunperinkin tietynlaisiksi
  # estimointia tehdessä tai kääntämällä merkit sitten että joku etäisyys alkuperäisistä estimaateista olisi mahdollisimman pieni tms?
  # Netsunajev (2013) tekee jotain tämän tyyppistä; mieti miten implementoit yleiseen tapaukseen bootstrappauksessa! Kiinnostuksen kohteena
  # olevan sokin merkki nyt ainakin pitää olla tietty, joten sehän käytännössä määrää sen miten päin ne merkit tulee. Esim rajoittamalla
  # diagonaalit W:ssä positiivisiksi jos muuta merkkirajoitetta ei kyseisellä sarakkeella jo valmiiksi ole (ja jos on niin sitä saraketta
  # ei muuteta merkin suhteen).

  # Rekursiivisella identifoinnilla tietysti sokkien järjestys on aina oikein.


  # NOTE which cumulative does not do anything yet! Maybe original + cumulative for all?
  # Depends on how the confidence bounds are created!
  # ALSO SCALING IS NOT IMPLEMENTED YET: scaling based in initial response

  # Return the results
  structure(list(all_irfs=point_est_IRF,
                 N=N,
                 ci=ci,
                 which_cumulative=which_cumulative,
                 seed=seed,
                 gsmvar=gsmvar),
            class="irf")
}
