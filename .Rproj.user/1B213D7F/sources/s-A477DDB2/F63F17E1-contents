#' @title Genetic algorithm for preliminary estimation of GMVAR model
#'
#' @description \code{GAfit} estimates specified GMVAR model using genetic algorithm. It's designed to find
#'   starting values for gradient based methods.
#'
#' @inheritParams loglikelihood_int
#' @inheritParams random_covmat
#' @param ngen a positive integer specifying the number of generations to be ran through in
#'   the genetic algorithm.
#' @param popsize a positive even integer specifying the population size in the genetic algorithm.
#'   Default is \code{10*n_params}.
#' @param smart_mu a positive integer specifying the generation after which the random mutations
#'   in the genetic algorithm are "smart". This means that mutating individuals will mostly mutate fairly
#'   close (or partially close) to the best fitting individual (which has the least regimes with time varying
#'   mixing weights practically at zero) so far.
#' @param initpop a list of parameter vectors from which the initial population of the genetic algorithm
#'   will be generated from. The parameter vectors should be...
#'   \describe{
#'     \item{\strong{For regular models:}}{
#'       Should be size \eqn{((M(pd^2+d+d(d+1)/2+1)-1)x1)} and have form
#'       \strong{\eqn{\theta}}\eqn{ = }(\strong{\eqn{\upsilon}}\eqn{_{1}},
#'       ...,\strong{\eqn{\upsilon}}\eqn{_{M}}, \eqn{\alpha_{1},...,\alpha_{M-1}}), where:
#'       \itemize{
#'         \item \strong{\eqn{\upsilon}}\eqn{_{m}} \eqn{ = (\phi_{m,0},}\strong{\eqn{\phi}}\eqn{_{m}}\eqn{,\sigma_{m})}
#'         \item \strong{\eqn{\phi}}\eqn{_{m}}\eqn{ = (vec(A_{m,1}),...,vec(A_{m,p})}
#'         \item and \eqn{\sigma_{m} = vech(\Omega_{m})}, m=1,...,M.
#'       }
#'     }
#'     \item{\strong{For constrained models:}}{
#'       Should be size \eqn{((M(d+d(d+1)/2+1)+q-1)x1)} and have form
#'       \strong{\eqn{\theta}}\eqn{ = (\phi_{1,0},...,\phi_{M,0},}\strong{\eqn{\psi}}
#'       \eqn{,\sigma_{1},...,\sigma_{M},\alpha_{1},...,\alpha_{M-1})}, where:
#'       \itemize{
#'         \item \strong{\eqn{\psi}} \eqn{(qx1)} satisfies (\strong{\eqn{\phi}}\eqn{_{1}}\eqn{,...,}
#'         \strong{\eqn{\phi}}\eqn{_{M}) =} \strong{\eqn{C \psi}}. Here \strong{\eqn{C}} is \eqn{(Mpd^2xq)}
#'         constraint matrix.
#'       }
#'     }
#'   }
#'   Above \eqn{\phi_{m,0}} is the intercept parameter, \eqn{A_{m,i}} denotes the \eqn{i}:th coefficient matrix of the \eqn{m}:th
#'   mixture component, \eqn{\Omega_{m}} denotes the error term covariance matrix of the \eqn{m}:th mixture component and
#'   \eqn{\alpha_{m}} is the mixing weight parameter.
#'   If \code{parametrization=="mean"}, just replace each \eqn{\phi_{m,0}} with regimewise mean \eqn{\mu_{m}}.
#'   \eqn{vec()} is vectorization operator that stacks columns of a given matrix into a vector. \eqn{vech()} stacks columns
#'   of a given matrix from the principal diagonal downwards (including elements on the diagonal) into a vector.
#'   The notations are in line with the cited article by \emph{Kalliovirta, Meitz and Saikkonen (2016)}.
#' @param mu_scale a size \eqn{(dx1)} vector defining \strong{means} of the normal distributions from which each
#'   mean parameter \eqn{\mu_{m}} is drawn from in random mutations. Default is \code{colMeans(data)}. Note that
#'   mean-parametrization is always used for optimization in \code{GAfit} - even when \code{parametrization=="intercept"}, but
#'   input (in \code{initpop}) and output (return value) parameter vectors may be intercept-parametrized.
#' @param mu_scale2 a size \eqn{(dx1)} strictly positive vector defining \strong{standard deviations} of the normal
#'   distributions from which each mean parameter \eqn{\mu_{m}} is drawn from in random mutations.
#'   Default is \code{2*sd(data[,i]), i=1,..,d}.
#' @param omega_scale a size \eqn{(dx1)} strictly positive vector specifying the scale and variability of the
#'   random covariance matrices in random mutations. The covariance matrices are drawn from (scaled) Wishart distribution.
#'   Expected values of the random covariance matrices are \code{diag(omega_scale)}. Standard deviations
#'   of the diagonal elements are \code{sqrt(2/d)*omega_scale[i]}
#'   and for non-diagonal elements they are \code{sqrt(1/d*omega_scale[i]*omega_scale[j])}.
#'   Note that for \code{d>4} this scale may need to be chosen carefully. Default in \code{GAfit} is
#'   \code{var(stats::ar(data[,i], order.max=10)$resid, na.rm=TRUE), i=1,...,d}.
#' @param ar_scale a positive real number adjusting how large AR parameter values are typically generated in some random
#'   mutations. See function \code{random_coefmats2} for details. This is ignored when estimating constrained models.
#' @param regime_force_scale a non-negative real number specifying how much should natural selection favour individuals
#'   with less regimes that have almost all mixing weights (practically) at zero. Set to zero for no favouring or large number
#'   for heavy favouring. Without any favouring the genetic algorithm gets more often stuck in an area of the parameter space where some
#'   regimes are wasted, but with too much favouring the best genes might never mix into the population and the algorithm might
#'   converge poorly. Default is \code{1} and it gives \eqn{2x} larger surviving probabilities for individuals with no wasted
#'   regimes compared to individuals with one wasted regime. Number \code{2} would give \eqn{3x} larger probabilities etc.
#' @param red_criteria a length 2 numeric vector specifying the criteria that is used to determine whether a regime is redundant or not.
#'   Any regime \code{m} which satisfies \code{sum(mixingWeights[,m] > red_criteria[1]) < red_criteria[2]*n_obs} will be considered "redundant".
#'   One should be careful when adjusting this argument.
#' @param to_return should the genetic algorithm return the best fitting individual which has "positive enough" mixing weights
#'   for as much regimes as possible (\code{"alt_ind"}) or the individual which has the highest log-likelihood in general (\code{"best_ind"}), but
#'   might have more wasted regimes? Default is \code{"alt_ind"}.
#' @param minval a real number defining the minimum value of the log-likelihood function that will be considered.
#'   Values smaller than this will be treated as they were \code{minval} and the corresponding individuals will
#'   never survive. The default is \code{-(10^(ceiling(log10(n_obs)) + d) - 1)}.
#' @details
#'    The genetic algorithm is mostly based on the description by \emph{Dorsey and Mayer (1995)}.
#'    It uses (slightly modified) individually adaptive crossover and mutation rates described by \emph{Patnaik and Srinivas (1994)}
#'    and employs (50\%) fitness inheritance discussed by \emph{Smith, Dike and Stegmann (1995)}.
#'
#'    By "redundant" or "unidentified" regimes we mean regimes that have the time varying mixing weights basically at zero for all t.
#'    The model would have the same log-likelihood value without redundant regimes and there is no purpose to have redundant regime in
#'    the model.
#' @return Returns estimated parameter vector which has the form described in \code{initpop}.
#' @references
#'  \itemize{
#'    \item Ansley C.F., Kohn R. 1986. A note on reparameterizing a vector autoregressive
#'          moving average model to enforce stationarity. \emph{Journal of statistical computation
#'          and simulation}, \strong{24}:2,  99-106.
#'    \item Dorsey R. E. and Mayer W. J. 1995. Genetic algorithms for estimation problems with multiple optima,
#'          nondifferentiability, and other irregular features. \emph{Journal of Business & Economic Statistics},
#'         \strong{13}, 53-66.
#'    \item Kalliovirta L., Meitz M. and Saikkonen P. 2016. Gaussian mixture vector autoregression.
#'          \emph{Journal of Econometrics}, \strong{192}, 485-498.
#'    \item Patnaik L.M. and Srinivas M. 1994. Adaptive Probabilities of Crossover and Mutation in Genetic Algorithms.
#'          \emph{Transactions on Systems, Man and Cybernetics} \strong{24}, 656-667.
#'    \item Smith R.E., Dike B.A., Stegmann S.A. 1995. Fitness inheritance in genetic algorithms.
#'          \emph{Proceedings of the 1995 ACM Symposium on Applied Computing}, 345-350.
#'  }
#'  @export


GAfit <- function(data, p, M, conditional=TRUE, parametrization=c("intercept", "mean"), constraints=NULL, ngen=200, popsize,
                  smart_mu=min(100, round(0.5*ngen)), initpop=NULL, mu_scale, mu_scale2, omega_scale, ar_scale=1,
                  regime_force_scale=1, red_criteria=c(0.05, 0.01), to_return=c("alt_ind", "best_ind"), minval) {

  # Required values and premilinary checks
  to_return <- match.arg(to_return)
  parametrization <- match.arg(parametrization)
  check_pMd(p=p, M=M)
  data <- check_data(data=data, p=p)
  d <- ncol(data)
  n_obs <- nrow(data)
  check_constraints(p=p, M=M, d=d, constraints=constraints)
  npars <- n_params(p=p, M=M, d=d, constraints=constraints)

  # Defaults and checks
  if(!all_pos_ints(c(ngen, smart_mu))) stop("Arguments ngen and smart_mu has to be positive integers")
  if(missing(popsize)) {
    popsize <- 10*npars
  } else if(popsize < 2 | popsize %% 2 != 0) {
    stop("The population size popsize must be even positive integer")
  }
  if(missing(minval)) {
    minval <- -(10^(ceiling(log10(n_obs)) + d) - 1)
  } else if(!is.numeric(minval)) {
    stop("Argument minval must be numeric")
  }
  if(missing(mu_scale)) {
    mu_scale <- colMeans(data)
  } else if(length(mu_scale) != d) {
    stop("Argument mu_scale must be numeric vector with length d")
  }
  if(missing(mu_scale2)) {
    mu_scale2 <- vapply(1:d, function(i1) 2*sd(data[,i1]), numeric(1))
  } else if(length(mu_scale2) != d | any(mu_scale2 <= 0)) {
    stop("Argument mu_scale2 must be strictly positive vector with length d")
  }
  if(missing(omega_scale)) {
    omega_scale <- vapply(1:d, function(i1) var(stats::ar(data[,i1], order.max=10)$resid, na.rm=TRUE), numeric(1))
  } else if(!(length(omega_scale) == d & all(omega_scale > 0))) {
    stop("omega_scale must be numeric vector with length d and positive elements")
  }
  if(length(ar_scale) != 1 | ar_scale <= 0) {
    stop("ar_scale must be positive and have length one")
  }
  if(length(regime_force_scale) != 1 | regime_force_scale < 0) {
    stop("regime_force_scale should be non-negative real number")
  }

  # Initial population
  if(is.null(initpop)) {
    if(is.null(constraints)) {
      G <- replicate(popsize, random_ind2(p=p, M=M, d=d, mu_scale=mu_scale, mu_scale2=mu_scale2, omega_scale=omega_scale, ar_scale=ar_scale))
    } else {
      G <- replicate(popsize, random_ind(p=p, M=M, d=d, constraints=constraints, mu_scale=mu_scale, mu_scale2=mu_scale2, omega_scale=omega_scale))
    }
  } else {
    stopifnot(is.list(initpop))
    for(i1 in 1:length(initpop)) {
      ind <- initpop[[i1]]
      tryCatch(check_parameters(p=p, M=M, d=d, params=ind, constraints=constraints),
               error=function(e) stop(paste("Problem with individual", i1, "in the initial population: "), e))
      if(parametrization=="intercept") {
        ind <- change_parametrization(p=p, M=M, d=d, params=ind, constraints=constraints, change_to="mean")
      }
      initpop[[i1]] <- sort_components(p=p, M=M, d=d, params=ind)
    }
    G <- replicate(popsize, initpop[[sample.int(length(initpop), size=1)]])
  }

  # Calculates the number of redundant regimes
  n_redundants <- function(M, mw) {
    sum(vapply(1:M, function(m) sum(mw[,m] > red_criteria[1]) < red_criteria[2]*n_obs, logical(1)))
  }

  # Initial setup
  generations <- array(dim=c(npars, popsize, ngen))
  logliks <- matrix(minval, nrow=ngen, ncol=popsize)
  redundants <- matrix(M, nrow=ngen, ncol=popsize) # Store the number of redundant regimes of each individual
  which_redundant_alt <- 1:M

  fill_lok_and_red <- function(i1, i2, lok_and_mw) {
    if(!is.list(lok_and_mw)) {
      logliks[i1, i2] <<- minval
      redundants[i1, i2] <<- M
    } else {
      logliks[i1, i2] <<- lok_and_mw$loglik
      redundants[i1, i2] <<- n_redundants(M=M, mw=lok_and_mw$mw) # Number of redundant regimes
    }
  }

  # Run through generations
  for(i1 in 1:ngen) {
    generations[, , i1] <- G

    # Calculate log-likelihoods and fitness inheritance. Use minval if values smaller than that.
    if(i1 == 1) {
      for(i2 in 1:popsize) {
        loks_and_mw <- loglikelihood_int(data, p, M, params=G[,i2], conditional=conditional, parametrization="mean",
                                         constraints=constraints, to_return="loglik_and_mw", check_params=TRUE, minval=minval)
        fill_lok_and_red(i1, i2, loks_and_mw)
      }
    } else {
      # Proportional fitness inheritance: individual has 50% change to inherit fitness if it's a result of crossover
      # Variable "I" tells the proportions of parent material.
      I2 <- rep(I, each=2)
      which_did_co <- which(1 - which_not_co == 1)
      if(length(which_did_co) > 0) {
        which_inherit <- sample(x=which_did_co, size=round(0.5*length(which_did_co)), replace=FALSE)
      } else {
        which_inherit <- numeric(0)
      }
      # survivor_liks holds the parent loglikelihood values: for odd number they are (index, index+1) and for even (index-1, index)
      if(length(which_inherit) > 0 & abs(max_lik - avg_lik) > abs(0.03*avg_lik)) { # No inheritance if massive mutations
        for(i2 in which_inherit) {
          if(i2 %% 2 == 0) {
            logliks[i1, i2] <- ((npars - I2[i2])/npars)*survivor_liks[i2-1] + (I2[i2]/npars)*survivor_liks[i2]
            redundants[i1, i2] <- max(c(survivor_redundants[i2-1], survivor_redundants[i2]))
          } else {
            logliks[i1, i2] <- (I2[i2]/npars)*survivor_liks[i2] + ((npars - I2[i2])/npars)*survivor_liks[i2+1]
            redundants[i1, i2] <- max(c(survivor_redundants[i2], survivor_redundants[i2+1]))
          }
        }
        which_no_inherit <- (1:popsize)[-which_inherit]
      } else {
        which_no_inherit <- 1:popsize
      }
      for(i2 in which_no_inherit) { # Calculate the rest log-likelihoods
          if((mutate[i2] == 0 & which_not_co[i2] == 1) | all(H[,i2] == G[,i2])) { # If nothing changed
            logliks[i1, i2] <- survivor_liks[i2]
            redundants[i1, i2] <- survivor_redundants[i2]
          } else {
            if(stat_mu == TRUE & mutate[i2] == 1) {
              loks_and_mw <- tryCatch(loglikelihood_int(data, p, M, params=G[,i2], conditional=conditional, parametrization="mean",
                                               constraints=constraints, to_return="loglik_and_mw",
                                               check_params=FALSE, minval=minval), error=function(e) minval)
            } else {
              loks_and_mw <- loglikelihood_int(data, p, M, params=G[,i2], conditional=conditional, parametrization="mean",
                                               constraints=constraints, to_return="loglik_and_mw", check_params=TRUE, minval=minval)
            }
            fill_lok_and_red(i1, i2, loks_and_mw)
         }
      }
    }
    logliks[i1, which(logliks[i1,] < minval)] <- minval
    redundants[i1, which(logliks[i1,] <= minval)] <- M

    # Calculate choosing probabilities
    if(length(unique(logliks[i1,])) == 1) {
      choosing_probs <- rep(1, popsize)
    } else {
      T_values <- logliks[i1,] + abs(min(logliks[i1,])) # Function T as described by Dorsey R. E. ja Mayer W. J., 1995
      T_values <- T_values/(1 + regime_force_scale*redundants[i1,]) # Favor individuals with least number of redundant regimes
      choosing_probs <- T_values/sum(T_values)
    }

    # Draw popsize individuals "with put back" and form the reproduction pool+
    survivors <- sample(1:popsize, size=popsize, replace=TRUE, prob=choosing_probs)
    H <- G[,survivors]

    # Calculate avg and max log-likelihood of the survivors
    survivor_liks <- logliks[i1, survivors]
    survivor_redundants <- redundants[i1, survivors]
    max_lik <- max(survivor_liks)
    avg_lik <- mean(survivor_liks)
    if(max_lik == avg_lik) avg_lik <- avg_lik + 0.5

    # Individually adaptive cross-over rates, Patnaik and Srinivas (1994) * "extra CO" 0.4
    indeces <- seq(1, popsize-1, by=2)
    parent_max <- vapply(indeces, function(i2) max(survivor_liks[i2], survivor_liks[i2+1]), numeric(1))
    co_rates <- vapply(1:length(indeces), function(i2) max(min((max_lik - parent_max[i2])/(max_lik - avg_lik), 1), 0.4), numeric(1))

    # Do the crossovers
    which_co <- rbinom(n=popsize/2, size=1, prob=co_rates)
    I <- round(runif(n=popsize/2, min=0.5 + 1e-16, max=npars - 0.5 - 1e-16))
    H2 <- vapply(1:(popsize/2), function(i2) {
            i3 <- indeces[i2]
            if(which_co[i2] == 1) {
              c(c(H[1:I[i2], i3], H[(I[i2]+1):npars, i3+1]), c(H[1:I[i2], i3+1], H[(I[i2]+1):npars, i3]))
            } else {
              c(H[,i3], H[,i3+1])
            }
          }, numeric(2*npars))
    H2 <- matrix(H2, nrow=npars, byrow=FALSE)

    # Get the best individual so far and check for reduntant regimes
    best_index0 <- which(logliks == max(logliks), arr.ind=TRUE)
    best_index <- best_index0[order(best_index0[,1], decreasing=FALSE)[1],] # First generation when the best loglik occured (because of fitness inheritance)
    best_ind <- generations[, best_index[2], best_index[1]]
    best_mw <- loglikelihood_int(data, p, M, params=best_ind, conditional=conditional, parametrization="mean",
                                 constraints=constraints, to_return="mw", check_params=FALSE, minval=minval)
    which_redundant <- which(vapply(1:M, function(i2) sum(best_mw[,i2] > red_criteria[1]) < red_criteria[2]*n_obs, logical(1))) # Which regimes are wasted

    # Keep track of "alternative best individual" that has less reduntant regimes than the current one (after the algorithm finds one)
    if(length(which_redundant) <= length(which_redundant_alt)) {
      alt_ind <- best_ind
      which_redundant_alt <- which_redundant
    }

    # Individually adaptive mutation rates, Patnaik and Srinivas (1994)
    which_not_co <- rep(1 - which_co, each=2)
    if(abs(max_lik - avg_lik) <= abs(0.03*avg_lik)) {
      mu_rates <- rep(0.7, popsize) # Massive mutations if converging
    } else {
      mu_rates <- 0.5*vapply(1:popsize, function(i2) min(which_not_co[i2], (max_lik - survivor_liks[i2])/(max_lik - avg_lik)), numeric(1))
    }

    # Do mutations and keep track if they are stationary for sure
    mutate <- rbinom(n=popsize, size=1, prob=mu_rates)
    which_mutate <- which(mutate == 1)
    if(i1 <= smart_mu & length(which_mutate) >= 1) { # Random mutations
      if(!is.null(constraints) | runif(1) > 0.5) { # Normal
        stat_mu <- FALSE
        H2[,which_mutate] <- vapply(1:length(which_mutate), function(x) random_ind(p=p, M=M, d=d, constraints=constraints,
                                                                                   mu_scale=mu_scale, mu_scale2=mu_scale2,
                                                                                   omega_scale=omega_scale), numeric(npars))
      } else { # Stationary by algorithm (slower), Ansley and Kohn (1986)
        stat_mu <- TRUE
        H2[,which_mutate] <- vapply(1:length(which_mutate), function(x) random_ind2(p=p, M=M, d=d, mu_scale=mu_scale,
                                                                                    mu_scale2=mu_scale2, omega_scale=omega_scale,
                                                                                    ar_scale=ar_scale), numeric(npars))
      }
    } else if(length(which_mutate) >= 1) { # Smart mutations
      stat_mu <- FALSE

      # If redundant regimes - smart mutate more
      if(length(which_redundant) >= 1) {
        mu_rates <- vapply(1:popsize, function(i2) which_not_co[i2]*max(0.1, mu_rates[i2]), numeric(1))
        mutate0 <- rbinom(n=popsize, size=1, prob=mu_rates)
        which_mutate0 <- which(mutate0 == 1)
        if(length(which_mutate0) > length(which_mutate)) {
          mutate <- mutate0
          which_mutate <- which_mutate0
        }
      }

      # Mutating accuracy
      accuracy <- abs(rnorm(length(which_mutate), mean=10, sd=15))

      if(!is.null(constraints) | length(which_redundant) <= length(which_redundant_alt) | runif(1) > 0.5) { # No regime combining with constrained models
        # Smart mutations close to alt_ind
        ind_to_use <- alt_ind
        rand_to_use <- which_redundant_alt
      } else {
        # Change redundant regime of best_ind to be non-redundant regime of alt_ind to form the individual used in smart mutations
        # Choose regime of best_ind to be changed
        which_to_change <- which_redundant[1]

        # Pick non_redundant regimes of best_ind and of alt_ind
        non_red_regs_best <- vapply((1:M)[-which_redundant], function(i2) pick_regime(p=p, M=M, d=d, params=best_ind, m=i2), numeric(p*d^2 + d + d*(d+1)/2))
        if(length(which_redundant_alt) == 0) {
          non_red_regs_alt <- vapply(1:M, function(i2) pick_regime(p=p, M=M, d=d, params=alt_ind, m=i2), numeric(p*d^2 + d + d*(d+1)/2))
        } else {
          non_red_regs_alt <- vapply((1:M)[-which_redundant_alt], function(i2) pick_regime(p=p, M=M, d=d, params=alt_ind, m=i2), numeric(p*d^2 + d + d*(d+1)/2))
        }

        # Calculate "distances" between the considered regimes
        dist_to_regime <- matrix(nrow=ncol(non_red_regs_best), ncol=ncol(non_red_regs_alt)) # Row for each non-red-reg-best and column for each non-red-reg-alt.
        for(i2 in 1:nrow(dist_to_regime)) {
          dist_to_regime[i2,] <- vapply(1:ncol(non_red_regs_alt), function(i3) regime_distance(regime_pars1=non_red_regs_best[,i2],
                                                                                               regime_pars2=non_red_regs_alt[,i3]), numeric(1))
        }

        # Which alt_ind regime, i.e. column should be used? Choose the one that with largest distance to avoid dublicating similar regimes
        reg_to_use <- non_red_regs_alt[,which(colMeans(dist_to_regime) == max(colMeans(dist_to_regime)))]

        # Change the chosen regime of best_ind to be the one chosen from alt_ind
        ind_to_use <- change_regime(p=p, M=M, d=d, params=best_ind, m=which_to_change, regime_pars=reg_to_use)

        # Should some regimes still be random?
        rand_to_use <- which_redundant[-which_to_change]
      }

      # Smart mutations
      H2[,which_mutate] <- vapply(1:length(which_mutate), function(i2) smart_ind(p, M, d, params=ind_to_use, constraints=constraints,
                                                                                 accuracy=accuracy[i2], which_random=rand_to_use,
                                                                                 mu_scale=mu_scale, mu_scale2=mu_scale2,
                                                                                 omega_scale=omega_scale, ar_scale=ar_scale), numeric(npars))
    }

    # Sort components by mixing weight parameters. No sorting if constraints are employed.
    if(is.null(constraints)) {
      H2 <- vapply(1:popsize, function(i2) sort_components(p=p, M=M, d=d, params=H2[,i2]), numeric(npars))
    }

    # Save the results and set up new generation
    G <- H2
  }

  if(to_return == "best_ind") {
    ret <- best_ind
  } else {
    ret <- alt_ind
  }

  # GA always optimizes with mean parametrization. Return intercept-parametrized estimate if parametrization=="intercept".
  if(parametrization=="mean") {
    return(ret)
  } else {
    return(change_parametrization(p=p, M=M, d=d, params=ret, constraints=constraints, change_to="intercept"))
  }
}

