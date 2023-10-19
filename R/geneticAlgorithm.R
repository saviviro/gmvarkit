#' @title Genetic algorithm for preliminary estimation of a GMVAR, StMVAR, or G-StMVAR model
#'
#' @description \code{GAfit} estimates the specified GMVAR, StMVAR, or G-StMVAR model using a genetic algorithm.
#'   It's designed to find starting values for gradient based methods.
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
#'     \item{\strong{For structural models:}}{
#'       Should have the form
#'       \strong{\eqn{\theta}}\eqn{ = (\phi_{1,0},...,\phi_{M,0},}\strong{\eqn{\phi}}\eqn{_{1},...,}\strong{\eqn{\phi}}\eqn{_{M},
#'       vec(W),}\strong{\eqn{\lambda}}\eqn{_{2},...,}\strong{\eqn{\lambda}}\eqn{_{M},\alpha_{1},...,\alpha_{M-1},}\strong{\eqn{\nu}}\eqn{)},
#'       where
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
#'         \item{\strong{If \eqn{\lambda_{mi}} are constrained:}}{Replace \strong{\eqn{\lambda}}\eqn{_{2},...,}\strong{\eqn{\lambda}}\eqn{_{M}}
#'          with \strong{\eqn{\gamma}} \eqn{(rx1)} that satisfies (\strong{\eqn{\lambda}}\eqn{_{2}}\eqn{,...,}
#'         \strong{\eqn{\lambda}}\eqn{_{M}) =} \strong{\eqn{C_{\lambda} \gamma}} where \eqn{C_{\lambda}} is a \eqn{(d(M-1) x r)}
#'          constraint matrix.}
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
#' @param conditional a logical argument specifying whether the conditional or exact log-likelihood function
#' @param mu_scale a size \eqn{(dx1)} vector defining \strong{means} of the normal distributions from which each
#'   mean parameter \eqn{\mu_{m}} is drawn from in random mutations. Default is \code{colMeans(data)}. Note that
#'   mean-parametrization is always used for optimization in \code{GAfit} - even when \code{parametrization=="intercept"}.
#'   However, input (in \code{initpop}) and output (return value) parameter vectors can be intercept-parametrized.
#' @param mu_scale2 a size \eqn{(dx1)} strictly positive vector defining \strong{standard deviations} of the normal
#'   distributions from which each mean parameter \eqn{\mu_{m}} is drawn from in random mutations.
#'   Default is \code{2*sd(data[,i]), i=1,..,d}.
#' @param omega_scale a size \eqn{(dx1)} strictly positive vector specifying the scale and variability of the
#'   random covariance matrices in random mutations. The covariance matrices are drawn from (scaled) Wishart
#'   distribution. Expected values of the random covariance matrices are \code{diag(omega_scale)}. Standard
#'   deviations of the diagonal elements are \code{sqrt(2/d)*omega_scale[i]}
#'   and for non-diagonal elements they are \code{sqrt(1/d*omega_scale[i]*omega_scale[j])}.
#'   Note that for \code{d>4} this scale may need to be chosen carefully. Default in \code{GAfit} is
#'   \code{var(stats::ar(data[,i], order.max=10)$resid, na.rm=TRUE), i=1,...,d}. This argument is ignored if
#'   structural model is considered.
#' @param W_scale a size \eqn{(dx1)} strictly positive vector partly specifying the scale and variability of the
#'   random covariance matrices in random mutations. The elements of the matrix \eqn{W} are drawn independently
#'   from such normal distributions that the expectation of the main \strong{diagonal} elements of the first
#'   regime's error term covariance matrix \eqn{\Omega_1 = WW'} is \code{W_scale}. The distribution of \eqn{\Omega_1}
#'   will be in some sense like a Wishart distribution but with the columns (elements) of \eqn{W} obeying the given
#'   constraints. The constraints are accounted for by setting the element to be always zero if it is subject to a zero
#'   constraint and for sign constraints the absolute value or negative the absolute value are taken, and then the
#'   variances of the elements of \eqn{W} are adjusted accordingly. This argument is ignored if reduced form model
#'   is considered.
#' @param lambda_scale a length \eqn{M - 1} vector specifying the \strong{standard deviation} of the mean zero normal
#'   distribution from which the eigenvalue \eqn{\lambda_{mi}} parameters are drawn from in random mutations.
#'   As the eigenvalues should always be positive, the absolute value is taken. The elements of \code{lambda_scale}
#'   should be strictly positive real numbers with the \eqn{m-1}th element giving the degrees of freedom for the \eqn{m}th
#'   regime. The expected value of the main \strong{diagonal} elements \eqn{ij} of the \eqn{m}th \eqn{(m>1)} error term covariance
#'   matrix will be \code{W_scale[i]*(d - n_i)^(-1)*sum(lambdas*ind_fun)} where the \eqn{(d x 1)} vector \code{lambdas} is
#'   drawn from the absolute value of the t-distribution, \code{n_i} is the number of zero constraints in the \eqn{i}th
#'   row of \eqn{W} and \code{ind_fun} is an indicator function that takes the value one iff the \eqn{ij}th element of
#'   \eqn{W} is not constrained to zero. Basically, larger lambdas (or smaller degrees of freedom) imply larger variance.
#'
#'   If the lambda parameters are \strong{constrained} with the \eqn{(d(M - 1) x r)} constraint matrix \eqn{C_lambda},
#'   then provide a length \eqn{r} vector specifying the standard deviation of the (absolute value of the) mean zero
#'   normal distribution each of the \eqn{\gamma} parameters are drawn from (the \eqn{\gamma} is a \eqn{(r x 1)} vector).
#'   The expected value of the main diagonal elements of the covariance matrices then depend on the constraints.
#'
#'   This argument is ignored if \eqn{M==1} or a reduced form model is considered. Default is \code{rep(3, times=M-1)}
#'   if lambdas are not constrained and \code{rep(3, times=r)} if lambdas are constrained.
#'
#'   As with omega_scale and W_scale, this argument should be adjusted carefully if specified by hand. \strong{NOTE}
#'   that if lambdas are constrained in some other way than restricting some of them to be identical, this parameter
#'   should be adjusted accordingly in order to the estimation succeed!
#' @param ar_scale a positive real number adjusting how large AR parameter values are typically proposed in construction
#'   of the initial population: larger value implies larger coefficients (in absolute value). After construction of the
#'   initial population, a new scale is drawn from \code{(0, 0.)} uniform distribution in each iteration.
#' @param upper_ar_scale the upper bound for \code{ar_scale} parameter (see above) in the random mutations. Setting
#'  this too high might lead to failure in proposing new parameters that are well enough inside the parameter space,
#'  and especially with large \code{p} one might want to try smaller upper bound (e.g., 0.5).
#' @param ar_scale2 a positive real number adjusting how large AR parameter values are typically proposed in some
#'   random mutations (if AR constraints are employed, in all random mutations): larger value implies \strong{smaller} coefficients
#'   (in absolute value). \strong{Values larger than 1 can be used if the AR coefficients are expected to be very small.
#'   If set smaller than 1, be careful as it might lead to failure in the creation of stationary parameter candidates}
#' @param regime_force_scale a non-negative real number specifying how much should natural selection favor individuals
#'   with less regimes that have almost all mixing weights (practically) at zero. Set to zero for no favoring or large
#'   number for heavy favoring. Without any favoring the genetic algorithm gets more often stuck in an area of the
#'   parameter space where some regimes are wasted, but with too much favouring the best genes might never mix into
#'   the population and the algorithm might converge poorly. Default is \code{1} and it gives \eqn{2x} larger surviving
#'   probability weights for individuals with no wasted regimes compared to individuals with one wasted regime.
#'   Number \code{2} would give \eqn{3x} larger probability weights etc.
#' @param red_criteria a length 2 numeric vector specifying the criteria that is used to determine whether a regime is
#'   redundant (or "wasted") or not.
#'   Any regime \code{m} which satisfies \code{sum(mixingWeights[,m] > red_criteria[1]) < red_criteria[2]*n_obs} will
#'   be considered "redundant". One should be careful when adjusting this argument (set \code{c(0, 0)} to fully disable
#'   the 'redundant regime' features from the algorithm).
#' @param pre_smart_mu_prob A number in \eqn{[0,1]} giving a probability of a "smart mutation" occuring randomly in each
#'   iteration before the iteration given by the argument \code{smart_mu}.
#' @param to_return should the genetic algorithm return the best fitting individual which has "positive enough" mixing
#'   weights for as many regimes as possible (\code{"alt_ind"}) or the individual which has the highest log-likelihood
#'   in general (\code{"best_ind"}) but might have more wasted regimes?
#' @param minval a real number defining the minimum value of the log-likelihood function that will be considered.
#'   Values smaller than this will be treated as they were \code{minval} and the corresponding individuals will
#'   never survive. The default is \code{-(10^(ceiling(log10(n_obs)) + d) - 1)}.
#' @param seed a single value, interpreted as an integer, or NULL, that sets seed for the random number generator in the beginning of
#'   the function call. If calling \code{GAfit} from \code{fitGSMVAR}, use the argument \code{seeds} instead of passing the argument
#'   \code{seed}.
#' @details
#'  The core of the genetic algorithm is mostly based on the description by \emph{Dorsey and Mayer (1995)}.
#'  It utilizes a slightly modified version of the individually adaptive crossover and mutation rates described
#'  by \emph{Patnaik and Srinivas (1994)} and employs (50\%) fitness inheritance discussed by
#'  \emph{Smith, Dike and Stegmann (1995)}.
#'
#'  By "redundant" or "wasted" regimes we mean regimes that have the time varying mixing weights practically at
#'  zero for almost all t. A model including redundant regimes would have about the same log-likelihood value without
#'  the redundant regimes and there is no purpose to have redundant regimes in a model.
#' @return Returns the estimated parameter vector which has the form described in \code{initpop}.
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
#'    \item Virolainen S. 2022. Structural Gaussian mixture vector autoregressive model with application to the asymmetric
#'      effects of monetary policy shocks. Unpublished working paper, available as arXiv:2007.04713.
#'    \item Virolainen S. 2022. Gaussian and Student's t mixture vector autoregressive model with application to the
#'      asymmetric effects of monetary policy shocks in the Euro area. Unpublished working
#'      paper, available as arXiv:2109.13648.
#'  }
#' @examples
#' \donttest{
#'  # Preliminary estimation of a G-StMVAR(1, 1, 1) model with 50 generations.
#'  GA_estimates <- GAfit(gdpdef, p=1, M=c(1, 1), model="G-StMVAR",
#'                        ngen=50, seed=1)
#'  GA_estimates
#' }
#' @export

GAfit <- function(data, p, M, model=c("GMVAR", "StMVAR", "G-StMVAR"), conditional=TRUE, parametrization=c("intercept", "mean"),
                  constraints=NULL, same_means=NULL, weight_constraints=NULL, structural_pars=NULL,
                  ngen=200, popsize, smart_mu=min(100, ceiling(0.5*ngen)), initpop=NULL,
                  mu_scale, mu_scale2, omega_scale, W_scale, lambda_scale,
                  ar_scale=0.2, upper_ar_scale=1, ar_scale2=1, regime_force_scale=1,
                  red_criteria=c(0.05, 0.01), pre_smart_mu_prob=0, to_return=c("alt_ind", "best_ind"),
                  minval, seed=NULL) {

  # Required values and preliminary checks
  set.seed(seed)
  model <- match.arg(model)
  to_return <- match.arg(to_return)
  parametrization <- match.arg(parametrization)
  check_pMd(p=p, M=M, model=model)
  data <- check_data(data=data, p=p)
  d <- ncol(data)
  n_obs <- nrow(data)
  check_same_means(parametrization=parametrization, same_means=same_means)
  check_constraints(p=p, M=M, d=d, constraints=constraints, same_means=same_means,
                    weight_constraints=weight_constraints, structural_pars=structural_pars)
  npars <- n_params(p=p, M=M, d=d, model=model, constraints=constraints, same_means=same_means,
                    weight_constraints=weight_constraints, structural_pars=structural_pars)
  M_orig <- M
  M <- sum(M)

  # For structural models, determine whether W constraints such that W and lambdas can be sorted for better estimation
  if(is.null(structural_pars) || !is.null(structural_pars$C_lambda) || !is.null(structural_pars$fixed_lambdas)) {
    sort_structural_pars <- FALSE
  } else {
    W <- structural_pars$W
    constrained_rows <- apply(W, MARGIN=1, FUN=function(x) !all(is.na(x)))
    if(sum(constrained_rows) > 1) { # Overidentifying constraints employed
      sort_structural_pars <- FALSE
    } else if(sum(constrained_rows) == 0) { # No constraints on W: can be sorted
      sort_structural_pars <- TRUE
    } else { # one constrained row: can be sorted if they all have the same sign constraint
      if(all(W[constrained_rows,] > 0) || all(W[constrained_rows,] < 0)) {
        sort_structural_pars <- TRUE
      } else {
        sort_structural_pars <- FALSE
      }
    }
  }

  # Defaults and checks
  if(!all_pos_ints(c(ngen, smart_mu))) stop("Arguments ngen and smart_mu have to be positive integers")
  if(missing(popsize)) {
    popsize <- 50*ceiling(sqrt(npars))
  } else if(popsize < 2 | popsize %% 2 != 0) {
    stop("The population size popsize must be even positive integer")
  }
  if(missing(minval)) {
    minval <- get_minval(data)
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
  if(missing(W_scale)) {
    W_scale <- omega_scale
  } else if(!(length(W_scale) == d & all(W_scale > 0))) {
    stop("W_scale must be numeric vector with length d and positive elements")
  }
  stopifnot(pre_smart_mu_prob >= 0 && pre_smart_mu_prob <= 1)
  if(is.null(structural_pars$C_lambda)) {
    n_lambs <- M - 1
  } else {
    n_lambs <- ncol(structural_pars$C_lambda)
  }
  if(missing(lambda_scale)) {
    lambda_scale <- rep(3, times=n_lambs) # numeric(0) if M == 1
  } else if(!(length(lambda_scale) == n_lambs & all(lambda_scale > 0.1))) {
    stop("lambda_scale must be numeric vector with length M-1 (r if lambdas are constrained) and elements larger than 0.1")
  }
  if(length(ar_scale) != 1 | ar_scale <= 0) {
    stop("ar_scale must be positive and have length one")
  }
  if(length(ar_scale2) != 1 | ar_scale2 <= 0) {
    stop("ar_scale2 must be positive and have length one")
  } else if(ar_scale2 > 1.5) {
    warning("Large ar_scale2 might lead to failure of the estimation process")
  }
  if(length(regime_force_scale) != 1 | regime_force_scale < 0) {
    stop("regime_force_scale should be non-negative real number")
  }

  # The initial population
  if(is.null(initpop)) {
    nattempts <- 20
    G <- numeric(0)
    for(i1 in 1:nattempts) {
      if(is.null(constraints)) {
        inds <- replicate(popsize, random_ind2(p=p, M=M_orig, d=d, model=model,
                                               same_means=same_means,
                                               weight_constraints=weight_constraints,
                                               structural_pars=structural_pars,
                                               mu_scale=mu_scale,
                                               mu_scale2=mu_scale2,
                                               omega_scale=omega_scale,
                                               ar_scale=ar_scale,
                                               W_scale=W_scale,
                                               lambda_scale=lambda_scale))
      } else {
        inds <- replicate(popsize, random_ind(p=p, M=M_orig, d=d, model=model,
                                              constraints=constraints,
                                              same_means=same_means,
                                              weight_constraints=weight_constraints,
                                              structural_pars=structural_pars,
                                              mu_scale=mu_scale,
                                              mu_scale2=mu_scale2,
                                              omega_scale=omega_scale,
                                              W_scale=W_scale,
                                              ar_scale2=ar_scale2,
                                              lambda_scale=lambda_scale))
      }
      ind_loks <- vapply(1:popsize, function(i2) loglikelihood_int(data=data, p=p, M=M_orig, params=inds[,i2], model=model,
                                                                   conditional=conditional, parametrization="mean",
                                                                   constraints=constraints, same_means=same_means,
                                                                   weight_constraints=weight_constraints,
                                                                   structural_pars=structural_pars,
                                                                   check_params=TRUE, to_return="loglik",
                                                                   minval=minval), numeric(1))
      G <- cbind(G, inds[, ind_loks > minval]) # Take good enough individuals
      if(ncol(G) >= popsize) {
        G <- G[, 1:popsize]
        break
      } else if(i1 == nattempts) {
        if(length(G) == 0) {
          stop(paste("Failed to create initial population with good enough individuals.",
                     "Scaling the individual series so that the AR coefficients (of a VAR model)",
                     "will not be very large (preferably less than one) may solve the problem.",
                     "If needed, another package may be used to fit linear VARs so see which",
                     "scalings produce relatively small AR coefficient estimates."))
        } else {
          G <- G[, sample.int(ncol(G), size=popsize, replace=TRUE)]
        }
      }
    }
  } else { # Initial population set by the user
    stopifnot(is.list(initpop))
    for(i1 in 1:length(initpop)) {
      ind <- initpop[[i1]]
      tryCatch(check_parameters(p=p, M=M_orig, d=d, params=ind, model=model, constraints=constraints,
                                parametrization=parametrization, same_means=same_means,
                                weight_constraints=weight_constraints, structural_pars=structural_pars),
               error=function(e) stop(paste("Problem with individual", i1, "in the initial population: "), e))
      if(parametrization == "intercept") { # This is never the case when !is.null(same_means)
        ind <- change_parametrization(p=p, M=M_orig, d=d, params=ind, model=model, constraints=constraints,
                                      weight_constraints=weight_constraints, structural_pars=structural_pars,
                                      change_to="mean")
      }
      if(is.null(constraints) && is.null(structural_pars$C_lambda) && is.null(structural_pars$fixed_lambdas) &&
         is.null(same_means) && is.null(weight_constraints)) {
        initpop[[i1]] <- sort_components(p=p, M=M_orig, d=d, params=ind, model=model, structural_pars=structural_pars)
      } else {
        initpop[[i1]] <- ind
      }
    }
    G <- replicate(popsize, initpop[[sample.int(length(initpop), size=1)]])
  }

  # Calculates the number of redundant regimes
  n_redundants <- function(M, mw) {
    sum(vapply(1:M, function(m) sum(mw[,m] > red_criteria[1]) < red_criteria[2]*n_obs, logical(1)))
  }

  # Initial setup
  generations <- array(NA, dim=c(npars, popsize, ngen))
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
        loks_and_mw <- loglikelihood_int(data=data, p=p, M=M_orig, params=G[,i2], model=model,
                                         conditional=conditional, parametrization="mean",
                                         constraints=constraints, same_means=same_means,
                                         weight_constraints=weight_constraints,
                                         structural_pars=structural_pars,
                                         to_return="loglik_and_mw",
                                         check_params=TRUE, minval=minval)
        fill_lok_and_red(i1, i2, loks_and_mw)
      }
    } else {
      # Proportional fitness inheritance: individual has 50% change to inherit fitness if it's a result of crossover.
      # Variable "I" tells the proportions of parent material.
      I2 <- rep(I, each=2)
      which_did_co <- which(1 - which_not_co == 1)
      if(length(which_did_co) > 0) {
        which_inherit <- sample(x=which_did_co, size=round(0.5*length(which_did_co)), replace=FALSE)
      } else {
        which_inherit <- numeric(0)
      }
      # survivor_liks holds the parent loglikelihood values: for odd number they are (index, index+1) and for even (index-1, index)
      if(length(which_inherit) > 0 & abs(max_lik - mean_lik) > abs(0.03*mean_lik)) { # No inheritance if massive mutations
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
              loks_and_mw <- tryCatch(loglikelihood_int(data=data, p=p, M=M_orig, params=G[,i2], model=model,
                                                        conditional=conditional, parametrization="mean",
                                                        constraints=constraints, same_means=same_means,
                                                        weight_constraints=weight_constraints,
                                                        structural_pars=structural_pars,
                                                        to_return="loglik_and_mw",
                                                        check_params=FALSE, minval=minval),
                                      error=function(e) minval)
            } else {
              loks_and_mw <- tryCatch(loglikelihood_int(data=data, p=p, M=M_orig, params=G[,i2], model=model,
                                                        conditional=conditional, parametrization="mean",
                                                        constraints=constraints, same_means=same_means,
                                                        weight_constraints=weight_constraints,
                                                        structural_pars=structural_pars,
                                                        to_return="loglik_and_mw",
                                                        check_params=TRUE, minval=minval),
                                      error=function(e) minval)
            }
            fill_lok_and_red(i1, i2, loks_and_mw)
         }
      }
    }

    # Take care of individuals that are not good enough + calculate the numbers redundant regimes
    if(anyNA(logliks[i1,])) {
      which_inds_na <- which(is.na(logliks[i1,]))
      logliks[i1, which_inds_na] <- minval
      warning(paste0("Removed NaNs from the logliks of inds ", which_inds_na, " in the GA round ", i1,
                     " with the seed ", seed, "."))
    }
    logliks[i1, which(logliks[i1,] < minval)] <- minval
    redundants[i1, which(logliks[i1,] <= minval)] <- M

    ## Selection and the reproduction pool ##
    if(length(unique(logliks[i1,])) == 1) {
      choosing_probs <- rep(1, popsize) # If all individuals are the same, the surviving probability weight is 1.
    } else {
      T_values <- logliks[i1,] + abs(min(logliks[i1,])) # Function T as described by Dorsey R. E. ja Mayer W. J., 1995
      T_values <- T_values/(1 + regime_force_scale*redundants[i1,]) # Favor individuals with least number of redundant regimes
      choosing_probs <- T_values/sum(T_values) # The surviving probability weights
    }

    # Draw popsize individuals with replacement and form the reproduction pool H.
    survivors <- sample(1:popsize, size=popsize, replace=TRUE, prob=choosing_probs)
    H <- G[,survivors]

    # Calculate mean and max log-likelihood of the survivors
    survivor_liks <- logliks[i1, survivors]
    survivor_redundants <- redundants[i1, survivors]
    max_lik <- max(survivor_liks)
    mean_lik <- mean(survivor_liks)
    if(max_lik == mean_lik) mean_lik <- mean_lik + 0.1 # We avoid dividing by zero when all the individuals are the same

    ## Cross-overs ##
    # Individually adaptive cross-over rates as described by Patnaik and Srinivas (1994) with the modification of
    # setting the crossover rate to be at least 0.4 for all individuals (so that the best genes mix in the population too).
    indeces <- seq(from=1, to=popsize - 1, by=2)
    parent_max <- vapply(indeces, function(i2) max(survivor_liks[i2], survivor_liks[i2+1]), numeric(1))
    co_rates <- vapply(1:length(indeces), function(i2) max(min((max_lik - parent_max[i2])/(max_lik - mean_lik), 1), 0.4), numeric(1))

    # Do the crossovers
    which_co <- rbinom(n=popsize/2, size=1, prob=co_rates)
    I <- round(runif(n=popsize/2, min=0.5 + 1e-16, max=npars - 0.5 - 1e-16)) # Break points
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
    # First gen when the best loglik occurred (because of fitness inheritance):
    best_index <- best_index0[order(best_index0[,1], decreasing=FALSE)[1],]
    best_ind <- generations[, best_index[2], best_index[1]]
    best_mw <- loglikelihood_int(data=data, p=p, M=M_orig, params=best_ind, model=model,
                                 conditional=conditional, parametrization="mean",
                                 constraints=constraints, same_means=same_means,
                                 weight_constraints=weight_constraints,
                                 structural_pars=structural_pars,
                                 to_return="mw",
                                 check_params=FALSE, minval=minval)
    # Which regimes are wasted:
    which_redundant <- which(vapply(1:M, function(i2) sum(best_mw[,i2] > red_criteria[1]) < red_criteria[2]*n_obs, logical(1)))

    # Keep track of "the alternative best individual" that has (weakly) less reduntant regimes than the current best one.
    if(length(which_redundant) <= length(which_redundant_alt)) {
      alt_ind <- best_ind
      which_redundant_alt <- which_redundant
    }


    ## Mutations ##
    which_not_co <- rep(1 - which_co, each=2)
    if(abs(max_lik - mean_lik) <= abs(0.03*mean_lik)) {
      mu_rates <- rep(0.7, popsize) # Massive mutations if converging
    } else {
      # Individually adaptive mutation rates, Patnaik and Srinivas (1994); we only mutate those who did not crossover.
      mu_rates <- 0.5*vapply(1:popsize, function(i2) min(which_not_co[i2], (max_lik - survivor_liks[i2])/(max_lik - mean_lik)), numeric(1))
    }

    # Do mutations and keep track if they are stationary for sure
    mutate <- rbinom(n=popsize, size=1, prob=mu_rates)
    which_mutate <- which(mutate == 1)
    pre_smart_mu <- runif(1, min=1e-6, max=1-1e-6) < pre_smart_mu_prob
    ar_scale <- runif(1, min=1e-6, max=upper_ar_scale - 1e-6) # Random AR scale
    if(i1 <= smart_mu & length(which_mutate) >= 1 & !pre_smart_mu) { # Random mutations
      if(!is.null(constraints) | runif(1, min=1e-6, max=1 - 1e-6) > 0.5) { # Regular, can be nonstationary
        stat_mu <- FALSE
        H2[,which_mutate] <- vapply(1:length(which_mutate), function(x) random_ind(p=p, M=M_orig, d=d, model=model,
                                                                                   constraints=constraints,
                                                                                   same_means=same_means,
                                                                                   weight_constraints=weight_constraints,
                                                                                   structural_pars=structural_pars,
                                                                                   mu_scale=mu_scale,
                                                                                   mu_scale2=mu_scale2,
                                                                                   omega_scale=omega_scale,
                                                                                   W_scale=W_scale,
                                                                                   lambda_scale=lambda_scale,
                                                                                   ar_scale2=ar_scale2), numeric(npars))
      } else { # For stationarity with algorithm (slower but can skip stationarity test), Ansley and Kohn (1986)
        stat_mu <- TRUE
        H2[,which_mutate] <- vapply(1:length(which_mutate), function(x) random_ind2(p=p, M=M_orig, d=d, model=model,
                                                                                    same_means=same_means,
                                                                                    weight_constraints=weight_constraints,
                                                                                    structural_pars=structural_pars,
                                                                                    mu_scale=mu_scale,
                                                                                    mu_scale2=mu_scale2,
                                                                                    omega_scale=omega_scale,
                                                                                    ar_scale=ar_scale,
                                                                                    W_scale=W_scale,
                                                                                    lambda_scale=lambda_scale), numeric(npars))
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

      ## 'Smart mutation': mutate close to a well fitting individual. We obviously don't mutate close to
      # redundant regimes but draw them at random ('rand_to_use' in what follows).
      if(!is.null(constraints) | !is.null(same_means) | !is.null(structural_pars) |
         length(which_redundant) <= length(which_redundant_alt) | runif(1) > 0.5) {
        # The first option for smart mutations: smart mutate to 'alt_ind' which is the best fitting individual
        # with the least redundant regimes.
        # Note that best_ind == alt_ind when length(which_redundant) <= length(which_redundant_alt).
        ind_to_use <- alt_ind
        rand_to_use <- which_redundant_alt
      } else {
        # Alternatively, if there exists an alternative individual with strictly less redundant regimes
        # than in the best_ind, a "regime combining procedure" might take place: take a redundant regime
        # of the best_ind and replace it with a nonredundant regime taken from alt_ind. Then do smart
        # mutation close to this new individual. For simplicity, regime combining is not considered for
        # models imposing linear constraints.

        # We want to take such nonredundant regime from alt_ind that is not similar to the nonredundant
        # regimes of best_ind. In order to choose such regime, we compare all nonredundant regimes of
        # best_ind to all nonredundant regimes of alt_ind. Then, we choose the nonredundant regime of
        # alt_ind which has the largest "distance" to the closest regime of best_ind.

        # In the case of StMVAR and G-StMVAR model, we remove the degrees of freedom parameters from
        # the comparison and treat them equally to the GMVAR type regimes. After doing the regime
        # combining procedure, we than add a random degrees of freedom parameter to the regime.
        # This way, overly large df parameters do not mess up the comparison, as StMVAR type regime
        # with overly degrees of freedom parameter is essentially a GMVAR type regime with a redundant
        # parameter. The trade-off is that a StMVAR type regime with small df parameter might not be
        # a good GMVAR type regime with the df parameter removed (or GMVAR regime might not be a
        # good StMVAR type regime with a df parameter added to it).

        # Choose regime of best_ind to be changed
        which_to_change <- which_redundant[1]

        # Pick the nonredundant regimes of best_ind and alt_ind
        non_red_regs_best <- vapply((1:M)[-which_redundant], function(i2) pick_regime(p=p, M=M_orig, d=d, params=best_ind,
                                                                                      model=model, m=i2, with_df=FALSE),
                                    numeric(p*d^2 + d + d*(d+1)/2))
        if(length(which_redundant_alt) == 0) { # Special case for technical reasons
          non_red_regs_alt <- vapply(1:M, function(i2) pick_regime(p=p, M=M_orig, d=d, params=alt_ind,
                                                                   model=model, m=i2, with_df=FALSE), numeric(p*d^2 + d + d*(d+1)/2))
        } else {
          non_red_regs_alt <- vapply((1:M)[-which_redundant_alt], function(i2) pick_regime(p=p, M=M_orig, d=d, params=alt_ind,
                                                                                           model=model, m=i2, with_df=FALSE),
                                     numeric(p*d^2 + d + d*(d+1)/2))
        }

        # Calculate the "distances" between the nonredundant regimes
        # Row for each non-red-reg-best and column for each non-red-reg-alt:
        dist_to_regime <- matrix(nrow=ncol(non_red_regs_best), ncol=ncol(non_red_regs_alt))
        for(i2 in 1:nrow(dist_to_regime)) {
          dist_to_regime[i2,] <- vapply(1:ncol(non_red_regs_alt), function(i3) regime_distance(regime_pars1=non_red_regs_best[,i2],
                                                                                               regime_pars2=non_red_regs_alt[,i3]), numeric(1))
        }

        # Which alt_ind regime, i.e. column should be used? Choose the one that with largest distance to the closest regime
        # to avoid duplicating similar regimes
        which_reg_to_use <- which(apply(dist_to_regime, 2, min) == max(apply(dist_to_regime, 2, min)))[1]

        # The obtain the regime of alt_ind that is used to replace a redundant regime in best_ind
        if(model == "GMVAR") { # GMVAR type regime, so no df parameter is included
          reg_to_use <- non_red_regs_alt[,which_reg_to_use]
        } else if(model == "StMVAR") { # Obtain the regime, including the degrees of freedom parameter if StMVAR type
          reg_to_use <-  pick_regime(p=p, M=M_orig, d=d, params=alt_ind, model=model, m=i2, with_df=TRUE)
        } else { # model == "G-StMVAR"
          if(which_to_change <= M_orig[1]) { # GMVAR type regime, so no df parameter is included
            reg_to_use <- non_red_regs_alt[,which_reg_to_use]
          } else { # StMVAR type regime, so df parameter is included
            if(which_reg_to_use <= M_orig[1]) { # We need to create df parameters because the replacing regime is GMVAR type
              reg_to_use <- c(non_red_regs_alt[,which_reg_to_use],  2 + rgamma(1, shape=0.3, rate=0.007))
            } else { # We obtain the replacing regime, which already constraint a df parameter
              reg_to_use <-  pick_regime(p=p, M=M_orig, d=d, params=alt_ind, model=model, m=i2, with_df=TRUE)
            }
          }
        }

        # Change the chosen regime of best_ind to be the one chosen from alt_ind
        ind_to_use <- change_regime(p=p, M=M_orig, d=d, params=best_ind, model=model, m=which_to_change, regime_pars=reg_to_use)

        # Should some regimes still be random?
        rand_to_use <- which_redundant[which_redundant != which_to_change]
      }

      # Smart mutations
      H2[,which_mutate] <- vapply(1:length(which_mutate), function(i2) smart_ind(p=p, M=M_orig, d=d,
                                                                                 params=ind_to_use,
                                                                                 model=model,
                                                                                 constraints=constraints,
                                                                                 same_means=same_means,
                                                                                 weight_constraints=weight_constraints,
                                                                                 structural_pars=structural_pars,
                                                                                 accuracy=accuracy[i2],
                                                                                 which_random=rand_to_use,
                                                                                 mu_scale=mu_scale,
                                                                                 mu_scale2=mu_scale2,
                                                                                 omega_scale=omega_scale,
                                                                                 ar_scale=ar_scale,
                                                                                 ar_scale2=ar_scale2,
                                                                                 W_scale=W_scale,
                                                                                 lambda_scale=lambda_scale), numeric(npars))
    }

    # If structural model is considered without overidentifying constraints on the W matrix and without constraints
    # on the lambda parameters, sort the columns of the W matrix by sorting the lambdas of the second regime to
    # increasing order. This gives statistical identification for structural models without overidentifying constraints.
    if(sort_structural_pars) {
      H2 <- vapply(1:popsize, function(i2) sort_W_and_lambdas(p=p, M=M_orig, d=d, params=H2[,i2], model=model), numeric(npars))
    }

    # Sort components according to the mixing weight parameters. No sorting if constraints are employed.
    if(is.null(constraints) && is.null(structural_pars$C_lambda) && is.null(structural_pars$fixed_lambdas) &&
       is.null(same_means) && is.null(weight_constraints)) {
      H2 <- vapply(1:popsize, function(i2) sort_components(p=p, M=M_orig, d=d, params=H2[,i2], model=model,
                                                           structural_pars=structural_pars), numeric(npars))
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
  if(parametrization == "mean") { # This is always the case with same_means
    return(ret)
  } else {
    return(change_parametrization(p=p, M=M_orig, d=d, params=ret, model=model, constraints=constraints,
                                  weight_constraints=weight_constraints, structural_pars=structural_pars,
                                  change_to="intercept"))
  }
}

