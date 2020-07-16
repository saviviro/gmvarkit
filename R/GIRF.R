
#' @title Estimate generalized impulse response function for a structural GMVAR model.
#'
#' @description \code{GIRF} estimate generalized impulse response function for
#'  a structural GMVAR model.
#'
#' @inheritParams simulateGMVAR
#' @inheritParams fitGMVAR
#' @param variables a numeric vector of length at most \eqn{d} (\code{=ncol(data)})
#'   and elements in \eqn{1,...,d} specifying the variables for which the GIRF
#'   should be estimated.
#' @param shock_size a vector with the same length as \code{variables} specifying
#'   the size of the structural shock for each variable. By default, the size of one standard
#'   deviation is used, calculated as the weighted average of the component process error term
#'   standard deviations with weights given by the mixing weight parameters.
#' @param N a positive integer specifying the horizon how far ahead should the generalized
#'   impulse responses be calculated?
#' @param R1 the number of repetitions used to estimate GIRF for each initial value?
#' @param R2 the number of initial values drawn from a stationary distribution of the process
#'   or of a specific regime? Ignored if the argument \code{init_value} is specified.
#' @param init_regimes a numeric vector of length at most \eqn{M} and elements in \eqn{1,...,M}
#'   specifying the regimes from which the initial values should be generated from. The initial
#'   values will be generated from a mixture distribution with the mixture components being the
#'   stationary distirbutions of the specific regimes and the (proportiional) mixing weights given
#'   by the mixing weight parameters of those regimes. Note that if \code{init_regimes=1:M}, the
#'   initial values are generated from the stationary distribution of the process and if
#'   \code{init_regimes=m}, the initial value are generated from the stationary distribution
#'   of the \eqn{m}th regime. Ignored if \code{init_value} is specified.
#' @param init_values a matrix or a multivariate class \code{'ts'} object with \eqn{d} columns
#'   and at least \eqn{p} rows specifying an initial value for the GIRF. The last \eqn{p} rows
#'   are taken to be the initial value assuming that the \strong{last} row is the most recent observation.
#' @param include_mixweights should the generalized impulse response be calculated for the mixing weights
#'   as well?
#' @details The model needs to be structural in order to use this function. A structural GMVAR model can
#'   be estimated by specifying the argument \code{structural_pars} in the function \code{fitGMVAR}.
#'   See the cited paper by Virolainen (2020) for details about the algorithm.
#' @inherit in_paramspace_int references
#' @examples
#'  \donttest{
#'  # To be filled in...
#'  }
#' @export

GIRF <- function(gmvar, variables, shock_size, N=10, R1=20, R2=20, init_regimes=1:M, init_values=NULL,
                 include_mixweights=TRUE, ncores=min(2, parallel::detectCores())) {
  d <- gmvar$model$d
  M <- gmvar$model$M
  stopifnot(!is.null(gmvar$model$structural_pars))
  stopifnot(length(variables) <= d && all(variables %in% 1:d) && length(unique(variables)) == length(variables))
  stopifnot(length(init_regimes) <= M && all(initregimes %in% 1:M) && length(unique(init_regimes)) == length(init_regimes))

  get_one_girf <- function() {
    NULL
    # call simulateGMVAR with init_values=init_values, nsimu=N
  }



}
