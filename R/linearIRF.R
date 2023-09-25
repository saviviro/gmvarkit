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
#' @param shock_size a non-zero scalar value specifying the common size for all scalar
#'   components of the structural shock.
#' @param N a positive integer specifying the horizon how far ahead should the
#'   linear impulse responses be calculated.
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
#' @param plot_res \code{TRUE} if the results should be plotted, \code{FALSE} if
#'   not.
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
#' @references FILL IN
#' @examples
#'  FILL IN
#' @export

linear_IRF <- function(gsmvar, which_shocks, shock_size=1, N=30, regime=1, which_cumulative=numeric(0),
                       scale=NULL, scale_type=c("instant", "peak"), scale_horizon=N,
                       ci=c(0.95, 0.80), plot_res=TRUE, seed, ...) {
  NULL
}
