
#' @title Vectorization operator.
#'
#' @description \code{vec} stacks colums of the given matrix to form a vector.
#'
#' @param A a size \eqn{(dxd)} square matrix to be vectorized.
#' @return a vector of size \eqn{(d^2x1)}.
#' @section Warning:
#'  No argument checks!

vec <- function(A) {
  as.vector(A)
}


#' @title Reverse vectorization operator.
#'
#' @description \code{unvec} reverse operator for \code{vec()}.
#'
#' @description \code{unvec} forms a square matrix from a vector of
#'  stacked colums, stacked by \code{vec()}.
#'
#' @param a a size \eqn{(d^2x1)} vector to be unvectorized into a \eqn{(dxd)} matrix.
#' @param d number of rows the square matrix to be formed.
#' @return a matrix of size \eqn{(dxd)}.
#' @section Warning:
#'  No argument checks!

unvec <- function(d, a) {
  matrix(a, nrow=d, byrow=FALSE)
}


#' @title Parsimonious vectorization operator for symmetric matrices.
#'
#' @description \code{vech} stacks colums of the given matrix from main diagonal
#'   downwards (including the main diagonal) to form a vector.
#'
#' @description \code{vech} stacks colums of the given matrix from
#'   the principal diagonal downwards (including elements on the diagonal) to form a vector.
#'
#' @param A a size \eqn{(dxd)} symmetric matrix to be vectorized parsimoniously.
#' @return a vector of size \eqn{(d(d+1)/2x1)}.
#' @section Warning:
#'  No argument checks!

vech <- function(A) {
  A[lower.tri(x=A, diag=TRUE), drop=TRUE]
}


#' @title Reverse operator of the parsimonious vectorization operator \code{vech()}.
#'
#' @description \code{unvech} creates a symmetric matrix from the given vector by
#'   copying the lower triangular part to be the upper triangular part as well.
#'
#' @param a a size \eqn{(d(d+1)/2x1)} vector to be unvectorized into a symmetric \eqn{(dxd)} matrix.
#' @param d number of rows the square matrix to be formed.
#' @return a symmetric matrix of size \eqn{(dxd)}.
#' @section Warning:
#'  no argument checks!

unvech <- function(d, a) {
  A <- matrix(nrow=d, ncol=d)
  upA <- upper.tri(A)
  A[!upA] <- a
  A[upA] <- t(A)[upA]
  A
}
