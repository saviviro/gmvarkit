
#' @title Vectorization operator
#'
#' @description \code{vec} stacks columns of the given matrix to form a vector.
#'
#' @param A a size \eqn{(dxd)} square matrix to be vectorized.
#' @return a vector of size \eqn{(d^2x1)}.
#' @section Warning:
#'  No argument checks!

vec <- function(A) {
  as.vector(A)
}


#' @title Vectorization operator that removes zeros
#'
#' @description \code{Wvec} stacks columns of the given matrix to form a vector
#'   and removes elements that are zeros.
#'
#' @param W a size \eqn{(dxd)} square matrix to be vectorized.
#' @return a vector of length \eqn{d^2 - n_zeros} where \eqn{n_zeros} is the
#'   number of zero entries in the matrix \code{W}.
#' @section Warning:
#'  No argument checks!

Wvec <- function(W) {
  W[W != 0]
}


#' @title Reverse vectorization operator that restores zeros
#'
#' @description \code{unWvec} forms a square matrix from a vector of
#'   stacked columns where zeros are removed according to structural
#'   parameter constaints.
#'
#' @inheritParams loglikelihood_int
#' @inheritParams unvec
#' @param Wvector a length \eqn{d^2 - n_zeros} vector where \eqn{n_zeros} is the
#'   number of zero entries in the matrix \code{W}.
#' @return a \eqn{(d x d)} matrix \eqn{W}.
#' @section Warning:
#'  No argument checks!

unWvec <- function(Wvector, d, structural_pars=NULL) {
  if(is.null(structural_pars)) stop("Structural parameters needed")
  W <- structural_pars$W
  n_zeros <- sum(W == 0, na.rm=TRUE)
  new_W <- numeric(d^2)
  new_W[W != 0] <- Wvector
  matrix(new_W, nrow=d, byrow=FALSE)
}


#' @title Reverse vectorization operator
#'
#' @description \code{unvec} forms a square matrix from a vector of
#'  stacked columns, stacked by \code{vec}.
#'
#' @param a a size \eqn{(d^2x1)} vector to be unvectorized into a \eqn{(dxd)} matrix.
#' @param d the number of rows in the square matrix to be formed.
#' @return a matrix of size \eqn{(dxd)}.
#' @section Warning:
#'  No argument checks!

unvec <- function(d, a) {
  matrix(a, nrow=d, byrow=FALSE)
}


#' @title Parsimonious vectorization operator for symmetric matrices
#'
#' @description \code{vech} stacks columns of the given matrix from
#'   the principal diagonal downwards (including elements on the diagonal) to form a vector.
#'
#' @param A a size \eqn{(dxd)} symmetric matrix to be vectorized parsimoniously.
#' @return a vector of size \eqn{(d(d+1)/2x1)}.
#' @section Warning:
#'  No argument checks!

vech <- function(A) {
  A[lower.tri(x=A, diag=TRUE), drop=TRUE]
}


#' @title Reverse operator of the parsimonious vectorization operator \code{vech}
#'
#' @description \code{unvech} creates a symmetric matrix from the given vector by
#'   copying the lower triangular part to be the upper triangular part as well.
#'
#' @param a a size \eqn{(d(d+1)/2x1)} vector to be unvectorized into a symmetric \eqn{(dxd)} matrix.
#' @param d number of rows the square matrix to be formed.
#' @return a symmetric matrix of size \eqn{(dxd)}.
#' @section Warning:
#'  No argument checks!

unvech <- function(d, a) {
  A <- matrix(nrow=d, ncol=d)
  upA <- upper.tri(A)
  A[!upA] <- a
  A[upA] <- t(A)[upA]
  A
}


#' @title Simultaneously diagonalize two covariance matrices
#'
#' @description \code{diag_Omegas} Simultaneously diagonalizes two covariance matrices using
#'   eigenvalue decomposition.
#'
#' @param Omega1 a positive definite \eqn{(dxd)} covariance matrix \eqn{(d>1)}
#' @param Omega2 another positive definite \eqn{(dxd)} covariance matrix
#' @details See the return value and Muirhead (1982), Theorem A9.9 for details.
#' @return Returns a length \eqn{d^2 + d} vector where the first \eqn{d^2} elements
#'   are \eqn{vec(W)} with the columns of \eqn{W} being (specific) eigenvectors of
#'   the matrix \eqn{\Omega_2\Omega_1^{-1}} and the rest \eqn{d} elements are the
#'   corresponding eigenvalues "lambdas". The result satisfies \eqn{WW' = Omega1} and
#'   \eqn{Wdiag(lambdas)W' = Omega2}.
#'
#'   If \code{Omega2} is not supplied, returns a vectorized symmetric (and pos. def.)
#'   square root matrix of \code{Omega1}.
#' @section Warning:
#'  No argument checks! Does not work with dimension \eqn{d=1}!
#' @references
#' \itemize{
#'   \item Muirhead R.J. 1982. Aspects of Multivariate Statistical Theory, \emph{Wiley}.
#' }
#' @examples
#' d <- 2
#' W0 <- matrix(1:(d^2), nrow=2)
#' lambdas0 <- 1:d
#' (Omg1 <- W0%*%t(W0))
#' (Omg2 <- W0%*%diag(lambdas0)%*%t(W0))
#' res <- diag_Omegas(Omg1, Omg2)
#' W <- matrix(res[1:(d^2)], nrow=d, byrow=FALSE)
#' tcrossprod(W) # == Omg1
#' lambdas <- res[(d^2 + 1):(d^2 + d)]
#' W%*%diag(lambdas)%*%t(W) # == Omg2
#' @export

diag_Omegas <- function(Omega1, Omega2) {
  eig1 <- eigen(Omega1, symmetric=TRUE)
  D <- diag(eig1$values) # Pos. def.
  H <- eig1$vectors # Orthogonal
  sqrt_omg1 <- H%*%sqrt(D)%*%t(H) # Symmetric and pos. def.
  if(missing(Omega2)) return(vec(sqrt_omg1))
  inv_sqrt_omg1 <- solve(sqrt_omg1)
  eig2 <- eigen(inv_sqrt_omg1%*%Omega2%*%inv_sqrt_omg1, symmetric=TRUE)
  lambdas <- eig2$values
  V <- eig2$vectors # Orthogonal
  W <- sqrt_omg1%*%V # all.equal(W%*%t(W), Omega1); all.equal(W%*%diag(lambdas)%*%t(W), Omega2)
  c(vec(W), lambdas)
}
