context("argument checking functions")
library(gmvarkit)

# The tests are very brief !!!

## A(M)(p)_(p)(M)(d)

# p=1, M=1, d=2
phi10_112 <- c(1.03, 2.36)
A11_112 <- matrix(c(1.25, 0.06, 0.04, 1.34), nrow=2, byrow=FALSE)
Omega1_112 <- matrix(c(0.93, -0.15, -0.15, 5.20), nrow=2, byrow=FALSE)

theta_112 <- c(phi10_112, vec(A11_112), vech(Omega1_112))

W_112 <- t(chol(Omega1_112))
theta_112s <- c(phi10_112, vec(A11_112), vec(W_112)) # SGMVAR
Omega1_112s <- tcrossprod(W_112)

theta_112sWC <- c(phi10_112, vec(A11_112), Wvec(W_112)) # SGMVAR W constrained

# p=1, M=2, d=2
phi10_122 <- c(1.03, 2.36)
A11_122 <- matrix(c(0.9, 0.3, -0.3, 0.9), nrow=2, byrow=FALSE)
Omega1_122 <- matrix(c(0.93, -0.15, -0.15, 5.20), nrow=2, byrow=FALSE)

phi20_122 <- c(1.79, 3.00)
A21_122 <- A11_122
Omega2_122 <- matrix(c(5.88, 3.56, 3.56, 9.80), nrow=2, byrow=FALSE)

alpha1_122 <- 0.37
upsilon1_122 <- c(phi10_122, vec(A11_122), vech(Omega1_122))
upsilon2_122 <- c(phi20_122, vec(A21_122), vech(Omega2_122))
theta_122 <- c(upsilon1_122, upsilon2_122, alpha1_122)

WL_122 <- diag_Omegas(Omega1_122, Omega2_122)
W_122 <- matrix(WL_122[1:(2^2)], nrow=2, byrow=FALSE)
lambdas_122 <- WL_122[(2^2 + 1):length(WL_122)]
theta_122s <- c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), vec(W_122), lambdas_122, alpha1_122) # SGMVAR

# p=2, M=2, d=2
phi10_222 <- c(1.03, 2.36)
A11_222 <- matrix(c(1.25, 0.06, 0.04, 1.34), nrow=2, byrow=FALSE)
A12_222 <- matrix(c(-0.29, -0.08, -0.05, -0.36), nrow=2, byrow=FALSE)
Omega1_222 <- matrix(c(0.93, -0.15, -0.15, 5.20), nrow=2, byrow=FALSE)

phi20_222 <- c(1.79, 3.00)
A21_222 <- matrix(c(1.20, 0.05, 0.05, 1.30), nrow=2, byrow=FALSE)
A22_222 <- matrix(c(-0.30, -0.10, -0.05, -0.40), nrow=2, byrow=FALSE)
Omega2_222 <- matrix(c(5.88, 3.56, 3.56, 9.80), nrow=2, byrow=FALSE)

alpha1_222 <- 0.37
upsilon1_222 <- c(phi10_222, vec(A11_222), vec(A12_222), vech(Omega1_222))
upsilon2_222 <- c(phi20_222, vec(A21_222), vec(A22_222), vech(Omega2_222))
theta_222 <- c(upsilon1_222, upsilon2_222, alpha1_222)

WL_222 <- diag_Omegas(Omega1_222, Omega2_222)
W_222 <- matrix(WL_222[1:(2^2)], nrow=2, byrow=FALSE)
lambdas_222 <- WL_222[(2^2 + 1):length(WL_222)]
theta_222s <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A21_222),
                vec(A22_222), vec(W_222), lambdas_222, alpha1_222) # SGMVAR

# p=3, M=3, d=2
phi10_332 <- c(1.03, 2.36)
A11_332 <- matrix(c(1.25, 0.06, 0.04, 1.34), nrow=2, byrow=FALSE)
A12_332 <- matrix(c(-0.29, -0.08, -0.05, -0.36), nrow=2, byrow=FALSE)
A13_332 <- matrix(c(1, 0.21, 0.12, 2), nrow=2, byrow=FALSE)
Omega1_332 <- matrix(c(0.93, -0.15, -0.15, 5.20), nrow=2, byrow=FALSE)

phi20_332 <- c(1.79, 3.00)
A21_332 <- matrix(c(1.20, 0.05, 0.05, 1.30), nrow=2, byrow=FALSE)
A22_332 <- matrix(c(-0.30, -0.10, -0.05, -0.40), nrow=2, byrow=FALSE)
A23_332 <- matrix(c(0.30, 0.10, 0.05, 0.40), nrow=2, byrow=FALSE)
Omega2_332 <- matrix(c(5.88, 3.56, 3.56, 9.80), nrow=2, byrow=FALSE)

phi30_332 <- c(1.79, 3.00)
A31_332 <- matrix(c(1.30, 0.03, 0.08, 1.33), nrow=2, byrow=FALSE)
A32_332 <- matrix(c(-0.50, -0.20, -0.01, -0.40), nrow=2, byrow=FALSE)
A33_332 <- matrix(c(0.50, 0.20, 0.01, 0.40), nrow=2, byrow=FALSE)
Omega3_332 <- matrix(c(5.00, 3.00, 3.00, 9.00), nrow=2, byrow=FALSE)

alpha1_332 <- 0.5
alpha2_332 <- 0.3
upsilon1_332 <- c(phi10_332, vec(A11_332), vec(A12_332), vec(A13_332), vech(Omega1_332))
upsilon2_332 <- c(phi20_332, vec(A21_332), vec(A22_332), vec(A23_332), vech(Omega2_332))
upsilon3_332 <- c(phi30_332, vec(A31_332), vec(A32_332), vec(A33_332), vech(Omega3_332))
theta_332 <- c(upsilon1_332, upsilon2_332, upsilon3_332, alpha1_332, alpha2_332)

W_332 <- matrix(c(-0.8924620, 0, 0.3653923, -2.1643472), nrow=2, byrow=FALSE)
lambdas2_332 <- c(7.16, 1.30)
lambdas3_332 <- c(6.10, 1.22)
theta_332s <- c(phi10_332, phi20_332, phi30_332, vec(A11_332), vec(A12_332), vec(A13_332),
                vec(A21_332), vec(A22_332), vec(A23_332), vec(A31_332), vec(A32_332), vec(A33_332),
                vec(W_332), lambdas2_332, lambdas3_332, alpha1_332, alpha2_332) # SGMVAR
Omega1_332s <- tcrossprod(W_332)
Omega2_332s <- W_332%*%tcrossprod(diag(lambdas2_332), W_332)
Omega3_332s <- W_332%*%tcrossprod(diag(lambdas3_332), W_332)
theta_332_froms <-  c(phi10_332, vec(A11_332), vec(A12_332), vec(A13_332), vech(Omega1_332s), # based on SGMVAR
                      phi20_332, vec(A21_332), vec(A22_332), vec(A23_332), vech(Omega2_332s),
                      phi30_332, vec(A31_332), vec(A32_332), vec(A33_332), vech(Omega3_332s),
                      alpha1_332, alpha2_332)

theta_332sWC <- c(phi10_332, phi20_332, phi30_332, vec(A11_332), vec(A12_332), vec(A13_332),
                  vec(A21_332), vec(A22_332), vec(A23_332), vec(A31_332), vec(A32_332), vec(A33_332),
                  Wvec(W_332), lambdas2_332, lambdas3_332, alpha1_332, alpha2_332) # SGMVAR W constrained

# p=1, M=2, d=3
phi10_123 <- c(1.1, 2.2, 3.3)
A11_123 <- matrix(c(0.1, 0.21, 0.31, 0.12, 0.2, 0.32, 0.13, 0.23, 0.1), nrow=3, byrow=FALSE)
Omega1_123 <- matrix(c(1, 0.22, 0.33, 0.22, 2, 0.44, 0.33, 0.44, 3), nrow=3, byrow=FALSE)

phi20_123 <- c(1.11, 2.22, 3.33)
A21_123 <- matrix(c(-0.1, -0.21, -0.31, -0.12, -0.2, -0.32, -0.13, -0.23, -2.1), nrow=3, byrow=FALSE)
Omega2_123 <- matrix(c(1.1, 0.222, 0.333, 0.222, 2.2, 0.444, 0.333, 0.444, 3.3), nrow=3, byrow=FALSE)

alpha1_123 <- 0.6
upsilon1_123 <- c(phi10_123, vec(A11_123), vech(Omega1_123))
upsilon2_123 <- c(phi20_123, vec(A21_123), vech(Omega2_123))
theta_123 <- c(upsilon1_123, upsilon2_123, alpha1_123)

WL_123 <- diag_Omegas(Omega1_123, Omega2_123)
W_123 <- matrix(WL_123[1:(3^2)], nrow=3, byrow=FALSE)
lambdas_123 <- WL_123[(3^2 + 1):length(WL_123)]
theta_123s <- c(phi10_123, phi20_123, vec(A11_123), vec(A21_123), vec(W_123), lambdas_123, alpha1_123) # SGMVAR

# p=2, M=1, d=3
phi10_213 <- c(1.1, 2.2, 3.3)
A11_213 <- matrix(c(1, 0.21, 0.31, 0.12, 1.3, 0.32, 0.13, 0.23, 1), nrow=3, byrow=FALSE)
A12_213 <- matrix(c(-0.1, -0.21, -0.31, -0.12, -0.2, -0.32, -0.13, -0.23, -0.3), nrow=3, byrow=FALSE)
Omega1_213 <- matrix(c(1, 0.22, 0.33, 0.22, 2, 0.44, 0.33, 0.44, 3), nrow=3, byrow=FALSE)

upsilon1_213 <- c(phi10_213, vec(A11_213), vec(A12_213), vech(Omega1_213))
theta_213 <- upsilon1_213

W_213 <- t(chol(Omega1_213))
theta_213s <- c(phi10_213, vec(A11_213), vec(A12_213), vec(W_213)) # SGMVAR

theta_213sWC <- c(phi10_213, vec(A11_213), vec(A12_213), Wvec(W_213)) # SGMVAR W constrained


# Constraining AR-parameters to be the same for all regimes
rbind_diags <- function(p, M, d) {
  I <- diag(p*d^2)
  Reduce(rbind, replicate(M, I, simplify=FALSE))
}


## A(M)(p)_(p)(M)(d)

# p=1, M=1, d=2
C_112 <- rbind_diags(p=1, M=1, d=2)
theta_112c <- c(phi10_112, vec(A11_112), vech(Omega1_112))
C_122 <- rbind_diags(p=1, M=2, d=2)
theta_122c <- c(phi10_122, phi20_122, vec(A11_122), vech(Omega1_122), vech(Omega2_122), alpha1_122)
theta_122c_expanded <- c(phi10_122, vec(A11_122), vech(Omega1_122), phi20_122, vec(A11_122), vech(Omega2_122), alpha1_122)

theta_112csWAR <- c(phi10_112, vec(A11_112), Wvec(W_112)) # SGMVAR W and AR

# p=1, M=2, d=2
C_222 <- rbind_diags(p=2, M=2, d=2)
theta_222c <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vech(Omega1_222), vech(Omega2_222), alpha1_222)
theta_222c_expanded <- c(phi10_222, vec(A11_222), vec(A12_222), vech(Omega1_222), phi20_222, vec(A11_222), vec(A12_222),
                         vech(Omega2_222), alpha1_222)

C_lambda_122 <- matrix(c(1, 1), nrow=2)
theta_122csL <- c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), vec(W_122), 0.5, alpha1_122) # SGMVAR lambdas
theta_122csL_expanded <- c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), vec(W_122), 0.5, 0.5, alpha1_122)
theta_122csLAR <- c(phi10_122, phi20_122, vec(A11_122), vec(W_122), 0.5, alpha1_122) # SGMVAR lambdas and AR
theta_122csLAR_expanded <- c(phi10_122, phi20_122, vec(A11_122), vec(A11_122), vec(W_122), 0.5, 0.5, alpha1_122)

C_lambda_222 <- matrix(c(1, 2), nrow=2)
theta_222csL <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A21_222),
                  vec(A22_222), vec(W_222), 0.2, alpha1_222) # SGMVAR lambdas
theta_222csL_expanded <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A21_222),
                           vec(A22_222), vec(W_222), 0.2, 2*0.2, alpha1_222)
theta_222csLAR <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(W_222), 0.2, alpha1_222) # SGMVAR lambdas and AR
theta_222csLAR_expanded <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222),
                             vec(W_222), 0.2, 2*0.2, alpha1_222) # SGMVAR lambdas and AR


# p=3, M=3, d=2
C_332 <- rbind_diags(p=3, M=3, d=2)
theta_332c <- c(phi10_332, phi20_332, phi30_332, vec(A11_332), vec(A12_332), vec(A13_332), vech(Omega1_332), vech(Omega2_332),
                vech(Omega3_332), alpha1_332, alpha2_332)
theta_332c_expanded <- c(phi10_332, vec(A11_332), vec(A12_332), vec(A13_332), vech(Omega1_332), phi20_332, vec(A11_332),
                         vec(A12_332), vec(A13_332), vech(Omega2_332), phi30_332, vec(A11_332), vec(A12_332), vec(A13_332),
                         vech(Omega3_332), alpha1_332, alpha2_332)

theta_332csWAR <- c(phi10_332, phi20_332, phi30_332, vec(A11_332), vec(A12_332), vec(A13_332), Wvec(W_332), lambdas2_332,
                    lambdas3_332, alpha1_332, alpha2_332) # SGMVAR W and AR
theta_332csWAR_expanded <- c(phi10_332, phi20_332, phi30_332, vec(A11_332), vec(A12_332), vec(A13_332), vec(A11_332),
                             vec(A12_332), vec(A13_332), vec(A11_332), vec(A12_332), vec(A13_332), vec(W_332), lambdas2_332,
                             lambdas3_332, alpha1_332, alpha2_332)

C_lambda_332 <- matrix(c(1, 1, 0, 0, 0, 0, 1, 1), nrow=4, byrow=FALSE)
theta_332csWL <- c(phi10_332, phi20_332, phi30_332, vec(A11_332), vec(A12_332), vec(A13_332), vec(A21_332), vec(A22_332), vec(A23_332),
                   vec(A31_332), vec(A32_332), vec(A33_332), Wvec(W_332), 1, 2, alpha1_332, alpha2_332) # SGMVAR W and L
theta_332csWL_expanded <- c(phi10_332, phi20_332, phi30_332, vec(A11_332), vec(A12_332), vec(A13_332), vec(A21_332), vec(A22_332), vec(A23_332),
                            vec(A31_332), vec(A32_332), vec(A33_332), vec(W_332), 1, 1, 2, 2, alpha1_332, alpha2_332)
theta_332csWLAR <- c(phi10_332, phi20_332, phi30_332, vec(A11_332), vec(A12_332), vec(A13_332), Wvec(W_332), 1, 2,
                     alpha1_332, alpha2_332) # SGMVAR W, L, and AR
theta_332csWLAR_expanded <- c(phi10_332, phi20_332, phi30_332, vec(A11_332), vec(A12_332), vec(A13_332), vec(A11_332), vec(A12_332), vec(A13_332),
                              vec(A11_332), vec(A12_332), vec(A13_332), vec(W_332), 1, 1, 2, 2, alpha1_332, alpha2_332)


# p=1, M=2, d=3
C_123 <- rbind_diags(p=1, M=2, d=3)
theta_123c <- c(phi10_123, phi20_123, vec(A11_123), vech(Omega1_123), vech(Omega2_123), alpha1_123)
theta_123c_expanded <- c(phi10_123, vec(A11_123), vech(Omega1_123), phi20_123, vec(A11_123), vech(Omega2_123), alpha1_123)

C_lambda_123 <- matrix(c(1, 1, 0, 0, 0, 1), nrow=3, byrow=FALSE)
theta_123csL <- c(phi10_123, phi20_123, vec(A11_123), vec(A21_123), vec(W_123), 1, 2, alpha1_123) # SGMVAR lambdas
theta_123csL_expanded <- c(phi10_123, phi20_123, vec(A11_123), vec(A21_123), vec(W_123), 1, 1, 2, alpha1_123)
theta_123csLAR <- c(phi10_123, phi20_123, vec(A11_123), vec(W_123), 1, 2, alpha1_123) # SGMVAR lambdas and AR
theta_123csLAR_expanded <- c(phi10_123, phi20_123, vec(A11_123), vec(A11_123), vec(W_123), 1, 1, 2, alpha1_123)

# p=2, M=1, d=3
C_213 <- rbind_diags(p=2, M=1, d=3)
theta_213c <- c(phi10_213, vec(A11_213), vec(A12_213), vech(Omega1_213))

C_213 <- rbind_diags(p=2, M=1, d=3)
theta_213c <- c(phi10_213, vec(A11_213), vec(A12_213), vech(Omega1_213))
theta_213csWAR <- c(phi10_213, vec(A11_213), vec(A12_213), Wvec(W_213))

# p=2, M=2, d=2, constraint AR-parameters to be the same for all regimes
# and constraint the of-diagonal elements of AR-matrices to be zero.
mat0 <- matrix(c(1, rep(0, 10), 1, rep(0, 8), 1, rep(0, 10), 1), nrow=2*2^2, byrow=FALSE)
C_222_2 <- rbind(mat0, mat0)
A21_222_c2 <- A11_222_c2 <- matrix(c(1.26, 0, 0, 1.34), nrow=2, byrow=FALSE)
A22_222_c2 <- A12_222_c2 <- matrix(c(-0.29, 0, 0, -0.36), nrow=2, byrow=FALSE)
phi10_222_c2 <- c(-0.11, 2.83)
phi20_222_c2 <- c(0.36, 3.19)
Omega1_222_c2 <- matrix(c(0.98, -0.33, -0.33, 5.24), nrow=2, byrow=FALSE)
Omega2_222_c2 <- matrix(c(5.60, 3.46, 3.46, 9.62), nrow=2, byrow=FALSE)
alpha1_222_c2 <- 0.35
theta_222_c2 <- c(phi10_222_c2, phi20_222_c2, 1.26, 1.34, -0.29, -0.36, vech(Omega1_222_c2),
                  vech(Omega2_222_c2), alpha1_222_c2)
theta_222_c2_expanded <- c(phi10_222_c2, vec(A11_222_c2), vec(A12_222_c2), vech(Omega1_222_c2),
                           phi20_222_c2, vec(A21_222_c2), vec(A22_222_c2), vech(Omega2_222_c2),
                           alpha1_222_c2)

WL_222c2 <- diag_Omegas(Omega1_222_c2, Omega2_222_c2)
W_222c2 <- matrix(WL_222c2[1:(2^2)], nrow=2, byrow=FALSE)
lambdas_222c2 <- WL_222c2[(2^2 + 1):length(WL_222c2)]
theta_222_c2s <- c(phi10_222_c2, phi20_222_c2, 1.26, 1.34, -0.29, -0.36, vec(W_222c2), lambdas_222c2, alpha1_222_c2)
theta_222_c2s_expanded <- c(phi10_222_c2, phi20_222_c2, vec(A11_222_c2), vec(A12_222_c2), vec(A11_222_c2), vec(A12_222_c2),
                            vec(W_222c2), lambdas_222c2, alpha1_222_c2)


test_that("is_stationary works correctly", {
  expect_false(is_stationary(p=1, M=1, d=2, params=theta_112))
  expect_true(is_stationary(p=2, M=2, d=2, params=theta_222))
  expect_false(is_stationary(p=3, M=3, d=2, params=theta_332))
  expect_false(is_stationary(p=1, M=2, d=3, params=theta_123))
  expect_false(is_stationary(p=2, M=1, d=3, params=theta_213))

  # SGMVAR
  expect_false(is_stationary(p=1, M=1, d=2, params=theta_112sWC, structural_pars=list(W=W_112)))
  expect_true(is_stationary(p=2, M=2, d=2, params=theta_222s, structural_pars=list(W=W_222)))
  expect_false(is_stationary(p=3, M=3, d=2, params=theta_332sWC, structural_pars=list(W=W_332)))
  expect_false(is_stationary(p=1, M=2, d=3, params=theta_123s, structural_pars=list(W=W_123)))
  expect_false(is_stationary(p=2, M=1, d=3, params=theta_213sWC, structural_pars=list(W=W_213)))
})

test_that("in_paramspace works correctly", {
  expect_false(in_paramspace(p=1, M=1, d=2, params=theta_112))
  expect_true(in_paramspace(p=2, M=2, d=2, params=theta_222))
  expect_false(in_paramspace(p=3, M=3, d=2, params=theta_332))
  expect_false(in_paramspace(p=1, M=2, d=3, params=theta_123))
  expect_false(in_paramspace(p=2, M=1, d=3, params=theta_213))

  expect_false(in_paramspace(p=1, M=1, d=2, params=theta_112c, constraints=C_112))
  expect_true(in_paramspace(p=2, M=2, d=2, params=theta_222c, constraints=C_222))
  expect_true(in_paramspace(p=2, M=2, d=2, params=theta_222_c2, constraints=C_222_2))

  # SGMVAR
  expect_false(in_paramspace(p=1, M=1, d=2, params=theta_112sWC, structural_pars=list(W=W_112)))
  expect_true(in_paramspace(p=2, M=2, d=2, params=theta_222s, structural_pars=list(W=W_222)))
  expect_false(in_paramspace(p=3, M=3, d=2, params=theta_332sWC, structural_pars=list(W=W_332)))
  expect_false(in_paramspace(p=1, M=2, d=3, params=theta_123s, structural_pars=list(W=W_123)))
  expect_false(in_paramspace(p=2, M=1, d=3, params=theta_213sWC, structural_pars=list(W=W_213)))

  expect_true(in_paramspace(p=2, M=2, d=2, params=theta_222csLAR, constraints=C_222, structural_pars=list(W=W_222, C_lambda=C_lambda_222)))
  expect_true(in_paramspace(p=2, M=2, d=2, params=theta_222_c2s, constraints=C_222_2, structural_pars=list(W=W_222c2)))
  theta_222_c3s <- theta_222_c2s
  theta_222_c3s[9] <- 1
  expect_false(in_paramspace(p=2, M=2, d=2, params=theta_222_c3s, constraints=C_222_2, structural_pars=list(W=W_222c2)))
  theta_222csLAR_2 <- theta_222csLAR
  theta_222csLAR_2[15] <- -1
  expect_false(in_paramspace(p=2, M=2, d=2, params=theta_222csLAR_2, constraints=C_222, structural_pars=list(W=W_222, C_lambda=C_lambda_222)))
  theta_222csLAR_3 <- theta_222csLAR
  theta_222csLAR_3[17] <- -0.1
  expect_false(in_paramspace(p=2, M=2, d=2, params=theta_222csLAR_3, constraints=C_222, structural_pars=list(W=W_222, C_lambda=C_lambda_222)))
  theta_112s_2 <- c(phi10_112, vec(A12_222), Wvec(W_112))
  expect_true(in_paramspace(p=1, M=1, d=2, params=theta_112s_2, structural_pars=list(W=W_112)))
})

test_that("check_parameters works correctly", {
  expect_error(check_parameters(p=1, M=1, d=2, params=theta_112))
  expect_error(check_parameters(p=2, M=3, d=2, params=theta_222))
  expect_error(check_parameters(p=3, M=3, d=2, params=theta_332))
  expect_error(check_parameters(p=1, M=2, d=3, params=theta_123))
  expect_error(check_parameters(p=2, M=1, d=3, params=theta_213))

  expect_error(check_parameters(p=1, M=1, d=2, params=theta_112c, constraints=C_112))
  expect_error(check_parameters(p=2, M=2, d=2, params=theta_222c, constraints=t(C_222)))
  expect_error(check_parameters(p=2, M=2, d=2, params=theta_222_c2, constraints=as.vector(C_222_2)))

  # SGMVAR
  expect_error(check_parameters(p=1, M=1, d=2, params=theta_112s, structural_pars=list(W=W_112)))
  expect_error(check_parameters(p=2, M=2, d=2, params=theta_222s, structural_pars=list(xW=W_222)))
  expect_error(check_parameters(p=2, M=2, d=2, params=theta_222s, structural_pars=list(W_222)))
  expect_error(check_parameters(p=2, M=2, d=2, params=theta_222s, structural_pars=list(W=W_123)))
  expect_error(check_parameters(p=2, M=2, d=2, params=c(theta_222s, 1), structural_pars=list(W=W_222)))
  expect_error(check_parameters(p=2, M=2, d=2, params=theta_222s, structural_pars=list(W=W_222, C_lambda=matrix(1))))

  expect_error(check_parameters(p=1, M=1, d=2, params=theta_112cs, constraints=C_112, structural_pars=list(W=W_112)))
  expect_error(check_parameters(p=2, M=2, d=2, params=theta_222c, constraints=C_222, structural_pars=list(W=W_123)))
  theta_222_c3s <- theta_222_c2s
  theta_222_c3s[9] <- 1
  expect_error(check_parameters(p=2, M=2, d=2, params=theta_222_c3s, constraints=C_222_2, structural_pars=list(W=W_222c2)))
  theta_222csLAR_2 <- theta_222csLAR
  theta_222csLAR_2[15] <- -1
  expect_error(check_parameters(p=2, M=2, d=2, params=theta_222csLAR_2, constraints=C_222, structural_pars=list(W=W_222, C_lambda=C_lambda_222)))
  theta_222csLAR_3 <- theta_222csLAR
  theta_222csLAR_3[17] <- -0.1
  expect_error(check_parameters(p=2, M=2, d=2, params=theta_222csLAR_3, constraints=C_222, structural_pars=list(W=W_222, C_lambda=C_lambda_222)))
})


test_that("check_constraints works correctly", {
  expect_error(check_constraints(p=1, M=1, d=2, constraints=cbind(C_112, C_112)))
  expect_error(check_constraints(p=1, M=1, d=2, constraints=rbind(C_112, C_112)))
  expect_error(check_constraints(p=1, M=2, d=2, constraints=cbind(C_122[,1:3], rep(0, 8))))
  expect_error(check_constraints(p=1, M=2, d=2, constraints=as.data.frame(C_122)))
  expect_error(check_constraints(p=2, M=2, d=2, constraints=C_222[-1,]))
  expect_error(check_constraints(p=2, M=2, d=2, constraints=C_222_2[-16,]))
  expect_error(check_constraints(p=3, M=3, d=2, constraints=t(C_332)))
  expect_error(check_constraints(p=1, M=2, d=3, constraints=cbind(C_123[,1:8], C_123[,8])))
  expect_error(check_constraints(p=2, M=1, d=3, constraints=as.vector(C_213)))

  # Structural parameter constraints
  expect_error(check_constraints(p=1, M=1, d=2, structural_pars=list(W=W_123)))
  expect_error(check_constraints(p=1, M=1, d=2, structural_pars=list(W=vec(W_112))))
  expect_error(check_constraints(p=1, M=1, d=2, structural_pars=list(W=W_112, C_lambda=numeric(0))))
  expect_error(check_constraints(p=3, M=3, d=2, structural_pars=list(A=1)))
  expect_error(check_constraints(p=3, M=3, d=2, structural_pars=list(W=W_332, C_lambda=matrix(1:3, nrow=3))))
  expect_error(check_constraints(p=1, M=2, d=2, structural_pars=list(W=W_122, C_lambda=matrix(1:6, nrow=2, ncol=3))))
  W_bad <- matrix(c(0, 0, 1, 2), nrow=2)
  expect_error(check_constraints(p=1, M=1, d=2, structural_pars=list(W=W_bad)))
  expect_error(check_constraints(p=1, M=1, d=2, structural_pars=list(W=t(W_bad))))
  C_lambda_bad <- matrix(c(1, -0.001), nrow=2)
  expect_error(check_constraints(p=1, M=2, d=2, structural_pars=list(W=W_122, C_lambda=C_lambda_bad)))
})



data_na <- eurusd; data_na[252, 2] <- NA
eurusd_mat <- t(t(eurusd))

test_that("check_data works correctly", {
  expect_equal(check_data(data=eurusd, p=10), eurusd)
  expect_equal(check_data(data=ts(eurusd_mat, start=c(1989, 1), deltat=1/12), p=1), eurusd)
  expect_equal(check_data(data=as.data.frame(eurusd_mat), p=2), eurusd_mat)

  expect_error(check_data(data=eurusd[,1, drop=FALSE], p=4))
  expect_error(check_data(data=as.data.frame(eurusd[,2]), p=4))
  expect_error(check_data(data=ts(eurusd[,2], start=c(1989, 1), deltat=1/12), p=3))
  expect_error(check_data(data=data_na, p=4))
  expect_error(check_data(data=eurusd, p=nrow(eurusd)))
})

test_that("check_pMd works correctly", {
  expect_error(check_pMd(p=1, M=1, d=1))
  expect_error(check_pMd(p=1, M=1.2, d=2))
  expect_error(check_pMd(p=0, M=1, d=2))
  expect_error(check_pMd(p=1, M=-1, d=2))
  expect_error(check_pMd(p=1.1, M=1, d=2))
  expect_error(check_pMd(p=1, M=1, d=2.2))
})

test_that("all_pos_ints works correctly", {
  expect_true(all_pos_ints(c(1, 2, 3)))
  expect_true(all_pos_ints(1))
  expect_true(all_pos_ints(list(1, 3, 100)))

  expect_false(all_pos_ints(c(1, 2, 0)))
  expect_false(all_pos_ints(-1))
  expect_false(all_pos_ints(0.1))
  expect_false(all_pos_ints(1.1))
  expect_false(all_pos_ints(list(1, 2, 3, 0.1)))
})

# p=2, M=1, d=2
phi10_212 <- c(1.03, 2.36)
A11_212 <- matrix(c(1.25, 0.06, 0.04, 1.34), nrow=2, byrow=FALSE)
A12_212 <- matrix(c(-0.29, -0.08, -0.05, -0.36), nrow=2, byrow=FALSE)
Omega1_212 <- matrix(c(0.93, -0.15, -0.15, 5.20), nrow=2, byrow=FALSE)

theta_212 <- c(phi10_212, vec(A11_212), vec(A12_212), vech(Omega1_212))

W_212 <- t(chol(Omega1_212))
theta_212sWC <- c(phi10_212, vec(A11_212), vec(A12_212), Wvec(W_212)) # SGMVAR, W constrained by Cholesky

C_212 <- rbind_diags(p=2, M=1, d=2)
theta_212c <- c(phi10_212, vec(A11_212), vec(A12_212), vech(Omega1_212))

theta_212csWAR <- c(phi10_212, vec(A11_212), vec(A11_212), Wvec(W_212)) # SGMVAR W and AR

test_that("n_params works correctly", {
  expect_equal(n_params(p=1, M=1, d=2), length(theta_112))
  expect_equal(n_params(p=1, M=2, d=2), length(theta_122))
  expect_equal(n_params(p=2, M=2, d=2), length(theta_222))
  expect_equal(n_params(p=3, M=3, d=2), length(theta_332))
  expect_equal(n_params(p=1, M=2, d=3), length(theta_123))
  expect_equal(n_params(p=2, M=1, d=3), length(theta_213))
  expect_equal(n_params(p=2, M=1, d=2), length(theta_212))

  expect_equal(n_params(p=1, M=1, d=2, constraints=C_112), length(theta_112c))
  expect_equal(n_params(p=1, M=2, d=2, constraints=C_122), length(theta_122c))
  expect_equal(n_params(p=2, M=2, d=2, constraints=C_222), length(theta_222c))
  expect_equal(n_params(p=3, M=3, d=2, constraints=C_332), length(theta_332c))
  expect_equal(n_params(p=1, M=2, d=3, constraints=C_123), length(theta_123c))
  expect_equal(n_params(p=2, M=1, d=3, constraints=C_213), length(theta_213c))
  expect_equal(n_params(p=2, M=1, d=2, constraints=C_212), length(theta_212c))

  # SGMVAR
  expect_equal(n_params(p=1, M=1, d=2, structural_pars=list(W=W_112)), length(theta_112sWC))
  expect_equal(n_params(p=1, M=2, d=2, structural_pars=list(W=W_122)), length(theta_122s))
  expect_equal(n_params(p=2, M=2, d=2, structural_pars=list(W=W_222, C_lambda=C_lambda_222)), length(theta_222csL))
  expect_equal(n_params(p=3, M=3, d=2, structural_pars=list(W=W_332)), length(theta_332sWC))
  expect_equal(n_params(p=1, M=2, d=3, structural_pars=list(W=W_123)), length(theta_123s))
  expect_equal(n_params(p=1, M=2, d=3, structural_pars=list(W=W_123, C_lambda=C_lambda_123)), length(theta_123csL))
  expect_equal(n_params(p=2, M=1, d=3, structural_pars=list(W=W_213)), length(theta_213sWC))
  expect_equal(n_params(p=2, M=1, d=2, structural_pars=list(W=W_212)), length(theta_212sWC))

  expect_equal(n_params(p=1, M=1, d=2, constraints=C_112, structural_pars=list(W=W_112)), length(theta_112csWAR))
  expect_equal(n_params(p=1, M=2, d=2, constraints=C_122, structural_pars=list(W=W_122, C_lambda=C_lambda_122)), length(theta_122csLAR))
  expect_equal(n_params(p=2, M=2, d=2, constraints=C_222, structural_pars=list(W=W_222, C_lambda=C_lambda_222)), length(theta_222csLAR))
  expect_equal(n_params(p=3, M=3, d=2, constraints=C_332, structural_pars=list(W=W_332, C_lambda=C_lambda_332)), length(theta_332csWLAR))
  expect_equal(n_params(p=1, M=2, d=3, constraints=C_123, structural_pars=list(W=W_123, C_lambda=C_lambda_123)), length(theta_123csLAR))
  expect_equal(n_params(p=2, M=1, d=3, constraints=C_213, structural_pars=list(W=W_213)), length(theta_213csWAR))
  expect_equal(n_params(p=2, M=1, d=2, constraints=C_212, structural_pars=list(W=W_212)), length(theta_212csWAR))

})


