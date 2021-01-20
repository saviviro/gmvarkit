context("Parameter reforms")
library(gmvarkit)

data0 <- t(t(eurusd))
data2 <- unname(data0)
data3 <- cbind(1:252, data2)
n_obs <- nrow(data3) # 252
dat2_1 <- reform_data(data2, p=1)
dat2_2 <- reform_data(data2, p=2)
dat2_3 <- reform_data(data2, p=3)
dat3_1 <- reform_data(data3, p=1)
dat3_2 <- reform_data(data3, p=2)

test_that("reform_data works correctly", {
  expect_equal(dat2_1, data2)

  expect_equal(nrow(dat2_2), n_obs - 2 + 1)
  expect_equal(dat2_2[1,], c(data2[2,], data2[1,]))
  expect_equal(dat2_2[23,], c(data2[24,], data2[23,]))
  expect_equal(dat2_2[nrow(dat2_2),], c(data2[252,], data2[251,]))

  expect_equal(nrow(dat2_3), n_obs - 3 + 1)
  expect_equal(dat2_3[1,], c(data2[3,], data2[2,], data2[1,]))
  expect_equal(dat2_3[100,], c(data2[102,], data2[101,], data2[100,]))
  expect_equal(dat2_3[nrow(dat2_3),], c(data2[252,], data2[251,], data2[250,]))

  expect_equal(dat3_1, data3)

  expect_equal(nrow(dat3_2), n_obs - 2 + 1)
  expect_equal(dat3_2[1,], c(data3[2,], data3[1,]))
  expect_equal(dat3_2[13,], c(data3[14,], data3[13,]))
  expect_equal(dat3_2[nrow(dat3_2),], c(data3[252,], data3[251,]))
})


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
A11_122 <- matrix(c(1, -0.06, -0.04, 1), nrow=2, byrow=FALSE)
Omega1_122 <- matrix(c(0.93, -0.15, -0.15, 5.20), nrow=2, byrow=FALSE)
phi20_122 <- c(1.79, 3.00)
A21_122 <- A11_122 + 0.02
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
Omega1_222s <- tcrossprod(W_222)
Omega2_222s <- W_222%*%tcrossprod(diag(lambdas_222), W_222)


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
A11_123 <- matrix(c(1, 0.21, 0.31, 0.12, 2, 0.32, 0.13, 0.23, 3), nrow=3, byrow=FALSE)
Omega1_123 <- matrix(c(1, 0.22, 0.33, 0.22, 2, 0.44, 0.33, 0.44, 3), nrow=3, byrow=FALSE)
phi20_123 <- c(1.11, 2.22, 3.33)
A21_123 <- matrix(c(-1, -0.21, -0.31, -0.12, -2, -0.32, -0.13, -0.23, -3), nrow=3, byrow=FALSE)
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
A11_213 <- matrix(c(1, 0.21, 0.31, 0.12, 2, 0.32, 0.13, 0.23, 3), nrow=3, byrow=FALSE)
A12_213 <- matrix(c(-1, -0.21, -0.31, -0.12, -2, -0.32, -0.13, -0.23, -3), nrow=3, byrow=FALSE)
Omega1_213 <- matrix(c(1, 0.22, 0.33, 0.22, 2, 0.44, 0.33, 0.44, 3), nrow=3, byrow=FALSE)
upsilon1_213 <- c(phi10_213, vec(A11_213), vec(A12_213), vech(Omega1_213))
theta_213 <- upsilon1_213

W_213 <- t(chol(Omega1_213))
theta_213s <- c(phi10_213, vec(A11_213), vec(A12_213), vec(W_213)) # SGMVAR

theta_213sWC <- c(phi10_213, vec(A11_213), vec(A12_213), Wvec(W_213)) # SGMVAR W constrained


## A(M)(p)_(p)(M)(d)
rbind_diags <- function(p, M, d) {
  I <- diag(p*d^2)
  Reduce(rbind, replicate(M, I, simplify=FALSE))
}

##  Constraining AR-parameters to be the same for all regimes

# p=1, M=1, d=2
C_112 <- rbind_diags(p=1, M=1, d=2)
theta_112c <- c(phi10_112, vec(A11_112), vech(Omega1_112))

theta_112csWAR <- c(phi10_112, vec(A11_112), Wvec(W_112)) # SGMVAR W and AR

# p=1, M=2, d=2
C_122 <- rbind_diags(p=1, M=2, d=2)
theta_122c <- c(phi10_122, phi20_122, vec(A11_122), vech(Omega1_122), vech(Omega2_122), alpha1_122)
theta_122c_expanded <- c(phi10_122, vec(A11_122), vech(Omega1_122), phi20_122, vec(A11_122), vech(Omega2_122), alpha1_122)

C_lambda_122 <- matrix(c(1, 1), nrow=2)
theta_122csL <- c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), vec(W_122), 0.5, alpha1_122) # SGMVAR lambdas
theta_122csL_expanded <- c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), vec(W_122), 0.5, 0.5, alpha1_122)
theta_122csLAR <- c(phi10_122, phi20_122, vec(A11_122), vec(W_122), 0.5, alpha1_122) # SGMVAR lambdas and AR
theta_122csLAR_expanded <- c(phi10_122, phi20_122, vec(A11_122), vec(A11_122), vec(W_122), 0.5, 0.5, alpha1_122)

# p=2, M=2, d=2
C_222 <- rbind_diags(p=2, M=2, d=2)
theta_222c <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vech(Omega1_222), vech(Omega2_222), alpha1_222)
theta_222c_expanded <- c(phi10_222, vec(A11_222), vec(A12_222), vech(Omega1_222), phi20_222, vec(A11_222), vec(A12_222),
                         vech(Omega2_222), alpha1_222)

C_lambda_222 <- matrix(c(1, 2), nrow=2)
theta_222csL <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A21_222),
                  vec(A22_222), vec(W_222), 0.2, alpha1_222) # SGMVAR lambdas
theta_222csL_expanded <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A21_222),
                           vec(A22_222), vec(W_222), 0.2, 2*0.2, alpha1_222)
theta_222csLAR <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(W_222), 0.2, alpha1_222) # SGMVAR lambdas and AR
theta_222csLAR_expanded <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222),
                             vec(W_222), 0.2, 2*0.2, alpha1_222)


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
theta_213csWAR <- c(phi10_213, vec(A11_213), vec(A12_213), Wvec(W_213))

# p=2, M=2, d=2, constrain AR-parameters to be the same for all regimes
# and constrain the of-diagonal elements of AR-matrices to be zero.
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
theta_222_c2s <- c(phi10_222_c2, phi20_222_c2, 1.26, 1.34, -0.29, -0.36, vec(W_222c2), lambdas_222c2, alpha1_222_c2) # SGMVAR AR
theta_222_c2s_expanded <- c(phi10_222_c2, phi20_222_c2, vec(A11_222_c2), vec(A12_222_c2), vec(A11_222_c2), vec(A12_222_c2),
                            vec(W_222c2), lambdas_222c2, alpha1_222_c2)

## Models with same_means

# p=1, M=1, d=2, same_means=list(1)
theta_112_int <- c(phi10_112, vec(A11_112), vech(Omega1_112))
theta_112_int_expanded <- theta_112_int

# p=1, M=2, d=2, same_means=list(1:2)
theta_122_int <- c(phi10_122, vec(A11_122), vec(A21_122), vech(Omega1_122), vech(Omega2_122), alpha1_122)
theta_122_int_expanded <- c(phi10_122, vec(A11_122), vech(Omega1_122), phi10_122, vec(A21_122), vech(Omega2_122), alpha1_122)

# p=1, M=2, d=2, same_means=list(1, 2)
theta_122_int2 <- c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), vech(Omega1_122), vech(Omega2_122), alpha1_122)
theta_122_int2_expanded <- c(phi10_122, vec(A11_122), vech(Omega1_122), phi20_122, vec(A21_122), vech(Omega2_122), alpha1_122)

# p=1, M=2, d=2, constraints=C_122, structural_pars=list(W=W_122, C_lambda=C_lambda_122), same_means=list(1:2)
theta_122csL_int <-  c(phi10_122, vec(A11_122), vec(W_122), 0.5, alpha1_122)
theta_122csL_int_expanded <- c(phi10_122, phi10_122, vec(A11_122), vec(A11_122), vec(W_122), 0.5, 0.5, alpha1_122)

# p=2, M=2, d=2, constraints=C_222, same_means=list(1:2)
theta_222c_int <- c(phi10_222, vec(A11_222), vec(A12_222), vech(Omega1_222), vech(Omega2_222), alpha1_222)
theta_222c_int_expanded <- c(phi10_222, vec(A11_222), vec(A12_222), vech(Omega1_222), phi10_222, vec(A11_222), vec(A12_222),
                             vech(Omega2_222), alpha1_222)

# p=2, M=2, d=2, constraints=C_222, structural_pars=list(W=W_222, C_lambda=C_lambda_222), same_means=list(1:2)
theta_222csLAR_int <- c(phi10_222, vec(A11_222), vec(A12_222), vec(W_222), 0.2, alpha1_222)
theta_222csLAR_int_expanded <-  c(phi10_222, phi10_222, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222),
                                  vec(W_222), 0.2, 2*0.2, alpha1_222)

# p=3, M=3, d=2, constraints=C_332, same_means=list(1, 2:3)
theta_332c_int <- c(phi10_332, phi20_332, vec(A11_332), vec(A12_332), vec(A13_332), vech(Omega1_332), vech(Omega2_332),
                    vech(Omega3_332), alpha1_332, alpha2_332)
theta_332c_int_expanded <- c(phi10_332, vec(A11_332), vec(A12_332), vec(A13_332), vech(Omega1_332), phi20_332, vec(A11_332),
                             vec(A12_332), vec(A13_332), vech(Omega2_332), phi20_332, vec(A11_332), vec(A12_332), vec(A13_332),
                             vech(Omega3_332), alpha1_332, alpha2_332)

# p=3, M=3, d=2, constraints=C_332, same_means=list(2, c(1, 3))
theta_332c_int2 <- c(phi10_332, phi20_332, vec(A11_332), vec(A12_332), vec(A13_332), vech(Omega1_332), vech(Omega2_332),
                    vech(Omega3_332), alpha1_332, alpha2_332)
theta_332c_int2_expanded <- c(phi20_332, vec(A11_332), vec(A12_332), vec(A13_332), vech(Omega1_332), phi10_332, vec(A11_332),
                              vec(A12_332), vec(A13_332), vech(Omega2_332), phi20_332, vec(A11_332), vec(A12_332), vec(A13_332),
                              vech(Omega3_332), alpha1_332, alpha2_332)

# p=3, M=3, d=2, constraints=C_332, same_means=list(1:3)
theta_332c_int3 <- c(phi10_332, vec(A11_332), vec(A12_332), vec(A13_332), vech(Omega1_332), vech(Omega2_332),
                     vech(Omega3_332), alpha1_332, alpha2_332)
theta_332c_int3_expanded <- c(phi10_332, vec(A11_332), vec(A12_332), vec(A13_332), vech(Omega1_332), phi10_332, vec(A11_332),
                              vec(A12_332), vec(A13_332), vech(Omega2_332), phi10_332, vec(A11_332), vec(A12_332), vec(A13_332),
                              vech(Omega3_332), alpha1_332, alpha2_332)

# p=3, M=3, d=2, constraints=C_332, structural_pars=list(W=W_332), same_means=list(c(1, 3), 2)
theta_332csWAR_int <- c(phi10_332, phi20_332, vec(A11_332), vec(A12_332), vec(A13_332), Wvec(W_332), lambdas2_332,
                        lambdas3_332, alpha1_332, alpha2_332)
theta_332csWAR_int_expanded <- c(phi10_332, phi20_332, phi10_332, vec(A11_332), vec(A12_332), vec(A13_332), vec(A11_332),
                                 vec(A12_332), vec(A13_332), vec(A11_332), vec(A12_332), vec(A13_332), vec(W_332), lambdas2_332,
                                 lambdas3_332, alpha1_332, alpha2_332)

# p=3, M=3, d=2, structural_pars=list(W=W_332, C_lambda=C_lambda_332), same_means=list(1:2, 3)
theta_332csWL_int <- c(phi10_332, phi30_332, vec(A11_332), vec(A12_332), vec(A13_332), vec(A21_332), vec(A22_332), vec(A23_332),
                       vec(A31_332), vec(A32_332), vec(A33_332), Wvec(W_332), 1, 2, alpha1_332, alpha2_332)
theta_332csWL_int_expanded <- c(phi10_332, phi10_332, phi30_332, vec(A11_332), vec(A12_332), vec(A13_332), vec(A21_332), vec(A22_332), vec(A23_332),
                               vec(A31_332), vec(A32_332), vec(A33_332), vec(W_332), 1, 1, 2, 2, alpha1_332, alpha2_332)

# p=3, M=3, d=2, constraints=C_332, structural_pars=list(W=W_332, C_lambda=C_lambda_332), same_means=list(2:3, 1)
theta_332csWLAR_int <- c(phi10_332, phi20_332, vec(A11_332), vec(A12_332), vec(A13_332), Wvec(W_332), 1, 2,
                         alpha1_332, alpha2_332)
theta_332csWLAR_int_expanded <- c(phi20_332, phi10_332, phi10_332, vec(A11_332), vec(A12_332), vec(A13_332), vec(A11_332), vec(A12_332), vec(A13_332),
                                  vec(A11_332), vec(A12_332), vec(A13_332), vec(W_332), 1, 1, 2, 2, alpha1_332, alpha2_332)

# p=1, M=2, d=3, same_means=list(1:2)
theta_123_int <- c(phi10_123, vec(A11_123), vec(A21_123), vech(Omega1_123), vech(Omega2_123), alpha1_123)
theta_123_int_expanded <- c(phi10_123, vec(A11_123), vech(Omega1_123), phi10_123, vec(A21_123), vech(Omega2_123), alpha1_123)

# p=1, M=2, d=3, constraints=C_123, same_means=list(1:2)
theta_123c_int <- c(phi10_123, vec(A11_123), vech(Omega1_123), vech(Omega2_123), alpha1_123)
theta_123c_int_expanded <- c(phi10_123, vec(A11_123), vech(Omega1_123), phi10_123, vec(A11_123), vech(Omega2_123), alpha1_123)

# p=1, M=2, d=3, structural_pars=list(W=W_123, C_lambda=C_lambda_123) same_means=list(1:2)
theta_123csL_int <- c(phi10_123, vec(A11_123), vec(A21_123), vec(W_123), 1, 2, alpha1_123)
theta_123csL_int_expanded <- c(phi10_123, phi10_123, vec(A11_123), vec(A21_123), vec(W_123), 1, 1, 2, alpha1_123)

# p=1, M=2, d=3, constraints=C_123, structural_pars=list(W=W_123, C_lambda=C_lambda_123) same_means=list(1:2)
theta_123csLAR_int <- c(phi10_123, vec(A11_123), vec(W_123), 1, 2, alpha1_123) # SGMVAR lambdas and AR
theta_123csLAR_int_expanded <- c(phi10_123, phi10_123, vec(A11_123), vec(A11_123), vec(W_123), 1, 1, 2, alpha1_123)

# p=2, M=1, d=3, constraints=C_213, same_means=list(1)
theta_213c_int <- c(phi10_213, vec(A11_213), vec(A12_213), vech(Omega1_213))
theta_213c_int_expanded <- theta_213c_int

# p=2, M=1, d=3, constraints=C_213, structural_pars=list(W=W_213), same_means=list(1)
theta_213csWAR_int <- c(phi10_213, vec(A11_213), vec(A12_213), Wvec(W_213))
theta_213csWAR_int_expanded <- c(phi10_213, vec(A11_213), vec(A12_213), vec(W_213))

test_that("reform_constrained_pars works correctly", {
  expect_equal(reform_constrained_pars(p=1, M=1, d=2, params=theta_112c, constraints=C_112), theta_112)
  expect_equal(reform_constrained_pars(p=1, M=2, d=2, params=theta_122c, constraints=C_122), theta_122c_expanded)
  expect_equal(reform_constrained_pars(p=2, M=2, d=2, params=theta_222, constraints=NULL), theta_222)
  expect_equal(reform_constrained_pars(p=2, M=2, d=2, params=theta_222c, constraints=C_222), theta_222c_expanded)
  expect_equal(reform_constrained_pars(p=3, M=3, d=2, params=theta_332c, constraints=C_332), theta_332c_expanded)
  expect_equal(reform_constrained_pars(p=1, M=2, d=3, params=theta_123c, constraints=C_123), theta_123c_expanded)
  expect_equal(reform_constrained_pars(p=2, M=1, d=3, params=theta_213c, constraints=C_213), theta_213)
  expect_equal(reform_constrained_pars(p=2, M=2, d=2, params=theta_222_c2, constraints=C_222_2), theta_222_c2_expanded)

  # SGMVAR only W constrained
  expect_equal(reform_constrained_pars(p=1, M=1, d=2, params=theta_112sWC, structural_pars=list(W=W_112)), theta_112s)
  expect_equal(reform_constrained_pars(p=3, M=3, d=2, params=theta_332sWC, structural_pars=list(W=W_332)), theta_332s)
  expect_equal(reform_constrained_pars(p=2, M=1, d=3, params=theta_213sWC, structural_pars=list(W=W_213)), theta_213s)

  # SGMVAR only lambdas constrained
  expect_equal(reform_constrained_pars(p=1, M=2, d=2, params=theta_122csL, structural_pars=list(W=W_122, C_lambda=C_lambda_122)), theta_122csL_expanded)
  expect_equal(reform_constrained_pars(p=2, M=2, d=2, params=theta_222csL, structural_pars=list(W=W_222, C_lambda=C_lambda_222)), theta_222csL_expanded)
  expect_equal(reform_constrained_pars(p=1, M=2, d=3, params=theta_123csL, structural_pars=list(W=W_123, C_lambda=C_lambda_123)), theta_123csL_expanded)

  # SGMVAR AR parameters only constrained
  expect_equal(reform_constrained_pars(p=2, M=2, d=2, params=theta_222_c2s, constraints=C_222_2, structural_pars=list(W=W_222c2)),
               theta_222_c2s_expanded)

  # SGMVAR W and AR parameters constrained
  expect_equal(reform_constrained_pars(p=1, M=1, d=2, params=theta_112csWAR, constraints=C_112, structural_pars=list(W=W_112)), theta_112s)
  expect_equal(reform_constrained_pars(p=3, M=3, d=2, params=theta_332csWAR, constraints=C_332, structural_pars=list(W=W_332)), theta_332csWAR_expanded)
  expect_equal(reform_constrained_pars(p=2, M=1, d=3, params=theta_213csWAR, constraints=C_213, structural_pars=list(W=W_213)), theta_213s)

  # SGMVAR lambdas and AR parameters
  expect_equal(reform_constrained_pars(p=1, M=2, d=2, params=theta_122csLAR, constraints=C_122, structural_pars=list(W=W_122, C_lambda=C_lambda_122)),
               theta_122csLAR_expanded)
  expect_equal(reform_constrained_pars(p=2, M=2, d=2, params=theta_222csLAR, constraints=C_222, structural_pars=list(W=W_222, C_lambda=C_lambda_222)),
               theta_222csLAR_expanded)
  expect_equal(reform_constrained_pars(p=1, M=2, d=3, params=theta_123csLAR, constraints=C_123, structural_pars=list(W=W_123, C_lambda=C_lambda_123)),
               theta_123csLAR_expanded)

  # SGMVAR W and Lambda parameters constrained
  expect_equal(reform_constrained_pars(p=3, M=3, d=2, params=theta_332csWL, structural_pars=list(W=W_332, C_lambda=C_lambda_332)),
               theta_332csWL_expanded)

  # SGMVAR W, Lambda, and AR parameters constrained
  expect_equal(reform_constrained_pars(p=3, M=3, d=2, params=theta_332csWLAR, constraints=C_332, structural_pars=list(W=W_332, C_lambda=C_lambda_332)),
               theta_332csWLAR_expanded)

  # Models with same_means
  expect_equal(reform_constrained_pars(p=1, M=1, d=2, params=theta_112_int, same_means=list(1)), theta_112_int_expanded)
  expect_equal(reform_constrained_pars(p=1, M=2, d=2, params=theta_122_int, same_means=list(1:2)), theta_122_int_expanded)
  expect_equal(reform_constrained_pars(p=1, M=2, d=2, params=theta_122_int2, same_means=list(1, 2)), theta_122_int2_expanded)
  expect_equal(reform_constrained_pars(p=1, M=2, d=2, params=theta_122csL_int, constraints=C_122, structural_pars=list(W=W_122, C_lambda=C_lambda_122),
                                       same_means=list(1:2)), theta_122csL_int_expanded)
  expect_equal(reform_constrained_pars(p=2, M=2, d=2, params=theta_222c_int, constraints=C_222, same_means=list(1:2)),
               theta_222c_int_expanded)
  expect_equal(reform_constrained_pars(p=2, M=2, d=2, params=theta_222csLAR_int, constraints=C_222, structural_pars=list(W=W_222, C_lambda=C_lambda_222),
                                       same_means=list(1:2)), theta_222csLAR_int_expanded)
  expect_equal(reform_constrained_pars(p=3, M=3, d=2, params=theta_332c_int, constraints=C_332, same_means=list(1, 2:3)),
               theta_332c_int_expanded)
  expect_equal(reform_constrained_pars(p=3, M=3, d=2, params=theta_332c_int2, constraints=C_332, same_means=list(2, c(1, 3))),
               theta_332c_int2_expanded)
  expect_equal(reform_constrained_pars(p=3, M=3, d=2, params=theta_332c_int3, constraints=C_332, same_means=list(1:3)),
               theta_332c_int3_expanded)
  expect_equal(reform_constrained_pars(p=3, M=3, d=2, params=theta_332csWAR_int, constraints=C_332, structural_pars=list(W=W_332), same_means=list(c(1, 3), 2)),
               theta_332csWAR_int_expanded)
  expect_equal(reform_constrained_pars(p=3, M=3, d=2, params=theta_332csWL_int, structural_pars=list(W=W_332, C_lambda=C_lambda_332), same_means=list(1:2, 3)),
               theta_332csWL_int_expanded)
  expect_equal(reform_constrained_pars(p=3, M=3, d=2, params=theta_332csWLAR_int, constraints=C_332, structural_pars=list(W=W_332, C_lambda=C_lambda_332), same_means=list(2:3, 1)),
               theta_332csWLAR_int_expanded)
  expect_equal(reform_constrained_pars(p=1, M=2, d=3, params=theta_123_int, same_means=list(1:2)),
               theta_123_int_expanded)
  expect_equal(reform_constrained_pars(p=1, M=2, d=3, params=theta_123c_int, constraints=C_123, same_means=list(1:2)),
               theta_123c_int_expanded)
  expect_equal(reform_constrained_pars(p=1, M=2, d=3, params=theta_123csL_int, structural_pars=list(W=W_123, C_lambda=C_lambda_123), same_means=list(1:2)),
               theta_123csL_int_expanded)
  expect_equal(reform_constrained_pars(p=1, M=2, d=3, params=theta_123csLAR_int, constraints=C_123, structural_pars=list(W=W_123, C_lambda=C_lambda_123), same_means=list(1:2)),
               theta_123csLAR_int_expanded)
  expect_equal(reform_constrained_pars(p=2, M=1, d=3, params=theta_213c_int, constraints=C_213, same_means=list(1)),
               theta_213c_int_expanded)
  expect_equal(reform_constrained_pars(p=2, M=1, d=3, params=theta_213csWAR_int, constraints=C_213, structural_pars=list(W=W_213), same_means=list(1)),
               theta_213csWAR_int_expanded)
})

test_that("reform_structural_pars works correctly", {
  expect_equal(reform_structural_pars(p=1, M=1, d=2, params=theta_112s, structural_pars=list(W=W_112)), theta_112, tolerance=1e-4)
  expect_equal(reform_structural_pars(p=1, M=2, d=2, params=theta_122s, structural_pars=list(W=W_122)), theta_122, tolerance=1e-4)
  expect_equal(reform_structural_pars(p=2, M=2, d=2, params=theta_222s, structural_pars=list(W=W_222)), theta_222, tolerance=1e-4)
  expect_equal(reform_structural_pars(p=3, M=3, d=2, params=theta_332s, structural_pars=list(W=W_332)), theta_332_froms, tolerance=1e-4)
  expect_equal(reform_structural_pars(p=1, M=2, d=3, params=theta_123s, structural_pars=list(W=W_123)), theta_123, tolerance=1e-4)
  expect_equal(reform_structural_pars(p=2, M=1, d=3, params=theta_213s, structural_pars=list(W=W_213)), theta_213, tolerance=1e-4)
})

allA_112 <- pick_allA(p=1, M=1, d=2, params=theta_112)
allA_122 <- pick_allA(p=1, M=2, d=2, params=theta_122)
allA_222 <- pick_allA(p=2, M=2, d=2, params=theta_222)
allA_332 <- pick_allA(p=3, M=3, d=2, params=theta_332)
allA_123 <- pick_allA(p=1, M=2, d=3, params=theta_123)
allA_213 <- pick_allA(p=2, M=1, d=3, params=theta_213)


lower_part <- function(p, d) {
  cbind(diag(nrow=d*(p-1)), matrix(0, nrow=d*(p - 1), ncol=d))
}

test_that("form_boldA works correctly", {
  expect_equal(form_boldA(p=1, M=1, d=2, all_A=allA_112)[, , 1], A11_112)

  expect_equal(form_boldA(p=1, M=2, d=2, all_A=allA_122)[, , 1], A11_122)
  expect_equal(form_boldA(p=1, M=2, d=2, all_A=allA_122)[, , 2], A21_122)

  expect_equal(form_boldA(p=2, M=2, d=2, all_A=allA_222)[, , 1], rbind(cbind(A11_222, A12_222), lower_part(p=2, d=2)))
  expect_equal(form_boldA(p=2, M=2, d=2, all_A=allA_222)[, , 2], rbind(cbind(A21_222, A22_222), lower_part(p=2, d=2)))

  expect_equal(form_boldA(p=3, M=3, d=2, all_A=allA_332)[, , 1], rbind(cbind(A11_332, A12_332, A13_332), lower_part(p=3, d=2)))
  expect_equal(form_boldA(p=3, M=3, d=2, all_A=allA_332)[, , 2], rbind(cbind(A21_332, A22_332, A23_332), lower_part(p=3, d=2)))
  expect_equal(form_boldA(p=3, M=3, d=2, all_A=allA_332)[, , 3], rbind(cbind(A31_332, A32_332, A33_332), lower_part(p=3, d=2)))

  expect_equal(form_boldA(p=1, M=2, d=3, all_A=allA_123)[, , 1], A11_123)
  expect_equal(form_boldA(p=1, M=2, d=3, all_A=allA_123)[, , 2], A21_123)

  expect_equal(form_boldA(p=2, M=1, d=3, all_A=allA_213)[, , 1], rbind(cbind(A11_213, A12_213), lower_part(p=2, d=3)))
})


## A(M)(p)_(p)(M)(d)

# p=3, M=3, d=2
alpha1_332_2 <- 0.1; alpha2_332_2 <- 0.4
theta_332_2 <- c(upsilon1_332, upsilon2_332, upsilon3_332, alpha1_332_2, alpha2_332_2)
alpha1_332_3 <- 0.5; alpha2_332_3 <- 0.2
theta_332_3 <- c(upsilon1_332, upsilon2_332, upsilon3_332, alpha1_332_3, alpha2_332_3)

theta_332s_2 <- c(phi10_332, phi20_332, phi30_332, vec(A11_332), vec(A12_332), vec(A13_332),
                  vec(A21_332), vec(A22_332), vec(A23_332), vec(A31_332), vec(A32_332), vec(A33_332),
                  Wvec(W_332), lambdas2_332, lambdas3_332, alpha1_332_2, alpha2_332_2) # SGMVAR

# p=3, M=4, d=2
upsilon1_342 <- upsilon1_332; upsilon2_342 <- upsilon2_332; upsilon3_342 <- upsilon3_332; upsilon4_342 <- upsilon3_342+0.07
alpha1_342 <- 0.1; alpha2_342 <- 0.2; alpha3_342 <- 0.3
theta_342 <- c(upsilon1_342, upsilon2_342, upsilon3_342, upsilon4_342, alpha1_342, alpha2_342, alpha3_342)
alpha1_342_2 <- 0.2; alpha2_342_2 <- 0.3; alpha3_342_2 <- 0.4
theta_342_2 <- c(upsilon1_342, upsilon2_342, upsilon3_342, upsilon4_342, alpha1_342_2, alpha2_342_2, alpha3_342_2)
alpha1_342_3 <- 0.3; alpha2_342_3 <- 0.4; alpha3_342_3 <- 0.1
theta_342_3 <- c(upsilon1_342, upsilon2_342, upsilon3_342, upsilon4_342, alpha1_342_3, alpha2_342_3, alpha3_342_3)

phi1_342 <- c(vec(A11_332), vec(A12_332), vec(A13_332))
phi2_342 <- c(vec(A21_332), vec(A22_332), vec(A23_332))
phi3_342 <- c(vec(A31_332), vec(A32_332), vec(A33_332))
phi4_342 <- phi3_342 + 0.1
phi30_342 <- phi30_332; phi20_342 <- phi20_332; phi10_342 <- phi10_332
phi40_342 <- phi30_342 + 0.1
lambdas2_342 <- lambdas2_332; lambdas3_342 <- lambdas3_332
lambdas4_342 <- lambdas3_342 + 0.2
W_342 <- W_332
theta_342s <- c(phi10_342, phi20_342, phi30_342, phi40_342, phi1_342, phi2_342, phi3_342, phi4_342,
                Wvec(W_342), lambdas2_342, lambdas3_342, lambdas4_342, alpha1_342, alpha2_342, alpha3_342) # SGMVAR
theta_342s_2 <- c(phi10_342, phi20_342, phi30_342, phi40_342, phi1_342, phi2_342, phi3_342, phi4_342,
                  Wvec(W_342), lambdas2_342, lambdas3_342, lambdas4_342, alpha1_342_2, alpha2_342_2, alpha3_342_2) # SGMVAR
theta_342s_3 <- c(phi10_342, phi20_342, phi30_342, phi40_342, phi1_342, phi2_342, phi3_342, phi4_342,
                  Wvec(W_342), lambdas2_342, lambdas3_342, lambdas4_342, alpha1_342_3, alpha2_342_3, alpha3_342_3) # SGMVAR

# p=1, M=2, d=3
alpha1_123_2 <- 0.3
theta_123_2 <- c(upsilon1_123, upsilon2_123, alpha1_123_2)

# p=1, M=3, d=3, structural_pars=list(W=W_133)
W_133 <- matrix(c(1, 2, 0, 4, 5, 6, 0, 8, 9), nrow=3)
phi10_133 <- c(0.388, 2.682, -0.388)
phi20_133 <- c(2.433, 5.961, 1.898)
phi30_133 <- c(-0.044, 3.139, 2.595)
A11_133 <- c(0.412, 0.03, 0.117, 0.005, -0.141, 0.038, -0.309, 0.252, 0.043)
A21_133 <- c(0.395, 0.088, -0.108, 0.111, -0.183, -0.227, 0.053, -0.087, 0.007)
A31_133 <- c(0.015, -0.114, -0.111, -0.027, 0.227, -0.275, 0.118, 0.064, 0.198)
W_133pars_with0 <- c(0.215, 0.302, 0, 0.384, 0.986, 1.421, 0, 1.296, 0.684)
lambdas2_133 <- c(1.277, 0.573, 1.225)
lambdas3_133 <- c(0.947, 1.241, 0.084)
alpha1_133_1 <- 0.5; alpha2_133_1 <- 0.1 # perm=c(1, 3, 2)
alpha1_133_2 <- 0.3; alpha2_133_2 <- 0.5 # perm=c(2, 1, 3)
alpha1_133_3 <- 0.1; alpha2_133_3 <- 0.2 # perm=c(3, 2, 1)
alpha1_133_4 <- 0.2; alpha2_133_4 <- 0.1 # perm=c(3, 1, 2)

theta_133s_1 <- c(phi10_133, phi20_133, phi30_133,
                  A11_133, A21_133, A31_133,
                  Wvec(W_133pars_with0), lambdas2_133, lambdas3_133,
                  alpha1_133_1, alpha2_133_1)
theta_133s_2 <- c(phi10_133, phi20_133, phi30_133,
                  A11_133, A21_133, A31_133,
                  Wvec(W_133pars_with0), lambdas2_133, lambdas3_133,
                  alpha1_133_2, alpha2_133_2)
theta_133s_3 <- c(phi10_133, phi20_133, phi30_133,
                  A11_133, A21_133, A31_133,
                  Wvec(W_133pars_with0), lambdas2_133, lambdas3_133,
                  alpha1_133_3, alpha2_133_3)
theta_133s_4 <- c(phi10_133, phi20_133, phi30_133,
                  A11_133, A21_133, A31_133,
                  Wvec(W_133pars_with0), lambdas2_133, lambdas3_133,
                  alpha1_133_4, alpha2_133_4)

set.seed(1); DAT3 <- matrix(round(rnorm(150), 3), ncol=3)

test_that("sort_components works correctly", {
  expect_equal(sort_components(p=1, M=1, d=2, params=theta_112), theta_112)
  expect_equal(sort_components(p=1, M=2, d=2, params=theta_122), c(upsilon2_122, upsilon1_122, 1-alpha1_122))
  expect_equal(sort_components(p=2, M=2, d=2, params=theta_222), c(upsilon2_222, upsilon1_222, 1-alpha1_222))

  expect_equal(sort_components(p=3, M=3, d=2, params=theta_332), theta_332)
  expect_equal(sort_components(p=3, M=3, d=2, params=theta_332_2), c(upsilon3_332, upsilon2_332, upsilon1_332, 1-alpha1_332_2-alpha2_332_2, alpha2_332_2))
  expect_equal(sort_components(p=3, M=3, d=2, params=theta_332_3), c(upsilon1_332, upsilon3_332, upsilon2_332, alpha1_332_3, 1-alpha1_332_3-alpha2_332_3))

  expect_equal(sort_components(p=3, M=4, d=2, params=theta_342), c(upsilon4_342, upsilon3_342, upsilon2_342, upsilon1_342, 1-alpha1_342-alpha2_342-alpha3_342, alpha3_342, alpha2_342))
  expect_equal(sort_components(p=3, M=4, d=2, params=theta_342_2), c(upsilon3_342, upsilon2_342, upsilon1_342, upsilon4_342, alpha3_342_2, alpha2_342_2, alpha1_342_2))
  expect_equal(sort_components(p=3, M=4, d=2, params=theta_342_3), c(upsilon2_342, upsilon1_342, upsilon4_342, upsilon3_342, alpha2_342_3, alpha1_342_3, 1-alpha1_342_3-alpha2_342_3-alpha3_342_3))

  expect_equal(sort_components(p=1, M=2, d=3, params=theta_123), theta_123)
  expect_equal(sort_components(p=1, M=2, d=3, params=theta_123_2),  c(upsilon2_123, upsilon1_123, 1-alpha1_123_2))
  expect_equal(sort_components(p=2, M=1, d=3, params=theta_213), theta_213)

  # SGMVAR (redecompose_Omegas is tested in matcal and should work)
  expect_equal(sort_components(p=1, M=2, d=2, params=theta_122s, structural_pars=list(W=W_122)),
               c(phi20_122, phi10_122, vec(A21_122),  vec(A11_122), redecompose_Omegas(M=2, d=2, W=W_122, lambdas=lambdas_122, perm=2:1),
                 1-alpha1_122), tolerance=1e-3)
  expect_equal(sort_components(p=3, M=3, d=2, params=theta_332s, structural_pars=list(W=W_332)), theta_332s)
  expect_equal(sort_components(p=3, M=3, d=2, params=theta_332s_2, structural_pars=list(W=W_332)),
               c(phi30_332, phi20_332, phi10_332, vec(A31_332), vec(A32_332), vec(A33_332),
                 vec(A21_332), vec(A22_332), vec(A23_332), vec(A11_332), vec(A12_332), vec(A13_332),
                 Wvec(redecompose_Omegas(M=3, d=2, W=W_332, lambdas=c(lambdas2_332, lambdas3_332), perm=c(3, 2, 1))),
                 1 - alpha1_332_2 - alpha2_332_2, alpha2_332_2), tolerance=1e-3)
  expect_equal(sort_components(p=3, M=4, d=2, params=theta_342s, structural_pars=list(W=W_342)), # perm=c(4, 3, 2, 1)
               c(phi40_342, phi30_342, phi20_342, phi10_342, phi4_342, phi3_342, phi2_342, phi1_342,
                 Wvec(redecompose_Omegas(M=4, d=2, W=W_342, lambdas=c(lambdas2_342, lambdas3_342, lambdas4_342), perm=c(4, 3, 2, 1))),
                 1 - alpha1_342 - alpha2_342 - alpha3_342, alpha3_342, alpha2_342), tolerance=1e-3)
  expect_equal(sort_components(p=3, M=4, d=2, params=theta_342s_2, structural_pars=list(W=W_342)), # perm=c(3, 2, 1, 4)
               c(phi30_342, phi20_342, phi10_342, phi40_342, phi3_342, phi2_342, phi1_342, phi4_342,
                 Wvec(redecompose_Omegas(M=4, d=2, W=W_342, lambdas=c(lambdas2_342, lambdas3_342, lambdas4_342), perm=c(3, 2, 1, 4))),
                 alpha3_342_2, alpha2_342_2, alpha1_342_2), tolerance=1e-3)
  expect_equal(sort_components(p=3, M=4, d=2, params=theta_342s_3, structural_pars=list(W=W_342)), # perm=c(2, 1, 4, 3)
               c(phi20_342, phi10_342, phi40_342, phi30_342, phi2_342, phi1_342, phi4_342, phi3_342,
                 Wvec(redecompose_Omegas(M=4, d=2, W=W_342, lambdas=c(lambdas2_342, lambdas3_342, lambdas4_342), perm=c(2, 1, 4, 3))),
                 alpha2_342_3, alpha1_342_3, 1 - alpha1_342_3 - alpha2_342_3 - alpha3_342_3), tolerance=1e-3)

  # p=1, M=3, d=3
  expect_equal(sort_components(p=1, M=3, d=3, params=theta_133s_1, structural_pars=list(W=W_133)), # # perm=c(1, 3, 2)
               c(phi10_133, phi30_133, phi20_133, A11_133, A31_133, A21_133,
                 Wvec(redecompose_Omegas(M=3, d=3, W=W_133pars_with0, lambdas=c(lambdas2_133, lambdas3_133), perm=c(1, 3, 2))),
                 alpha1_133_1, 1-alpha1_133_1-alpha2_133_1), tolerance=1e-3)
  expect_equal(loglikelihood(data=DAT3, p=1, M=3, params=sort_components(p=1, M=3, d=3, params=theta_133s_1, structural_pars=list(W=W_133)),
                             structural_pars=list(W=W_133), conditional=FALSE),
               loglikelihood(data=DAT3, p=1, M=3, params=theta_133s_1, structural_pars=list(W=W_133), conditional=FALSE), tolerance=1e-3)
  expect_equal(sort_components(p=1, M=3, d=3, params=theta_133s_2, structural_pars=list(W=W_133)), # # perm=c(2, 1, 3)
               c(phi20_133, phi10_133, phi30_133, A21_133, A11_133, A31_133,
                 Wvec(redecompose_Omegas(M=3, d=3, W=W_133pars_with0, lambdas=c(lambdas2_133, lambdas3_133), perm=c(2, 1, 3))),
                 alpha2_133_2, alpha1_133_2), tolerance=1e-3)
  expect_equal(loglikelihood(data=DAT3, p=1, M=3, params=sort_components(p=1, M=3, d=3, params=theta_133s_2, structural_pars=list(W=W_133)),
                             structural_pars=list(W=W_133), conditional=FALSE),
               loglikelihood(data=DAT3, p=1, M=3, params=theta_133s_2, structural_pars=list(W=W_133), conditional=FALSE), tolerance=1e-3)
  expect_equal(sort_components(p=1, M=3, d=3, params=theta_133s_3, structural_pars=list(W=W_133)), # # perm=c(3, 2, 1)
               c(phi30_133, phi20_133, phi10_133, A31_133, A21_133, A11_133,
                 Wvec(redecompose_Omegas(M=3, d=3, W=W_133pars_with0, lambdas=c(lambdas2_133, lambdas3_133), perm=c(3, 2, 1))),
                 1-alpha2_133_3-alpha1_133_3, alpha2_133_3), tolerance=1e-3)
  expect_equal(loglikelihood(data=DAT3, p=1, M=3, params=sort_components(p=1, M=3, d=3, params=theta_133s_3, structural_pars=list(W=W_133)),
                             structural_pars=list(W=W_133), conditional=FALSE),
               loglikelihood(data=DAT3, p=1, M=3, params=theta_133s_3, structural_pars=list(W=W_133), conditional=FALSE), tolerance=1e-3)
  expect_equal(sort_components(p=1, M=3, d=3, params=theta_133s_4, structural_pars=list(W=W_133)), # # perm=c(3, 1, 2)
               c(phi30_133, phi10_133, phi20_133, A31_133, A11_133, A21_133,
                 Wvec(redecompose_Omegas(M=3, d=3, W=W_133pars_with0, lambdas=c(lambdas2_133, lambdas3_133), perm=c(3, 1, 2))),
                 1-alpha2_133_4-alpha1_133_4, alpha1_133_4), tolerance=1e-3)
  expect_equal(loglikelihood(data=DAT3, p=1, M=3, params=sort_components(p=1, M=3, d=3, params=theta_133s_4, structural_pars=list(W=W_133)),
                             structural_pars=list(W=W_133), conditional=FALSE),
               loglikelihood(data=DAT3, p=1, M=3, params=theta_133s_4, structural_pars=list(W=W_133), conditional=FALSE), tolerance=1e-3)
})


calc_mu <- function(p, M, d, params, constraints=NULL, structural_pars=NULL) {
  params <- reform_constrained_pars(p, M, d, params, constraints=constraints, same_means=NULL,
                                    structural_pars=structural_pars)
  all_A <- pick_allA(p=p, M=M, d=d, params=params, structural_pars=structural_pars)
  all_phi0 <- pick_phi0(p=p, M=M, d=d, params=params, structural_pars=structural_pars)
  vapply(1:M, function(m) solve(diag(d) - rowSums(all_A[, , , m, drop=FALSE], dims=2), all_phi0[,m]), numeric(d))
}

theta_112_mu <- change_parametrization(p=1, M=1, d=2, params=theta_112, change_to="mean")
theta_122_mu <- change_parametrization(p=1, M=2, d=2, params=theta_122, change_to="mean")
theta_222_mu <- change_parametrization(p=2, M=2, d=2, params=theta_222, change_to="mean")
theta_332_mu <- change_parametrization(p=3, M=3, d=2, params=theta_332, change_to="mean")
theta_342_mu <- change_parametrization(p=3, M=4, d=2, params=theta_342, change_to="mean")
theta_123_mu <- change_parametrization(p=1, M=2, d=3, params=theta_123, change_to="mean")
theta_213_mu <- change_parametrization(p=2, M=1, d=3, params=theta_213, change_to="mean")

theta_112c_mu <- change_parametrization(p=1, M=1, d=2, params=theta_112c, constraints=C_112, change_to="mean")
theta_222c_mu <- change_parametrization(p=2, M=2, d=2, params=theta_222c, constraints=C_222, change_to="mean")
theta_222c_mu2 <- change_parametrization(p=2, M=2, d=2, params=theta_222_c2, constraints=C_222_2, change_to="mean")
theta_123c_mu <- change_parametrization(p=1, M=2, d=3, params=theta_123c, constraints=C_123, change_to="mean")

# SGMVAR
theta_112s_mu <- change_parametrization(p=1, M=1, d=2, params=theta_112s, structural_pars=list(W=W_112), change_to="mean")
theta_222s_mu <- change_parametrization(p=2, M=2, d=2, params=theta_222s, structural_pars=list(W=W_222), change_to="mean")
theta_123s_mu <- change_parametrization(p=1, M=2, d=3, params=theta_123s, structural_pars=list(W=W_123), change_to="mean")

theta_112csWAR_mu <- change_parametrization(p=1, M=1, d=2, params=theta_112csWAR, structural_pars=list(W=W_112), constraints=C_112, change_to="mean")
theta_222csLAR_mu <- change_parametrization(p=2, M=2, d=2, params=theta_222csLAR, structural_pars=list(W=W_222), constraints=C_222, change_to="mean")
theta_123csL_mu <- change_parametrization(p=1, M=2, d=3, params=theta_123csL, structural_pars=list(W=W_123), constraints=NULL, change_to="mean")

test_that("change_parametrization works correctly", {
  expect_equal(pick_phi0(p=1, M=1, d=2, params=theta_112_mu), calc_mu(p=1, M=1, d=2, params=theta_112))
  expect_equal(change_parametrization(p=1, M=1, d=2, params=theta_112_mu, change_to="intercept"), theta_112)

  expect_equal(pick_phi0(p=1, M=2, d=2, params=theta_122_mu), calc_mu(p=1, M=2, d=2, params=theta_122))
  expect_equal(change_parametrization(p=1, M=2, d=2, params=theta_122_mu, change_to="intercept"), theta_122)

  expect_equal(pick_phi0(p=2, M=2, d=2, params=theta_222_mu), calc_mu(p=2, M=2, d=2, params=theta_222))
  expect_equal(change_parametrization(p=2, M=2, d=2, params=theta_222_mu, change_to="intercept"), theta_222)

  expect_equal(pick_phi0(p=3, M=3, d=2, params=theta_332_mu), calc_mu(p=3, M=3, d=2, params=theta_332))
  expect_equal(change_parametrization(p=3, M=3, d=2, params=theta_332_mu, change_to="intercept"), theta_332)

  expect_equal(pick_phi0(p=3, M=4, d=2, params=theta_342_mu), calc_mu(p=3, M=4, d=2, params=theta_342))
  expect_equal(change_parametrization(p=3, M=4, d=2, params=theta_342_mu, change_to="intercept"), theta_342)

  expect_equal(pick_phi0(p=1, M=2, d=3, params=theta_123_mu), calc_mu(p=1, M=2, d=3, params=theta_123))
  expect_equal(change_parametrization(p=1, M=2, d=3, params=theta_123_mu, change_to="intercept"), theta_123)

  expect_equal(pick_phi0(p=2, M=1, d=3, params=theta_213_mu), calc_mu(p=2, M=1, d=3, params=theta_213))
  expect_equal(change_parametrization(p=2, M=1, d=3, params=theta_213_mu, change_to="intercept"), theta_213)

  expect_equal(matrix(theta_112c_mu[1:(1*2)], nrow=2, byrow=FALSE), calc_mu(p=1, M=1, d=2, constraints=C_112, params=theta_112c))
  expect_equal(change_parametrization(p=1, M=1, d=2, params=theta_112c_mu, constraints=C_112, change_to="intercept"), theta_112c)

  expect_equal(matrix(theta_222c_mu[1:(2*2)], nrow=2, byrow=FALSE), calc_mu(p=2, M=2, d=2, constraints=C_222, params=theta_222c))
  expect_equal(change_parametrization(p=2, M=2, d=2, params=theta_222c_mu, constraints=C_222, change_to="intercept"), theta_222c)

  expect_equal(matrix(theta_222c_mu2[1:(2*2)], nrow=2, byrow=FALSE), calc_mu(p=2, M=2, d=2, constraints=C_222_2, params=theta_222_c2))
  expect_equal(change_parametrization(p=2, M=2, d=2, params=theta_222c_mu2, constraints=C_222_2, change_to="intercept"), theta_222_c2)

  expect_equal(matrix(theta_123c_mu[1:(2*3)], nrow=3, byrow=FALSE), calc_mu(p=1, M=2, d=3, constraints=C_123, params=theta_123c))
  expect_equal(change_parametrization(p=1, M=2, d=3, params=theta_123c_mu, constraints=C_123, change_to="intercept"), theta_123c)

  # SGMVAR
  expect_equal(pick_phi0(p=1, M=1, d=2, params=theta_112s_mu, structural_pars=list(W=W_112)),
               calc_mu(p=1, M=1, d=2, params=theta_112s, structural_pars=list(W=W_112)))
  expect_equal(change_parametrization(p=1, M=1, d=2, params=theta_112s_mu, structural_pars=list(W=W_112), change_to="intercept"), theta_112s)

  expect_equal(pick_phi0(p=2, M=2, d=2, params=theta_222s_mu, structural_pars=list(W=W_222)),
               calc_mu(p=2, M=2, d=2, params=theta_222s, structural_pars=list(W=W_222)))
  expect_equal(change_parametrization(p=2, M=2, d=2, params=theta_222s_mu, structural_pars=list(W=W_222), change_to="intercept"), theta_222s)

  expect_equal(pick_phi0(p=1, M=2, d=3, params=theta_123s_mu, structural_pars=list(W=W_123)),
               calc_mu(p=1, M=2, d=3, params=theta_123s, structural_pars=list(W=W_123)))
  expect_equal(change_parametrization(p=1, M=2, d=3, params=theta_123s_mu, change_to="intercept", structural_pars=list(W=W_123)), theta_123s)

  expect_equal(matrix(theta_112csWAR_mu[1:(1*2)], nrow=2, byrow=FALSE),
               calc_mu(p=1, M=1, d=2, constraints=C_112, structural_pars=list(W=W_112), params=theta_112csWAR))
  expect_equal(change_parametrization(p=1, M=1, d=2, params=theta_112csWAR_mu, constraints=C_112, structural_pars=list(W=W_112), change_to="intercept"),
               theta_112csWAR)

  expect_equal(matrix(theta_222csLAR_mu[1:(2*2)], nrow=2, byrow=FALSE),
               calc_mu(p=2, M=2, d=2, constraints=C_222, structural_pars=list(W=W_222), params=theta_222csLAR))
  expect_equal(change_parametrization(p=2, M=2, d=2, params=theta_222csLAR_mu, constraints=C_222, structural_pars=list(W=W_222), change_to="intercept"),
               theta_222csLAR)

  expect_equal(matrix(theta_123csL_mu[1:(2*3)], nrow=3, byrow=FALSE),
               calc_mu(p=1, M=2, d=3, constraints=NULL, structural_pars=list(W=W_123), params=theta_123csL))
  expect_equal(change_parametrization(p=1, M=2, d=3, params=theta_123csL_mu, constraints=NULL, structural_pars=list(W=W_123), change_to="intercept"),
               theta_123csL)
})


theta_122_cr <- c(upsilon2_122, upsilon2_122, alpha1_122)
theta_222_cr <- c(upsilon1_222, upsilon1_222, alpha1_222)
theta_332_cr <- c(upsilon1_332, upsilon2_332, upsilon2_332, alpha1_332, alpha2_332)
theta_123_cr1 <- c(upsilon1_123, upsilon1_123, alpha1_123)
theta_123_cr2 <- c(upsilon2_123, upsilon2_123, alpha1_123)

theta_222s <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A21_222),
                vec(A22_222), vec(W_222), lambdas_222, alpha1_222)

## A(M)(p)_(p)(M)(d)
theta_122_crs <- c(phi10_122, phi10_122, vec(A11_122), vec(A11_122), Wvec(W_122), lambdas_122, alpha1_122)
theta_122_crs2 <- c(phi20_122, phi20_122, vec(A21_122), vec(A21_122), Wvec(W_122), lambdas_122, alpha1_122)
theta_222_crs <- c(phi10_222, phi10_222, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222), Wvec(W_222), lambdas_222, alpha1_222)
theta_332_crs <- c(phi10_332, phi20_332, phi20_332, vec(A11_332), vec(A12_332), vec(A13_332),
                   vec(A21_332), vec(A22_332), vec(A23_332), vec(A21_332), vec(A22_332), vec(A23_332),
                   Wvec(W_332), lambdas2_332, lambdas3_332, alpha1_332, alpha2_332)
theta_123_crs <- c(phi10_123, phi10_123, vec(A11_123), vec(A11_123), Wvec(W_123), lambdas_123, alpha1_123)
theta_123_crs2 <- c(phi20_123, phi20_123, vec(A21_123), vec(A21_123), Wvec(W_123), lambdas_123, alpha1_123)



test_that("change_regime works correctly", {
  expect_equal(change_regime(p=1, M=2, d=2, params=theta_122, m=1, regime_pars=upsilon2_122), theta_122_cr)
  expect_equal(change_regime(p=2, M=2, d=2, params=theta_222, m=2, regime_pars=upsilon1_222), theta_222_cr)
  expect_equal(change_regime(p=3, M=3, d=2, params=theta_332, m=3, regime_pars=upsilon2_332), theta_332_cr)
  expect_equal(change_regime(p=1, M=2, d=3, params=theta_123, m=2, regime_pars=upsilon1_123), theta_123_cr1)
  expect_equal(change_regime(p=1, M=2, d=3, params=theta_123, m=1, regime_pars=upsilon2_123), theta_123_cr2)

  # SGMVAR
  rpars122_m1 <-c(phi20_122, vec(A21_122))
  rpars122_m2 <-c(phi10_122, vec(A11_122))
  expect_equal(change_regime(p=1, M=2, d=2, params=theta_122s, m=2, regime_pars=rpars122_m2,
                             structural_pars=list(W=W_122)), theta_122_crs, tol=1-6)
  expect_equal(change_regime(p=1, M=2, d=2, params=theta_122s, m=1, regime_pars=rpars122_m1,
                             structural_pars=list(W=W_122)), theta_122_crs2, tol=1e-6)

  rpars222_m2 <- c(phi10_222, vec(A11_222), vec(A12_222))
  expect_equal(change_regime(p=2, M=2, d=2, params=theta_222s, m=2, regime_pars=rpars222_m2,
                             structural_pars=list(W=W_222)), theta_222_crs, tol=1-6)

  rpars332_m3 <- c(phi20_332, vec(A21_332), vec(A22_332), vec(A23_332))
  expect_equal(change_regime(p=3, M=3, d=2, params=theta_332sWC, m=3, regime_pars=rpars332_m3,
                             structural_pars=list(W=W_332)), theta_332_crs, tol=1-6)

  rpars123_m2 <- c(phi10_123, vec(A11_123))
  expect_equal(change_regime(p=1, M=2, d=3, params=theta_123s, m=2, regime_pars=rpars123_m2,
                             structural_pars=list(W=W_123)), theta_123_crs, tol=1-6)
  rpars123_m1 <- c(phi20_123, vec(A21_123))
  expect_equal(change_regime(p=1, M=2, d=3, params=theta_123s, m=1, regime_pars=rpars123_m1,
                             structural_pars=list(W=W_123)), theta_123_crs2, tol=1-6)
})


params122s <- c(1.03, 2.36, 1.79, 3, 1, -0.06, -0.04, 1, 1.02, -0.04, -0.02, 1.02, -0.89, -0.72, 0.37, -2.16, 7.16, 1.3, 0.37)
params222s <- c(1.03, 2.36, 1.79, 3, 1.25, 0.06, 0.04, 1.34, -0.29, -0.08, -0.05, -0.36, 1.2, 0.05, 0.05, 1.3, -0.3, -0.1,
                -0.05, -0.4, -0.89, -0.72, 0.37, -2.16, 7.16, 1.3, 0.37)
params222s_2 <- c(1.03, 2.36, 1.79, 3, 1.25, 0.06, 0.04, 1.34, -0.29, -0.08, -0.05, -0.36, 1.2, 0.05, 0.05, 1.3, -0.3, -0.1,
                  -0.05, -0.4, -0.89, -0.72, 0.37, -2.16, 1.16, 1.3, 0.37)
params123s <- c(1.1, 2.2, 3.3, 1.11, 2.22, 3.33, 1, 0.21, 0.31, 0.12, 2, 0.32, 0.13, 0.23, 3, -1, -0.21, -0.31, -0.12, -2,
                -0.32, -0.13, -0.23, -3, 0.47, 0.39, -1.23, 0.58, -1, 0.18, -0.67, -0.92, -1.2, 1.13, 1.12, 1.08, 0.6)
params123s_2 <- c(1.1, 2.2, 3.3, 1.11, 2.22, 3.33, 1, 0.21, 0.31, 0.12, 2, 0.32, 0.13, 0.23, 3, -1, -0.21, -0.31, -0.12, -2,
                  -0.32, -0.13, -0.23, -3, 0.47, 0.39, -1.23, 0.58, -1, 0.18, -0.67, -0.92, -1.2, 1.08, 1.13, 1.12, 0.6)
params123s_3 <- c(1.1, 2.2, 3.3, 1.11, 2.22, 3.33, 1, 0.21, 0.31, 0.12, 2, 0.32, 0.13, 0.23, 3, -1, -0.21, -0.31, -0.12, -2,
                  -0.32, -0.13, -0.23, -3, 0.47, 0.39, -1.23, 0.58, -1, 0.18, -0.67, -0.92, -1.2, 1.12, 1.08, 1.13, 0.6)
params222cs <- c(-0.11, 2.83, 0.36, 3.19, 1.26, 1.34, -0.29, -0.36, -0.88, -0.76, 0.46, -2.16, 6.97, 1.2, 0.35)
params123cs <- c(1.1, 2.2, 3.3, 1.11, 2.22, 3.33, 1, 0.21, 0.31, 0.12, 2, 0.32, 0.13, 0.23, 3, -1, -0.21, -0.31, -0.12, -2,
                 -0.32, -0.13, -0.23, -3, 0.47, 0.39, -1.23, 0.58, -1, 0.18, -0.67, -0.92, -1.2, 3, 1, 2, 0.6)
params122csm <- c(1.03, 2.36, 1, -0.06, -0.04, 1, -0.89, -0.72, 0.37, -2.16, 2, 1, 0.37)
params222csm <- c(1.03, 2.36, 1.25, 0.06, 0.04, 1.34, -0.29, -0.08, -0.05, -0.36, -0.89, -0.72, 0.37, -2.16, 2, 1, 0.37)
params123csm <- c(1.1, 2.2, 3.3, 1, 0.21, 0.31, 0.12, 2, 0.32, 0.13, 0.23, 3, -1, -0.21, -0.31, -0.12, -2, -0.32, -0.13,
                  -0.23, -3, 0.47, 0.39, -1.23, 0.58, -1, 0.18, -0.67, -0.92, -1.2, 3, 2, 1, 0.6)

test_that("sort_W_and_lambdas works correctly", {
  expect_equal(sort_W_and_lambdas(p=1, M=2, d=2, params=params122s),
                c(1.03, 2.36, 1.79, 3, 1, -0.06, -0.04, 1, 1.02, -0.04, -0.02, 1.02, 0.37, -2.16, -0.89, -0.72, 1.3, 7.16, 0.37))
  expect_equal(sort_W_and_lambdas(p=2, M=2, d=2, params=params222s),
               c(1.03, 2.36, 1.79, 3, 1.25, 0.06, 0.04, 1.34, -0.29, -0.08, -0.05, -0.36, 1.2, 0.05, 0.05, 1.3, -0.3, -0.1,
                 -0.05, -0.4, 0.37, -2.16, -0.89, -0.72, 1.3, 7.16, 0.37))
  expect_equal(sort_W_and_lambdas(p=2, M=2, d=2, params=params222s_2), params222s_2)
  expect_equal(sort_W_and_lambdas(p=1, M=2, d=3, params=params123s),
               c(1.1, 2.2, 3.3, 1.11, 2.22, 3.33, 1, 0.21, 0.31, 0.12, 2, 0.32, 0.13, 0.23, 3, -1, -0.21, -0.31, -0.12, -2,
                 -0.32, -0.13, -0.23, -3, -0.67, -0.92, -1.2, 0.58, -1, 0.18, 0.47, 0.39, -1.23, 1.08, 1.12, 1.13, 0.6))
  expect_equal(sort_W_and_lambdas(p=1, M=2, d=3, params=params123s_2),
               c(1.1, 2.2, 3.3, 1.11, 2.22, 3.33, 1, 0.21, 0.31, 0.12, 2, 0.32, 0.13, 0.23, 3, -1, -0.21, -0.31, -0.12, -2,
                 -0.32, -0.13, -0.23, -3, 0.47, 0.39, -1.23, -0.67, -0.92, -1.2, 0.58, -1, 0.18, 1.08, 1.12, 1.13, 0.6))
  expect_equal(sort_W_and_lambdas(p=1, M=2, d=3, params=params123s_3),
               c(1.1, 2.2, 3.3, 1.11, 2.22, 3.33, 1, 0.21, 0.31, 0.12, 2, 0.32, 0.13, 0.23, 3, -1, -0.21, -0.31, -0.12, -2,
                 -0.32, -0.13, -0.23, -3, 0.58, -1, 0.18, 0.47, 0.39, -1.23, -0.67, -0.92, -1.2, 1.08, 1.12, 1.13, 0.6))
  expect_equal(sort_W_and_lambdas(p=2, M=2, d=2, params=params222cs),
               c(-0.11, 2.83, 0.36, 3.19, 1.26, 1.34, -0.29, -0.36, 0.46, -2.16, -0.88, -0.76, 1.2, 6.97, 0.35))
  expect_equal(sort_W_and_lambdas(p=1, M=2, d=3, params=params123cs),
               c(1.1, 2.2, 3.3, 1.11, 2.22, 3.33, 1, 0.21, 0.31, 0.12, 2, 0.32, 0.13, 0.23, 3, -1, -0.21, -0.31, -0.12, -2,
                 -0.32, -0.13, -0.23, -3, 0.58, -1, 0.18, -0.67, -0.92, -1.2, 0.47, 0.39, -1.23, 1, 2, 3, 0.6))
  expect_equal(sort_W_and_lambdas(p=1, M=2, d=2, params=params122csm),
               c(1.03, 2.36, 1, -0.06, -0.04, 1, 0.37, -2.16, -0.89, -0.72, 1, 2, 0.37))
  expect_equal(sort_W_and_lambdas(p=2, M=2, d=2, params=params222csm),
               c(1.03, 2.36, 1.25, 0.06, 0.04, 1.34, -0.29, -0.08, -0.05, -0.36, 0.37, -2.16, -0.89, -0.72, 1, 2, 0.37))
  expect_equal(sort_W_and_lambdas(p=1, M=2, d=3, params=params123csm),
               c(1.1, 2.2, 3.3, 1, 0.21, 0.31, 0.12, 2, 0.32, 0.13, 0.23, 3, -1, -0.21, -0.31, -0.12, -2, -0.32, -0.13,
                 -0.23, -3, -0.67, -0.92, -1.2, 0.58, -1, 0.18, 0.47, 0.39, -1.23, 1, 2, 3, 0.6))
})


test_that("sort_and_standardize_alphas works correctly", {
  expect_equal(sort_and_standardize_alphas(alphas=c(0.3, 0.1, 0.6), constraints=1), c(0.3, 0.1, 0.6))
  expect_equal(sort_and_standardize_alphas(alphas=c(0.3, 0.2, 0.5), same_means=1), c(0.3, 0.2, 0.5))
  expect_equal(sort_and_standardize_alphas(alphas=c(0.3, 0.2, 0.5), structural_pars=list(W=1, C_lambda=1)), c(0.3, 0.2, 0.5))
  expect_equal(sort_and_standardize_alphas(alphas=c(3, 1, 6), constraints=1), c(0.3, 0.1, 0.6))

  expect_equal(sort_and_standardize_alphas(alphas=c(0.9, 0.1)), c(0.9, 0.1))
  expect_equal(sort_and_standardize_alphas(alphas=c(0.4, 0.6)), c(0.6, 0.4))
  expect_equal(sort_and_standardize_alphas(alphas=c(4, 6)), c(0.6, 0.4))
  expect_equal(sort_and_standardize_alphas(alphas=c(0.3, 0.1, 0.6)), c(0.6, 0.3, 0.1))

  expect_equal(sort_and_standardize_alphas(alphas=c(8, 2), structural_pars=list(1)), c(0.8, 0.2))
  expect_equal(sort_and_standardize_alphas(alphas=c(0.3, 0.7), structural_pars=list(1)), c(0.7, 0.3))
  expect_equal(sort_and_standardize_alphas(alphas=c(0.3, 0.1, 0.6), structural_pars=list(1)), c(0.6, 0.3, 0.1))
  expect_equal(sort_and_standardize_alphas(alphas=c(3, 2, 5), structural_pars=list(1)), c(0.5, 0.3, 0.2))
})


