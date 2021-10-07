context("Pick parameters")
library(gmvarkit)


## A(M)(p)_(p)(M)(d)

# p=1, M=1, d=2
phi10_112 <- c(1.03, 2.36)
A11_112 <- matrix(c(1.25, 0.06, 0.04, 1.34), nrow=2, byrow=FALSE)
Omega1_112 <- matrix(c(0.93, -0.15, -0.15, 5.20), nrow=2, byrow=FALSE)

theta_112 <- upsilon1_112 <- c(phi10_112, vec(A11_112), vech(Omega1_112))

W_112 <- t(chol(Omega1_112))
theta_112s <- c(phi10_112, vec(A11_112), vec(W_112)) # SGMVAR
Omega1_112s <- tcrossprod(W_112)

# p=1, M=2, d=2
phi10_122 <- c(1.03, 2.36)
A11_122 <- matrix(c(1, -0.06, -0.04, 1), nrow=2, byrow=FALSE)
Omega1_122 <- matrix(c(0.93, -0.15, -0.15, 5.20), nrow=2, byrow=FALSE)

phi20_122 <- c(1.79, 3.00)
A21_122 <- A11_122
Omega2_122 <- matrix(c(5.88, 3.56, 3.56, 9.80), nrow=2, byrow=FALSE)

alpha1_122 <- 0.37
upsilon1_122 <- c(phi10_122, vec(A11_122), vech(Omega1_122))
upsilon2_122 <- c(phi20_122, vec(A21_122), vech(Omega2_122))
theta_122 <- c(upsilon1_122, upsilon2_122, alpha1_122)

W_122 <- matrix(c(-0.89, -0.71, 0.36, -2.10), nrow=2, ncol=2, byrow=FALSE)
lambdas_122 <- c(7.20, 1.30)
theta_122s <- c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), vec(W_122), lambdas_122, alpha1_122) # SGMVAR
Omega1_122s <- tcrossprod(W_122)
Omega2_122s <- W_122%*%tcrossprod(diag(lambdas_122), W_122)

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

W_222 <- matrix(c(-0.88, -0.72, 0.37, -2.20), nrow=2, ncol=2, byrow=FALSE)
lambdas_222 <- c(7.10, 1.40)
theta_222s <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A21_222),
                vec(A22_222), vec(W_222), lambdas_222, alpha1_222) # SGMVAR
Omega1_222s <- tcrossprod(W_222)
Omega2_222s <- W_222%*%tcrossprod(diag(lambdas_222), W_222)

theta_222t <- c(theta_222, 20, 25) # StMVAR
theta_222gs <- c(theta_222, 25) # G-StMVAR

theta_222ts <- c(theta_222s, 10, 15) # SStMVAR
theta_222gss <- c(theta_222s, 15) # SG-StMVAR


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

W_332 <- matrix(c(-0.8924620, -0.7180539, 0.3653923, -2.1643472), nrow=2, byrow=FALSE)
lambdas2_332 <- c(7.16, 1.30)
lambdas3_332 <- c(6.10, 1.22)
theta_332s <- c(phi10_332, phi20_332, phi30_332, vec(A11_332), vec(A12_332), vec(A13_332),
                vec(A21_332), vec(A22_332), vec(A23_332), vec(A31_332), vec(A32_332), vec(A33_332),
                vec(W_332), lambdas2_332, lambdas3_332, alpha1_332, alpha2_332) # SGMVAR
Omega1_332s <- tcrossprod(W_332)
Omega2_332s <- W_332%*%tcrossprod(diag(lambdas2_332), W_332)
Omega3_332s <- W_332%*%tcrossprod(diag(lambdas3_332), W_332)

theta_332t <- c(theta_332, 10, 20, 30) # StMVAR
theta_332gs <- c(theta_332, 20, 30) # G-StMVAR, M1=1, M2=2

theta_332ts <- c(theta_332s, 10, 20, 30) # SStMVAR
theta_332gss <- c(theta_332s, 30) # SG-StMVAR, M1=2, M2=1

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

W_123 <- matrix(c(0.4694488, 0.3943046, -1.2316702, 0.5795625, -0.9982353, 0.1833377,
                  -0.6661270, -0.9208964, -1.2039003), nrow=3, byrow=FALSE)
lambdas_123 <- c(1.1221914, 1.1166332, 1.0763873)
theta_123s <- c(phi10_123, phi20_123, vec(A11_123), vec(A21_123), vec(W_123), lambdas_123, alpha1_123) # SGMVAR
Omega1_123s <- tcrossprod(W_123)
Omega2_123s <- W_123%*%tcrossprod(diag(lambdas_123), W_123)

theta_123t <- c(theta_123, 20, 25) # StMVAR
theta_123gs <- c(theta_123, 25) # G-StMVAR

theta_123ts <- c(theta_123s, 20, 25) # StMVAR
theta_123gss <- c(theta_123s, 25) # G-StMVAR

# p=2, M=1, d=3
phi10_213 <- c(1.1, 2.2, 3.3)
A11_213 <- matrix(c(1, 0.21, 0.31, 0.12, 2, 0.32, 0.13, 0.23, 3), nrow=3, byrow=FALSE)
A12_213 <- matrix(c(-1, -0.21, -0.31, -0.12, -2, -0.32, -0.13, -0.23, -3), nrow=3, byrow=FALSE)
Omega1_213 <- matrix(c(1, 0.22, 0.33, 0.22, 2, 0.44, 0.33, 0.44, 3), nrow=3, byrow=FALSE)

upsilon1_213 <- c(phi10_213, vec(A11_213), vec(A12_213), vech(Omega1_213))
theta_213 <- upsilon1_213

W_213 <- t(chol(Omega1_213))
theta_213s <- c(phi10_213, vec(A11_213), vec(A12_213), vec(W_213)) # SGMVAR
Omega1_213s <- tcrossprod(W_213)


## A(M)(p)_(p)(M)(d)
rbind_diags <- function(p, M, d) {
  I <- diag(p*d^2)
  Reduce(rbind, replicate(M, I, simplify=FALSE))
}

# Constraining AR-parameters to be the same for all regimes

# p=1, M=1, d=2
C_112 <- rbind_diags(p=1, M=1, d=2)
theta_112c <- c(phi10_112, vec(A11_112), vech(Omega1_112))

theta_112csWAR <- c(phi10_112, vec(A11_112), Wvec(W_112)) # SGMVAR W and AR

theta_112tc <- c(theta_112c, 10) # StMVAR
theta_112tcsWAR <- c(theta_112csWAR, 10) # StMVAR

# p=2, M=2, d=2
C_222 <- rbind_diags(p=2, M=2, d=2)
theta_222c <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vech(Omega1_222), vech(Omega2_222), alpha1_222)

C_lambda_222 <- matrix(c(1, 2), nrow=2)
theta_222csL <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A21_222),
                  vec(A22_222), vec(W_222), 0.2, alpha1_222) # SGMVAR lambdas
theta_222csLAR <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(W_222), 0.2, alpha1_222) # SGMVAR lambdas and AR

theta_222gscsL <- c(theta_222csL, 20) # SG-StMVAR, M1=1, M2=1
theta_222tcsLAR <- c(theta_222csLAR, 10, 20) # SStMVAR

# p=1, M=2, d=3
C_123 <- rbind_diags(p=1, M=2, d=3)
theta_123c <- c(phi10_123, phi20_123, vec(A11_123), vech(Omega1_123), vech(Omega2_123), alpha1_123)

C_lambda_123 <- matrix(c(1, 1, 0, 0, 0, 1), nrow=3, byrow=FALSE)
theta_123csL <- c(phi10_123, phi20_123, vec(A11_123), vec(A21_123), vec(W_123), 1, 2, alpha1_123) # SGMVAR lambdas
theta_123csLAR <- c(phi10_123, phi20_123, vec(A11_123), vec(W_123), 1, 2, alpha1_123) # SGMVAR lambdas and AR

theta_123tc <- c(theta_123c, 10, 20) # StMVAR
theta_123tcsL <- c(theta_123csL, 10, 20) # StMVAR

## Models with same_means

# p=1, M=1, d=2, same_means=list(1)
theta_112_int <- c(phi10_112, vec(A11_112), vech(Omega1_112))
theta_112t_int <- c(theta_112_int, 10)

# p=2, M=2, d=2, constraints=C_222, same_means=list(1:2)
theta_222c_int <- c(phi10_222, vec(A11_222), vec(A12_222), vech(Omega1_222), vech(Omega2_222), alpha1_222)
theta_222tc_int <- c(theta_222c_int, 10, 20) # StMVAR

# p=3, M=3, d=2, constraints=C_332, same_means=list(1, 2:3)
C_332 <- rbind_diags(p=3, M=3, d=2)
theta_332c_int <- c(phi10_332, phi20_332, vec(A11_332), vec(A12_332), vec(A13_332), vech(Omega1_332), vech(Omega2_332),
                    vech(Omega3_332), alpha1_332, alpha2_332)
theta_332gsc_int <- c(theta_332c_int, 20, 30) # G-StMVAR, M1=1, M2=2

# p=1, M=2, d=3, constraints=C_123, structural_pars=list(W=W_123, C_lambda=C_lambda_123) same_means=list(1:2)
theta_123csLAR_int <- c(phi10_123, vec(A11_123), vec(W_123), 1, 2, alpha1_123) # SGMVAR lambdas and AR
theta_123tcsLAR_int <- c(theta_123csLAR_int, 10, 20) # SStMVAR



test_that("pick_Ami works correctly", {
  expect_equal(pick_Ami(p=1, M=1, d=2, params=theta_112, m=1, i=1, unvec=TRUE), A11_112)

  expect_equal(pick_Ami(p=1, M=2, d=2, params=theta_122, m=1, i=1, unvec=TRUE), A11_122)
  expect_equal(pick_Ami(p=1, M=2, d=2, params=theta_122, m=2, i=1, unvec=TRUE), A21_122)

  expect_equal(pick_Ami(p=2, M=2, d=2, params=theta_222, m=1, i=1, unvec=TRUE), A11_222)
  expect_equal(pick_Ami(p=2, M=2, d=2, params=theta_222, m=2, i=2, unvec=TRUE), A22_222)

  expect_equal(pick_Ami(p=2, M=2, d=2, params=theta_222t, m=1, i=1, unvec=TRUE), A11_222) # StMVAR
  expect_equal(pick_Ami(p=2, M=2, d=2, params=theta_222t, m=2, i=2, unvec=TRUE), A22_222) # StMVAR

  expect_equal(pick_Ami(p=2, M=c(1, 1), d=2, params=theta_222gs, m=1, i=1, unvec=TRUE), A11_222) # G-StMVAR
  expect_equal(pick_Ami(p=2, M=c(1, 1), d=2, params=theta_222gs, m=2, i=2, unvec=TRUE), A22_222) # G-StMVAR

  expect_equal(pick_Ami(p=3, M=3, d=2, params=theta_332, m=1, i=2, unvec=TRUE), A12_332)
  expect_equal(pick_Ami(p=3, M=3, d=2, params=theta_332, m=2, i=3, unvec=TRUE), A23_332)
  expect_equal(pick_Ami(p=3, M=3, d=2, params=theta_332, m=3, i=1, unvec=TRUE), A31_332)

  expect_equal(pick_Ami(p=1, M=2, d=3, params=theta_123, m=1, i=1, unvec=TRUE), A11_123)
  expect_equal(pick_Ami(p=1, M=2, d=3, params=theta_123, m=2, i=1, unvec=TRUE), A21_123)

  expect_equal(pick_Ami(p=2, M=1, d=3, params=theta_213, m=1, i=1, unvec=TRUE), A11_213)
  expect_equal(pick_Ami(p=2, M=1, d=3, params=theta_213, m=1, i=2, unvec=TRUE), A12_213)

  # Structural models
  expect_equal(pick_Ami(p=1, M=1, d=2, params=theta_112s, m=1, i=1, structural_pars=list(W=W_112), unvec=TRUE), A11_112)

  expect_equal(pick_Ami(p=1, M=2, d=2, params=theta_122s, m=1, i=1, structural_pars=list(W=W_122), unvec=TRUE), A11_122)
  expect_equal(pick_Ami(p=1, M=2, d=2, params=theta_122s, m=2, i=1, structural_pars=list(W=W_122), unvec=TRUE), A21_122)

  expect_equal(pick_Ami(p=2, M=2, d=2, params=theta_222s, m=1, i=1, structural_pars=list(W=W_222), unvec=TRUE), A11_222)
  expect_equal(pick_Ami(p=2, M=2, d=2, params=theta_222s, m=2, i=2, structural_pars=list(W=W_222), unvec=TRUE), A22_222)

  expect_equal(pick_Ami(p=3, M=3, d=2, params=theta_332s, m=1, i=2, structural_pars=list(W=W_332), unvec=TRUE), A12_332)
  expect_equal(pick_Ami(p=3, M=3, d=2, params=theta_332s, m=2, i=3, structural_pars=list(W=W_332), unvec=TRUE), A23_332)
  expect_equal(pick_Ami(p=3, M=3, d=2, params=theta_332s, m=3, i=1, structural_pars=list(W=W_332), unvec=TRUE), A31_332)

  expect_equal(pick_Ami(p=1, M=2, d=3, params=theta_123s, m=1, i=1, structural_pars=list(W=W_123), unvec=TRUE), A11_123)
  expect_equal(pick_Ami(p=1, M=2, d=3, params=theta_123s, m=2, i=1, structural_pars=list(W=W_123), unvec=TRUE), A21_123)

  expect_equal(pick_Ami(p=1, M=2, d=3, params=theta_123ts, m=1, i=1, structural_pars=list(W=W_123), unvec=TRUE), A11_123) # SStMVAR
  expect_equal(pick_Ami(p=1, M=2, d=3, params=theta_123ts, m=2, i=1, structural_pars=list(W=W_123), unvec=TRUE), A21_123) # SStMVAR

  expect_equal(pick_Ami(p=1, M=c(1, 1), d=3, params=theta_123gss, m=1, i=1, structural_pars=list(W=W_123), unvec=TRUE), A11_123) # SG-StMVAR
  expect_equal(pick_Ami(p=1, M=c(1, 1), d=3, params=theta_123gss, m=2, i=1, structural_pars=list(W=W_123), unvec=TRUE), A21_123) # SG-StMVAR

  expect_equal(pick_Ami(p=2, M=1, d=3, params=theta_213s, m=1, i=1, structural_pars=list(W=W_213), unvec=TRUE), A11_213)
  expect_equal(pick_Ami(p=2, M=1, d=3, params=theta_213s, m=1, i=2, structural_pars=list(W=W_213), unvec=TRUE), A12_213)
})


test_that("pick_Am works correctly", {
  expect_equal(pick_Am(p=1, M=1, d=2, params=theta_112, m=1)[, , 1], A11_112)

  expect_equal(pick_Am(p=1, M=2, d=2, params=theta_122, m=1)[, , 1], A11_122)
  expect_equal(pick_Am(p=1, M=2, d=2, params=theta_122, m=2)[, , 1], A21_122)

  expect_equal(pick_Am(p=2, M=2, d=2, params=theta_222, m=1)[, , 2], A12_222)
  expect_equal(pick_Am(p=2, M=2, d=2, params=theta_222, m=2)[, , 1], A21_222)

  expect_equal(pick_Am(p=3, M=3, d=2, params=theta_332, m=1)[, , 3], A13_332)
  expect_equal(pick_Am(p=3, M=3, d=2, params=theta_332, m=2)[, , 2], A22_332)
  expect_equal(pick_Am(p=3, M=3, d=2, params=theta_332, m=3)[, , 1], A31_332)

  expect_equal(pick_Am(p=1, M=2, d=3, params=theta_123, m=1)[, , 1], A11_123)
  expect_equal(pick_Am(p=1, M=2, d=3, params=theta_123, m=2)[, , 1], A21_123)

  expect_equal(pick_Am(p=1, M=2, d=3, params=theta_123t, m=1)[, , 1], A11_123) # StMVAR
  expect_equal(pick_Am(p=1, M=2, d=3, params=theta_123t, m=2)[, , 1], A21_123) # StMVAR

  expect_equal(pick_Am(p=1, M=c(1, 1), d=3, params=theta_123gs, m=1)[, , 1], A11_123) # G-StMVAR
  expect_equal(pick_Am(p=1, M=c(1, 1), d=3, params=theta_123gs, m=2)[, , 1], A21_123) # G-StMVAR

  expect_equal(pick_Am(p=2, M=1, d=3, params=theta_213, m=1)[, , 1], A11_213)
  expect_equal(pick_Am(p=2, M=1, d=3, params=theta_213, m=1)[, , 2], A12_213)

  # Structural
  expect_equal(pick_Am(p=1, M=1, d=2, params=theta_112s, m=1, structural_pars=list(W=W_112))[, , 1], A11_112)

  expect_equal(pick_Am(p=1, M=2, d=2, params=theta_122s, m=1, structural_pars=list(W=W_122))[, , 1], A11_122)
  expect_equal(pick_Am(p=1, M=2, d=2, params=theta_122s, m=2, structural_pars=list(W=W_122))[, , 1], A21_122)

  expect_equal(pick_Am(p=2, M=2, d=2, params=theta_222s, m=1, structural_pars=list(W=W_222))[, , 2], A12_222)
  expect_equal(pick_Am(p=2, M=2, d=2, params=theta_222s, m=2, structural_pars=list(W=W_222))[, , 1], A21_222)

  expect_equal(pick_Am(p=2, M=2, d=2, params=theta_222ts, m=1, structural_pars=list(W=W_222))[, , 2], A12_222) # SStMVAR
  expect_equal(pick_Am(p=2, M=2, d=2, params=theta_222ts, m=2, structural_pars=list(W=W_222))[, , 1], A21_222) # SStMVAR

  expect_equal(pick_Am(p=2, M=c(1, 1), d=2, params=theta_222gss, m=1, structural_pars=list(W=W_222))[, , 2], A12_222) # SG-StMVAR
  expect_equal(pick_Am(p=2, M=c(1, 1), d=2, params=theta_222gss, m=2, structural_pars=list(W=W_222))[, , 1], A21_222) # SG-StMVAR

  expect_equal(pick_Am(p=3, M=3, d=2, params=theta_332s, m=1, structural_pars=list(W=W_332))[, , 3], A13_332)
  expect_equal(pick_Am(p=3, M=3, d=2, params=theta_332s, m=2, structural_pars=list(W=W_332))[, , 2], A22_332)
  expect_equal(pick_Am(p=3, M=3, d=2, params=theta_332s, m=3, structural_pars=list(W=W_332))[, , 1], A31_332)

  expect_equal(pick_Am(p=1, M=2, d=3, params=theta_123s, m=1, structural_pars=list(W=W_123))[, , 1], A11_123)
  expect_equal(pick_Am(p=1, M=2, d=3, params=theta_123s, m=2, structural_pars=list(W=W_123))[, , 1], A21_123)

  expect_equal(pick_Am(p=2, M=1, d=3, params=theta_213s, m=1, structural_pars=list(W=W_213))[, , 1], A11_213)
  expect_equal(pick_Am(p=2, M=1, d=3, params=theta_213s, m=1, structural_pars=list(W=W_213))[, , 2], A12_213)
})


test_that("pick_allA works correctly", {
  expect_equal(pick_allA(p=1, M=1, d=2, params=theta_112)[, , 1, 1], A11_112)

  expect_equal(pick_allA(p=1, M=2, d=2, params=theta_122)[, , 1, 1], A11_122)
  expect_equal(pick_allA(p=1, M=2, d=2, params=theta_122)[, , 1, 2], A21_122)

  expect_equal(pick_allA(p=2, M=2, d=2, params=theta_222)[, , 1, 1], A11_222)
  expect_equal(pick_allA(p=2, M=2, d=2, params=theta_222)[, , 2, 1], A12_222)
  expect_equal(pick_allA(p=2, M=2, d=2, params=theta_222)[, , 2, 2], A22_222)

  expect_equal(pick_allA(p=2, M=2, d=2, params=theta_222t)[, , 1, 1], A11_222) # StMVAR
  expect_equal(pick_allA(p=2, M=2, d=2, params=theta_222t)[, , 2, 1], A12_222) # StMVAR
  expect_equal(pick_allA(p=2, M=2, d=2, params=theta_222t)[, , 2, 2], A22_222) # StMVAR

  expect_equal(pick_allA(p=3, M=3, d=2, params=theta_332)[, , 3, 1], A13_332)
  expect_equal(pick_allA(p=3, M=3, d=2, params=theta_332)[, , 1, 2], A21_332)
  expect_equal(pick_allA(p=3, M=3, d=2, params=theta_332)[, , 2, 3], A32_332)
  expect_equal(pick_allA(p=3, M=3, d=2, params=theta_332)[, , 3, 3], A33_332)

  expect_equal(pick_allA(p=1, M=2, d=3, params=theta_123)[, , 1, 1], A11_123)
  expect_equal(pick_allA(p=1, M=2, d=3, params=theta_123)[, , 1, 2], A21_123)

  expect_equal(pick_allA(p=1, M=c(1, 1), d=3, params=theta_123gs)[, , 1, 1], A11_123) # G-StMVAR
  expect_equal(pick_allA(p=1, M=c(1, 1), d=3, params=theta_123gs)[, , 1, 2], A21_123) # G-StMVAR

  expect_equal(pick_allA(p=2, M=1, d=3, params=theta_213)[, , 1, 1], A11_213)
  expect_equal(pick_allA(p=2, M=1, d=3, params=theta_213)[, , 2, 1], A12_213)

  # Structural
  expect_equal(pick_allA(p=1, M=1, d=2, params=theta_112s, structural_pars=list(W=W_112))[, , 1, 1], A11_112)

  expect_equal(pick_allA(p=1, M=2, d=2, params=theta_122s, structural_pars=list(W=W_122))[, , 1, 1], A11_122)
  expect_equal(pick_allA(p=1, M=2, d=2, params=theta_122s, structural_pars=list(W=W_122))[, , 1, 2], A21_122)

  expect_equal(pick_allA(p=2, M=2, d=2, params=theta_222s, structural_pars=list(W=W_222))[, , 1, 1], A11_222)
  expect_equal(pick_allA(p=2, M=2, d=2, params=theta_222s, structural_pars=list(W=W_222))[, , 2, 1], A12_222)
  expect_equal(pick_allA(p=2, M=2, d=2, params=theta_222s, structural_pars=list(W=W_222))[, , 2, 2], A22_222)

  expect_equal(pick_allA(p=2, M=c(1, 1), d=2, params=theta_222gss, structural_pars=list(W=W_222))[, , 1, 1], A11_222) # SG-StMVAR
  expect_equal(pick_allA(p=2, M=c(1, 1), d=2, params=theta_222gss, structural_pars=list(W=W_222))[, , 2, 1], A12_222) # SG-StMVAR
  expect_equal(pick_allA(p=2, M=c(1, 1), d=2, params=theta_222gss, structural_pars=list(W=W_222))[, , 2, 2], A22_222) # SG-StMVAR

  expect_equal(pick_allA(p=3, M=3, d=2, params=theta_332s, structural_pars=list(W=W_332))[, , 3, 1], A13_332)
  expect_equal(pick_allA(p=3, M=3, d=2, params=theta_332s, structural_pars=list(W=W_332))[, , 1, 2], A21_332)
  expect_equal(pick_allA(p=3, M=3, d=2, params=theta_332s, structural_pars=list(W=W_332))[, , 2, 3], A32_332)
  expect_equal(pick_allA(p=3, M=3, d=2, params=theta_332s, structural_pars=list(W=W_332))[, , 3, 3], A33_332)

  expect_equal(pick_allA(p=1, M=2, d=3, params=theta_123s, structural_pars=list(W=W_123))[, , 1, 1], A11_123)
  expect_equal(pick_allA(p=1, M=2, d=3, params=theta_123s, structural_pars=list(W=W_123))[, , 1, 2], A21_123)

  expect_equal(pick_allA(p=1, M=c(1, 1), d=3, params=theta_123ts, structural_pars=list(W=W_123))[, , 1, 1], A11_123) # SStMVAR
  expect_equal(pick_allA(p=1, M=c(1, 1), d=3, params=theta_123ts, structural_pars=list(W=W_123))[, , 1, 2], A21_123) # SStMVAR

  expect_equal(pick_allA(p=2, M=1, d=3, params=theta_213s, structural_pars=list(W=W_213))[, , 1, 1], A11_213)
  expect_equal(pick_allA(p=2, M=1, d=3, params=theta_213s, structural_pars=list(W=W_213))[, , 2, 1], A12_213)
})


test_that("pick_phi0 works correctly", {
  expect_equal(pick_phi0(p=1, M=1, d=2, params=theta_112)[, 1], phi10_112)

  expect_equal(pick_phi0(p=1, M=2, d=2, params=theta_122)[, 1], phi10_122)
  expect_equal(pick_phi0(p=1, M=2, d=2, params=theta_122)[, 2], phi20_122)

  expect_equal(pick_phi0(p=2, M=2, d=2, params=theta_222)[, 1], phi10_222)
  expect_equal(pick_phi0(p=2, M=2, d=2, params=theta_222)[, 2], phi20_222)

  expect_equal(pick_phi0(p=2, M=c(1, 1), d=2, params=theta_222gs)[, 1], phi10_222) # G-StMVAR
  expect_equal(pick_phi0(p=2, M=c(1, 1), d=2, params=theta_222gs)[, 2], phi20_222) # G-StMVAR

  expect_equal(pick_phi0(p=3, M=3, d=2, params=theta_332)[, 1], phi10_332)
  expect_equal(pick_phi0(p=3, M=3, d=2, params=theta_332)[, 2], phi20_332)
  expect_equal(pick_phi0(p=3, M=3, d=2, params=theta_332)[, 3], phi30_332)

  expect_equal(pick_phi0(p=1, M=2, d=3, params=theta_123)[, 1], phi10_123)
  expect_equal(pick_phi0(p=1, M=2, d=3, params=theta_123)[, 2], phi20_123)

  expect_equal(pick_phi0(p=1, M=2, d=3, params=theta_123t)[, 1], phi10_123) # StMVAR
  expect_equal(pick_phi0(p=1, M=2, d=3, params=theta_123t)[, 2], phi20_123) # StMVAR

  expect_equal(pick_phi0(p=2, M=1, d=3, params=theta_213)[, 1], phi10_213)

  # Structural
  expect_equal(pick_phi0(p=1, M=1, d=2, params=theta_112s, structural_pars=list(W=W_112))[, 1], phi10_112)

  expect_equal(pick_phi0(p=1, M=2, d=2, params=theta_122s, structural_pars=list(W=W_122))[, 1], phi10_122)
  expect_equal(pick_phi0(p=1, M=2, d=2, params=theta_122s, structural_pars=list(W=W_122))[, 2], phi20_122)

  expect_equal(pick_phi0(p=2, M=2, d=2, params=theta_222s, structural_pars=list(W=W_222))[, 1], phi10_222)
  expect_equal(pick_phi0(p=2, M=2, d=2, params=theta_222s, structural_pars=list(W=W_222))[, 2], phi20_222)

  expect_equal(pick_phi0(p=2, M=2, d=2, params=theta_222ts, structural_pars=list(W=W_222))[, 1], phi10_222) # StMVAR
  expect_equal(pick_phi0(p=2, M=2, d=2, params=theta_222ts, structural_pars=list(W=W_222))[, 2], phi20_222) # StMVAR

  expect_equal(pick_phi0(p=3, M=3, d=2, params=theta_332s, structural_pars=list(W=W_332))[, 1], phi10_332)
  expect_equal(pick_phi0(p=3, M=3, d=2, params=theta_332s, structural_pars=list(W=W_332))[, 2], phi20_332)
  expect_equal(pick_phi0(p=3, M=3, d=2, params=theta_332s, structural_pars=list(W=W_332))[, 3], phi30_332)

  expect_equal(pick_phi0(p=1, M=2, d=3, params=theta_123s, structural_pars=list(W=W_123))[, 1], phi10_123)
  expect_equal(pick_phi0(p=1, M=2, d=3, params=theta_123s, structural_pars=list(W=W_123))[, 2], phi20_123)

  expect_equal(pick_phi0(p=1, M=c(1, 1), d=3, params=theta_123gss, structural_pars=list(W=W_123))[, 1], phi10_123) # G-StMVAR
  expect_equal(pick_phi0(p=1, M=c(1, 1), d=3, params=theta_123gss, structural_pars=list(W=W_123))[, 2], phi20_123) # G-StMVAR

  expect_equal(pick_phi0(p=2, M=1, d=3, params=theta_213s, structural_pars=list(W=W_213))[, 1], phi10_213)
})


test_that("pick_all_phi0_A works correctly", {
  expect_equal(pick_all_phi0_A(p=1, M=1, d=2, params=theta_112), as.matrix(c(phi10_112, vec(A11_112))))
  expect_equal(pick_all_phi0_A(p=1, M=2, d=2, params=theta_122), cbind(c(phi10_122, vec(A11_122)),
                                                                       c(phi20_122, vec(A21_122))))
  expect_equal(pick_all_phi0_A(p=2, M=2, d=2, params=theta_222), cbind(c(phi10_222, vec(A11_222), vec(A12_222)),
                                                                       c(phi20_222, vec(A21_222), vec(A22_222))))
  expect_equal(pick_all_phi0_A(p=2, M=2, d=2, params=theta_222t), cbind(c(phi10_222, vec(A11_222), vec(A12_222)),
                                                                        c(phi20_222, vec(A21_222), vec(A22_222)))) # StMVAR
  expect_equal(pick_all_phi0_A(p=2, M=c(1, 1), d=2, params=theta_222gs), cbind(c(phi10_222, vec(A11_222), vec(A12_222)),
                                                                               c(phi20_222, vec(A21_222), vec(A22_222)))) # G-StMVAR

  expect_equal(pick_all_phi0_A(p=3, M=3, d=2, params=theta_332), cbind(c(phi10_332, vec(A11_332), vec(A12_332), vec(A13_332)),
                                                                       c(phi20_332, vec(A21_332), vec(A22_332), vec(A23_332)),
                                                                       c(phi30_332, vec(A31_332), vec(A32_332), vec(A33_332))))
  expect_equal(pick_all_phi0_A(p=1, M=2, d=3, params=theta_123), cbind(c(phi10_123, vec(A11_123)),
                                                                       c(phi20_123, vec(A21_123))))
  expect_equal(pick_all_phi0_A(p=1, M=2, d=3, params=theta_123t), cbind(c(phi10_123, vec(A11_123)),
                                                                        c(phi20_123, vec(A21_123)))) # StMVAR
  expect_equal(pick_all_phi0_A(p=1, M=c(1, 1), d=3, params=theta_123gs), cbind(c(phi10_123, vec(A11_123)),
                                                                               c(phi20_123, vec(A21_123)))) # G-SMVAR

  expect_equal(pick_all_phi0_A(p=2, M=1, d=3, params=theta_213), as.matrix(c(phi10_213, vec(A11_213), vec(A12_213))))

  # Structural
  expect_equal(pick_all_phi0_A(p=1, M=1, d=2, params=theta_112s, structural_pars=list(W=W_112)),
               as.matrix(c(phi10_112, vec(A11_112))))
  expect_equal(pick_all_phi0_A(p=1, M=2, d=2, params=theta_122s, structural_pars=list(W=W_122)),
               cbind(c(phi10_122, vec(A11_122)), c(phi20_122, vec(A21_122))))
  expect_equal(pick_all_phi0_A(p=2, M=2, d=2, params=theta_222s, structural_pars=list(W=W_222)),
               cbind(c(phi10_222, vec(A11_222), vec(A12_222)), c(phi20_222, vec(A21_222), vec(A22_222))))
  expect_equal(pick_all_phi0_A(p=2, M=2, d=2, params=theta_222ts, structural_pars=list(W=W_222)),
               cbind(c(phi10_222, vec(A11_222), vec(A12_222)), c(phi20_222, vec(A21_222), vec(A22_222)))) # SStMVAR
  expect_equal(pick_all_phi0_A(p=2, M=c(1, 1), d=2, params=theta_222gss, structural_pars=list(W=W_222)),
               cbind(c(phi10_222, vec(A11_222), vec(A12_222)), c(phi20_222, vec(A21_222), vec(A22_222)))) # SG-StMVAR
  expect_equal(pick_all_phi0_A(p=3, M=3, d=2, params=theta_332s, structural_pars=list(W=W_332)),
               cbind(c(phi10_332, vec(A11_332), vec(A12_332), vec(A13_332)),
                     c(phi20_332, vec(A21_332), vec(A22_332), vec(A23_332)),
                     c(phi30_332, vec(A31_332), vec(A32_332), vec(A33_332))))
  expect_equal(pick_all_phi0_A(p=1, M=2, d=3, params=theta_123s, structural_pars=list(W=W_123)),
               cbind(c(phi10_123, vec(A11_123)), c(phi20_123, vec(A21_123))))
  expect_equal(pick_all_phi0_A(p=2, M=1, d=3, params=theta_213s, structural_pars=list(W=W_213)),
               as.matrix(c(phi10_213, vec(A11_213), vec(A12_213))))
})


test_that("pick_Omegas works correctly", {
  expect_equal(pick_Omegas(p=1, M=1, d=2, params=theta_112)[, , 1], Omega1_112)

  expect_equal(pick_Omegas(p=1, M=2, d=2, params=theta_122)[, , 1], Omega1_122)
  expect_equal(pick_Omegas(p=1, M=2, d=2, params=theta_122)[, , 2], Omega2_122)

  expect_equal(pick_Omegas(p=2, M=2, d=2, params=theta_222)[, , 1], Omega1_222)
  expect_equal(pick_Omegas(p=2, M=2, d=2, params=theta_222)[, , 2], Omega2_222)

  expect_equal(pick_Omegas(p=2, M=2, d=2, params=theta_222t)[, , 1], Omega1_222) # StMVAR
  expect_equal(pick_Omegas(p=2, M=2, d=2, params=theta_222t)[, , 2], Omega2_222) # StMVAR

  expect_equal(pick_Omegas(p=3, M=3, d=2, params=theta_332)[, , 1], Omega1_332)
  expect_equal(pick_Omegas(p=3, M=3, d=2, params=theta_332)[, , 2], Omega2_332)
  expect_equal(pick_Omegas(p=3, M=3, d=2, params=theta_332)[, , 3], Omega3_332)

  expect_equal(pick_Omegas(p=1, M=2, d=3, params=theta_123)[, , 1], Omega1_123)
  expect_equal(pick_Omegas(p=1, M=2, d=3, params=theta_123)[, , 2], Omega2_123)

  expect_equal(pick_Omegas(p=1, M=c(1, 1), d=3, params=theta_123gs)[, , 1], Omega1_123) # G-StMVAR
  expect_equal(pick_Omegas(p=1, M=c(1, 1), d=3, params=theta_123gs)[, , 2], Omega2_123) # G-StMVAR

  expect_equal(pick_Omegas(p=2, M=1, d=3, params=theta_213)[, , 1], Omega1_213)

  # Structural
  expect_equal(pick_Omegas(p=1, M=1, d=2, params=theta_112s, structural_pars=list(W=W_112))[, , 1], Omega1_112s)

  expect_equal(pick_Omegas(p=1, M=2, d=2, params=theta_122s, structural_pars=list(W=W_122))[, , 1], Omega1_122s)
  expect_equal(pick_Omegas(p=1, M=2, d=2, params=theta_122s, structural_pars=list(W=W_122))[, , 2], Omega2_122s)

  expect_equal(pick_Omegas(p=2, M=2, d=2, params=theta_222s, structural_pars=list(W=W_222))[, , 1], Omega1_222s)
  expect_equal(pick_Omegas(p=2, M=2, d=2, params=theta_222s, structural_pars=list(W=W_222))[, , 2], Omega2_222s)

  expect_equal(pick_Omegas(p=2, M=c(1, 1), d=2, params=theta_222gss, structural_pars=list(W=W_222))[, , 1], Omega1_222s) # SG-StMVAR
  expect_equal(pick_Omegas(p=2, M=c(1, 1), d=2, params=theta_222gss, structural_pars=list(W=W_222))[, , 2], Omega2_222s) # SG-StMVAR

  expect_equal(pick_Omegas(p=3, M=3, d=2, params=theta_332s, structural_pars=list(W=W_332))[, , 1], Omega1_332s)
  expect_equal(pick_Omegas(p=3, M=3, d=2, params=theta_332s, structural_pars=list(W=W_332))[, , 2], Omega2_332s)
  expect_equal(pick_Omegas(p=3, M=3, d=2, params=theta_332s, structural_pars=list(W=W_332))[, , 3], Omega3_332s)

  expect_equal(pick_Omegas(p=1, M=2, d=3, params=theta_123s, structural_pars=list(W=W_123))[, , 1], Omega1_123s)
  expect_equal(pick_Omegas(p=1, M=2, d=3, params=theta_123s, structural_pars=list(W=W_123))[, , 2], Omega2_123s)

  expect_equal(pick_Omegas(p=1, M=2, d=3, params=theta_123ts, structural_pars=list(W=W_123))[, , 1], Omega1_123s) # StMVAR
  expect_equal(pick_Omegas(p=1, M=2, d=3, params=theta_123ts, structural_pars=list(W=W_123))[, , 2], Omega2_123s) # StMVAR

  expect_equal(pick_Omegas(p=2, M=1, d=3, params=theta_213s, structural_pars=list(W=W_213))[, , 1], Omega1_213s)
})


test_that("pick_alphas works correctly", {
  expect_equal(pick_alphas(p=1, M=1, d=2, params=theta_112), 1)
  expect_equal(pick_alphas(p=1, M=2, d=2, params=theta_122), c(alpha1_122, 1-alpha1_122))
  expect_equal(pick_alphas(p=2, M=2, d=2, params=theta_222), c(alpha1_222, 1-alpha1_222))
  expect_equal(pick_alphas(p=2, M=2, d=2, params=theta_222t, model="StMVAR"), c(alpha1_222, 1-alpha1_222))
  expect_equal(pick_alphas(p=2, M=c(1, 1), d=2, params=theta_222gs, model="G-StMVAR"), c(alpha1_222, 1-alpha1_222))
  expect_equal(pick_alphas(p=3, M=3, d=2, params=theta_332), c(alpha1_332, alpha2_332, 1-alpha1_332-alpha2_332))
  expect_equal(pick_alphas(p=3, M=3, d=2, params=theta_332t, model="StMVAR"), c(alpha1_332, alpha2_332, 1-alpha1_332-alpha2_332))
  expect_equal(pick_alphas(p=3, M=c(1, 2), d=2, params=theta_332gs, model="G-StMVAR"), c(alpha1_332, alpha2_332, 1-alpha1_332-alpha2_332))
  expect_equal(pick_alphas(p=1, M=2, d=3, params=theta_123), c(alpha1_123, 1-alpha1_123))
  expect_equal(pick_alphas(p=1, M=2, d=3, params=theta_123t, model="StMVAR"), c(alpha1_123, 1-alpha1_123))
  expect_equal(pick_alphas(p=1, M=c(1, 1), d=3, params=theta_123gs, model="G-StMVAR"), c(alpha1_123, 1-alpha1_123))
  expect_equal(pick_alphas(p=2, M=1, d=3, params=theta_213), 1)

  # Structural
  expect_equal(pick_alphas(p=1, M=1, d=2, params=theta_112s), 1)
  expect_equal(pick_alphas(p=1, M=2, d=2, params=theta_122s), c(alpha1_122, 1-alpha1_122))
  expect_equal(pick_alphas(p=2, M=2, d=2, params=theta_222s), c(alpha1_222, 1-alpha1_222))
  expect_equal(pick_alphas(p=2, M=2, d=2, params=theta_222ts, model="StMVAR"), c(alpha1_222, 1-alpha1_222))
  expect_equal(pick_alphas(p=2, M=c(1, 1), d=2, params=theta_222gs, model="G-StMVAR"), c(alpha1_222, 1-alpha1_222))
  expect_equal(pick_alphas(p=3, M=3, d=2, params=theta_332s), c(alpha1_332, alpha2_332, 1-alpha1_332-alpha2_332))
  expect_equal(pick_alphas(p=3, M=3, d=2, params=theta_332ts, model="StMVAR"), c(alpha1_332, alpha2_332, 1-alpha1_332-alpha2_332))
  expect_equal(pick_alphas(p=3, M=c(2, 1), d=2, params=theta_332gss, model="G-StMVAR"), c(alpha1_332, alpha2_332, 1-alpha1_332-alpha2_332))
  expect_equal(pick_alphas(p=1, M=2, d=3, params=theta_123s), c(alpha1_123, 1-alpha1_123))
  expect_equal(pick_alphas(p=2, M=1, d=3, params=theta_213s), 1)
})

test_that("pick_df works correctly", {
  expect_equal(pick_df(M=2, params=theta_222t, model="StMVAR"), c(20, 25))
  expect_equal(pick_df(M=c(1, 1), params=theta_222gs, model="G-StMVAR"), c(25))
  expect_equal(pick_df(M=3, params=theta_332t, model="StMVAR"), c(10, 20, 30))
  expect_equal(pick_df(M=c(1, 2), params=theta_332gs, model="G-StMVAR"), c(20, 30))
  expect_equal(pick_df(M=2, params=theta_123), numeric(0))
  expect_equal(pick_df(M=2, params=theta_123t, model="StMVAR"), c(20, 25))
  expect_equal(pick_df(M=c(1, 1), params=theta_123gs, model="G-StMVAR"), c(25))

  # Structural
  expect_equal(pick_df(M=1, params=theta_112s), numeric(0))
  expect_equal(pick_df(M=2, params=theta_222ts, model="StMVAR"), c(10, 15))
  expect_equal(pick_df(M=c(1, 1), params=theta_222gs, model="G-StMVAR"), c(25))
  expect_equal(pick_df(M=3, params=theta_332ts, model="StMVAR"), c(10, 20, 30))
  expect_equal(pick_df(M=c(2, 1), params=theta_332gss, model="G-StMVAR"), c(30))
})

test_that("pick_W works correctly", {
  expect_null(pick_W(p=1, M=1, d=2, params=theta_112))

  # Structural
  expect_equal(pick_W(p=1, M=1, d=2, params=theta_112s, structural_pars=list(W_112)), W_112)
  expect_equal(pick_W(p=1, M=2, d=2, params=theta_122s, structural_pars=list(W_122)), W_122)
  expect_equal(pick_W(p=2, M=2, d=2, params=theta_222s, structural_pars=list(W_222)), W_222)
  expect_equal(pick_W(p=2, M=2, d=2, params=theta_222ts, structural_pars=list(W_222)), W_222) # SStMVAR
  expect_equal(pick_W(p=2, M=c(1, 1), d=2, params=theta_222gss, structural_pars=list(W_222)), W_222) # SG-StMVAR
  expect_equal(pick_W(p=3, M=3, d=2, params=theta_332s, structural_pars=list(W_332)), W_332)
  expect_equal(pick_W(p=3, M=3, d=2, params=theta_332ts, structural_pars=list(W_332)), W_332) # SStMVAR
  expect_equal(pick_W(p=3, M=c(2, 1), d=2, params=theta_332gss, structural_pars=list(W_332)), W_332) # SG-StMVAR
  expect_equal(pick_W(p=1, M=2, d=3, params=theta_123s, structural_pars=list(W_123)), W_123)
  expect_equal(pick_W(p=1, M=2, d=3, params=theta_123ts, structural_pars=list(W_123)), W_123) # SStMVAR
  expect_equal(pick_W(p=1, M=c(1, 1), d=3, params=theta_123gss, structural_pars=list(W_123)), W_123) # SG-StMVAR
  expect_equal(pick_W(p=2, M=1, d=3, params=theta_213s, structural_pars=list(W_213)), W_213)
})

test_that("pick_lambdas works correctly", {
  expect_equal(pick_lambdas(p=1, M=1, d=2, params=theta_112), numeric(0))

  # Structural
  expect_equal(pick_lambdas(p=1, M=1, d=2, params=theta_112s, structural_pars=list(W_112)), numeric(0))
  expect_equal(pick_lambdas(p=1, M=2, d=2, params=theta_122s, structural_pars=list(W_122)), lambdas_122)
  expect_equal(pick_lambdas(p=2, M=2, d=2, params=theta_222s, structural_pars=list(W_222)), lambdas_222)
  expect_equal(pick_lambdas(p=2, M=2, d=2, params=theta_222ts, structural_pars=list(W_222)), lambdas_222) # SStMVAR
  expect_equal(pick_lambdas(p=2, M=c(1, 1), d=2, params=theta_222gss, structural_pars=list(W_222)), lambdas_222) # GS-StMVAR
  expect_equal(pick_lambdas(p=3, M=3, d=2, params=theta_332s, structural_pars=list(W_332)), c(lambdas2_332, lambdas3_332))
  expect_equal(pick_lambdas(p=3, M=3, d=2, params=theta_332ts, structural_pars=list(W_332)), c(lambdas2_332, lambdas3_332)) # SStMVAR
  expect_equal(pick_lambdas(p=3, M=c(2, 1), d=2, params=theta_332gss, structural_pars=list(W_332)), c(lambdas2_332, lambdas3_332)) # SG-StMVAR
  expect_equal(pick_lambdas(p=1, M=2, d=3, params=theta_123s, structural_pars=list(W_123)), lambdas_123)
  expect_equal(pick_lambdas(p=1, M=2, d=3, params=theta_123ts, structural_pars=list(W_123)), lambdas_123) # SStMVAR
  expect_equal(pick_lambdas(p=1, M=c(1, 1), d=3, params=theta_123s, structural_pars=list(W_123)), lambdas_123) # SG-StMVAR
  expect_equal(pick_lambdas(p=2, M=1, d=3, params=theta_213s, structural_pars=list(W_213)), numeric(0))
})


## A(M)(p)_(p)(M)(d)
theta_213sWC <- c(phi10_213, vec(A11_213), vec(A12_213), Wvec(W_213)) # SGMVAR
theta_112sWC <- c(phi10_112, vec(A11_112), Wvec(W_112)) # SGMVAR
theta_213tsWC <- c(theta_213sWC, 10) # StMVAR
theta_112tsWC <- c(theta_112sWC, 15) # StMVAR

test_that("pick_regime works correctly", {
  expect_equal(pick_regime(p=1, M=1, d=2, params=theta_112, m=1), upsilon1_112)
  expect_equal(pick_regime(p=1, M=2, d=2, params=theta_122, m=1), upsilon1_122)
  expect_equal(pick_regime(p=1, M=2, d=2, params=theta_122, m=2), upsilon2_122)
  expect_equal(pick_regime(p=2, M=2, d=2, params=theta_222, m=1), upsilon1_222)
  expect_equal(pick_regime(p=2, M=2, d=2, params=theta_222, m=2), upsilon2_222)
  expect_equal(pick_regime(p=2, M=2, d=2, params=theta_222t, m=1, model="StMVAR", with_df=TRUE), c(upsilon1_222, 20))
  expect_equal(pick_regime(p=2, M=2, d=2, params=theta_222t, m=2, model="StMVAR", with_df=TRUE), c(upsilon2_222, 25))
  expect_equal(pick_regime(p=2, M=2, d=2, params=theta_222t, m=2, model="StMVAR", with_df=FALSE), upsilon2_222)
  expect_equal(pick_regime(p=2, M=c(1, 1), d=2, params=theta_222gs, m=1, model="G-StMVAR", with_df=TRUE), c(upsilon1_222))
  expect_equal(pick_regime(p=2, M=c(1, 1), d=2, params=theta_222gs, m=2, model="G-StMVAR", with_df=TRUE), c(upsilon2_222, 25))
  expect_equal(pick_regime(p=2, M=c(1, 1), d=2, params=theta_222gs, m=2, model="G-StMVAR", with_df=FALSE), c(upsilon2_222))
  expect_equal(pick_regime(p=3, M=3, d=2, params=theta_332, m=1), upsilon1_332)
  expect_equal(pick_regime(p=3, M=3, d=2, params=theta_332, m=2), upsilon2_332)
  expect_equal(pick_regime(p=3, M=3, d=2, params=theta_332, m=3), upsilon3_332)
  expect_equal(pick_regime(p=3, M=3, d=2, params=theta_332t, m=1, model="StMVAR"), c(upsilon1_332, 10))
  expect_equal(pick_regime(p=3, M=3, d=2, params=theta_332t, m=2, model="StMVAR"), c(upsilon2_332, 20))
  expect_equal(pick_regime(p=3, M=3, d=2, params=theta_332t, m=3, model="StMVAR"), c(upsilon3_332, 30))
  expect_equal(pick_regime(p=3, M=c(1, 2), d=2, params=theta_332gs, m=1, model="G-StMVAR"), upsilon1_332)
  expect_equal(pick_regime(p=3, M=c(1, 2), d=2, params=theta_332gs, m=2, model="G-StMVAR"), c(upsilon2_332, 20))
  expect_equal(pick_regime(p=3, M=c(1, 2), d=2, params=theta_332gs, m=3, model="G-StMVAR"), c(upsilon3_332, 30))
  expect_equal(pick_regime(p=1, M=2, d=3, params=theta_123, m=1), upsilon1_123)
  expect_equal(pick_regime(p=1, M=2, d=3, params=theta_123, m=2), upsilon2_123)
  expect_equal(pick_regime(p=1, M=2, d=3, params=theta_123t, m=1, model="StMVAR"), c(upsilon1_123, 20))
  expect_equal(pick_regime(p=1, M=2, d=3, params=theta_123t, m=2, model="StMVAR"), c(upsilon2_123, 25))
  expect_equal(pick_regime(p=1, M=c(1, 1), d=3, params=theta_123t, m=1, model="G-StMVAR"), c(upsilon1_123))
  expect_equal(pick_regime(p=1, M=c(1, 1), d=3, params=theta_123t, m=2, model="G-StMVAR"), c(upsilon2_123, 25))
  expect_equal(pick_regime(p=2, M=1, d=3, params=theta_213, m=1), upsilon1_213)

  # Structural
  expect_equal(pick_regime(p=1, M=1, d=2, params=theta_112sWC, m=1, structural_pars=list(W=W_112)),
               c(phi10_112, vec(A11_112)))
  expect_equal(pick_regime(p=1, M=1, d=2, params=theta_112tsWC, m=1, model="StMVAR", structural_pars=list(W=W_112)),
               c(phi10_112, vec(A11_112), 15))
  expect_equal(pick_regime(p=1, M=1, d=2, params=theta_112tsWC, m=1, model="StMVAR", structural_pars=list(W=W_112), with_df=FALSE),
               c(phi10_112, vec(A11_112)))
  expect_equal(pick_regime(p=1, M=2, d=2, params=theta_122s, m=1, structural_pars=list(W=W_122)),
               c(phi10_122, vec(A11_122)))
  expect_equal(pick_regime(p=1, M=2, d=2, params=theta_122s, m=2, structural_pars=list(W=W_122)),
               c(phi20_122, vec(A21_122)))
  expect_equal(pick_regime(p=2, M=2, d=2, params=theta_222s, m=1, structural_pars=list(W=W_222)),
               c(phi10_222, vec(A11_222), vec(A12_222)))
  expect_equal(pick_regime(p=2, M=2, d=2, params=theta_222s, m=2, structural_pars=list(W=W_222)),
               c(phi20_222, vec(A21_222), vec(A22_222)))
  expect_equal(pick_regime(p=2, M=2, d=2, params=theta_222ts, m=1, model="StMVAR", structural_pars=list(W=W_222)),
               c(phi10_222, vec(A11_222), vec(A12_222), 10))
  expect_equal(pick_regime(p=2, M=2, d=2, params=theta_222ts, m=2, model="StMVAR", structural_pars=list(W=W_222)),
               c(phi20_222, vec(A21_222), vec(A22_222), 15))
  expect_equal(pick_regime(p=2, M=c(1, 1), d=2, params=theta_222gss, m=1, model="G-StMVAR", structural_pars=list(W=W_222)),
               c(phi10_222, vec(A11_222), vec(A12_222)))
  expect_equal(pick_regime(p=2, M=c(1, 1), d=2, params=theta_222gss, m=2, model="G-StMVAR", structural_pars=list(W=W_222)),
               c(phi20_222, vec(A21_222), vec(A22_222), 15))
  expect_equal(pick_regime(p=3, M=3, d=2, params=theta_332s, m=1, structural_pars=list(W=W_332)),
               c(phi10_332, vec(A11_332), vec(A12_332), vec(A13_332)))
  expect_equal(pick_regime(p=3, M=3, d=2, params=theta_332s, m=2, structural_pars=list(W=W_332)),
               c(phi20_332, vec(A21_332), vec(A22_332), vec(A23_332)))
  expect_equal(pick_regime(p=3, M=3, d=2, params=theta_332s, m=3, structural_pars=list(W=W_332)),
               c(phi30_332, vec(A31_332), vec(A32_332), vec(A33_332)))
  expect_equal(pick_regime(p=3, M=3, d=2, params=theta_332ts, m=1, model="StMVAR", structural_pars=list(W=W_332)),
               c(phi10_332, vec(A11_332), vec(A12_332), vec(A13_332), 10))
  expect_equal(pick_regime(p=3, M=3, d=2, params=theta_332ts, m=2, model="StMVAR", structural_pars=list(W=W_332)),
               c(phi20_332, vec(A21_332), vec(A22_332), vec(A23_332), 20))
  expect_equal(pick_regime(p=3, M=3, d=2, params=theta_332ts, m=3, model="StMVAR", structural_pars=list(W=W_332)),
               c(phi30_332, vec(A31_332), vec(A32_332), vec(A33_332), 30))
  expect_equal(pick_regime(p=3, M=c(2, 1), d=2, params=theta_332gss, m=1, model="G-StMVAR", structural_pars=list(W=W_332)),
               c(phi10_332, vec(A11_332), vec(A12_332), vec(A13_332)))
  expect_equal(pick_regime(p=3, M=c(2, 1), d=2, params=theta_332gss, m=2, model="G-StMVAR", structural_pars=list(W=W_332)),
               c(phi20_332, vec(A21_332), vec(A22_332), vec(A23_332)))
  expect_equal(pick_regime(p=3, M=c(2, 1), d=2, params=theta_332gss, m=3, model="G-StMVAR", structural_pars=list(W=W_332)),
               c(phi30_332, vec(A31_332), vec(A32_332), vec(A33_332), 30))
  expect_equal(pick_regime(p=1, M=2, d=3, params=theta_123s, m=1, structural_pars=list(W=W_123)),
               c(phi10_123, vec(A11_123)))
  expect_equal(pick_regime(p=1, M=2, d=3, params=theta_123s, m=2, structural_pars=list(W=W_123)),
               c(phi20_123, vec(A21_123)))
  expect_equal(pick_regime(p=2, M=1, d=3, params=theta_213sWC, m=1, structural_pars=list(W=W_213)),
               c(phi10_213, vec(A11_213), vec(A12_213)))
  expect_equal(pick_regime(p=2, M=1, d=3, params=theta_213tsWC, m=1, model="StMVAR", structural_pars=list(W=W_213)),
               c(phi10_213, vec(A11_213), vec(A12_213), 10))

  # Constrained model
  expect_equal(pick_regime(p=1, M=1, d=2, params=theta_112c, model="GMVAR", m=1, constraints=C_112),
               c(phi10_112, vech(Omega1_112)))
  expect_equal(pick_regime(p=1, M=1, d=2, params=theta_112tc, model="StMVAR", m=1, constraints=C_112),
               c(phi10_112, vech(Omega1_112), 10))
  expect_equal(pick_regime(p=1, M=1, d=2, params=theta_112tcsWAR, model="StMVAR", m=1, constraints=C_112, structural_pars=list(W=W_112)),
               c(phi10_112, 10))
  expect_equal(pick_regime(p=2, M=2, d=2, params=theta_222c, model="GMVAR", m=1, constraints=C_222),
               c(phi10_222, vech(Omega1_222)))
  expect_equal(pick_regime(p=2, M=2, d=2, params=theta_222csL, model="GMVAR", m=1, structural_pars=list(W=W_222, C_lambda=C_lambda_222)),
               c(phi10_222, vec(A11_222), vec(A12_222)))
  expect_equal(pick_regime(p=2, M=c(1, 1), d=2, params=theta_222gscsL, model="G-StMVAR", m=2, structural_pars=list(W=W_222, C_lambda=C_lambda_222)),
               c(phi20_222, vec(A21_222), vec(A22_222), 20))
  expect_equal(pick_regime(p=2, M=2, d=2, params=theta_222csLAR, model="GMVAR", m=1, constraints=C_222,
                           structural_pars=list(W=W_222, C_lambda=C_lambda_222)), c(phi10_222))
  expect_equal(pick_regime(p=2, M=2, d=2, params=theta_222tcsLAR, model="StMVAR", m=2, constraints=C_222,
                           structural_pars=list(W=W_222, C_lambda=C_lambda_222)), c(phi20_222, 20))
  expect_equal(pick_regime(p=2, M=2, d=2, params=theta_222tcsLAR, model="StMVAR", m=2, constraints=C_222,
                           structural_pars=list(W=W_222, C_lambda=C_lambda_222), with_df=FALSE), c(phi20_222))
  expect_equal(pick_regime(p=1, M=2, d=3, params=theta_123c, model="GMVAR", m=2, constraints=C_123),
               c(phi20_123, vech(Omega2_123)))
  expect_equal(pick_regime(p=1, M=2, d=3, params=theta_123tc, model="StMVAR", m=2, constraints=C_123),
               c(phi20_123, vech(Omega2_123), 20))
  expect_equal(pick_regime(p=1, M=2, d=3, params=theta_123tc, model="StMVAR", m=2, constraints=C_123, with_df=FALSE),
               c(phi20_123, vech(Omega2_123)))
  expect_equal(pick_regime(p=1, M=2, d=3, params=theta_123tc, model="StMVAR", m=1, constraints=C_123, with_df=FALSE),
               c(phi10_123, vech(Omega1_123)))
  expect_equal(pick_regime(p=1, M=2, d=3, params=theta_123tcsL, model="StMVAR", m=1, constraints=C_123,
                           structural_pars=list(W=W_123, C_lambda=C_lambda_123)), c(phi10_123, 10))
  expect_equal(pick_regime(p=2, M=2, d=2, params=theta_222c_int, model="GMVAR", m=1, constraints=C_222, same_means=list(1:2)),
               c(vech(Omega1_222)))
  expect_equal(pick_regime(p=2, M=2, d=2, params=theta_222tc_int, model="StMVAR", m=2, constraints=C_222, same_means=list(1:2)),
               c(vech(Omega2_222), 20))
  expect_equal(pick_regime(p=3, M=3, d=2, params=theta_332c_int, model="GMVAR", m=3, constraints=C_332, same_means=list(1, 2:3)),
               c(vech(Omega3_332)))
  expect_equal(pick_regime(p=3, M=c(1, 2), d=2, params=theta_332gsc_int, model="G-StMVAR", m=1, constraints=C_332, same_means=list(1, 2:3)),
               c(vech(Omega1_332)))
  expect_equal(pick_regime(p=3, M=c(1, 2), d=2, params=theta_332gsc_int, model="G-StMVAR", m=2, constraints=C_332, same_means=list(1, 2:3)),
               c(vech(Omega2_332), 20))
  expect_equal(pick_regime(p=3, M=c(1, 2), d=2, params=theta_332gsc_int, model="G-StMVAR", m=3, constraints=C_332, same_means=list(1, 2:3)),
               c(vech(Omega3_332), 30))
  expect_equal(pick_regime(p=1, M=2, d=3, params=theta_123csLAR_int, model="GMVAR", m=2, constraints=C_123,
                           structural_pars=list(W=W_123, C_lambda=C_lambda_123), same_means=list(1:2)), numeric(0))
  expect_equal(pick_regime(p=1, M=2, d=3, params=theta_123tcsLAR_int, model="StMVAR", m=1, constraints=C_123,
                           structural_pars=list(W=W_123, C_lambda=C_lambda_123), same_means=list(1:2)), 10)
  expect_equal(pick_regime(p=1, M=2, d=3, params=theta_123tcsLAR_int, model="StMVAR", m=2, constraints=C_123,
                           structural_pars=list(W=W_123, C_lambda=C_lambda_123), same_means=list(1:2), with_df=FALSE), numeric(0))

})


params222 <- c(-11.904, 154.684, 1.314, 0.145, 0.094, 1.292, -0.389,
 -0.070, -0.109, -0.281, 0.920, -0.025, 4.839, 11.633, 124.983, 1.248,
  0.077, -0.040, 1.266, -0.272, -0.074, 0.034, -0.313, 5.855, 3.570,
  9.838, 0.740)
mod222 <- GSMVAR(d=2, p=2, M=2, params=params222, parametrization="mean")

mod222t <- GSMVAR(d=2, p=2, M=2, params=c(params222, 10, 20), model="StMVAR", parametrization="mean") # StMVAR
mod222gs <- GSMVAR(d=2, p=2, M=c(1, 1), params=c(params222, 20), model="G-StMVAR", parametrization="mean") # G-StMVAR


## A(M)(p)_(p)(M)(d)
rbind_diags <- function(p, M, d) {
  I <- diag(p*d^2)
  Reduce(rbind, replicate(M, I, simplify=FALSE))
}

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
mod222c <- GSMVAR(d=2, p=2, M=2, params=theta_222_c2, constraints=C_222_2)

# SGSMVAR models
C_112 <- rbind_diags(p=1, M=1, d=2)
A11_112 <- matrix(c(0.25, 0.06, 0.04, 0.34), nrow=2, byrow=FALSE) # Re-define as stationary
theta_112csWAR <- c(phi10_112, vec(A11_112), Wvec(W_112)) # SGMVAR W and AR
mod112csWAR <- GSMVAR(p=1, M=1, d=2, params=theta_112csWAR, structural_pars=list(W=W_112))

mod112tcsWAR <- GSMVAR(p=1, M=1, d=2, params=c(theta_112csWAR, 10), model="StMVAR", structural_pars=list(W=W_112)) # SStMVAR


C_222 <- rbind_diags(p=2, M=2, d=2)
C_lambda_222 <- matrix(c(1, 2), nrow=2)
theta_222csLAR <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(W_222), 0.2, alpha1_222) # SGMVAR lambdas and AR
mod222csLAR <- GSMVAR(p=2, M=2, d=2, params=theta_222csLAR, constraints=C_222,
                     structural_pars=list(W=W_222, C_lambda=C_lambda_222))

mod222gscsLAR <- GSMVAR(p=2, M=c(1, 1), d=2, params=c(theta_222csLAR, 20), model="G-StMVAR", constraints=C_222,
                        structural_pars=list(W=W_222, C_lambda=C_lambda_222)) # SG-StMVAR


test_that("get_boldA_eigens works correctly", {
  expect_equal(get_boldA_eigens(mod222)[,1], c(0.9917467, 0.9112338, 0.4566127, 0.2464068), tolerance=1e-5)
  expect_equal(get_boldA_eigens(mod222t)[2,], c(0.9112338, 0.9285837), tolerance=1e-5)
  expect_equal(get_boldA_eigens(mod222gs)[3,], c(0.4566127, 0.3125171), tolerance=1e-5)
  expect_equal(get_boldA_eigens(mod222c)[,2], c(0.9681610, 0.9569557, 0.3718390, 0.3030443), tolerance=1e-5)

  # SGSMVAR
  expect_equal(get_boldA_eigens(mod112csWAR)[,1], c(0.3615207, 0.2284793), tolerance=1e-5)
  expect_equal(get_boldA_eigens(mod112tcsWAR)[,1], c(0.3615207, 0.2284793), tolerance=1e-5)
  expect_equal(get_boldA_eigens(mod222csLAR)[,2], c(0.9819784, 0.9215689, 0.4260533, 0.2603994), tolerance=1e-5)
  expect_equal(get_boldA_eigens(mod222gscsLAR)[4,], c(.2603994, 0.2603994), tolerance=1e-5)
})

test_that("get_omega_eigens works correctly", {
  expect_equal(get_omega_eigens(mod222)[,1], c(4.8391595, 0.9198405), tolerance=1e-5)
  expect_equal(get_omega_eigens(mod222t)[2,], c(0.9198405, 3.7585944), tolerance=1e-5)
  expect_equal(get_omega_eigens(mod222gs)[1,], c(4.839159, 11.934406), tolerance=1e-5)
  expect_equal(get_omega_eigens(mod222c)[,2], c(11.611462, 3.608538), tolerance=1e-5)

  # SGSMVAR
  expect_equal(get_omega_eigens(mod112csWAR)[,1], c(5.2052628, 0.9247372), tolerance=1e-5)
  expect_equal(get_omega_eigens(mod112tcsWAR)[,1], c(5.2052628, 0.9247372), tolerance=1e-5)
  expect_equal(get_omega_eigens(mod222csLAR)[,2], c(2.0610439, 0.1882761), tolerance=1e-5)
  expect_equal(get_omega_eigens(mod222gscsLAR)[,1], c(5.3657061, 0.9039939), tolerance=1e-5)
})





















