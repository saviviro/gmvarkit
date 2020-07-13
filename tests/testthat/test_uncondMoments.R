context("unconditional moments")
library(gmvarkit)


## A(M)(p)_(p)(M)(d)

# p=1, M=1, d=2
phi10_112 <- c(1.03, 2.36)
A11_112 <- matrix(c(0.25, 0.06, 0.04, 0.34), nrow=2, byrow=FALSE)
Omega1_112 <- matrix(c(0.93, -0.15, -0.15, 5.20), nrow=2, byrow=FALSE)

theta_112 <- upsilon1_112 <- c(phi10_112, vec(A11_112), vech(Omega1_112))

W_112 <- t(chol(Omega1_112))
theta_112sWC <- c(phi10_112, vec(A11_112), Wvec(W_112)) # SGMVAR, W constrained by Cholesky

# p=1, M=2, d=2
phi10_122 <- c(1.03, 2.36)
A11_122 <- matrix(c(0.1, -0.06, -0.04, 0.1), nrow=2, byrow=FALSE)
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
A11_332 <- matrix(c(0.25, 0.06, 0.04, 0.34), nrow=2, byrow=FALSE)
A12_332 <- matrix(c(-0.29, -0.08, -0.05, -0.36), nrow=2, byrow=FALSE)
A13_332 <- matrix(c(0.1, 0.21, 0.12, 0.2), nrow=2, byrow=FALSE)
Omega1_332 <- matrix(c(0.93, -0.15, -0.15, 5.20), nrow=2, byrow=FALSE)

phi20_332 <- c(1.79, 3.00)
A21_332 <- matrix(c(0.20, 0.05, 0.05, 0.30), nrow=2, byrow=FALSE)
A22_332 <- matrix(c(-0.30, -0.10, -0.05, -0.40), nrow=2, byrow=FALSE)
A23_332 <- matrix(c(0.30, 0.10, 0.05, 0.40), nrow=2, byrow=FALSE)
Omega2_332 <- matrix(c(5.88, 3.56, 3.56, 9.80), nrow=2, byrow=FALSE)

phi30_332 <- c(1.79, 3.00)
A31_332 <- matrix(c(0.30, 0.03, 0.08, 0.33), nrow=2, byrow=FALSE)
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
A11_123 <- matrix(c(0.1, 0.21, 0.31, 0.12, 0.2, 0.32, 0.13, 0.23, 0.3), nrow=3, byrow=FALSE)
Omega1_123 <- matrix(c(1, 0.22, 0.33, 0.22, 2, 0.44, 0.33, 0.44, 3), nrow=3, byrow=FALSE)

phi20_123 <- c(1.11, 2.22, 3.33)
A21_123 <- matrix(c(-0.1, -0.21, -0.31, -0.12, -0.2, -0.32, -0.13, -0.23, -0.3), nrow=3, byrow=FALSE)
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
A11_213 <- matrix(c(0.1, 0.21, 0.31, 0.12, 0.2, 0.32, 0.13, 0.23, 0.3), nrow=3, byrow=FALSE)
A12_213 <- matrix(c(-0.1, -0.21, -0.31, -0.12, -0.2, -0.32, -0.13, -0.23, -0.3), nrow=3, byrow=FALSE)
Omega1_213 <- matrix(c(1, 0.22, 0.33, 0.22, 2, 0.44, 0.33, 0.44, 3), nrow=3, byrow=FALSE)

upsilon1_213 <- c(phi10_213, vec(A11_213), vec(A12_213), vech(Omega1_213))
theta_213 <- upsilon1_213

W_213 <- t(chol(Omega1_213))
theta_213s <- c(phi10_213, vec(A11_213), vec(A12_213), vec(W_213)) # SGMVAR

theta_213sWC <- c(phi10_213, vec(A11_213), vec(A12_213), Wvec(W_213)) # SGMVAR W constrained

calc_mu <- function(p, M, d, params, constraints=NULL, structural_pars=NULL) {
  params <- reform_constrained_pars(p=p, M=M, d=d, params=params, constraints=constraints, structural_pars=structural_pars)
  all_A <- pick_allA(p=p, M=M, d=d, params=params, structural_pars=structural_pars)
  all_phi0 <- pick_phi0(p=p, M=M, d=d, params=params, structural_pars=structural_pars)
  vapply(1:M, function(m) solve(diag(d) - rowSums(all_A[, , , m, drop=FALSE], dims=2), all_phi0[,m]), numeric(d))
}

theta_112_mu <- change_parametrization(p=1, M=1, d=2, params=theta_112, change_to="mean")
theta_122_mu <- change_parametrization(p=1, M=2, d=2, params=theta_122, change_to="mean")
theta_222_mu <- change_parametrization(p=2, M=2, d=2, params=theta_222, change_to="mean")
theta_332_mu <- change_parametrization(p=3, M=3, d=2, params=theta_332, change_to="mean")
theta_123_mu <- change_parametrization(p=1, M=2, d=3, params=theta_123, change_to="mean")
theta_213_mu <- change_parametrization(p=2, M=1, d=3, params=theta_213, change_to="mean")

theta_112sWC_mu <- change_parametrization(p=1, M=1, d=2, params=theta_112sWC, structural_pars=list(W=W_112), change_to="mean")
theta_222s_mu <- change_parametrization(p=2, M=2, d=2, params=theta_222s, structural_pars=list(W=W_222), change_to="mean")

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

# p=2, M=2, d=2
C_222 <- rbind_diags(p=2, M=2, d=2)
theta_222c <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vech(Omega1_222), vech(Omega2_222), alpha1_222)

C_lambda_222 <- matrix(c(1, 2), nrow=2)
theta_222csL <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A21_222),
                  vec(A22_222), vec(W_222), 0.2, alpha1_222) # SGMVAR lambdas
theta_222csLAR <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(W_222), 0.2, alpha1_222) # SGMVAR lambdas and AR

# p=1, M=2, d=3
C_123 <- rbind_diags(p=1, M=2, d=3)
theta_123c <- c(phi10_123, phi20_123, vec(A11_123), vech(Omega1_123), vech(Omega2_123), alpha1_123)

C_lambda_123 <- matrix(c(1, 1, 0, 0, 0, 1), nrow=3, byrow=FALSE)
theta_123csL <- c(phi10_123, phi20_123, vec(A11_123), vec(A21_123), vec(W_123), 1, 2, alpha1_123) # SGMVAR lambdas
theta_123csLAR <- c(phi10_123, phi20_123, vec(A11_123), vec(W_123), 1, 2, alpha1_123) # SGMVAR lambdas and AR

# p=2, M=2, d=2, constraint AR-parameters to be the same for all regimes
# and constraint the of-diagonal elements of AR-matrices to be zero.
mat0 <- matrix(c(1, rep(0, 10), 1, rep(0, 8), 1, rep(0, 10), 1), nrow=2*2^2, byrow=FALSE)
C_222_2 <- rbind(mat0, mat0)
A21_222_c2 <- A11_222_c2 <- matrix(c(0.26, 0, 0, 0.34), nrow=2, byrow=FALSE)
A22_222_c2 <- A12_222_c2 <- matrix(c(-0.29, 0, 0, -0.36), nrow=2, byrow=FALSE)
phi10_222_c2 <- c(-0.11, 2.83)
phi20_222_c2 <- c(0.36, 3.19)
Omega1_222_c2 <- matrix(c(0.98, -0.33, -0.33, 5.24), nrow=2, byrow=FALSE)
Omega2_222_c2 <- matrix(c(5.60, 3.46, 3.46, 9.62), nrow=2, byrow=FALSE)
alpha1_222_c2 <- 0.35
theta_222_c2 <- c(phi10_222_c2, phi20_222_c2, 0.26, 0.34, -0.29, -0.36, vech(Omega1_222_c2),
                  vech(Omega2_222_c2), alpha1_222_c2)

WL_222c2 <- diag_Omegas(Omega1_222_c2, Omega2_222_c2)
W_222c2 <- matrix(WL_222c2[1:(2^2)], nrow=2, byrow=FALSE)
lambdas_222c2 <- WL_222c2[(2^2 + 1):length(WL_222c2)]
theta_222_c2s <- c(phi10_222_c2, phi20_222_c2, 1.26, 1.34, -0.29, -0.36, vec(W_222c2), lambdas_222c2, alpha1_222_c2) # SGMVAR AR

theta_112c_mu <- change_parametrization(p=1, M=1, d=2, params=theta_112c, constraints=C_112, change_to="mean")
theta_222c_mu <- change_parametrization(p=2, M=2, d=2, params=theta_222c, constraints=C_222, change_to="mean")
theta_222c_mu2 <- change_parametrization(p=2, M=2, d=2, params=theta_222_c2, constraints=C_222_2, change_to="mean")
theta_123c_mu <- change_parametrization(p=1, M=2, d=3, params=theta_123c, constraints=C_123, change_to="mean")

theta_112c_mu_exp <- reform_constrained_pars(p=1, M=1, d=2, params=theta_112c_mu, constraints=C_112)
theta_222c_mu_exp <- reform_constrained_pars(p=2, M=2, d=2, params=theta_222c_mu, constraints=C_222)
theta_222c_mu2_exp <- reform_constrained_pars(p=2, M=2, d=2, params=theta_222c_mu2, constraints=C_222_2)
theta_123c_mu_exp <- reform_constrained_pars(p=1, M=2, d=3, params=theta_123c_mu, constraints=C_123)

theta_112cs_mu <- change_parametrization(p=1, M=1, d=2, params=theta_112csWAR, constraints=C_112,
                                         structural_pars=list(W=W_112), change_to="mean")
theta_222cs_mu <- change_parametrization(p=2, M=2, d=2, params=theta_222csLAR, constraints=C_222,
                                         structural_pars=list(W=W_222, C_lambda=C_lambda_222), change_to="mean")
theta_112cs_mu_exp <- reform_constrained_pars(p=1, M=1, d=2, params=theta_112cs_mu, constraints=C_112, structural_pars=list(W=W_112))
theta_222cs_mu_exp <- reform_constrained_pars(p=2, M=2, d=2, params=theta_222cs_mu, constraints=C_222,
                                              structural_pars=list(W=W_222, C_lambda_222))


test_that("get_regime_means_int works correctly", {
  expect_equal(get_regime_means_int(p=1, M=1, d=2, params=theta_112, parametrization="intercept", constraints=NULL),
               pick_phi0(p=1, M=1, d=2, params=theta_112_mu))
  expect_equal(get_regime_means_int(p=1, M=1, d=2, params=theta_112_mu, parametrization="mean", constraints=NULL),
               calc_mu(p=1, M=1, d=2, params=theta_112))
  expect_equal(get_regime_means_int(p=1, M=2, d=2, params=theta_122, parametrization="intercept", constraints=NULL),
               pick_phi0(p=1, M=2, d=2, params=theta_122_mu))
  expect_equal(get_regime_means_int(p=1, M=2, d=2, params=theta_122_mu, parametrization="mean", constraints=NULL),
              calc_mu(p=1, M=2, d=2, params=theta_122))
  expect_equal(get_regime_means_int(p=2, M=2, d=2, params=theta_222, parametrization="intercept", constraints=NULL),
               pick_phi0(p=2, M=2, d=2, params=theta_222_mu))
  expect_equal(get_regime_means_int(p=2, M=2, d=2, params=theta_222_mu, parametrization="mean", constraints=NULL),
               calc_mu(p=2, M=2, d=2, params=theta_222))
  expect_equal(get_regime_means_int(p=3, M=3, d=2, params=theta_332, parametrization="intercept", constraints=NULL),
               pick_phi0(p=3, M=3, d=2, params=theta_332_mu))
  expect_equal(get_regime_means_int(p=3, M=3, d=2, params=theta_332_mu, parametrization="mean", constraints=NULL),
               calc_mu(p=3, M=3, d=2, params=theta_332))
  expect_equal(get_regime_means_int(p=1, M=2, d=3, params=theta_123, parametrization="intercept", constraints=NULL),
               pick_phi0(p=1, M=2, d=3, params=theta_123_mu))
  expect_equal(get_regime_means_int(p=1, M=2, d=3, params=theta_123_mu, parametrization="mean", constraints=NULL),
               calc_mu(p=1, M=2, d=3, params=theta_123))
  expect_equal(get_regime_means_int(p=2, M=1, d=3, params=theta_213, parametrization="intercept", constraints=NULL),
               pick_phi0(p=2, M=1, d=3, params=theta_213_mu))
  expect_equal(get_regime_means_int(p=2, M=1, d=3, params=theta_213_mu, parametrization="mean", constraints=NULL),
               calc_mu(p=2, M=1, d=3, params=theta_213))

  expect_equal(get_regime_means_int(p=1, M=1, d=2, params=theta_112c, parametrization="intercept", constraints=C_112),
               pick_phi0(p=1, M=1, d=2, params=theta_112c_mu_exp))
  expect_equal(get_regime_means_int(p=1, M=1, d=2, params=theta_112c_mu, parametrization="mean", constraints=C_112),
               calc_mu(p=1, M=1, d=2, params=theta_112c, constraints=C_112))
  expect_equal(get_regime_means_int(p=2, M=2, d=2, params=theta_222c, parametrization="intercept", constraints=C_222),
               pick_phi0(p=2, M=2, d=2, params=theta_222c_mu_exp))
  expect_equal(get_regime_means_int(p=2, M=2, d=2, params=theta_222c_mu, parametrization="mean", constraints=C_222),
               calc_mu(p=2, M=2, d=2, params=theta_222c, constraints=C_222))
  expect_equal(get_regime_means_int(p=2, M=2, d=2, params=theta_222_c2, parametrization="intercept", constraints=C_222_2),
               pick_phi0(p=2, M=2, d=2, params=theta_222c_mu2_exp))
  expect_equal(get_regime_means_int(p=2, M=2, d=2, params=theta_222c_mu2, parametrization="mean", constraints=C_222_2),
               calc_mu(p=2, M=2, d=2, params=theta_222_c2, constraints=C_222_2))
  expect_equal(get_regime_means_int(p=1, M=2, d=3, params=theta_123c, parametrization="intercept", constraints=C_123),
               pick_phi0(p=1, M=2, d=3, params=theta_123c_mu_exp))
  expect_equal(get_regime_means_int(p=1, M=2, d=3, params=theta_123c_mu, parametrization="mean", constraints=C_123),
               calc_mu(p=1, M=2, d=3, params=theta_123c, constraints=C_123))

  # SGMVAR
  expect_equal(get_regime_means_int(p=1, M=1, d=2, params=theta_112sWC, parametrization="intercept", constraints=NULL,
                                    structural_pars=list(W=W_112)),
               pick_phi0(p=1, M=1, d=2, params=theta_112sWC_mu, structural_pars=list(W=W_112)))
  expect_equal(get_regime_means_int(p=1, M=1, d=2, params=theta_112sWC_mu, parametrization="mean", constraints=NULL,
                                    structural_pars=list(W=W_112)),
               calc_mu(p=1, M=1, d=2, params=theta_112sWC, structural_pars=list(W=W_112)))
  expect_equal(get_regime_means_int(p=2, M=2, d=2, params=theta_222s, parametrization="intercept", constraints=NULL,
                                    structural_pars=list(W=W_222)),
               pick_phi0(p=2, M=2, d=2, params=theta_222s_mu, structural_pars=list(W=W_222)))
  expect_equal(get_regime_means_int(p=2, M=2, d=2, params=theta_222s_mu, parametrization="mean", constraints=NULL,
                                    structural_pars=list(W=W_222)),
               calc_mu(p=2, M=2, d=2, params=theta_222s, structural_pars=list(W=W_222)))

  expect_equal(get_regime_means_int(p=1, M=1, d=2, params=theta_112csWAR, parametrization="intercept", constraints=C_112,
                                    structural_pars=list(W=W_112)),
               pick_phi0(p=1, M=1, d=2, params=theta_112cs_mu_exp, structural_pars=list(W=W_112)))
  expect_equal(get_regime_means_int(p=1, M=1, d=2, params=theta_112cs_mu, parametrization="mean", constraints=C_112,
                                    structural_pars=list(W=W_112)),
               calc_mu(p=1, M=1, d=2, params=theta_112csWAR, constraints=C_112, structural_pars=list(W=W_112)))
  expect_equal(get_regime_means_int(p=2, M=2, d=2, params=theta_222csLAR, parametrization="intercept", constraints=C_222,
                                    structural_pars=list(W=W_222, C_lambda=C_lambda_222)),
               pick_phi0(p=2, M=2, d=2, params=theta_222cs_mu_exp, structural_pars=list(W=W_222, C_lambda=C_lambda_222)))
  expect_equal(get_regime_means_int(p=2, M=2, d=2, params=theta_222cs_mu, parametrization="mean", constraints=C_222,
                                    structural_pars=list(W=W_222, C_lambda=C_lambda_222)),
               calc_mu(p=2, M=2, d=2, params=theta_222csLAR, constraints=C_222, structural_pars=list(W=W_222, C_lambda=C_lambda_222)))
})


test_that("get_regime_autocovs_int works correctly", {
  expect_equal(get_regime_autocovs_int(p=1, M=1, d=2, params=theta_112, constraints=NULL)[, 2, 2, 1],
               c(0.2201706, 1.9959185), tolerance=1e-6)
  expect_equal(get_regime_autocovs_int(p=1, M=2, d=2, params=theta_122, constraints=NULL)[, 1, 1, 2],
               c(5.926843, 3.528684), tolerance=1e-5)
  expect_equal(get_regime_autocovs_int(p=2, M=2, d=2, params=theta_222, constraints=NULL)[, 1, 3, 2],
               c(36.50580, 13.77247), tolerance=1e-5)
  expect_equal(get_regime_autocovs_int(p=3, M=3, d=2, params=theta_332, constraints=NULL)[2, , 4, 3],
               c(2.89655, 4.28663), tolerance=1e-6)
  expect_equal(get_regime_autocovs_int(p=1, M=2, d=3, params=theta_123, constraints=NULL)[1, , 2, 2],
               c(-0.3194244, -0.5931677, -0.8712061), tolerance=1e-6)
  expect_equal(get_regime_autocovs_int(p=2, M=1, d=3, params=theta_213, constraints=NULL)[, 3, 3, 1],
               c(-0.582898, -1.038395, -1.413783), tolerance=1e-6)
  expect_equal(get_regime_autocovs_int(p=2, M=2, d=2, params=theta_222c, constraints=C_222)[, 1, 3, 2],
               c(102.509805, 3.535469), tolerance=1e-6)

  # SGMVAR
  expect_equal(get_regime_autocovs_int(p=1, M=1, d=2, params=theta_112sWC, constraints=NULL, structural_pars=list(W=W_112))[, 2, 2, 1],
               c(0.2201706, 1.9959185), tolerance=1e-6)
  expect_equal(get_regime_autocovs_int(p=1, M=2, d=2, params=theta_122s, constraints=NULL, structural_pars=list(W=W_122))[, 1, 1, 2],
               c(5.926843, 3.528684), tolerance=1e-5)
  expect_equal(get_regime_autocovs_int(p=2, M=2, d=2, params=theta_222s, constraints=NULL, structural_pars=list(W=W_222))[, 1, 3, 2],
               c(36.50580, 13.77247), tolerance=1e-5)
  expect_equal(get_regime_autocovs_int(p=3, M=3, d=2, params=theta_332sWC, constraints=NULL, structural_pars=list(W=W_332))[2, , 4, 3],
               c(1.396341, 2.174643), tolerance=1e-5)
  expect_equal(get_regime_autocovs_int(p=1, M=2, d=3, params=theta_123s, constraints=NULL, structural_pars=list(W=W_123))[1, , 2, 2],
               c(-0.3194244, -0.5931677, -0.8712061), tolerance=1e-6)
  expect_equal(get_regime_autocovs_int(p=2, M=1, d=3, params=theta_213sWC, constraints=NULL, structural_pars=list(W=W_213))[, 3, 3, 1],
               c(-0.582898, -1.038395, -1.413783), tolerance=1e-6)
  expect_equal(get_regime_autocovs_int(p=1, M=2, d=3, params=theta_123csLAR, constraints=C_123, structural_pars=list(W=W_123, C_lambda=C_lambda_123))[, 1, 2, 1],
               c(0.2987136, 0.5610925, 0.8093760), tolerance=1e-6)
})


test_that("uncond_moments_int works correctly", {
  expect_equal(uncond_moments_int(p=1, M=1, d=2, params=theta_112, parametrization="intercept", constraints=NULL)$uncond_mean,
               c(1.571661, 3.718636), tolerance=1e-6)
  expect_equal(uncond_moments_int(p=1, M=1, d=2, params=theta_112, parametrization="intercept", constraints=NULL)$autocors[, 1, 1],
               c(1.00000000, -0.02484575), tolerance=1e-6)
  expect_equal(uncond_moments_int(p=1, M=2, d=2, params=theta_122, parametrization="intercept", constraints=NULL)$uncond_mean,
               c(1.544567, 2.967251), tolerance=1e-6)
  expect_equal(uncond_moments_int(p=1, M=2, d=2, params=theta_122, parametrization="intercept", constraints=NULL)$autocors[, 2, 2],
               c(0.002322773, 0.095293066), tolerance=1e-6)
  expect_equal(uncond_moments_int(p=2, M=2, d=2, params=theta_222, parametrization="intercept", constraints=NULL)$uncond_mean,
               c(9.4270, 58.7715), tolerance=1e-4)
  expect_equal(uncond_moments_int(p=2, M=2, d=2, params=theta_222, parametrization="intercept", constraints=NULL)$autocors[1, , 3],
               c(0.9623361, -0.8532178), tolerance=1e-6)
  expect_equal(uncond_moments_int(p=3, M=3, d=2, params=theta_332, parametrization="intercept", constraints=NULL)$uncond_mean,
               c(2.108140, 3.872402), tolerance=1e-6)
  expect_equal(uncond_moments_int(p=3, M=3, d=2, params=theta_332, parametrization="intercept", constraints=NULL)$autocors[2, , 2],
               c(0.08213133, 0.17093255), tolerance=1e-6)
  expect_equal(uncond_moments_int(p=1, M=2, d=3, params=theta_123, parametrization="intercept", constraints=NULL)$uncond_mean,
               c(2.263278, 4.278151, 6.265532), tolerance=1e-6)
  expect_equal(uncond_moments_int(p=1, M=2, d=3, params=theta_123, parametrization="intercept", constraints=NULL)$autocors[, 3, 1],
               c(0.7763858, 0.8178659, 1.0000000), tolerance=1e-6)
  expect_equal(uncond_moments_int(p=2, M=1, d=3, params=theta_213, parametrization="intercept", constraints=NULL)$uncond_mean,
               c(1.1, 2.2, 3.3), tolerance=1e-1)
  expect_equal(uncond_moments_int(p=2, M=1, d=3, params=theta_213, parametrization="intercept", constraints=NULL)$autocors[3, , 3],
               c(-0.2324331, -0.2750250, -0.2901025), tolerance=1e-6)
  expect_equal(uncond_moments_int(p=2, M=2, d=2, params=theta_222c, parametrization="intercept", constraints=C_222)$uncond_mean,
               c(4.24, 133.92), tolerance=1e-2)
  expect_equal(uncond_moments_int(p=2, M=2, d=2, params=theta_222c, parametrization="intercept", constraints=C_222)$autocors[, 2, 1],
               c(0.1983422, 1.0000000), tolerance=1e-6)

  # SGMVAR
  expect_equal(uncond_moments_int(p=1, M=1, d=2, params=theta_112sWC, parametrization="intercept", constraints=NULL, structural_pars=list(W=W_112))$uncond_mean,
               c(1.571661, 3.718636), tolerance=1e-6)
  expect_equal(uncond_moments_int(p=1, M=1, d=2, params=theta_112sWC, parametrization="intercept", constraints=NULL, structural_pars=list(W=W_112))$autocors[, 1, 1],
               c(1.00000000, -0.02484575), tolerance=1e-6)
  expect_equal(uncond_moments_int(p=1, M=2, d=2, params=theta_122s, parametrization="intercept", constraints=NULL, structural_pars=list(W=W_122))$uncond_mean,
               c(1.544567, 2.967251), tolerance=1e-6)
  expect_equal(uncond_moments_int(p=1, M=2, d=2, params=theta_122s, parametrization="intercept", constraints=NULL, structural_pars=list(W=W_122))$autocors[, 2, 2],
               c(0.002322773, 0.095293066), tolerance=1e-6)
  expect_equal(uncond_moments_int(p=2, M=2, d=2, params=theta_222s, parametrization="intercept", constraints=NULL, structural_pars=list(W=W_222))$uncond_mean,
               c(9.4270, 58.7715), tolerance=1e-4)
  expect_equal(uncond_moments_int(p=2, M=2, d=2, params=theta_222s, parametrization="intercept", constraints=NULL, structural_pars=list(W=W_222))$autocors[1, , 3],
               c(0.9623361, -0.8532178), tolerance=1e-6)
  expect_equal(uncond_moments_int(p=3, M=3, d=2, params=theta_332s, parametrization="intercept", constraints=NULL, structural_pars=list(W=W_332))$uncond_mean,
               c(2.108140, 3.872402), tolerance=1e-6)
  expect_equal(uncond_moments_int(p=3, M=3, d=2, params=theta_332s, parametrization="intercept", constraints=NULL, structural_pars=list(W=W_332))$autocors[2, , 2],
               c(0.4336740, 0.4434898), tolerance=1e-6)
  expect_equal(uncond_moments_int(p=1, M=2, d=3, params=theta_123s, parametrization="intercept", constraints=NULL, structural_pars=list(W=W_123))$uncond_mean,
               c(2.263278, 4.278151, 6.265532), tolerance=1e-6)
  expect_equal(uncond_moments_int(p=1, M=2, d=3, params=theta_123s, parametrization="intercept", constraints=NULL, structural_pars=list(W=W_123))$autocors[, 3, 1],
               c(0.7763858, 0.8178659, 1.0000000), tolerance=1e-6)
  expect_equal(uncond_moments_int(p=2, M=1, d=3, params=theta_213sWC, parametrization="intercept", constraints=NULL, structural_pars=list(W=W_213))$uncond_mean,
               c(1.1, 2.2, 3.3), tolerance=1e-1)
  expect_equal(uncond_moments_int(p=2, M=1, d=3, params=theta_213sWC, parametrization="intercept", constraints=NULL, structural_pars=list(W=W_213))$autocors[3, , 3],
               c(-0.2324331, -0.2750250, -0.2901025), tolerance=1e-6)
  expect_equal(uncond_moments_int(p=2, M=2, d=2, params=theta_222csLAR, parametrization="intercept", constraints=C_222, structural_pars=list(W=W_222, C_lambda=C_lambda_222))$uncond_mean,
               c(4.24, 133.92), tolerance=1e-2)
  expect_equal(uncond_moments_int(p=2, M=2, d=2, params=theta_222csLAR, parametrization="intercept", constraints=C_222, structural_pars=list(W=W_222, C_lambda=C_lambda_222))$autocors[, 2, 1],
               c(0.2277718, 1.0000000), tolerance=1e-6)
})


test_that("non_int uncond moment functions work", {
  mod122 <- GMVAR(p=1, M=2, d=2, params=theta_122)
  mod112csWAR <- GMVAR(p=1, M=1, d=2, params=theta_112csWAR, structural_pars=list(W=W_112))
  mod222csLAR <- GMVAR(p=2, M=2, d=2, params=theta_222csLAR, constraints=C_222,
                       structural_pars=list(W=W_222, C_lambda=C_lambda_222))

  unc122 <- uncond_moments(mod122)
  unc112csWAR <- uncond_moments(mod112csWAR)
  unc222csLAR <- uncond_moments(mod222csLAR)

  expect_equal(unc122, uncond_moments_int(p=1, M=2, d=2, params=theta_122), tolerance=1e-6)
  expect_equal(unc112csWAR, uncond_moments_int(p=1, M=1, d=2, params=theta_112csWAR,
                                          structural_pars=list(W=W_112)), tolerance=1e-6)
  expect_equal(unc222csLAR$autocovs[1, 2, ], c(27.91380, 27.59758, 27.16074), tolerance=1e-4)
  expect_equal(unc222csLAR$autocors[2, , 1], c(0.2277718, 1.0000000), tolerance=1e-4)
  expect_equal(unc222csLAR$uncond_mean, c(4.24, 133.92), tolerance=1e-4)

  reg_means122 <- get_regime_means(mod122)
  reg_means112csWAR <- get_regime_means(mod112csWAR)
  reg_means222csLAR <- get_regime_means(mod222csLAR)

  expect_equal(reg_means122[2, ], c(2.553492, 3.210253), tolerance=1e-5)
  expect_equal(reg_means112csWAR[, 1], c(1.571661, 3.718636), tolerance=1e-5)
  expect_equal(reg_means222csLAR[, 2], c(9.666667, 140.333333), tolerance=1e-5)

  reg_autocovs122 <- get_regime_autocovs(mod122)
  reg_autocovs112csWAR <- get_regime_autocovs(mod112csWAR)
  reg_autocovs222csLAR <- get_regime_autocovs(mod222csLAR)

  expect_equal(reg_autocovs122[2 ,2 ,1 , ], c(5.258146, 9.877770), tolerance=1e-5)
  expect_equal(reg_autocovs112csWAR[1, ,2 , 1], c(0.2477767, 0.2201706), tolerance=1e-5)
  expect_equal(reg_autocovs222csLAR[1, , 2, 2], c(9.302809, -21.716886), tolerance=1e-5)
})



































