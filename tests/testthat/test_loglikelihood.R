context("log-likelihood")
library(gmvarkit)

# These tests are very brief!

## A(M)(p)_(p)(M)(d), only d=2 in the example time series
data <- cbind(10*eurusd[,1], 100*eurusd[,2])

# p=1, M=1, d=2
phi10_112 <- c(1.03, 2.36)
A11_112 <- matrix(c(0.8, 0.06, 0.04, 0.9), nrow=2, byrow=FALSE)
Omega1_112 <- matrix(c(0.93, -0.15, -0.15, 5.20), nrow=2, byrow=FALSE)

theta_112 <- c(phi10_112, vec(A11_112), vech(Omega1_112))

W_112 <- t(chol(Omega1_112))
theta_112sWC <- c(phi10_112, vec(A11_112), Wvec(W_112)) # SGMVAR, W constrained by Cholesky

# p=2, M=1, d=2
phi10_212 <- c(1.03, 2.36)
A11_212 <- matrix(c(1.25, 0.06, 0.04, 1.34), nrow=2, byrow=FALSE)
A12_212 <- matrix(c(-0.29, -0.08, -0.05, -0.36), nrow=2, byrow=FALSE)
Omega1_212 <- matrix(c(0.93, -0.15, -0.15, 5.20), nrow=2, byrow=FALSE)

theta_212 <- c(phi10_212, vec(A11_212), vec(A12_212), vech(Omega1_212))

W_212 <- t(chol(Omega1_212))
theta_212sWC <- c(phi10_212, vec(A11_212), vec(A12_212), Wvec(W_212)) # SGMVAR, W constrained by Cholesky

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
A21_222 <- A11_222
A22_222 <- A12_222
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

# p=1, M=4, d=2
theta_142 <- c(22.743598, 98.461525, -0.006778, 0.212983, -0.12094, -0.061364,
               1.313218, -3.305944, 14.459664, 44.482964, 82.615109, -0.344603,
               0.082359, -0.168216, 0.341988, 4.012248, -3.555077, 6.449586,
               14.672977, 100.177417, 0.345826, -0.338611, -0.18687, -0.221931,
               12.086199, 3.677137, 1.331716, 17.668096, 129.042416, 0.628735,
               -0.026376, 0.185899, -0.199485, 0.470336, 0.980442, 7.146605,
               0.427396, 0.417413, 0.142379)


theta_112_mu <- change_parametrization(p=1, M=1, d=2, params=theta_112, change_to="mean")
theta_212_mu <- change_parametrization(p=2, M=1, d=2, params=theta_212, change_to="mean")
theta_122_mu <- change_parametrization(p=1, M=2, d=2, params=theta_122, change_to="mean")
theta_222_mu <- change_parametrization(p=2, M=2, d=2, params=theta_222, change_to="mean")

theta_112sWC_mu <- change_parametrization(p=1, M=1, d=2, params=theta_112sWC, structural_pars=list(W=W_112), change_to="mean")
theta_222s_mu <- change_parametrization(p=2, M=2, d=2, params=theta_222s, structural_pars=list(W=W_222), change_to="mean")


test_that("loglikelihood_int works correctly", {
  expect_equal(loglikelihood_int(data=data, p=1, M=1, params=theta_112, conditional=FALSE), -7013.364, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=data, p=2, M=1, params=theta_212, conditional=FALSE), -1428.165, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=data, p=1, M=2, params=theta_122, conditional=FALSE), -30712.15, tolerance=1e-2)
  expect_equal(loglikelihood_int(data=data, p=2, M=2, params=theta_222, conditional=FALSE), -1095.952, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=data, p=2, M=2, params=theta_222, conditional=TRUE), -1084.186, tolerance=1e-3)

  expect_equal(loglikelihood_int(data=data, p=1, M=1, params=theta_112_mu, conditional=FALSE, parametrization="mean"), -7013.364, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=data, p=2, M=1, params=theta_212_mu, conditional=FALSE, parametrization="mean"), -1428.165, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=data, p=1, M=2, params=theta_122_mu, conditional=FALSE, parametrization="mean"), -30712.15, tolerance=1e-2)
  expect_equal(loglikelihood_int(data=data, p=2, M=2, params=theta_222_mu, conditional=FALSE, parametrization="mean"), -1095.952, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=data, p=2, M=2, params=theta_222_mu, conditional=TRUE, parametrization="mean"), -1084.186, tolerance=1e-3)

  expect_equal(loglikelihood_int(data=data, p=1, M=4, params=theta_142, conditional=FALSE, parametrization="intercept"), -20948.16, tolerance=1e-2)
  expect_equal(loglikelihood_int(data=data, p=1, M=4, params=theta_142, conditional=TRUE, parametrization="intercept"), -20947.31, tolerance=1e-2)
  expect_equal(loglikelihood_int(data=data, p=1, M=4, params=theta_142, conditional=TRUE, parametrization="mean"), -30299.49, tolerance=1e-2)

  # SGMVAR
  expect_equal(loglikelihood_int(data=data, p=1, M=1, params=theta_112sWC, structural_pars=list(W=W_112), conditional=FALSE),
               -7013.364, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=data, p=2, M=1, params=theta_212sWC, structural_pars=list(W=W_212), conditional=FALSE),
               -1428.165, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=data, p=1, M=2, params=theta_122s, structural_pars=list(W=W_122), conditional=FALSE),
               -30712.15, tolerance=1e-2)
  expect_equal(loglikelihood_int(data=data, p=2, M=2, params=theta_222s, structural_pars=list(W=W_222), conditional=FALSE),
               -1095.952, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=data, p=2, M=2, params=theta_222s, structural_pars=list(W=W_222), conditional=TRUE),
               -1084.186, tolerance=1e-3)

  expect_equal(loglikelihood_int(data=data, p=1, M=1, params=theta_112sWC_mu, structural_pars=list(W=W_112), conditional=FALSE, parametrization="mean"),
               -7013.364, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=data, p=2, M=2, params=theta_222s_mu, structural_pars=list(W=W_222), conditional=FALSE, parametrization="mean"),
               -1095.952, tolerance=1e-3)
})


# Constraining AR-parameters to be the same for all regimes
rbind_diags <- function(p, M, d) {
  I <- diag(p*d^2)
  Reduce(rbind, replicate(M, I, simplify=FALSE))
}

# p=1, M=1, d=2
C_112 <- rbind_diags(p=1, M=1, d=2)
theta_112c <- c(phi10_112, vec(A11_112), vech(Omega1_112))

theta_112csWAR <- c(phi10_112, vec(A11_112), Wvec(W_112)) # SGMVAR W and AR

# p=2, M=1, d=2
C_212 <- rbind_diags(p=2, M=1, d=2)
theta_212c <- c(phi10_212, vec(A11_212), vec(A12_212), vech(Omega1_212))

theta_212csWAR <- c(phi10_212, vec(A11_212), vec(A12_212), Wvec(W_212)) # SGMVAR W and AR

# p=1, M=2, d=2
C_122 <- rbind_diags(p=1, M=2, d=2)
theta_122c <- c(phi10_122, phi20_122, vec(A11_122), vech(Omega1_122), vech(Omega2_122), alpha1_122)

C_lambda_122 <- matrix(c(7, 1), nrow=2)
theta_122csLAR <- c(phi10_122, phi20_122, vec(A11_122), vec(W_122), 1, alpha1_122) # SGMVAR lambdas and AR

# p=2, M=2, d=2
C_222 <- rbind_diags(p=2, M=2, d=2)
theta_222c <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vech(Omega1_222), vech(Omega2_222), alpha1_222)

C_lambda_222 <- matrix(c(1, 2), nrow=2)
theta_222csLAR <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(W_222), 0.2, alpha1_222) # SGMVAR lambdas and AR

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

WL_222c2 <- diag_Omegas(Omega1_222_c2, Omega2_222_c2)
W_222c2 <- matrix(WL_222c2[1:(2^2)], nrow=2, byrow=FALSE)
lambdas_222c2 <- WL_222c2[(2^2 + 1):length(WL_222c2)]
theta_222_c2s <- c(phi10_222_c2, phi20_222_c2, 1.26, 1.34, -0.29, -0.36, vec(W_222c2), lambdas_222c2, alpha1_222_c2) # SGMVAR AR

test_that("loglikelihood_int works correctly for constrained models", {
  expect_equal(loglikelihood_int(data=data, p=1, M=1, params=theta_112c, conditional=FALSE, constraints=C_112), -7013.364, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=data, p=2, M=1, params=theta_212c, conditional=FALSE, constraints=C_212), -1428.165, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=data, p=1, M=2, params=theta_122c, conditional=TRUE, constraints=C_122), -30639.31, tolerance=1e-2)
  expect_equal(loglikelihood_int(data=data, p=2, M=2, params=theta_222c, conditional=FALSE, constraints=C_222), -1095.952, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=data, p=2, M=2, params=theta_222_c2, conditional=TRUE, constraints=C_222_2), -1096.98, tolerance=1e-2)

  # SGMVAR
  expect_equal(loglikelihood_int(data=data, p=1, M=1, params=theta_112csWAR, conditional=FALSE, constraints=C_112,
                                 structural_pars=list(W=W_112)), -7013.364, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=data, p=2, M=1, params=theta_212csWAR, conditional=FALSE, constraints=C_212,
                                 structural_pars=list(W=W_212)), -1428.165, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=data, p=1, M=2, params=theta_122csLAR, conditional=TRUE, constraints=C_122,
                                 structural_pars=list(W=W_122, C_lambda=C_lambda_122)), -33364.45, tolerance=1e-2)
  expect_equal(loglikelihood_int(data=data, p=2, M=2, params=theta_222csLAR, conditional=FALSE, constraints=C_222,
                                 structural_pars=list(W=W_222, C_lambda=C_lambda_222)), -1647.087, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=data, p=2, M=2, params=theta_222_c2s, conditional=TRUE, constraints=C_222_2,
                                 structural_pars=list(W=W_222)), -1096.98, tolerance=1e-2)
})


test_that("user loglikelihood works correctly", {
  expect_equal(loglikelihood(data=data, p=1, M=1, params=theta_112, parametrization="intercept", conditional=FALSE), -7013.364, tolerance=1e-3)
  expect_equal(loglikelihood(data, p=2, M=2, params=theta_222, parametrization="intercept", conditional=TRUE), -1084.186, tolerance=1e-3)
  expect_equal(loglikelihood(data=data, p=1, M=2, params=theta_122_mu, conditional=FALSE, parametrization="mean"), -30712.15, tolerance=1e-2)
  expect_equal(loglikelihood(data=data, p=2, M=2, params=theta_222c, conditional=FALSE, parametrization="intercept", constraints=C_222), -1095.952, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=data, p=2, M=2, params=theta_222_c2, conditional=TRUE, parametrization="intercept", constraints=C_222_2), -1096.98, tolerance=1e-2)

  # SGMVAR
  expect_equal(loglikelihood(data=data, p=1, M=1, params=theta_112sWC, parametrization="intercept", conditional=FALSE,
                             structural_pars=list(W=W_112)), -7013.364, tolerance=1e-3)
  expect_equal(loglikelihood(data, p=2, M=2, params=theta_222s, parametrization="intercept", conditional=TRUE,
                             structural_pars=list(W=W_222)), -1084.186, tolerance=1e-3)
  expect_equal(loglikelihood(data=data, p=1, M=1, params=theta_112sWC_mu, conditional=FALSE, parametrization="mean",
                             structural_pars=list(W=W_112)), -7013.364, tolerance=1e-2)
  expect_equal(loglikelihood(data=data, p=2, M=2, params=theta_222csLAR, conditional=FALSE, parametrization="intercept", constraints=C_222,
                             structural_pars=list(W=W_222, C_lambda=C_lambda_222)), -1647.087, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=data, p=2, M=2, params=theta_222_c2s, conditional=TRUE, parametrization="intercept",
                                 constraints=C_222_2, structural_pars=list(W=W_222)), -1096.98, tolerance=1e-2)

})


test_that("cond_moments works correctly", {
  expect_equal(cond_moments(data=data, p=1, M=1, params=theta_112, parametrization="intercept", to_return="total_ccovs")[, 2, 20], c(-0.15, 5.20), tolerance=1e-2)
  expect_equal(cond_moments(data=data, p=1, M=1, params=theta_112, parametrization="intercept", to_return="regime_cmeans")[1, , 1], c(6.172968, 105.106160), tolerance=1e-5)
  expect_equal(cond_moments(data=data, p=1, M=1, params=theta_112, parametrization="intercept", to_return="total_cmeans")[100, ], c(1.611256, 105.138260), tolerance=1e-5)
  expect_equal(cond_moments(data, p=2, M=2, params=theta_222, parametrization="intercept", to_return="total_ccovs")[1, , 200], c(1.2645665, 0.1038158), tolerance=1e-5)
  expect_equal(cond_moments(data, p=2, M=2, params=theta_222, parametrization="intercept", to_return="total_cmeans")[1, ], c(2.752923, 112.571724), tolerance=1e-5)
  expect_equal(cond_moments(data=data, p=1, M=2, params=theta_122_mu, parametrization="mean", to_return="total_ccovs")[, 2, 100], c(3.56, 9.80), tolerance=1e-2)
  expect_equal(cond_moments(data=data, p=2, M=2, params=theta_222c, parametrization="intercept", constraints=C_222, to_return="total_cmeans")[13, ], c(24.78782, 122.40034), tolerance=1e-5)
  expect_equal(cond_moments(data=data, p=2, M=2, params=theta_222_c2, parametrization="intercept", constraints=C_222_2, to_return="total_ccovs")[ , 2, 13], c(3.459987, 9.619985), tolerance=1e-5)

  # SGMVAR
  expect_equal(cond_moments(data=data, p=1, M=1, params=theta_112sWC, parametrization="intercept",
                            structural_pars=list(W=W_112), to_return="total_ccovs")[, 2, 20], c(-0.15, 5.20), tolerance=1e-2)
  expect_equal(cond_moments(data=data, p=1, M=1, params=theta_112sWC, parametrization="intercept",
                            structural_pars=list(W=W_112), to_return="regime_cmeans")[1, , 1], c(6.172968, 105.106160), tolerance=1e-5)
  expect_equal(cond_moments(data=data, p=1, M=1, params=theta_112sWC, parametrization="intercept",
                            structural_pars=list(W=W_112), to_return="total_cmeans")[100, ], c(1.611256, 105.138260), tolerance=1e-5)
  expect_equal(cond_moments(data, p=2, M=2, params=theta_222s, parametrization="intercept",
                            structural_pars=list(W=W_222), to_return="total_ccovs")[1, , 200],c(1.2645665, 0.1038158), tolerance=1e-5)
  expect_equal(cond_moments(data, p=2, M=2, params=theta_222s, parametrization="intercept",
                            structural_pars=list(W=W_222), to_return="total_cmeans")[1, ], c(2.752923, 112.571724), tolerance=1e-5)
  expect_equal(cond_moments(data=data, p=2, M=2, params=theta_222s_mu, parametrization="mean",
                            structural_pars=list(W=W_222), to_return="total_ccovs")[, 2, 100], c(0.01131515, 5.39283069), tolerance=1e-5)
  expect_equal(cond_moments(data=data, p=2, M=2, params=theta_222csLAR, parametrization="intercept",
                            constraints=C_222, structural_pars=list(W=W_222, C_lambda=C_lambda_222), to_return="total_cmeans")[13, ],
               c(24.68823, 122.31647), tolerance=1e-5)
  expect_equal(cond_moments(data=data, p=2, M=2, params=theta_222_c2s, parametrization="intercept",
                            constraints=C_222_2, structural_pars=list(W=W_222), to_return="total_ccovs")[ , 2, 13],
               c(3.459987, 9.619985), tolerance=1e-5)

})























