context("log-likelihood")
library(gmvarkit)

# These tests are very brief!

## A(M)(p)_(p)(M)(d), only d=2 in the example time series
data <- cbind(10*eurusd[,1], 100*eurusd[,2])
set.seed(1); data2 <- cbind(data, round(rnorm(nrow(data), 3)))

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
theta_122s_mu <- change_parametrization(p=1, M=1, d=2, params=theta_122s, structural_pars=list(W=W_122), change_to="mean")
theta_222s_mu <- change_parametrization(p=2, M=2, d=2, params=theta_222s, structural_pars=list(W=W_222), change_to="mean")

# p=4, M=1, d=3, data2
theta_413 <- c(0.955, 1.968, 5.832, 0.384, 0.895, 1.519, 1.175, 0.067, -0.68, 1.517, -0.654, -0.237, -1.343, -1.188, -0.403,
              0.98, -0.285, -1.787, 0.546, -0.545, -1.633, 0.277, -0.33, 2.037, 0.563, 1.126, -1.264, 0.179, 2.345, -1.126,
              -0.152, -0.82, 0.428, -0.332, 1.776, -0.976, -0.256, 0.374, -0.82, 1.964, -0.42, 0.798, 3.931, -0.436, 0.404)

# p=4, M=2, d=3, data2
theta_423 <- c(3.402, 1.922, 5.069, -0.414, -1.492, 1.286, -0.351, 1.074, -0.613, 0.504, 0.386, -0.193, 1.751, 0.754, -0.118,
               -1.15, -1.44, -0.301, -0.048, 0.229, -0.575, 0.444, 1.361, -0.128, 1.267, 0.48, 1.709, -0.7, -0.247, 0.295,
               -0.553, -0.332, 0.243, 0.478, 1.026, 0.157, -0.591, -0.648, -1.037, 0.388, -0.464, -0.408, 0.689, 1.137, 4.44,
               1.717, 3.82, 4.153, -1.21, -0.002, -0.061, 0.898, -0.349, 0.269, -0.228, -0.55, -0.38, 0.031, 1.835, -0.561,
               1.542, 0.19, -0.257, -0.212, -0.014, -0.244, 0.556, -0.07, 0.048, 0.929, -1.442, 0.508, -0.383, -0.122, -0.515,
               0.3, -1.578, -0.12, 0.33, -1.092, 0.21, -0.371, -0.634, 0.613, 0.26, -0.536, -0.009, 1.414, -0.158, 0.109, 0.822)

# p=7, M=2, d=2
theta_722 <- c(0.93, 1.505, -0.74, 1.106, -0.304, 2.488, 0.107, -2.008, -2.701, -2.937, 2.696, -0.594, 2.472, 2.991, 1.15, 2.901,
               -0.63, -5.265, -2.446, -0.921, 4.506, 5.536, -0.109, -0.302, -3.8, -1.724, 0.023, -0.137, 0.664, -0.214, 1.875,
               -0.795, 0.501, -0.001, 0.664, 3.122, 2.272, -0.125, -1.659, -2.756, 1.051, -1.534, -2.175, -0.832, -11.509, 2.103,
               1.778, 4.395, 20.155, -0.55, 3.99, -4.606, -16.731, -2.187, -9.805, 2.456, 8.306, 3.297, 11.725, -0.494, -1.857,
               -1.317, -4.49, 1.542, -0.15, 0.291, 0.601)

# p=6, M=3, d=2
theta_632 <- c(-0.476, 1.711, -0.686, 0.786, 0.663, -0.668, 1.299, 0.204, 0.198, -0.276, 1.582, -0.379, -1.283, 0.614, 0.425, 1.522,
               -0.951, 1.032, -0.231, 1.743, 0.167, -1.887, -0.309, 0.15, 0.089, -1.959, 0.907, 0.985, 2.304, 1.082, 1.835, 0.24,
               -0.373, 1.39, 0.618, 0.787, -1.068, 0.508, 1.358, 1.416, 0.443, -1.82, -0.103, -0.729, 0.674, 0.236, 0.218, -1.159,
               0.131, 0.508, 0.332, 0.53, 0.285, -0.35, 0.354, 0.023, -0.106, 3.524, 0.419, -0.002, -1.14, -0.544, 1.281, 2.777,
               -1.086, -2.576, -1.891, -3.575, 1.624, 0.173, 2.726, 4.498, -0.005, -0.975, -3.034, -3.234, 1.349, 1.182, 2.611, 1.811,
               -0.743, 0.358, -1.082, 0.037, 1.51, -0.558, 1.022, 0.461, 0.353)

# p=8, M=1, d=2
theta_812 <- c(2.54, 1.341, -0.316, -0.405, -0.181, 0.744, 0.505, -0.795, -0.221, 0.109, 1.382, 1.453, -0.659, -1.004, 0.44, 0.314,
               0.437, 1.082, -2.219, -0.177, 1.155, -0.081, -0.663, 0.478, -0.869, -0.114, 0.881, -0.098, -0.868, 0.448, 0.308, -0.535,
               0.122, 0.165, 0.051, 0.086, 3.732)

# Same means
mu_112 <- pick_phi0(p=1, M=1, d=2, params=theta_112_mu)
mu_212 <- pick_phi0(p=2, M=1, d=2, params=theta_212_mu)
mu_122 <- pick_phi0(p=1, M=2, d=2, params=theta_122_mu)
mu_222 <- pick_phi0(p=2, M=2, d=2, params=theta_222_mu)
mu_112sWC <- pick_phi0(p=1, M=1, d=2, params=theta_112sWC_mu, structural_pars=list(W=W_112))
mu_122s <- pick_phi0(p=1, M=2, d=2, params=theta_122s_mu, structural_pars=list(W=W_122))
mu_222s <- pick_phi0(p=2, M=2, d=2, params=theta_222s_mu, structural_pars=list(W=W_222))

theta_112_int <- c(mu_112, vec(A11_112), vech(Omega1_112)) # same_means=list(1)
theta_212_int <- c(mu_212, vec(A11_212), vec(A12_212), vech(Omega1_212)) # same_means=list(1)
theta_122_int <- c(mu_122[,1], vec(A11_122), vec(A21_122), vech(Omega1_122), vech(Omega2_122), alpha1_122) # same_means=list(1:2)
theta_122_int2 <- c(vec(mu_122), vec(A11_122), vec(A21_122), vech(Omega1_122), vech(Omega2_122), alpha1_122) # same_means=list(1, 2)
theta_222_int <- c(mu_222[,1], vec(A11_222), vec(A12_222), vec(A21_222), vec(A22_222), vech(Omega1_222), vech(Omega2_222), alpha1_222) # same_means=list(1:2)
theta_112sWC_int <- c(mu_112sWC, vec(A11_112), Wvec(W_112)) # structural_pars=list(W=W_112), same_means=list(1)
theta_122s_int <- c(mu_122s[,1], vec(A11_122), vec(A21_122), vec(W_122), lambdas_122, alpha1_122) # structural_pars=list(W=W_122), same_means=list(1:2)
theta_222s_int <- c(mu_222s[,1], vec(A11_222), vec(A12_222), vec(A21_222),
                    vec(A22_222), vec(W_222), lambdas_222, alpha1_222)  # structural_pars=list(W=W_222), same_means=list(1:2)

test_that("loglikelihood_int works correctly", {
  # Regular
  expect_equal(loglikelihood_int(data=data, p=1, M=1, params=theta_112, conditional=FALSE), -7013.364, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=data, p=2, M=1, params=theta_212, conditional=FALSE), -1428.165, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=data, p=1, M=2, params=theta_122, conditional=FALSE), -30712.15, tolerance=1e-2)
  expect_equal(loglikelihood_int(data=data, p=2, M=2, params=theta_222, conditional=FALSE), -1095.952, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=data, p=2, M=2, params=theta_222, conditional=TRUE), -1084.186, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=500*data, p=2, M=2, params=theta_222, conditional=FALSE), -4467536, tolerance=1)

  expect_equal(loglikelihood_int(data=data, p=1, M=1, params=theta_112_mu, conditional=FALSE, parametrization="mean"), -7013.364, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=data, p=2, M=1, params=theta_212_mu, conditional=FALSE, parametrization="mean"), -1428.165, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=data, p=1, M=2, params=theta_122_mu, conditional=FALSE, parametrization="mean"), -30712.15, tolerance=1e-2)
  expect_equal(loglikelihood_int(data=data, p=2, M=2, params=theta_222_mu, conditional=FALSE, parametrization="mean"), -1095.952, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=data, p=2, M=2, params=theta_222_mu, conditional=TRUE, parametrization="mean"), -1084.186, tolerance=1e-3)

  expect_equal(loglikelihood_int(data=data, p=1, M=4, params=theta_142, conditional=FALSE, parametrization="intercept"), -20948.16, tolerance=1e-2)
  expect_equal(loglikelihood_int(data=data, p=1, M=4, params=theta_142, conditional=TRUE, parametrization="intercept"), -20947.31, tolerance=1e-2)
  expect_equal(loglikelihood_int(data=data, p=1, M=4, params=theta_142, conditional=TRUE, parametrization="mean"), -30299.49, tolerance=1e-2)

  # p*d >= 12
  expect_equal(loglikelihood_int(data=data2, p=4, M=1, params=theta_413, conditional=FALSE), -223585.4, tolerance=1e-1)
  expect_equal(loglikelihood_int(data=data2, p=4, M=2, params=theta_423, conditional=FALSE), -176134.4, tolerance=1e-1)
  expect_equal(loglikelihood_int(data=data, p=7, M=2, params=theta_722, conditional=FALSE), -119265.1, tolerance=1e-1)
  expect_equal(loglikelihood_int(data=data, p=6, M=3, params=theta_632, conditional=FALSE), -181090, tolerance=1)
  expect_equal(loglikelihood_int(data=data, p=8, M=1, params=theta_812, conditional=FALSE), -188731, tolerance=1)

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

  # same_means
  expect_equal(loglikelihood_int(data=data, p=1, M=1, params=theta_112_int, parametrization="mean", same_means=list(1), conditional=FALSE),
               -7013.364, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=data, p=2, M=1, params=theta_212_int, parametrization="mean", same_means=list(1), conditional=FALSE),
               -1428.165, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=data, p=1, M=2, params=theta_122_int, parametrization="mean", same_means=list(1:2), conditional=FALSE),
               -31768.84, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=data, p=1, M=2, params=theta_122_int2, parametrization="mean", same_means=list(1, 2), conditional=FALSE),
               -30712.15, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=data, p=2, M=2, params=theta_222_int, parametrization="mean", same_means=list(1:2), conditional=FALSE),
               -1107.581, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=data, p=1, M=1, params=theta_112sWC_int, parametrization="mean", structural_pars=list(W=W_112), same_means=list(1), conditional=FALSE),
               -7013.364, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=data, p=1, M=2, params=theta_122s_int, parametrization="mean", structural_pars=list(W=W_122), same_means=list(1:2), conditional=FALSE),
               -33079.93, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=data, p=2, M=2, params=theta_222s_int, parametrization="mean", structural_pars=list(W=W_222), same_means=list(1:2), conditional=FALSE),
               -1107.581, tolerance=1e-3)
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

# Same means
theta_112c_mu <- change_parametrization(p=1, M=1, d=2, params=theta_112c, constraints=C_112, change_to="mean")
theta_122c_mu <- change_parametrization(p=1, M=2, d=2, params=theta_122c, constraints=C_122, change_to="mean")
theta_222c_mu <- change_parametrization(p=2, M=2, d=2, params=theta_222c, constraints=C_222, change_to="mean")
theta_112csWAR_mu <- change_parametrization(p=1, M=1, d=2, params=theta_112csWAR, constraints=C_112, structural_pars=list(W=W_112), change_to="mean")
theta_122csLAR_mu <- change_parametrization(p=1, M=2, d=2, params=theta_122csLAR, constraints=C_112, structural_pars=list(W=W_122, C_lambda=C_lambda_122), change_to="mean")
theta_222c2s_mu <- change_parametrization(p=2, M=2, d=2, params=theta_222_c2s, constraints=C_222_2, structural_pars=list(W=W_222), change_to="mean")

mu_112c <- theta_112c_mu[1:2]
mu_122c <-  theta_122c_mu[1:4]
mu_222c <-  theta_222c_mu[1:4]
mu_112csWAR <- theta_112csWAR_mu[1:2]
mu_122csLAR <- theta_122csLAR_mu[1:4]
mu_222c2s <- theta_222c2s_mu[1:4]

theta_112c_int <- c(mu_112c, vec(A11_112), vech(Omega1_112)) # same_means=list(1)
theta_122c_int <- c(mu_122c, vec(A11_122), vech(Omega1_122), vech(Omega2_122), alpha1_122) # same_means=list(1, 2)
theta_222c_int <- c(mu_222c[1:2], vec(A11_222), vec(A12_222), vech(Omega1_222), vech(Omega2_222), alpha1_222) # same_means=list(1:2)

theta_112csWAR_int <- c(mu_112c, vec(A11_112), Wvec(W_112)) # same_means=list(1), structural_pars=list(W=W_112)
theta_122csLAR_int <- c(mu_122csLAR[1:2], vec(A11_122), vec(W_122), 1, alpha1_122) # same_means=list(1:2), structural_pars=list(W=W_122, C_lambda=C_lambda_122)
theta_222c2s_int <- c(mu_222c2s[1:2], 1.26, 1.34, -0.29, -0.36, vec(W_222c2), lambdas_222c2, alpha1_222_c2) # constraints=C_222_2, same_means=list(1:2), structural_pars=list(W=W_222)


test_that("loglikelihood_int works correctly for constrained models", {
   # Regular
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
  expect_equal(loglikelihood_int(data=data, p=2, M=2, params=theta_222c2s_mu, conditional=TRUE, constraints=C_222_2,
                                 structural_pars=list(W=W_222), parametrization="mean"), -1096.98, tolerance=1e-2)

  # Same means
  expect_equal(loglikelihood_int(data=data, p=1, M=1, params=theta_112c_int, same_means=list(1), conditional=FALSE, constraints=C_112, parametrization="mean"),
               -7013.364, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=data, p=1, M=2, params=theta_122c_int, same_means=list(1, 2), conditional=FALSE, constraints=C_122, parametrization="mean"),
               -30712.15, tolerance=1e-2)
  expect_equal(loglikelihood_int(data=data, p=2, M=2, params=theta_222c_int, same_means=list(1:2), conditional=FALSE, constraints=C_222, parametrization="mean"),
               -1107.581, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=data, p=1, M=1, params=theta_112csWAR_int, same_means=list(1), conditional=FALSE, constraints=C_112, structural_pars=list(W=W_112), parametrization="mean"),
               -7013.364, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=data, p=1, M=2, params=theta_122csLAR_int, same_means=list(1:2), conditional=FALSE, constraints=C_122,
                                 structural_pars=list(W=W_122, C_lambda=C_lambda_122), parametrization="mean"),
               -34521.1, tolerance=1e-1)
  expect_equal(loglikelihood_int(data=200*data, p=1, M=2, params=theta_122csLAR_int, same_means=list(1:2), conditional=FALSE, constraints=C_122,
                                 structural_pars=list(W=W_122, C_lambda=C_lambda_122), parametrization="mean"),
               -3583346, tolerance=1)
  expect_equal(loglikelihood_int(data=data, p=2, M=2, params=theta_222c2s_int, same_means=list(1:2), conditional=FALSE, constraints=C_222_2,
                                 structural_pars=list(W=W_222), parametrization="mean"),
               -1121.216, tolerance=1e-3)
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

  # Same means
  expect_equal(loglikelihood(data=data, p=1, M=1, params=theta_112c_int, same_means=list(1), conditional=FALSE, constraints=C_112, parametrization="mean"),
               -7013.364, tolerance=1e-3)
  expect_equal(loglikelihood(data=data, p=1, M=2, params=theta_122c_int, same_means=list(1, 2), conditional=FALSE, constraints=C_122, parametrization="mean"),
               -30712.15, tolerance=1e-2)
  expect_equal(loglikelihood(data=data, p=2, M=2, params=theta_222c_int, same_means=list(1:2), conditional=FALSE, constraints=C_222, parametrization="mean"),
               -1107.581, tolerance=1e-3)
  expect_equal(loglikelihood(data=data, p=1, M=1, params=theta_112csWAR_int, same_means=list(1), conditional=FALSE, constraints=C_112, structural_pars=list(W=W_112), parametrization="mean"),
               -7013.364, tolerance=1e-3)
  expect_equal(loglikelihood(data=data, p=1, M=2, params=theta_122csLAR_int, same_means=list(1:2), conditional=FALSE, constraints=C_122,
                                 structural_pars=list(W=W_122, C_lambda=C_lambda_122), parametrization="mean"),
               -34521.1, tolerance=1e-1)
  expect_equal(loglikelihood(data=200*data, p=1, M=2, params=theta_122csLAR_int, same_means=list(1:2), conditional=FALSE, constraints=C_122,
                                 structural_pars=list(W=W_122, C_lambda=C_lambda_122), parametrization="mean"),
               -3583346, tolerance=1)
  expect_equal(loglikelihood(data=data, p=2, M=2, params=theta_222c2s_int, same_means=list(1:2), conditional=FALSE, constraints=C_222_2,
                                 structural_pars=list(W=W_222), parametrization="mean"),
               -1121.216, tolerance=1e-3)
  expect_error(loglikelihood(data=data, p=2, M=2, params=theta_222c2s_int, same_means=list(1:2), conditional=FALSE, constraints=C_222_2,
                             structural_pars=list(W=W_222), parametrization="intercept"))

  expect_equal(loglikelihood(data=data, p=2, M=2, params=theta_222_int, parametrization="mean", same_means=list(1:2), conditional=FALSE),
               -1107.581, tolerance=1e-3)
  expect_equal(loglikelihood(data=data, p=2, M=2, params=theta_222s_int, parametrization="mean", structural_pars=list(W=W_222), same_means=list(1:2), conditional=FALSE),
               -1107.581, tolerance=1e-3)

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

  # Same means
  expect_equal(cond_moments(data=data, p=1, M=1, params=theta_112c_int, same_means=list(1), constraints=C_112, parametrization="mean",
                            to_return="regime_cmeans")[2, , 1],
               c(7.42388, 103.95566), tolerance=1e-5)
  expect_equal(cond_moments(data=data, p=1, M=2, params=theta_122c_int, same_means=list(1, 2), constraints=C_122, parametrization="mean", to_return="total_ccovs")[2, , 200],
               c(3.56, 9.80), tolerance=1e-2)
  expect_equal(cond_moments(data=data, p=2, M=2, params=theta_222c_int, same_means=list(1:2), constraints=C_222, parametrization="mean", to_return="total_cmeans")[1,],
               c(2.5687, 112.4166), tolerance=1e-3)
  expect_equal(cond_moments(data=data, p=1, M=1, params=theta_112csWAR_int, same_means=list(1), constraints=C_112, structural_pars=list(W=W_112), parametrization="mean",
                            to_return="total_cmeans")[1, ],
               c(6.172968, 105.106160), tolerance=1e-4)
  expect_equal(cond_moments(data=data, p=1, M=2, params=theta_122csLAR_int, same_means=list(1:2), constraints=C_122,
                             structural_pars=list(W=W_122, C_lambda=C_lambda_122), parametrization="mean", to_return="total_ccovs")[1, , 100],
               c(5.708931, 3.695015), tolerance=1e-4)
  expect_equal(cond_moments(data=200*data, p=1, M=2, params=theta_122csLAR_int, same_means=list(1:2), constraints=C_122,
                             structural_pars=list(W=W_122, C_lambda=C_lambda_122), parametrization="mean", to_return="regime_cmeans")[1, , 2],
               c(-6715.682, 20586.296), tolerance=1e-3)
  expect_equal(cond_moments(data=data, p=2, M=2, params=theta_222c2s_int, same_means=list(1:2), constraints=C_222_2,
                             structural_pars=list(W=W_222), parametrization="mean", to_return="total_ccovs")[, 2, 13],
               c(3.459943, 9.619934), tolerance=1e-4)
  expect_equal(cond_moments(data=data, p=2, M=2, params=theta_222_int, parametrization="mean", same_means=list(1:2), to_return="total_cmeans")[13,],
               c(24.02782, 121.76034), tolerance=1e-4)
  expect_equal(cond_moments(data=data, p=2, M=2, params=theta_222s_int, parametrization="mean", structural_pars=list(W=W_222), same_means=list(1:2),
                            to_return="total_ccovs")[2, , 1],
               c(1.392173, 7.112128), tolerance=1e-4)
})























