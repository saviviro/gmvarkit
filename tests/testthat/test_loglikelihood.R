context("log-likelihood")
library(gmvarkit)

set.seed(1); data2 <- cbind(gdpdef, round(rnorm(nrow(gdpdef)), 3))

# p=1, M=1, d=2
phi10_112 <- c(0.65, 0.7)
A11_112 <- matrix(c(0.29, 0.02, -0.14, 0.9), nrow=2, byrow=FALSE)
Omega1_112 <- matrix(c(0.60, 0.01, 0.01, 0.07), nrow=2, byrow=FALSE)

theta_112 <- c(phi10_112, vec(A11_112), vech(Omega1_112))

W_112 <- t(chol(Omega1_112))
theta_112sWC <- c(phi10_112, vec(A11_112), Wvec(W_112)) # SGMVAR, W constrained by Cholesky

# p=2, M=1, d=2
phi10_212 <- c(0.53, 0.03)
A11_212 <- matrix(c(0.23, 0.02, -0.17, 0.66), nrow=2, byrow=FALSE)
A12_212 <- matrix(c(0.18, 0.02, 0.04, 0.26), nrow=2, byrow=FALSE)
Omega1_212 <- matrix(c(0.58, 0.01, 0.01, 0.06), nrow=2, byrow=FALSE)

theta_212 <- c(phi10_212, vec(A11_212), vec(A12_212), vech(Omega1_212))

W_212 <- t(chol(Omega1_212))
theta_212sWC <- c(phi10_212, vec(A11_212), vec(A12_212), Wvec(W_212)) # SGMVAR, W constrained by Cholesky

# p=1, M=2, d=2
phi10_122 <- c(0.55, 0.11)
A11_122 <- matrix(c(0.34, 0.05, -0.01, 0.72), nrow=2, byrow=FALSE)
Omega1_122 <- matrix(c(0.58, 0.01, 0.01, 0.06), nrow=2, byrow=FALSE)

phi20_122 <- c(0.17, 0.25)
A21_122 <- A11_122
Omega2_122 <- matrix(c(0.50, -0.01, -0.01, 0.20), nrow=2, byrow=FALSE)

alpha1_122 <- 0.60
upsilon1_122 <- c(phi10_122, vec(A11_122), vech(Omega1_122))
upsilon2_122 <- c(phi20_122, vec(A21_122), vech(Omega2_122))
theta_122 <- c(upsilon1_122, upsilon2_122, alpha1_122)

WL_122 <- diag_Omegas(Omega1_122, Omega2_122)
W_122 <- matrix(WL_122[1:(2^2)], nrow=2, byrow=FALSE)
lambdas_122 <- WL_122[(2^2 + 1):length(WL_122)]
theta_122s <- c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), vec(W_122), lambdas_122, alpha1_122) # SGMVAR


# p=2, M=2, d=2
phi10_222 <- c(0.36, 0.12)
A11_222 <- matrix(c(0.22, 0.06, -0.15, 0.39), nrow=2, byrow=FALSE)
A12_222 <- matrix(c(0.41, -0.01, 0.08, 0.3), nrow=2, byrow=FALSE)
Omega1_222 <- matrix(c(0.21, 0.01, 0.01, 0.03), nrow=2, byrow=FALSE)

phi20_222 <- c(0.48, 0.07)
A21_222 <- matrix(c(0.22, 0.02, -0.12, 0.72), nrow=2, byrow=FALSE)
A22_222 <- matrix(c(0.09, 0.03, 0.04, 0.19), nrow=2, byrow=FALSE)
Omega2_222 <- matrix(c(1.10, 0.01, 0.01, 0.11), nrow=2, byrow=FALSE)

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
theta_142 <- c(0.680381, 0.099736, 0.502918, 0.080781, -0.627348, 0.674579,
               0.37666, 0.020433, 0.016747, 0.285857, 0.205766, 0.392568,
               0.047474, 0.317407, 0.683117, 0.415324, -0.059849, 0.079795,
               1.927008, 0.687905, 0.036475, -0.014841, -0.793236, 0.638711,
               1.281068, 0.017391, 0.135752, 1.716725, 0.668851, -0.184389,
               -0.155109, -1.746672, -0.348001, 0.43192, -0.028927, 0.002053,
               0.428343, 0.347302, 0.152716)

theta_112_mu <- change_parametrization(p=1, M=1, d=2, params=theta_112, change_to="mean")
theta_212_mu <- change_parametrization(p=2, M=1, d=2, params=theta_212, change_to="mean")
theta_122_mu <- change_parametrization(p=1, M=2, d=2, params=theta_122, change_to="mean")
theta_222_mu <- change_parametrization(p=2, M=2, d=2, params=theta_222, change_to="mean")
theta_112sWC_mu <- change_parametrization(p=1, M=1, d=2, params=theta_112sWC, structural_pars=list(W=W_112), change_to="mean")
theta_122s_mu <- change_parametrization(p=1, M=1, d=2, params=theta_122s, structural_pars=list(W=W_122), change_to="mean")
theta_222s_mu <- change_parametrization(p=2, M=2, d=2, params=theta_222s, structural_pars=list(W=W_222), change_to="mean")

# p=4, M=1, d=3, data2 (note: very bad fit)
theta_413 <- c(0.955, 1.968, 5.832, 0.384, 0.895, 1.519, 1.175, 0.067, -0.68, 1.517, -0.654, -0.237, -1.343, -1.188, -0.403,
              0.98, -0.285, -1.787, 0.546, -0.545, -1.633, 0.277, -0.33, 2.037, 0.563, 1.126, -1.264, 0.179, 2.345, -1.126,
              -0.152, -0.82, 0.428, -0.332, 1.776, -0.976, -0.256, 0.374, -0.82, 1.964, -0.42, 0.798, 3.931, -0.436, 0.404)

# p=3, M=2, d=3, data2
theta_323 <- c(0.352, 0.085, 0.19, 0.166, 0.059, -0.35, -0.003, 0.351, 0.019, -0.003, 0.003, -0.015, 0.465, -0.021, -0.051,
               0.281, 0.128, -0.452, -0.011, 0.004, -0.048, -0.059, 0.036, 0.207, -0.284, 0.243, 0.315, 0.023, 0.001, 0.055,
               0.205, 0, 0.086, 0.031, -0.01, 0.957, 0.657, 0.047, 0.383, 0.289, -0.002, -0.069, -0.156, 0.7, 0.152, 0.071,
               -0.031, -0.055, 0.093, 0.044, 0.041, 0.058, 0.079, -0.658, 0.032, 0.054, -0.036, -0.057, 0.038, -0.107,
               -0.078, 0.138, 0.445, -0.166, 0.049, -0.237, 1.143, 0.02, -0.271, 0.095, 0.018, 0.688, 0.719)

# p=7, M=2, d=2
theta_722 <- c(0.556, -0.027, 0.267, 0.021, -0.15, 0.622, 0.078, 0.026, 0.158, 0.04, -0.027, 0.015, -0.242, 0.099, 0.05,
               0.057, -0.235, 0.27, -0.069, -0.002, 0.249, -0.135, -0.014, 0.032, 0.054, 0.058, -0.09, 0.081, 0.132,
               -0.053, 0.843, 0.021, 0.082, 0.537, 0.061, 0.167, 0.051, -0.112, 0.308, 0.372, 0.039, 0.211, 0.227,
               -0.139, 0.066, -0.371, 0.208, 0.142, -0.046, -0.071, 0.076, -0.164, 0.012, -0.163, 0.03, 0.165, 0.006,
               -0.496, 0.234, 0.271, -0.089, 0.225, -0.259, 0.218, 0.002, 0.021, 0.721)

# p=6, M=3, d=2
theta_632 <- c(0.946496, 0.080268, 0.030026, 0.454984, 0.153281, -0.597408, -0.343207, 0.27003, -0.009369, 0.398815,
               -0.557508, -0.606283, 0.570008, 0.729047, -0.636236, -0.324584, 0.784139, 0.157365, -0.64247, 0.100334,
               0.851107, -0.089478, -0.469049, 0.404556, 0.17415, -0.2745, 0.02781, 0.039447, 0.059343, 0.585345, -0.001013,
               0.325745, 0.052822, 0.020458, 0.464934, 0.25861, 0.021979, 0.005002, 0.028551, -0.161162, 0.033817, -0.407718,
               0.239431, 0.11308, -0.013225, 0.035824, 0.067477, -0.038552, -0.005739, 0.011806, 0.101888, 0.088671, -0.018675,
               -0.150412, 0.006073, 0.337413, 0.001703, 0.028568, 2.701142, 0.370117, 0.104281, -0.019346, -0.509732, 0.564074,
               -0.009026, -0.006529, 0.238403, 0.049195, -0.069648, 0.011972, -0.486556, 0.14413, -0.030511, 0.092299, -0.388031,
               0.253254, -0.125366, 0.010989, 0.221103, -0.201323, -0.086355, 0.064991, -0.174453, -0.091006, 1.080211, -0.009982,
               0.101493, 0.990494, 0.008251)

# p=8, M=1, d=2 (note: very bad fit)
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
  expect_equal(loglikelihood_int(data=gdpdef, p=1, M=1, params=theta_112, conditional=FALSE), -1065.289, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=gdpdef, p=2, M=1, params=theta_212, conditional=FALSE), -290.4883, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=gdpdef, p=1, M=2, params=theta_122, conditional=FALSE), -309.2328, tolerance=1e-2)
  expect_equal(loglikelihood_int(data=gdpdef, p=2, M=2, params=theta_222, conditional=FALSE), -241.6249, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=gdpdef, p=2, M=2, params=theta_222, conditional=TRUE), -236.5213, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=500*gdpdef, p=2, M=2, params=theta_222, conditional=FALSE), -960378.7, tolerance=1)

  expect_equal(loglikelihood_int(data=gdpdef, p=1, M=1, params=theta_112_mu, conditional=FALSE, parametrization="mean"), -1065.289, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=gdpdef, p=2, M=1, params=theta_212_mu, conditional=FALSE, parametrization="mean"), -290.4883, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=gdpdef, p=1, M=2, params=theta_122_mu, conditional=FALSE, parametrization="mean"), -309.2328, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=gdpdef, p=2, M=2, params=theta_222_mu, conditional=FALSE, parametrization="mean"), -241.6249, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=gdpdef, p=2, M=2, params=theta_222_mu, conditional=TRUE, parametrization="mean"), -236.5213, tolerance=1e-3)

  expect_equal(loglikelihood_int(data=gdpdef, p=1, M=4, params=theta_142, conditional=FALSE, parametrization="intercept"), -226.3188, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=gdpdef, p=1, M=4, params=theta_142, conditional=TRUE, parametrization="intercept"), -223.825, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=gdpdef, p=1, M=4, params=theta_142, conditional=TRUE, parametrization="mean"), -355.9515, tolerance=1e-3)

  # p*d >= 12
  expect_equal(loglikelihood_int(data=data2, p=4, M=1, params=theta_413, conditional=FALSE), -73586.24, tolerance=1e-1)
  expect_equal(loglikelihood_int(data=data2, p=3, M=2, params=theta_323, conditional=FALSE), -536.6736, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=gdpdef, p=7, M=2, params=theta_722, conditional=FALSE), -202.0608, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=gdpdef, p=6, M=3, params=theta_632, conditional=FALSE), -179.4021, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=gdpdef, p=8, M=1, params=theta_812, conditional=FALSE), -15542.64, tolerance=1e-1)

  # SGMVAR
  expect_equal(loglikelihood_int(data=gdpdef, p=1, M=1, params=theta_112sWC, structural_pars=list(W=W_112), conditional=FALSE),
               -1065.289, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=gdpdef, p=2, M=1, params=theta_212sWC, structural_pars=list(W=W_212), conditional=FALSE),
               -290.4883, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=gdpdef, p=1, M=2, params=theta_122s, structural_pars=list(W=W_122), conditional=FALSE),
               -309.232, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=gdpdef, p=2, M=2, params=theta_222s, structural_pars=list(W=W_222), conditional=FALSE),
               -241.6249, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=gdpdef, p=2, M=2, params=theta_222s, structural_pars=list(W=W_222), conditional=TRUE),
               -236.5213, tolerance=1e-3)

  expect_equal(loglikelihood_int(data=gdpdef, p=1, M=1, params=theta_112sWC_mu, structural_pars=list(W=W_112), conditional=FALSE, parametrization="mean"),
               -1065.289, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=gdpdef, p=2, M=2, params=theta_222s_mu, structural_pars=list(W=W_222), conditional=FALSE, parametrization="mean"),
               -241.6249, tolerance=1e-3)

  # same_means
  expect_equal(loglikelihood_int(data=gdpdef, p=1, M=1, params=theta_112_int, parametrization="mean", same_means=list(1), conditional=FALSE),
               -1065.289, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=gdpdef, p=2, M=1, params=theta_212_int, parametrization="mean", same_means=list(1), conditional=FALSE),
               -290.4883, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=gdpdef, p=1, M=2, params=theta_122_int, parametrization="mean", same_means=list(1:2), conditional=FALSE),
               -317.2695, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=gdpdef, p=1, M=2, params=theta_122_int2, parametrization="mean", same_means=list(1, 2), conditional=FALSE),
               -309.2328, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=gdpdef, p=2, M=2, params=theta_222_int, parametrization="mean", same_means=list(1:2), conditional=FALSE),
               -245.8273, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=gdpdef, p=1, M=1, params=theta_112sWC_int, parametrization="mean", structural_pars=list(W=W_112), same_means=list(1), conditional=FALSE),
               -1065.289, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=gdpdef, p=1, M=2, params=theta_122s_int, parametrization="mean", structural_pars=list(W=W_122), same_means=list(1:2), conditional=FALSE),
               -328.9261, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=gdpdef, p=2, M=2, params=theta_222s_int, parametrization="mean", structural_pars=list(W=W_222), same_means=list(1:2), conditional=FALSE),
               -245.8273, tolerance=1e-3)
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
A21_222_c2 <- A11_222_c2 <- matrix(c(0.20, 0, 0, 60), nrow=2, byrow=FALSE)
A22_222_c2 <- A12_222_c2 <- matrix(c(0.24, 0, 0, 0.24), nrow=2, byrow=FALSE)
phi10_222_c2 <- c(0.34, 0.18)
phi20_222_c2 <- c(0.47, 0.10)
Omega1_222_c2 <- matrix(c(1.08, -0.02, -0.02, 0.11), nrow=2, byrow=FALSE)
Omega2_222_c2 <- matrix(c(0.22, 0.01, 0.01, 0.03), nrow=2, byrow=FALSE)
alpha1_222_c2 <- 0.35
theta_222_c2 <- c(phi10_222_c2, phi20_222_c2, 0.26, 0.03, -0.01, -0.50, vech(Omega1_222_c2),
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
  expect_equal(loglikelihood_int(data=gdpdef, p=1, M=1, params=theta_112c, conditional=FALSE, constraints=C_112), -1065.289, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=gdpdef, p=2, M=1, params=theta_212c, conditional=FALSE, constraints=C_212), -290.4883, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=gdpdef, p=1, M=2, params=theta_122c, conditional=TRUE, constraints=C_122), -306.7991, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=gdpdef, p=2, M=2, params=theta_222c, conditional=FALSE, constraints=C_222), -295.6763, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=gdpdef, p=2, M=2, params=theta_222_c2, conditional=TRUE, constraints=C_222_2), -2168.156, tolerance=1e-2)

  # SGMVAR
  expect_equal(loglikelihood_int(data=gdpdef, p=1, M=1, params=theta_112csWAR, conditional=FALSE, constraints=C_112,
                                 structural_pars=list(W=W_112)), -1065.289, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=gdpdef, p=2, M=1, params=theta_212csWAR, conditional=FALSE, constraints=C_212,
                                 structural_pars=list(W=W_212)), -290.4883, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=gdpdef, p=1, M=2, params=theta_122csLAR, conditional=TRUE, constraints=C_122,
                                 structural_pars=list(W=W_122, C_lambda=C_lambda_122)), -310.5903, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=gdpdef, p=2, M=2, params=theta_222csLAR, conditional=FALSE, constraints=C_222,
                                 structural_pars=list(W=W_222, C_lambda=C_lambda_222)), -602.0334, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=gdpdef, p=2, M=2, params=theta_222_c2s, conditional=TRUE, constraints=C_222_2,
                                 structural_pars=list(W=W_222c2)), -472.9432, tolerance=1e-2)
  expect_equal(loglikelihood_int(data=gdpdef, p=2, M=2, params=theta_222c2s_mu, conditional=TRUE, constraints=C_222_2,
                                 structural_pars=list(W=W_222c2), parametrization="mean"), -472.9432, tolerance=1e-2)

  # Same means
  expect_equal(loglikelihood_int(data=gdpdef, p=1, M=1, params=theta_112c_int, same_means=list(1), conditional=FALSE, constraints=C_112, parametrization="mean"),
               -1065.289, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=gdpdef, p=1, M=2, params=theta_122c_int, same_means=list(1, 2), conditional=FALSE, constraints=C_122, parametrization="mean"),
               -309.2328, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=gdpdef, p=2, M=2, params=theta_222c_int, same_means=list(1:2), conditional=FALSE, constraints=C_222, parametrization="mean"),
               -284.3073, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=gdpdef, p=1, M=1, params=theta_112csWAR_int, same_means=list(1), conditional=FALSE, constraints=C_112,
                                 structural_pars=list(W=W_112), parametrization="mean"),
               -1065.2894, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=gdpdef, p=1, M=2, params=theta_122csLAR_int, same_means=list(1:2), conditional=FALSE, constraints=C_122,
                                 structural_pars=list(W=W_122, C_lambda=C_lambda_122), parametrization="mean"),
               -318.24, tolerance=1e-2)
  expect_equal(loglikelihood_int(data=200*gdpdef, p=1, M=2, params=theta_122csLAR_int, same_means=list(1:2), conditional=FALSE, constraints=C_122,
                                 structural_pars=list(W=W_122, C_lambda=C_lambda_122), parametrization="mean"),
               -280885.9, tolerance=1e-1)
  expect_equal(loglikelihood_int(data=gdpdef, p=2, M=2, params=theta_222c2s_int, same_means=list(1:2), conditional=FALSE, constraints=C_222_2,
                                 structural_pars=list(W=W_222c2), parametrization="mean"),
               -490.3359, tolerance=1e-3)
})


test_that("user loglikelihood works correctly", {
  expect_equal(loglikelihood(data=gdpdef, p=1, M=1, params=theta_112, parametrization="intercept", conditional=FALSE), -1065.289, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222, parametrization="intercept", conditional=TRUE), -236.5213, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=1, M=2, params=theta_122_mu, conditional=FALSE, parametrization="mean"), -309.2328, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222c, conditional=FALSE, parametrization="intercept", constraints=C_222), -295.6763, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=gdpdef, p=2, M=2, params=theta_222_c2, conditional=TRUE, parametrization="intercept", constraints=C_222_2), -2168.156, tolerance=1e-3)

  # SGMVAR
  expect_equal(loglikelihood(data=gdpdef, p=1, M=1, params=theta_112sWC, parametrization="intercept", conditional=FALSE,
                             structural_pars=list(W=W_112)), -1065.289, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222s, parametrization="intercept", conditional=TRUE,
                             structural_pars=list(W=W_222)), -236.5213, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=1, M=1, params=theta_112sWC_mu, conditional=FALSE, parametrization="mean",
                             structural_pars=list(W=W_112)), -1065.289, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222csLAR, conditional=FALSE, parametrization="intercept", constraints=C_222,
                             structural_pars=list(W=W_222, C_lambda=C_lambda_222)), -602.0334, tolerance=1e-3)
  expect_equal(loglikelihood_int(data=gdpdef, p=2, M=2, params=theta_222_c2s, conditional=TRUE, parametrization="intercept",
                                 constraints=C_222_2, structural_pars=list(W=W_222c2)), -472.9432, tolerance=1e-3)

  # Same means
  expect_equal(loglikelihood(data=gdpdef, p=1, M=1, params=theta_112c_int, same_means=list(1), conditional=FALSE, constraints=C_112, parametrization="mean"),
               -1065.289, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=1, M=2, params=theta_122c_int, same_means=list(1, 2), conditional=FALSE, constraints=C_122, parametrization="mean"),
               -309.2328, tolerance=1e-2)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222c_int, same_means=list(1:2), conditional=FALSE, constraints=C_222, parametrization="mean"),
               -284.3073, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=1, M=1, params=theta_112csWAR_int, same_means=list(1), conditional=FALSE, constraints=C_112,
                             structural_pars=list(W=W_112), parametrization="mean"),
               -1065.289, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=1, M=2, params=theta_122csLAR_int, same_means=list(1:2), conditional=FALSE, constraints=C_122,
                                 structural_pars=list(W=W_122, C_lambda=C_lambda_122), parametrization="mean"),
               -318.24, tolerance=1e-2)
  expect_equal(loglikelihood(data=200*gdpdef, p=1, M=2, params=theta_122csLAR_int, same_means=list(1:2), conditional=FALSE, constraints=C_122,
                                 structural_pars=list(W=W_122, C_lambda=C_lambda_122), parametrization="mean"),
               -280885.9, tolerance=1e-1)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222c2s_int, same_means=list(1:2), conditional=FALSE, constraints=C_222_2,
                                 structural_pars=list(W=W_222c2), parametrization="mean"),
               -490.3359, tolerance=1e-3)
  expect_error(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222c2s_int, same_means=list(1:2), conditional=FALSE, constraints=C_222_2,
                             structural_pars=list(W=W_222c2), parametrization="intercept"))

  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222_int, parametrization="mean", same_means=list(1:2), conditional=FALSE),
               -245.8273, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222s_int, parametrization="mean", structural_pars=list(W=W_222),
                             same_means=list(1:2), conditional=FALSE),
               -245.8273, tolerance=1e-3)

})


test_that("cond_moments works correctly", {
  expect_equal(cond_moments(data=gdpdef, p=1, M=1, params=theta_112, parametrization="intercept", to_return="total_ccovs")[, 2, 20], c(0.01, 0.07), tolerance=1e-2)
  expect_equal(cond_moments(data=gdpdef, p=1, M=1, params=theta_112, parametrization="intercept", to_return="regime_cmeans")[1, , 1], c(1.174027, 0.948235), tolerance=1e-5)
  expect_equal(cond_moments(data=gdpdef, p=1, M=1, params=theta_112, parametrization="intercept", to_return="total_cmeans")[100, ], c(1.148798, 1.424473), tolerance=1e-5)
  expect_equal(cond_moments(data=gdpdef, p=2, M=2, params=theta_222, parametrization="intercept", to_return="total_ccovs")[1, , 200], c(1.099995625, 0.009998022), tolerance=1e-5)
  expect_equal(cond_moments(data=gdpdef, p=2, M=2, params=theta_222, parametrization="intercept", to_return="total_cmeans")[1, ], c(1.2107632, 0.3326596), tolerance=1e-5)
  expect_equal(cond_moments(data=gdpdef, p=1, M=2, params=theta_122_mu, parametrization="mean", to_return="total_ccovs")[, 2, 100], c(0.005306351, 0.070671908), tolerance=1e-5)
  expect_equal(cond_moments(data=gdpdef, p=2, M=2, params=theta_222c, parametrization="intercept", constraints=C_222, to_return="total_cmeans")[13, ], c(1.3920590, 0.3398279), tolerance=1e-5)
  expect_equal(cond_moments(data=gdpdef, p=2, M=2, params=theta_222_c2, parametrization="intercept", constraints=C_222_2, to_return="total_ccovs")[ , 2, 13], c(-0.008241019, 0.073319271), tolerance=1e-5)

  # SGMVAR
  expect_equal(cond_moments(data=gdpdef, p=1, M=1, params=theta_112sWC, parametrization="intercept",
                            structural_pars=list(W=W_112), to_return="total_ccovs")[, 2, 20], c(0.01, 0.07), tolerance=1e-2)
  expect_equal(cond_moments(data=gdpdef, p=1, M=1, params=theta_112sWC, parametrization="intercept",
                            structural_pars=list(W=W_112), to_return="regime_cmeans")[1, , 1], c(1.174027, 0.948235), tolerance=1e-5)
  expect_equal(cond_moments(data=gdpdef, p=1, M=1, params=theta_112sWC, parametrization="intercept",
                            structural_pars=list(W=W_112), to_return="total_cmeans")[100, ], c( 1.148798, 1.424473), tolerance=1e-5)
  expect_equal(cond_moments(data=gdpdef, p=2, M=2, params=theta_222s, parametrization="intercept",
                            structural_pars=list(W=W_222), to_return="total_ccovs")[1, , 200], c(1.099995625, 0.009998022), tolerance=1e-5)
  expect_equal(cond_moments(data=gdpdef, p=2, M=2, params=theta_222s, parametrization="intercept",
                            structural_pars=list(W=W_222), to_return="total_cmeans")[1, ], c(1.2107632, 0.3326596), tolerance=1e-5)
  expect_equal(cond_moments(data=gdpdef, p=2, M=2, params=theta_222s_mu, parametrization="mean",
                            structural_pars=list(W=W_222), to_return="total_ccovs")[, 2, 100], c(-0.01783773, 0.08487145), tolerance=1e-5)
  expect_equal(cond_moments(data=gdpdef, p=2, M=2, params=theta_222csLAR, parametrization="intercept",
                            constraints=C_222, structural_pars=list(W=W_222, C_lambda=C_lambda_222), to_return="total_cmeans")[13, ],
               c(1.3460926, 0.3589806), tolerance=1e-5)
  expect_equal(cond_moments(data=gdpdef, p=2, M=2, params=theta_222_c2s, parametrization="intercept",
                            constraints=C_222_2, structural_pars=list(W=W_222c2), to_return="total_ccovs")[ , 2, 13],
               c(-0.01999998, 0.10999991), tolerance=1e-5)

  # Same means
  expect_equal(cond_moments(data=gdpdef, p=1, M=1, params=theta_112c_int, same_means=list(1), constraints=C_112, parametrization="mean",
                            to_return="regime_cmeans")[2, , 1],
               c(1.2825737, 0.8828394), tolerance=1e-5)
  expect_equal(cond_moments(data=gdpdef, p=1, M=2, params=theta_122c_int, same_means=list(1, 2), constraints=C_122, parametrization="mean",
                            to_return="total_ccovs")[2, , 200],
               c(-0.01446378, 0.14487362), tolerance=1e-5)
  expect_equal(cond_moments(data=gdpdef, p=2, M=2, params=theta_222c_int, same_means=list(1:2), constraints=C_222, parametrization="mean",
                            to_return="total_cmeans")[1,],
               c(1.6388497, 0.3657483), tolerance=1e-3)
  expect_equal(cond_moments(data=gdpdef, p=1, M=1, params=theta_112csWAR_int, same_means=list(1), constraints=C_112, structural_pars=list(W=W_112), parametrization="mean",
                            to_return="total_cmeans")[1, ],
               c(1.174027, 0.948235), tolerance=1e-4)
  expect_equal(cond_moments(data=gdpdef, p=1, M=2, params=theta_122csLAR_int, same_means=list(1:2), constraints=C_122,
                             structural_pars=list(W=W_122, C_lambda=C_lambda_122), parametrization="mean", to_return="total_ccovs")[1, , 100],
               c(0.5812121401, 0.0003112651), tolerance=1e-5)
  expect_equal(cond_moments(data=200*gdpdef, p=1, M=2, params=theta_122csLAR_int, same_means=list(1:2), constraints=C_122,
                             structural_pars=list(W=W_122, C_lambda=C_lambda_122), parametrization="mean", to_return="regime_cmeans")[1, , 2],
               c(130.61308, 52.88054), tolerance=1e-4)
  expect_equal(cond_moments(data=gdpdef, p=2, M=2, params=theta_222c2s_int, same_means=list(1:2), constraints=C_222_2,
                             structural_pars=list(W=W_222c2), parametrization="mean", to_return="total_ccovs")[, 2, 13],
               c(-0.02, 0.11), tolerance=1e-2)
  expect_equal(cond_moments(data=gdpdef, p=2, M=2, params=theta_222_int, parametrization="mean", same_means=list(1:2), to_return="total_cmeans")[13,],
               c(1.1300052, 0.3267675), tolerance=1e-4)
  expect_equal(cond_moments(data=gdpdef, p=2, M=2, params=theta_222s_int, parametrization="mean", structural_pars=list(W=W_222), same_means=list(1:2),
                            to_return="total_ccovs")[2, , 1],
               c(0.01160004, 0.10666026), tolerance=1e-4)
})


