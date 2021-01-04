context("quantile residuals")
library(gmvarkit)

# NOTE that some test use random elements that might change if random_ind2 or simulateGMVAR is modified

data <- cbind(10*eurusd[,1], 100*eurusd[,2])

## A(M)(p)_(p)(M)(d)

# p=1, M=1, d=2, parametrization="mean"
phi10_112 <- c(1.07, 127.71)
A11_112 <- matrix(c(0.99, 0.00, -0.01, 0.99), nrow=2)
Omega1_112 <- matrix(c(4.05, 2.22, 2.22, 2.22), nrow=2)
theta_112 <- c(phi10_112, vec(A11_112), vech(Omega1_112))
mod_112 <- GMVAR(data, p=1, M=1, d=2, params=theta_112, conditional=TRUE, parametrization="mean", constraints=NULL)

W_112 <- t(chol(Omega1_112))
theta_112sWC <- c(phi10_112, vec(A11_112), Wvec(W_112)) # SGMVAR, W constrained by Cholesky
mod_112s <- GMVAR(data, p=1, M=1, d=2, params=theta_112sWC, conditional=TRUE, parametrization="mean", constraints=NULL,
                  structural_pars=list(W=W_112))

# p=2, M=2, d=2, no constraints, GMVAR-paper
phi10_222 <- c(1.03, 2.36)
A11_222 <- matrix(c(1.25, 0.06, 0.04, 1.34), nrow=2, byrow=FALSE)
A12_222 <- matrix(c(-0.29, -0.08, -0.05, -0.36), nrow=2, byrow=FALSE)
Omega1_222 <- matrix(c(0.93, -0.15, -0.15, 5.20), nrow=2, byrow=FALSE)

phi20_222 <- c(1.79, 3.00)
A21_222 <- A11_222; A22_222 <- A12_222
Omega2_222 <- matrix(c(5.88, 3.56, 3.56, 9.80), nrow=2, byrow=FALSE)

alpha1_222 <- 0.37

upsilon1_222 <- c(phi10_222, vec(A11_222), vec(A12_222), vech(Omega1_222))
upsilon2_222 <- c(phi20_222, vec(A21_222), vec(A22_222), vech(Omega2_222))
theta_222 <- c(upsilon1_222, upsilon2_222, alpha1_222)
mod_222 <- GMVAR(data, p=2, M=2, d=2, params=theta_222, conditional=TRUE, parametrization="intercept", constraints=NULL)

WL_222 <- diag_Omegas(Omega1_222, Omega2_222)
W_222 <- matrix(WL_222[1:(2^2)], nrow=2, byrow=FALSE)
lambdas_222 <- WL_222[(2^2 + 1):length(WL_222)]
theta_222s <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A21_222),
                vec(A22_222), vec(W_222), lambdas_222, alpha1_222) # SGMVAR
mod_222s <- GMVAR(data, p=2, M=2, d=2, params=theta_222s, conditional=TRUE, parametrization="intercept", constraints=NULL,
                 structural_pars=list(W=W_222))

# p=1, M=4, d=2
theta_142 <- c(22.743598, 98.461525, -0.006778, 0.212983, -0.12094, -0.061364,
               1.313218, -3.305944, 14.459664, 44.482964, 82.615109, -0.344603,
               0.082359, -0.168216, 0.341988, 4.012248, -3.555077, 6.449586,
               14.672977, 100.177417, 0.345826, -0.338611, -0.18687, -0.221931,
               12.086199, 3.677137, 1.331716, 17.668096, 129.042416, 0.628735,
               -0.026376, 0.185899, -0.199485, 0.470336, 0.980442, 7.146605,
               0.427396, 0.417413, 0.142379)
mod_142_int <- GMVAR(data, p=1, M=4, d=2, params=theta_142, conditional=TRUE, parametrization="intercept")
mod_142_mean <- GMVAR(data, p=1, M=4, d=2, params=theta_142, conditional=TRUE, parametrization="mean")


# p=2, M=2, d=2, SGMVAR AR params constrained to be the same in both regimes
rbind_diags <- function(p, M, d) {
  I <- diag(p*d^2)
  Reduce(rbind, replicate(M, I, simplify=FALSE))
}
C_222 <- rbind_diags(p=2, M=2, d=2)
C_lambda_222 <- matrix(c(7, 1), nrow=2)
theta_222cs <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(W_222), 1, alpha1_222)
mod_222csLAR <- GMVAR(data, p=2, M=2, d=2, params=theta_222cs, conditional=TRUE, parametrization="intercept",
                   constraints=C_222, structural_pars=list(W=W_222, C_lambda=C_lambda_222))

# p=2, M=2, d=2, AR parameters same for both regimes, non-diagonals zero, intercept
theta_222c <- c(0.3552775, 3.1929675, -0.1143198, 2.8294743, 1.2633425, 1.3375150, -0.2919742, -0.3624010,
                5.5971764, 3.4559442, 9.6221422, 0.9820759, -0.3267521, 5.2358855, 0.6501600)
mat0 <- matrix(c(1, rep(0, 10), 1, rep(0, 8), 1, rep(0, 10), 1), nrow=2*2^2, byrow=FALSE)
C_222c <- rbind(mat0, mat0)
mod_222c <- GMVAR(data, p=2, M=2, d=2, params=theta_222c, conditional=TRUE, parametrization="intercept", constraints=C_222c)


# p=3, M=2, d=3, no constraints, rand_ind and simulated data
set.seed(42)
theta_323 <- random_ind2(p=3, M=2, d=3, mu_scale=c(-10, 0, 5), mu_scale2=1:3, omega_scale=1:3, ar_scale=1)
mod_323 <- GMVAR(p=3, M=2, d=3, params=theta_323, conditional=FALSE, parametrization="mean", constraints=NULL)
sim_323 <- simulateGMVAR(mod_323, nsimu=500)$sample
mod_323 <- add_data(data=sim_323, gmvar=mod_323)

# p=2, M=2, d=2, constraints=C_222, structural_pars=list(W=W_222, C_lambda=C_lambda_222), same_means=list(1:2)
theta_222csLAR_int <- c(phi10_222, vec(A11_222), vec(A12_222), vec(W_222), 0.2, alpha1_222)
mod_222csLAR_int <- GMVAR(data, p=2, M=2, params=theta_222csLAR_int, conditional=TRUE, parametrization="mean",
                          constraints=C_222, structural_pars=list(W=W_222, C_lambda=C_lambda_222), same_means=list(1:2))

# Get quantile residuals
res_112 <- quantile_residuals(mod_112)
res_222 <- quantile_residuals(mod_222)
res_222c <- quantile_residuals(mod_222c)
res_323 <- quantile_residuals(mod_323)

res_142_int <- quantile_residuals(mod_142_int)
res_142_mean <- quantile_residuals(mod_142_mean)

res_112s <- quantile_residuals(mod_112s)
res_222s <- quantile_residuals(mod_222s)
res_222csLAR <- quantile_residuals(mod_222csLAR)

res222csLAR_int <- quantile_residuals(mod_222csLAR_int)


test_that("quantile_residuals works correctly", {
  expect_equal(res_112[1,], c(0.7421619, -2.3382223), tolerance=1e-6)
  expect_equal(res_112[2,], c(0.1764158, -1.3046112), tolerance=1e-6)
  expect_equal(res_112[100,], c(-0.2447183, 0.6045296), tolerance=1e-6)
  expect_equal(res_112[251,], c(-0.6408632, -2.0843120), tolerance=1e-6)

  expect_equal(res_222[1:2, 1], c(0.1694553, 1.1810660), tolerance=1e-6)
  expect_equal(res_222[100:102, 1], c(1.2474611, -0.2844403, 0.3318623), tolerance=1e-6)
  expect_equal(res_222[212:215, 1], c(0.1436343, 0.3750201, 0.9942028, -0.5531501), tolerance=1e-6)
  expect_equal(res_222[1:2, 2], c(-0.277522646, -0.000144752), tolerance=1e-6)
  expect_equal(res_222[100:102, 2], c(-0.628521, -1.309668, -1.053083), tolerance=1e-6)
  expect_equal(res_222[212:215, 2], c(-0.3754761, 1.3751472, 1.0785548, -1.4534341), tolerance=1e-6)

  expect_equal(res_222c[1,], c(0.09074664, -0.19603850), tolerance=1e-6)
  expect_equal(res_222c[100,], c(1.2506933, -0.4634526), tolerance=1e-6)

  expect_equal(res_323[1,], c(-0.04254252, 0.37966523, -1.43671034), tolerance=1e-6)
  expect_equal(res_323[13,], c(-1.1160317, -0.1883617, -0.4064072), tolerance=1e-6)
  expect_equal(res_323[150,], c(0.09401363, 0.44800577, 1.11962898), tolerance=1e-6)
  expect_equal(res_323[497,], c(0.3658767, 0.4336302, -0.6696751), tolerance=1e-6)

  expect_equal(res_142_int[122:125, 1], c(-0.1932272, -0.6044207, -1.1737746, -1.1299033), tolerance=1e-6)
  expect_equal(res_142_int[12:15, 2], c(0.3681731, 1.7120364, 0.7318820, -0.1006260), tolerance=1e-6)
  expect_equal(res_142_mean[132:135, 1], c(-4.506558, -4.203856, -4.220869, -3.854244), tolerance=1e-6)
  expect_equal(res_142_int[2:5, 2], c(1.612173, 2.645782, 4.811656, 8.014016), tolerance=1e-6)


  # SGMVAR
  expect_equal(res_112s[3,], c(0.5815407, -1.0931660), tol=1e-6)
  expect_equal(res_112s, res_112, tol=1e-6)
  expect_equal(res_222s, res_222, tol=1e-6)
  expect_equal(res_222csLAR[1,], c(0.1658934, -0.2827984), tol=1e-6)
  expect_equal(res_222csLAR[250,], c(-0.7458269, -0.8943336), tol=1e-6)

  # Same_means
  expect_equal(res222csLAR_int[c(1, 5, 100, 200)], c(1.2904864, 2.6818951, 2.3194764, 0.1348595), tolerance=1e-6)
})

test_that("quantile_residuals_int works correctly", {
  qr112 <- quantile_residuals_int(data, p=1, M=1, params=theta_112, conditional=TRUE, parametrization="mean", constraints=NULL)
  qr222csLAR <- quantile_residuals_int(data, p=2, M=2, params=theta_222cs, conditional=TRUE, parametrization="intercept",
                                       constraints=C_222, structural_pars=list(W=W_222, C_lambda=C_lambda_222))
  qr222csLAR_int <- quantile_residuals_int(data, p=2, M=2, params=theta_222csLAR_int, conditional=TRUE, parametrization="mean",
                                           constraints=C_222, structural_pars=list(W=W_222, C_lambda=C_lambda_222), same_means=list(1:2))
  expect_equal(qr112, res_112, tol=1e-6)
  expect_equal(qr222csLAR, res_222csLAR, tol=1e-6)
  expect_equal(qr222csLAR_int[1:3], c(1.290486, 2.472298, 5.088865), tolerance=1e-6)
})








