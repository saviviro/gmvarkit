context("quantile residuals")
library(gmvarkit)

# NOTE that some test use random elements that might change if random_ind2 or simulateGMVAR is modified

data <- cbind(10*eurusd[,1], 100*eurusd[,2])

## A(M)(p)_(p)(M)(d)

# p=1, M=1, d=2, parametrization="mean"
phi10_112 <- c(1.07, 127.71)
A11_112 <- matrix(c(0.99, 0.00, -0.01, 0.99), nrow=2)
Omega1_112 <- matrix(c(4.05, 2.22, 8.87, 2.22), nrow=2)
theta_112 <- c(phi10_112, vec(A11_112), vech(Omega1_112))
mod_112 <- GMVAR(data, p=1, M=1, d=2, params=theta_112, conditional=TRUE, parametrization="mean", constraints=NULL)

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

# p=2, M=2, d=2, AR paramars same, non-diagonals zero, intercept
theta_222c <- c(0.3552775, 3.1929675, -0.1143198, 2.8294743, 1.2633425, 1.3375150, -0.2919742, -0.3624010,
                5.5971764, 3.4559442, 9.6221422, 0.9820759, -0.3267521, 5.2358855, 0.6501600)
mat0 <- matrix(c(1, rep(0, 10), 1, rep(0, 8), 1, rep(0, 10), 1), nrow=2*2^2, byrow=FALSE) # Laske paperilla monimutkaisemmat rajoitteet
C_222c <- rbind(mat0, mat0)
mod_222c <- GMVAR(data, p=2, M=2, d=2, params=theta_222c, conditional=TRUE, parametrization="intercept", constraints=C_222c)


# p=3, M=2, d=3, no constraints, rand_ind and simulated data
set.seed(42)
theta_323 <- random_ind2(p=3, M=2, d=3, mu_scale=c(-10, 0, 5), mu_scale2=1:3, omega_scale=1:3, ar_scale=1)
mod_323 <- GMVAR(p=3, M=2, d=3, params=theta_323, conditional=FALSE, parametrization="mean", constraints=NULL)
sim_323 <- simulateGMVAR(mod_323, nsimu=500)$sample
mod_323 <- add_data(data=sim_323, gmvar=mod_323)

res_112 <- quantile_residuals(mod_112)
res_222 <- quantile_residuals(mod_222)
res_222c <- quantile_residuals(mod_222c)
res_323 <- quantile_residuals(mod_323)

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

  expect_equal(res_323[1,], c(-0.5599706, 0.7861472, -1.1307067), tolerance=1e-6)
  expect_equal(res_323[13,], c(-1.0996495, -0.4214492, -0.2435677), tolerance=1e-6)
  expect_equal(res_323[150,], c(-0.6546153, 0.9539610, 0.3529515), tolerance=1e-6)
  expect_equal(res_323[497,], c(0.6395263, 0.1111915, -0.5907638), tolerance=1e-6)
})
