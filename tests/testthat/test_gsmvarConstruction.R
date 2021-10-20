context("GSMVAR construction")
library(gmvarkit)

## A(M)(p)_(p)(M)(d)
rbind_diags <- function(p, M, d) {
  I <- diag(p*d^2)
  Reduce(rbind, replicate(M, I, simplify=FALSE))
}

# p=1, M=2, d=2
params122 <- c(0.623, -0.129, 0.959, 0.089, -0.006, 1.006, 1.746, 0.804, 5.804, 3.245, 7.913,
               0.952, -0.037, -0.019, 0.943, 6.926, 3.982, 12.135, 0.789) # p=1, M=2, d=2

params122t <- c(params122, 10, 20) # StMVAR
params122gs <- c(params122, 20) # G-StMVAR

# p=2, M=2, d=2, constrained
C_mat <- rbind(diag(2*2^2), diag(2*2^2))
params222c <- c(1.031, 2.356, 1.786, 3.000, 1.250, 0.060, 0.036, 1.335, -0.290, -0.083, -0.047,
                -0.356, 0.934, -0.152, 5.201, 5.883, 3.560, 9.799, 0.368)  # p=2, M=2, d=2, AR parameters restricted same for both regimes

params222gsc <- c(params222c, 20) # G-StMVAR

# p=2, M=1, d=2, SGMVAR, W constrained
phi10_212 <- c(1.03, 2.36)
A11_212 <- matrix(c(1.25, 0.06, 0.04, 1.34), nrow=2, byrow=FALSE)
A12_212 <- matrix(c(-0.29, -0.08, -0.05, -0.36), nrow=2, byrow=FALSE)
Omega1_212 <- matrix(c(0.93, -0.15, -0.15, 5.20), nrow=2, byrow=FALSE)

W_212 <- t(chol(Omega1_212))
theta_212csW <- c(phi10_212, vec(A11_212), vec(A12_212), Wvec(W_212))

theta_212tcsW <- c(theta_212csW, 10) # SStMVAR

# p=1, M=2, d=2, SGMVAR, AR parameters and lambdas constrained
phi10_122 <- c(1.03, 2.36)
A11_122 <- matrix(c(0.9, 0.3, -0.3, 0.9), nrow=2, byrow=FALSE)
Omega1_122 <- matrix(c(0.93, -0.15, -0.15, 5.20), nrow=2, byrow=FALSE)

phi20_122 <- c(1.79, 3.00)
A21_122 <- A11_122
Omega2_122 <- matrix(c(5.88, 3.56, 3.56, 9.80), nrow=2, byrow=FALSE)

alpha1_122 <- 0.37

C_122 <- rbind_diags(p=1, M=2, d=2)
WL_122 <- diag_Omegas(Omega1_122, Omega2_122)
W_122 <- matrix(WL_122[1:(2^2)], nrow=2, byrow=FALSE)
C_lambda_122 <- matrix(c(7, 1), nrow=2)
theta_122csLAR <- c(phi10_122, phi20_122, vec(A11_122), vec(W_122), 1, alpha1_122)

theta_122gscsLAR <- c(theta_122csLAR, 10) # SG-StMVAR

# p=2, M=2, d=2, constraints=C_222, same_means=list(1:2)
C_222 <- rbind_diags(p=2, M=2, d=2)
phi10_222 <- c(1.03, 2.36)
A11_222 <- matrix(c(1.25, 0.06, 0.04, 1.34), nrow=2, byrow=FALSE)
A12_222 <- matrix(c(-0.29, -0.08, -0.05, -0.36), nrow=2, byrow=FALSE)
Omega1_222 <- matrix(c(0.93, -0.15, -0.15, 5.20), nrow=2, byrow=FALSE)
Omega2_222 <- matrix(c(5.88, 3.56, 3.56, 9.80), nrow=2, byrow=FALSE)
alpha1_222 <- 0.37

theta_222c_int <- c(-5, 123, vec(A11_222), vec(A12_222), vech(Omega1_222), vech(Omega2_222), alpha1_222)

theta_222gsc_int <- c(theta_222c_int, 20) # G-StMVAR

# p=2, M=2, d=2, constraints=C_222, structural_pars=list(W=W_222, C_lambda=C_lambda_222), same_means=list(1:2)
WL_222 <- diag_Omegas(Omega1_222, Omega2_222)
W_222 <- matrix(WL_222[1:(2^2)], nrow=2, byrow=FALSE)
C_lambda_222 <- matrix(c(1, 2), nrow=2)
theta_222csLAR_int <- c(phi10_222, vec(A11_222), vec(A12_222), vec(W_222), 0.2, alpha1_222)

theta_222tcsLAR_int <- c(theta_222csLAR_int, 10, 20) # SStMVAR

test_that("GSMVAR works correctly", {
  mod122 <- GSMVAR(gdpdef, p=1, M=2, params=params122)
  mod122t <- GSMVAR(gdpdef, p=1, M=2, params=params122t, model="StMVAR")
  mod122gs <- GSMVAR(gdpdef, p=1, M=c(1, 1), params=params122gs, model="G-StMVAR")
  mod222c <- GSMVAR(gdpdef, p=2, M=2, params=params222c, constraints=C_mat)
  mod222gsc <- GSMVAR(gdpdef, p=2, M=c(1, 1), params=params222gsc, model="G-StMVAR", constraints=C_mat)
  expect_equal(mod122$params, params122)
  expect_equal(mod122t$params, params122t)
  expect_equal(mod122gs$params, params122gs)
  expect_equal(mod222c$params, params222c)
  expect_equal(mod222gsc$params, params222gsc)

  # SGSMVAR
  mod212csW <- GSMVAR(gdpdef, p=2, M=1, params=theta_212csW, structural_pars=list(W=W_212))
  mod212tcsW <- GSMVAR(gdpdef, p=2, M=1, params=theta_212tcsW, model="StMVAR", structural_pars=list(W=W_212))
  mod122csLAR <- GSMVAR(gdpdef, p=1, M=2, params=theta_122csLAR, constraints=C_122,
                       structural_pars=list(W=W_122, C_lambda=C_lambda_122))
  mod122gscsLAR <- GSMVAR(gdpdef, p=1, M=c(1, 1), params=theta_122gscsLAR, model="G-StMVAR", constraints=C_122,
                        structural_pars=list(W=W_122, C_lambda=C_lambda_122))
  expect_equal(mod212csW$params, theta_212csW)
  expect_equal(mod212tcsW$params, theta_212tcsW)
  expect_equal(mod122csLAR$params, theta_122csLAR)
  expect_equal(mod122gscsLAR$params, theta_122gscsLAR)

  # Same_means
  mod222c_int <- GSMVAR(gdpdef, p=2, M=2, params=theta_222c_int, parametrization="mean", constraints=C_222, same_means=list(1:2))
  mod222gsc_int <- GSMVAR(gdpdef, p=2, M=c(1, 1), params=theta_222gsc_int, model="G-StMVAR", parametrization="mean", constraints=C_222, same_means=list(1:2))
  mod222csLAR_int <- GSMVAR(gdpdef, p=2, M=2, params=theta_222csLAR_int, parametrization="mean", constraints=C_222,
                           structural_pars=list(W=W_222, C_lambda=C_lambda_222), same_means=list(1:2))
  mod222tcsLAR_int <- GSMVAR(gdpdef, p=2, M=2, params=theta_222tcsLAR_int, model="StMVAR", parametrization="mean", constraints=C_222,
                            structural_pars=list(W=W_222, C_lambda=C_lambda_222), same_means=list(1:2))
  expect_equal(mod222c_int$params, theta_222c_int)
  expect_equal(mod222gsc_int$params, theta_222gsc_int)
  expect_equal(mod222tcsLAR_int$params, theta_222tcsLAR_int)
})


test_that("add_data works correctly", {
  mod122 <- GSMVAR(p=1, M=2, d=2, params=params122)
  mod122_2 <- add_data(data=gdpdef, gsmvar=mod122)
  expect_equal(mod122_2$data, gdpdef)

  # SGSMVAR
  mod122csLAR <- GSMVAR(p=1, M=2, d=2, params=theta_122csLAR, constraints=C_122,
                       structural_pars=list(W=W_122, C_lambda=C_lambda_122))
  mod122csLAR_2 <- add_data(data=gdpdef, gsmvar=mod122csLAR)
  expect_equal(mod122csLAR_2$data, gdpdef)

  mod122gscsLAR <- GSMVAR(p=1, M=c(1, 1), d=2, params=theta_122gscsLAR,
                          model="G-StMVAR", constraints=C_122,
                        structural_pars=list(W=W_122, C_lambda=C_lambda_122))
  mod122gscsLAR_2 <- add_data(data=gdpdef, gsmvar=mod122gscsLAR)
  expect_equal(mod122gscsLAR_2$data, gdpdef)

  # Same means
  mod222c_int <- GSMVAR(p=2, M=2, d=2, params=theta_222c_int, parametrization="mean", constraints=C_222, same_means=list(1:2))
  mod222c_int2 <- add_data(data=gdpdef, gsmvar=mod222c_int)
  expect_equal(mod222c_int2$data, gdpdef)

  mod222tc_int <- GSMVAR(p=2, M=2, d=2, params=c(theta_222c_int, 3, 4), model="StMVAR",
                         parametrization="mean", constraints=C_222, same_means=list(1:2))
  mod222tc_int2 <- add_data(data=gdpdef, gsmvar=mod222tc_int)
  expect_equal(mod222tc_int2$data, gdpdef)
})


test_that("swap_parametrization works correctly", {
  mod122 <- GSMVAR(gdpdef, p=1, M=2, params=params122)
  mod222c <- GSMVAR(gdpdef, p=2, M=2, params=params222c, constraints=C_mat)
  mod222tc <- GSMVAR(gdpdef, p=2, M=2, params=c(params222c, 10, 20), model="StMVAR", constraints=C_mat)
  mod122_2 <- swap_parametrization(mod122)
  mod222c_2 <- swap_parametrization(mod222c)
  mod222tc_2 <- swap_parametrization(mod222tc)
  expect_equal(mod122_2$params, c(-10.291667, 174.159722, 0.959, 0.089, -0.006, 1.006, 1.746, 0.804,
                                  5.804, 17.028037, 127.771274, 0.952, -0.037, -0.019, 0.943, 6.926,
                                  3.982, 12.135, 0.789), tolerance=1e-4)
  expect_equal(mod222c_2$params, c(-7.265758, 120.148211, 7.67632, 134.449744, 1.25, 0.06, 0.036, 1.335,
                                   -0.29, -0.083, -0.047, -0.356, 0.934, -0.152, 5.201, 5.883, 3.56,
                                   9.799, 0.368), tolerance=1e-4)
  expect_equal(mod222tc_2$params, c(-7.265758, 120.148211, 7.67632, 134.449744, 1.25, 0.06, 0.036, 1.335,
                                   -0.29, -0.083, -0.047, -0.356, 0.934, -0.152, 5.201, 5.883, 3.56,
                                   9.799, 0.368, 10, 20), tolerance=1e-4)

  # SGSMVAR
  mod212csW <- GSMVAR(gdpdef, p=2, M=1, params=theta_212csW, structural_pars=list(W=W_212))
  mod122csLAR <- GSMVAR(gdpdef, p=1, M=2, params=theta_122csLAR, constraints=C_122,
                        structural_pars=list(W=W_122, C_lambda=C_lambda_122))
  mod122gscsLAR <- GSMVAR(gdpdef, p=1, M=c(1, 1), params=c(theta_122csLAR, 3), model="G-StMVAR", constraints=C_122,
                          structural_pars=list(W=W_122, C_lambda=C_lambda_122))
  mod212csW_2 <- swap_parametrization(mod212csW)
  mod122csLAR_2 <- swap_parametrization(mod122csLAR)
  mod122gscsLAR_2 <- swap_parametrization(mod122gscsLAR)
  expect_equal(mod212csW_2$params, c(-5, 123, 1.25, 0.06, 0.04, 1.34, -0.29, -0.08, -0.05, -0.36, 0.96437,
                                     -0.15554, 2.27504), tolerance=1e-4)
  expect_equal(mod122csLAR_2$params, c(-6.05, 5.45, -7.21, 8.37, 0.9, 0.3, -0.3, 0.9, -0.89246, -0.71805,
                                       0.36539, -2.16435, 1, 0.37), tolerance=1e-4)
  expect_equal(mod122gscsLAR_2$params, c(-6.05, 5.45, -7.21, 8.37, 0.9, 0.3, -0.3, 0.9, -0.89246, -0.71805,
                                          0.36539, -2.16435, 1, 0.37, 3), tolerance=1e-4)
})


test_that("alt_gsmvar does not throw errors", {
  gsmvar11 <- suppressMessages(fitGSMVAR(gdpdef[1:20,], p=1, M=1, ncalls=2, ncores=1, maxit=1, seeds=1:2, print_res=FALSE, ngen=1))
  tmp <- alt_gsmvar(gsmvar11, which_round=1, calc_cond_moments=FALSE, calc_std_errors=FALSE)
  expect_true(TRUE) # There needs to be a test inside each test_that
})

# GMVAR(1,2), d=2 model:
params122 <- c(0.623, -0.129, 0.959, 0.089, -0.006, 1.006, 1.746,
               0.804, 5.804, 3.245, 7.913, 0.952, -0.037, -0.019,
               0.943, 6.926, 3.982, 12.135, 0.789)
mod122 <- GSMVAR(gdpdef, p=1, M=2, params=params122)
mod122s <- gsmvar_to_sgsmvar(mod122)

mod122t <- GSMVAR(gdpdef, p=1, M=2, params=c(params122, 10, 20), model="StMVAR")
mod122ts <- gsmvar_to_sgsmvar(mod122t)

# GMVAR(2,2), d=2 model with AR-parameters restricted to be
# the same for both regimes:
C_mat <- rbind(diag(2*2^2), diag(2*2^2))
params222c <- c(1.031, 2.356, 1.786, 3.000, 1.250, 0.060, 0.036,
                1.335, -0.290, -0.083, -0.047, -0.356, 0.934, -0.152,
                5.201, 5.883, 3.560, 9.799, 0.368)
mod222c <- GSMVAR(p=2, M=2, d=2, params=params222c, constraints=C_mat)
mod222cs <- gsmvar_to_sgsmvar(mod222c)

mod222gsc <- GSMVAR(p=2, M=c(1, 1), d=2, params=c(params222c, 10), model="G-StMVAR", constraints=C_mat)
mod222gscs <- gsmvar_to_sgsmvar(mod222gsc)


# p=2, M=2, d=2, constraints=C_222, same_means=list(1:2)
mod222c_int <- GSMVAR(gdpdef, p=2, M=2, params=theta_222c_int, parametrization="mean", constraints=C_222, same_means=list(1:2))
mod222cs_int <- gsmvar_to_sgsmvar(mod222c_int)



test_that("gsmvar_to_sgsmvar works correctly", {
  expect_equal(mod122s$params, c(0.623, -0.129, 3.245, 7.913, 0.959, 0.089, -0.006, 1.006, 0.952, -0.037, -0.019, 0.943,
                                 1.31216, 0.87892, -0.15571, 2.2431, 3.99732, 1.79808, 0.789), tolerance=1e-3)
  expect_equal(mod122ts$params, c(0.623, -0.129, 3.245, 7.913, 0.959, 0.089, -0.006, 1.006, 0.952, -0.037, -0.019, 0.943,
                                  1.31216, 0.87892, -0.15571, 2.2431, 3.99732, 1.79808, 0.789, 10, 20), tolerance=1e-3)
  expect_equal(mod222cs$params, c(1.031, 2.356, 1.786, 3, 1.25, 0.06, 0.036, 1.335, -0.29, -0.083, -0.047, -0.356, 0.89382,
                                  0.71976, -0.36753, 2.16401, 7.14351, 1.30222, 0.368), tolerance=1e-3)
  expect_equal(mod222gscs$params, c(1.031, 2.356, 1.786, 3, 1.25, 0.06, 0.036, 1.335, -0.29, -0.083, -0.047, -0.356, 0.89382,
                                  0.71976, -0.36753, 2.16401, 7.14351, 1.30222, 0.368, 10), tolerance=1e-3)
  expect_equal(mod222cs_int$params,
               c(-5.0000000, 123.0000000, 1.2500000, 0.0600000, 0.0400000, 1.3400000, -0.2900000, -0.0800000, -0.0500000,
                 -0.3600000, 0.8924620, 0.7180539, -0.3653923, 2.1643472, 7.1638990, 1.3035363, 0.3700000), tolerance=1e-6)
})


mod222cs2 <- reorder_W_columns(mod222cs, perm=2:1)

params112 <- c(0.899784972, 1.707574085, 0.989440729, 0.002636335, -0.006957649, 0.986606661,
               1.961053552, 0.452934215, 0.452934215, 2.944004374)
W112 <- matrix(NA, nrow=2, ncol=2)
diag(W112) <- c(1, 1)
mod112 <- GSMVAR(gdpdef, p=1, M=1, params=params112, structural_pars=list(W=W112))
mod112_2 <- reorder_W_columns(mod112, perm=2:1)

mod112t <- GSMVAR(gdpdef, p=1, M=1, params=c(params112, 10), model="StMVAR", structural_pars=list(W=W112))
mod112t_2 <- reorder_W_columns(mod112t, perm=2:1)


# p=2, M=2, d=2, constraints=C_222, structural_pars=list(W=W_222), same_means=list(1:2)
theta_222csAR_int <- c(phi10_222, vec(A11_222), vec(A12_222), vec(W_222), 0.2, 0.4, alpha1_222)
mod222csAR_int <- GSMVAR(gdpdef, p=2, M=2, params=theta_222csAR_int, parametrization="mean",
                         constraints=C_222, structural_pars=list(W=W_222), same_means=list(1:2))
mod222csAR_int2 <- reorder_W_columns(mod222csAR_int, perm=2:1)

mod222gscsAR_int <- GSMVAR(gdpdef, p=2, M=c(1, 1), params=c(theta_222csAR_int, 20), model="G-StMVAR",
                           parametrization="mean", constraints=C_222, structural_pars=list(W=W_222),
                           same_means=list(1:2))
mod222gscsAR_int2 <- reorder_W_columns(mod222gscsAR_int, perm=2:1)

mod222tcsAR_int <- GSMVAR(gdpdef, p=2, M=2, params=c(theta_222csAR_int, 10, 20), model="StMVAR",
                          parametrization="mean", constraints=C_222, structural_pars=list(W=W_222),
                          same_means=list(1:2))
mod222tcsAR_int2 <- reorder_W_columns(mod222tcsAR_int, perm=2:1)


test_that("reorder_W_columns works correctly", {
  expect_equal(mod222cs2$params, c(1.031, 2.356, 1.786, 3, 1.25, 0.06, 0.036, 1.335, -0.29, -0.083, -0.047, -0.356, -0.36753,
                                   2.16401, 0.89382, 0.71976, 1.30222, 7.14351, 0.368), tolerance=1e-3)
  expect_equal(mod112_2$params, c(0.89978, 1.70757, 0.98944, 0.00264, -0.00696, 0.98661, 0.45293, 2.944, 1.96105, 0.45293), tolerance=1e-3)
  expect_equal(mod112t_2$params, c(0.89978, 1.70757, 0.98944, 0.00264, -0.00696, 0.98661, 0.45293, 2.944, 1.96105, 0.45293, 10), tolerance=1e-3)
  expect_equal(mod222csAR_int2$params,
               c(1.0300000, 2.3600000, 1.2500000, 0.0600000, 0.0400000, 1.3400000, -0.2900000, -0.0800000, -0.0500000, -0.3600000, 0.3653923,
                 -2.1643472, -0.8924620, -0.7180539, 0.4000000, 0.2000000, 0.3700000), tolerance=1e-6)
  expect_equal(mod222gscsAR_int2$params,
               c(1.0300000, 2.3600000, 1.2500000, 0.0600000, 0.0400000, 1.3400000, -0.2900000, -0.0800000, -0.0500000, -0.3600000, 0.3653923,
                 -2.1643472, -0.8924620, -0.7180539, 0.4000000, 0.2000000, 0.3700000, 20), tolerance=1e-6)
  expect_equal(mod222tcsAR_int2$params,
               c(1.0300000, 2.3600000, 1.2500000, 0.0600000, 0.0400000, 1.3400000, -0.2900000, -0.0800000, -0.0500000, -0.3600000, 0.3653923,
                 -2.1643472, -0.8924620, -0.7180539, 0.4000000, 0.2000000, 0.3700000, 10, 20), tolerance=1e-6)
})

mod222cs3 <- swap_W_signs(mod222cs, which_to_swap=1:2)
mod112_3 <- swap_W_signs(mod112, which_to_swap=2)
mod112t_3 <- swap_W_signs(mod112t, which_to_swap=2)
mod222csAR_int3 <- swap_W_signs(mod222csAR_int, which_to_swap=1)
mod222gscsAR_int3 <- swap_W_signs(mod222gscsAR_int, which_to_swap=1)
mod222tcsAR_int3 <- swap_W_signs(mod222tcsAR_int, which_to_swap=1)


test_that("swap_W_signs works correctly", {
  expect_equal(mod222cs3$params, c(1.031, 2.356, 1.786, 3, 1.25, 0.06, 0.036, 1.335, -0.29, -0.083, -0.047, -0.356, -0.89382,
                                   -0.71976, 0.36753, -2.16401, 7.14351, 1.30222, 0.368), tolerance=1e-3)
  expect_equal(mod112_3$params, c(0.89978, 1.70757, 0.98944, 0.00264, -0.00696, 0.98661, 1.96105, 0.45293, -0.45293, -2.944),
               tolerance=1e-3)
  expect_equal(mod112t_3$params, c(0.89978, 1.70757, 0.98944, 0.00264, -0.00696, 0.98661, 1.96105, 0.45293, -0.45293, -2.944, 10),
               tolerance=1e-3)
  expect_equal(mod222csAR_int3$params,
               c(1.0300000, 2.3600000, 1.2500000, 0.0600000, 0.0400000, 1.3400000, -0.2900000, -0.0800000, -0.0500000, -0.3600000,
                 0.8924620, 0.7180539, 0.3653923, -2.1643472, 0.2000000, 0.4000000, 0.3700000), tolerance=1e-6)
  expect_equal(mod222gscsAR_int3$params,
               c(1.0300000, 2.3600000, 1.2500000, 0.0600000, 0.0400000, 1.3400000, -0.2900000, -0.0800000, -0.0500000, -0.3600000,
                 0.8924620, 0.7180539, 0.3653923, -2.1643472, 0.2000000, 0.4000000, 0.3700000, 20), tolerance=1e-6)
  expect_equal(mod222tcsAR_int3$params,
               c(1.0300000, 2.3600000, 1.2500000, 0.0600000, 0.0400000, 1.3400000, -0.2900000, -0.0800000, -0.0500000, -0.3600000,
                 0.8924620, 0.7180539, 0.3653923, -2.1643472, 0.2000000, 0.4000000, 0.3700000, 10, 20), tolerance=1e-6)
})


test_that("update_numtols works correctly", {
  mod112t_4 <- update_numtols(mod112t, stat_tol=1e-5, posdef_tol=1e-6, df_tol=1e-7)
  mod222csLAR_4 <- update_numtols(mod222csAR_int, stat_tol=1e-5, posdef_tol=1e-6, df_tol=1e-7)
  mod222gscsLAR_4 <- update_numtols(mod222gscsAR_int, stat_tol=1e-5, posdef_tol=1e-6, df_tol=1e-7)
  new_tols <- list(stat_tol=1e-5, posdef_tol=1e-6, df_tol=1e-7)

  expect_equal(mod112t_4$num_tols, new_tols, tolerance=1e-8)
  expect_equal(mod222csLAR_4$num_tols, new_tols, tolerance=1e-8)
  expect_equal(mod222gscsLAR_4$num_tols, new_tols, tolerance=1e-8)
})


# StMVAR(1, 2), d=2 model:
params12 <- c(0.5453, 0.1157, 0.331, 0.0537, -0.0422, 0.7089, 0.4181, 0.0018,
              0.0413, 1.6004, 0.4843, 0.1256, -0.0311, -0.6139, 0.7221, 1.2123, -0.0357,
              0.1381, 0.8337)
mod12t <- suppressWarnings(GSMVAR(gdpdef, p=1, M=2, params=c(params12, 7, 90000), model="StMVAR"))
mod12gs <- suppressWarnings(GSMVAR(gdpdef, p=1, M=c(1, 1), params=c(params12, 90), model="G-StMVAR"))

test_that("sgmvar_to_gstmvar works correctly", {
  new_mod12t_1 <- suppressMessages(stmvar_to_gstmvar(mod12t, maxdf=100, maxit=1))
  new_mod12t_2 <- suppressMessages(stmvar_to_gstmvar(mod12t, maxdf=5, maxit=1))
  new_mod12t_3 <- suppressWarnings(suppressMessages(stmvar_to_gstmvar(mod12t, maxdf=90001, maxit=1)))
  new_mod12gs_1 <- suppressWarnings(suppressMessages(stmvar_to_gstmvar(mod12gs, maxdf=100, maxit=1)))
  new_mod12gs_2 <- suppressWarnings(suppressMessages(stmvar_to_gstmvar(mod12gs, maxdf=50, maxit=1)))

  expect_equal(new_mod12t_1$model$model, "G-StMVAR")
  expect_equal(new_mod12t_1$model$M, c(1, 1))
  expect_equal(new_mod12t_2$model$model, "GMVAR")
  expect_equal(new_mod12t_2$model$M, 2)
  expect_equal(new_mod12t_3$model$model, "StMVAR")
  expect_equal(new_mod12t_3$model$M, 2)
  expect_equal(new_mod12gs_1$model$model, "G-StMVAR")
  expect_equal(new_mod12gs_1$model$M, c(1, 1))
  expect_equal(new_mod12gs_2$model$model, "GMVAR")
  expect_equal(new_mod12gs_2$model$M, 2)
})
