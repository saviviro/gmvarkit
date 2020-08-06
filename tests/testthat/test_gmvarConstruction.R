context("GMVAR construction")
library(gmvarkit)

data <- cbind(10*eurusd[,1], 100*eurusd[,2])


## A(M)(p)_(p)(M)(d)
rbind_diags <- function(p, M, d) {
  I <- diag(p*d^2)
  Reduce(rbind, replicate(M, I, simplify=FALSE))
}

# p=1, M=2, d=2
params122 <- c(0.623, -0.129, 0.959, 0.089, -0.006, 1.006, 1.746, 0.804, 5.804, 3.245, 7.913,
               0.952, -0.037, -0.019, 0.943, 6.926, 3.982, 12.135, 0.789) # p=1, M=2, d=2

# p=2, M=2, d=2, constrained
C_mat <- rbind(diag(2*2^2), diag(2*2^2))
params222c <- c(1.031, 2.356, 1.786, 3.000, 1.250, 0.060, 0.036, 1.335, -0.290, -0.083, -0.047,
                -0.356, 0.934, -0.152, 5.201, 5.883, 3.560, 9.799, 0.368) # p=2, M=2, d=2, AR parameters restricted same for both regimes

# p=2, M=1, d=2, SGMVAR, W constrained
phi10_212 <- c(1.03, 2.36)
A11_212 <- matrix(c(1.25, 0.06, 0.04, 1.34), nrow=2, byrow=FALSE)
A12_212 <- matrix(c(-0.29, -0.08, -0.05, -0.36), nrow=2, byrow=FALSE)
Omega1_212 <- matrix(c(0.93, -0.15, -0.15, 5.20), nrow=2, byrow=FALSE)

W_212 <- t(chol(Omega1_212))
theta_212csW <- c(phi10_212, vec(A11_212), vec(A12_212), Wvec(W_212))

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


test_that("GMVAR works correctly", {
  mod122 <- GMVAR(data, p=1, M=2, params=params122)
  mod222c <- GMVAR(data, p=2, M=2, params=params222c, constraints=C_mat)
  expect_equal(mod122$params, params122)
  expect_equal(mod222c$params, params222c)

  # SGMVAR
  mod212csW <- GMVAR(data, p=2, M=1, params=theta_212csW, structural_pars=list(W=W_212))
  mod122csLAR <- GMVAR(data, p=1, M=2, params=theta_122csLAR, constraints=C_122,
                       structural_pars=list(W=W_122, C_lambda=C_lambda_122))
  expect_equal(mod212csW$params, theta_212csW)
  expect_equal(mod122csLAR$params, theta_122csLAR)
})


test_that("add_data works correctly", {
  mod122 <- GMVAR(p=1, M=2, d=2, params=params122)
  mod122_2 <- add_data(data=data, gmvar=mod122)
  expect_equal(mod122_2$data, data)

  # SGMVAR
  mod122csLAR <- GMVAR(p=1, M=2, d=2, params=theta_122csLAR, constraints=C_122,
                       structural_pars=list(W=W_122, C_lambda=C_lambda_122))
  mod122csLAR_2 <- add_data(data=data, gmvar=mod122csLAR)
  expect_equal(mod122csLAR_2$data, data)
})


test_that("swap_parametrization works correctly", {
  mod122 <- GMVAR(data, p=1, M=2, params=params122)
  mod222c <- GMVAR(data, p=2, M=2, params=params222c, constraints=C_mat)
  mod122_2 <- swap_parametrization(mod122)
  mod222c_2 <- swap_parametrization(mod222c)
  expect_equal(mod122_2$params, c(-10.291667, 174.159722, 0.959, 0.089, -0.006, 1.006, 1.746, 0.804,
                                  5.804, 17.028037, 127.771274, 0.952, -0.037, -0.019, 0.943, 6.926,
                                  3.982, 12.135, 0.789), tolerance=1e-4)
  expect_equal(mod222c_2$params, c(-7.265758, 120.148211, 7.67632, 134.449744, 1.25, 0.06, 0.036, 1.335,
                                   -0.29, -0.083, -0.047, -0.356, 0.934, -0.152, 5.201, 5.883, 3.56,
                                   9.799, 0.368), tolerance=1e-4)

  # SGMVAR
  mod212csW <- GMVAR(data, p=2, M=1, params=theta_212csW, structural_pars=list(W=W_212))
  mod122csLAR <- GMVAR(data, p=1, M=2, params=theta_122csLAR, constraints=C_122,
                       structural_pars=list(W=W_122, C_lambda=C_lambda_122))
  mod212csW_2 <- swap_parametrization(mod212csW)
  mod122csLAR_2 <- swap_parametrization(mod122csLAR)
  expect_equal(mod212csW_2$params, c(-5, 123, 1.25, 0.06, 0.04, 1.34, -0.29, -0.08, -0.05, -0.36, 0.96437,
                                     -0.15554, 2.27504), tolerance=1e-4)
  expect_equal(mod122csLAR_2$params, c(-6.05, 5.45, -7.21, 8.37, 0.9, 0.3, -0.3, 0.9, -0.89246, -0.71805,
                                       0.36539, -2.16435, 1, 0.37), tolerance=1e-4)
})


test_that("alt_gmvar does not throw errors", {
  gmvar11 <- suppressMessages(fitGMVAR(data[1:20,], p=1, M=1, ncalls=2, ncores=1, maxit=1, seeds=1:2, print_res=FALSE, ngen=1))
  tmp <- alt_gmvar(gmvar11, which_round=1, calc_cond_moments=FALSE, calc_std_errors=FALSE)
  expect_true(TRUE)
})

# GMVAR(1,2), d=2 model:
params122 <- c(0.623, -0.129, 0.959, 0.089, -0.006, 1.006, 1.746,
               0.804, 5.804, 3.245, 7.913, 0.952, -0.037, -0.019,
               0.943, 6.926, 3.982, 12.135, 0.789)
mod122 <- GMVAR(data, p=1, M=2, params=params122)
mod122s <- gmvar_to_sgmvar(mod122)

# GMVAR(2,2), d=2 model with AR-parameters restricted to be
# the same for both regimes:
C_mat <- rbind(diag(2*2^2), diag(2*2^2))
params222c <- c(1.031, 2.356, 1.786, 3.000, 1.250, 0.060, 0.036,
                1.335, -0.290, -0.083, -0.047, -0.356, 0.934, -0.152,
                5.201, 5.883, 3.560, 9.799, 0.368)
mod222c <- GMVAR(p=2, M=2, d=2, params=params222c, constraints=C_mat)
mod222cs <- gmvar_to_sgmvar(mod222c)

test_that("gmvar_to_sgmvar works correctly", {
  expect_equal(mod122s$params, c(0.623, -0.129, 3.245, 7.913, 0.959, 0.089, -0.006, 1.006, 0.952, -0.037, -0.019, 0.943,
                                 1.31216, 0.87892, -0.15571, 2.2431, 3.99732, 1.79808, 0.789), tolerance=1e-3)
  expect_equal(mod222cs$params, c(1.031, 2.356, 1.786, 3, 1.25, 0.06, 0.036, 1.335, -0.29, -0.083, -0.047, -0.356, 0.89382,
                                  0.71976, -0.36753, 2.16401, 7.14351, 1.30222, 0.368), tolerance=1e-3)
})


mod222cs2 <- reorder_W_columns(mod222cs, perm=2:1)

params112 <- c(0.899784972, 1.707574085, 0.989440729, 0.002636335, -0.006957649, 0.986606661,
               1.961053552, 0.452934215, 0.452934215, 2.944004374)
W112 <- matrix(NA, nrow=2, ncol=2)
diag(W112) <- c(1, 1)
mod112 <- GMVAR(data, p=1, M=1, params=params112, structural_pars=list(W=W112))
mod112_2 <- reorder_W_columns(mod112, perm=2:1)

test_that("reorder_W_columns works correctly", {
  expect_equal(mod222cs2$params, c(1.031, 2.356, 1.786, 3, 1.25, 0.06, 0.036, 1.335, -0.29, -0.083, -0.047, -0.356, -0.36753,
                                   2.16401, 0.89382, 0.71976, 1.30222, 7.14351, 0.368), tolerance=1e-3)
  expect_equal(mod112_2$params, c(0.89978, 1.70757, 0.98944, 0.00264, -0.00696, 0.98661, 0.45293, 2.944, 1.96105, 0.45293), tolerance=1e-3)

})

mod222cs3 <- swap_W_signs(mod222cs, which_to_swap=1:2)
mod112_3 <- swap_W_signs(mod112, which_to_swap=2)

test_that("swap_W_signs works correctly", {
  expect_equal(mod222cs3$params, c(1.031, 2.356, 1.786, 3, 1.25, 0.06, 0.036, 1.335, -0.29, -0.083, -0.047, -0.356, -0.89382,
                                   -0.71976, 0.36753, -2.16401, 7.14351, 1.30222, 0.368), tolerance=1e-3)
  expect_equal(mod112_3$params, c(0.89978, 1.70757, 0.98944, 0.00264, -0.00696, 0.98661, 1.96105, 0.45293, -0.45293, -2.944),
               tolerance=1e-3)
})










































































