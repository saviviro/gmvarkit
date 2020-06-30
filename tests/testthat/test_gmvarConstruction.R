context("GMVAR construction")
library(gmvarkit)

data <- cbind(10*eurusd[,1], 100*eurusd[,2])

params122 <- c(0.623, -0.129, 0.959, 0.089, -0.006, 1.006, 1.746, 0.804, 5.804, 3.245, 7.913,
               0.952, -0.037, -0.019, 0.943, 6.926, 3.982, 12.135, 0.789) # p=1, M=2, d=2
C_mat <- rbind(diag(2*2^2), diag(2*2^2))
params222c <- c(1.031, 2.356, 1.786, 3.000, 1.250, 0.060, 0.036, 1.335, -0.290, -0.083, -0.047,
                -0.356, 0.934, -0.152, 5.201, 5.883, 3.560, 9.799, 0.368) # p=2, M=2, d=2, AR parameters restricted same for both regimes

test_that("GMVAR works correctly", {
  mod122 <- GMVAR(data, p=1, M=2, params=params122)
  mod222c <- GMVAR(data, p=2, M=2, params=params222c, constraints=C_mat)
  expect_equal(mod122$params, params122)
  expect_equal(mod222c$params, params222c)
})


test_that("add_data works correctly", {
  mod122 <- GMVAR(p=1, M=2, d=2, params=params122)
  mod122_2 <- add_data(data, mod122)
  expect_equal(mod122_2$data, data)
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
})


test_that("alt_gmvar does not throw errors", {
  gmvar11 <- suppressMessages(fitGMVAR(data[1:20,], p=1, M=1, ncalls=2, ncores=1, maxit=1, seeds=1:2, print_res=FALSE, ngen=1))
  tmp <- alt_gmvar(gmvar11, which_round=1, calc_cond_moments=FALSE, calc_std_errors=FALSE)
  expect_true(TRUE)
})
