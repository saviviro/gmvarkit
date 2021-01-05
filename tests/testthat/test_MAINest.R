context("Functions in MAINest")
library(gmvarkit)

data <- cbind(10*eurusd[,1], 100*eurusd[,2])

test_that("fitGMVAR does not throw errors", {
   tmp <- suppressMessages(fitGMVAR(data, p=1, M=1, ncalls=1, ncores=1, maxit=1, seeds=1, print_res=FALSE, ngen=2))
   W_112 <- matrix(c(1, -1, NA, 1), nrow=2)
   tmp <- suppressMessages(fitGMVAR(data, p=1, M=1, structural_pars=list(W=W_112),
                                    ncalls=1, ncores=1, maxit=1, seeds=1, print_res=FALSE, ngen=2))
   tmp <- suppressMessages(fitGMVAR(data, p=1, M=2, conditional=FALSE, structural_pars=list(W=W_112), same_means=list(1:2),
                                    parametrization="mean", ncalls=1, ncores=1, maxit=1, seeds=1, print_res=FALSE, ngen=1))
   expect_true(TRUE)
})


test_that("get_minval works correctly", {
   expect_equal(get_minval(data), -99999)
   expect_equal(get_minval(cbind(rep(0, 1001), rep(0, 1001), rep(0, 1001))), -9999999)
})


params_122 <- c(1.273, 0.077, 0.954, 0.143, -0.013, 1.007, 2.437, 1.328, 5.331, 2.612,
                5.871, 0.96, -0.04, -0.016, 0.958, 7.527, 1.767, 10.702, 0.794)

C_mat <- rbind(diag(2*2^2), diag(2*2^2))
params_222c <- c(-0.2394, 1.62068, 0.80599, 2.52145, 1.11216, 0.0043, -0.00799, 1.3045,
                 -0.17324, -0.04317, 0.00787, -0.32015, 2.05112, 0.50578, 8.2455, 5.31921,
                 3.53916, 11.71524, 0.55429)

test_that("iterate_more works correctly", {
   gmvar122 <- suppressMessages(iterate_more(GMVAR(data, p=1, M=2, d=2, params=params_122), maxit=2))
   expect_equal(gmvar122$params, c(1.273004, 0.076998, 0.953971, 0.142986, -0.012552, 1.006759, 2.436997, 1.327998, 5.331001,
                                   2.612, 5.871, 0.959994, -0.039999, -0.015939, 0.95798, 7.526999, 1.767002, 10.702, 0.793997),
                tolerance=1e-5)

   gmvar222c <- suppressMessages(iterate_more(GMVAR(data, p=2, M=2, d=2, params=params_222c, constraints=C_mat), maxit=2))
   expect_equal(gmvar222c$params, c(-0.239397, 1.62068, 0.805992, 2.521451, 1.11222, 0.004325, -0.007545, 1.304624, -0.173193,
                                    -0.043143, 0.008296, -0.320026, 2.051118, 0.505781, 8.245499, 5.319211, 3.53916, 11.71524,
                                    0.554291), tolerance=1e-5)
})
