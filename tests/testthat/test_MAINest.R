context("Functions in MAINest")
library(gmvarkit)

test_that("fitGSMVAR does not throw errors", {
   tmp <- suppressMessages(fitGSMVAR(gdpdef, p=1, M=1, ncalls=1, ncores=1, maxit=1, seeds=1, print_res=FALSE, ngen=2))
   W_112 <- matrix(c(1, -1, NA, 1), nrow=2)
   tmp <- suppressMessages(fitGSMVAR(gdpdef, p=1, M=1, structural_pars=list(W=W_112),
                                    ncalls=1, ncores=1, maxit=1, seeds=1, print_res=FALSE, ngen=2))
   tmp <- suppressMessages(fitGSMVAR(gdpdef, p=1, M=2, conditional=FALSE, structural_pars=list(W=W_112), same_means=list(1:2),
                                    parametrization="mean", ncalls=1, ncores=1, maxit=1, seeds=1, print_res=FALSE, ngen=1))
   expect_true(TRUE)
})


test_that("get_minval works correctly", {
   expect_equal(get_minval(gdpdef), -99999)
   expect_equal(get_minval(cbind(rep(0, 1001), rep(0, 1001), rep(0, 1001))), -9999999)
})


params_122 <- c(1.273, 0.077, 0.954, 0.143, -0.013, 1.007, 2.437, 1.328, 5.331, 2.612,
                5.871, 0.96, -0.04, -0.016, 0.958, 7.527, 1.767, 10.702, 0.794)

C_mat <- rbind(diag(2*2^2), diag(2*2^2))
params_222c <- c(-0.2394, 1.62068, 0.80599, 2.52145, 1.11216, 0.0043, -0.00799, 1.3045,
                 -0.17324, -0.04317, 0.00787, -0.32015, 2.05112, 0.50578, 8.2455, 5.31921,
                 3.53916, 11.71524, 0.55429)

test_that("iterate_more works correctly", {
   gsmvar122 <- suppressMessages(iterate_more(GSMVAR(gdpdef, p=1, M=2, d=2, params=params_122), maxit=2))
   expect_equal(gsmvar122$params, c(1.05403, 0.11766, 0.70969, 0.18783, -0.18745, 1.03718, 2.44473,
                                   1.33305, 5.29359, 2.612, 5.871, 0.96, -0.04, -0.016, 0.958, 7.527,
                                   1.767, 10.702, 0.794),
                tolerance=1e-3)

   gsmvar222c <- suppressMessages(iterate_more(GSMVAR(gdpdef, p=2, M=2, d=2, params=params_222c, constraints=C_mat), maxit=2))
   expect_equal(gsmvar222c$params, c(-0.23898, 1.62018, 0.80578, 2.5218, 1.10873, 0.00901, -0.00753, 1.30718, -0.17653,
                                    -0.03837, 0.00851, -0.31741, 2.05107, 0.50571, 8.24552, 5.3192, 3.53919, 11.71521,
                                    0.55438), tolerance=1e-4)
})

# G-StMVAR(1, 1, 1)
sg_stmvar12_pars <- c(1.59756, 0.48318, 0.54543, 0.11584, 0.12584, -0.03109, -0.61263, 0.7226, 0.33097, 0.05368, -0.04244,
                      0.70867, 0.9545, 0.15627, -0.54858, 0.33714, 0.36483, 0.28422, 0.1666, 7.56742)
sg_stmvar12 <- GSMVAR(gdpdef, p=1, M=c(1, 1), d=2, params=sg_stmvar12_pars, model="G-StMVAR",
                      structural_pars=list(W=matrix(nrow=2, ncol=2)))
new_W1 <- matrix(c(1, 0, -1, NA), nrow=2)

# Cholesky VAR(1)
var1pars <- c(0.649528, 0.066507, 0.288525, 0.021766, -0.144026, 0.897103, 0.601791, -0.002944, 0.067224)
var1 <- GSMVAR(gdpdef, p=1, M=1, d=2, params=var1pars, model="GMVAR")
new_W2 <- matrix(c(1, NA, 0, 1), nrow=2)

test_that("estimate_ works correctly", {
  new_fit12gss <- estimate_sgsmvar(sg_stmvar12, new_W=new_W1, ncalls=1, ncores=1, seeds=1)
  expect_equal(new_fit12gss$params, c(1.68044468, 0.50581789, 0.53848797, 0.11334361, 0.11786532,
                                      -0.03078259, -0.65462904, 0.71336338, 0.33229816, 0.05399712,
                                      -0.03006147, 0.71400402, 1.10621116, -0.01436285, 0.37288158,
                                      0.34530781, 0.30271439, 0.14914748, 7.46245286), tolerance=1e-3)

  new_fit11s <- estimate_sgsmvar(var1, new_W=new_W2, ncalls=1, ncores=1, seeds=1)
  expect_equal(new_fit11s$params, c(0.649528, 0.066507, 0.288525, 0.021766, -0.144026, 0.897103,
                                    0.775748, -0.003775, 0.259252), tolerance=1e-3)
})
