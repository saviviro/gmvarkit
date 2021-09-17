context("simulateGMVAR")
library(gmvarkit)

# NOTE that some elements of these tests use random elements obtained from simulation algorithms

## A(M)(p)_(p)(M)(d)

# p=1, M=1, d=2, parametrization="mean"
phi10_112 <- c(0.75, 0.8)
A11_112 <- matrix(c(0.29, 0.02, -0.14, 0.9), nrow=2, byrow=FALSE)
Omega1_112 <- matrix(c(0.60, 0.01, 0.01, 0.07), nrow=2, byrow=FALSE)
theta_112 <- c(phi10_112, vec(A11_112), vech(Omega1_112))
mod_112 <- GMVAR(gdpdef, p=1, M=1, d=2, params=theta_112, conditional=TRUE, parametrization="mean", constraints=NULL)

# p=2, M=2, d=2, no constraints, structural
phi10_222 <- c(0.36, 0.12)
A11_222 <- matrix(c(0.22, 0.06, -0.15, 0.39), nrow=2, byrow=FALSE)
A12_222 <- matrix(c(0.41, -0.01, 0.08, 0.3), nrow=2, byrow=FALSE)
Omega1_222 <- matrix(c(0.21, 0.01, 0.01, 0.03), nrow=2, byrow=FALSE)

phi20_222 <- c(0.48, 0.07)
A21_222 <- matrix(c(0.22, 0.02, -0.12, 0.72), nrow=2, byrow=FALSE)
A22_222 <- matrix(c(0.09, 0.03, 0.04, 0.19), nrow=2, byrow=FALSE)
Omega2_222 <- matrix(c(1.10, 0.01, 0.01, 0.11), nrow=2, byrow=FALSE)

alpha1_222 <- 0.37

WL_222 <- diag_Omegas(Omega1_222, Omega2_222)
W_222 <- matrix(WL_222[1:(2^2)], nrow=2, byrow=FALSE)
lambdas_222 <- WL_222[(2^2 + 1):length(WL_222)]
theta_222s <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A21_222),
                vec(A22_222), vec(W_222), lambdas_222, alpha1_222) # SGMVAR
mod_222s <- GMVAR(gdpdef, p=2, M=2, d=2, params=theta_222s, conditional=TRUE, parametrization="intercept", constraints=NULL,
                  structural_pars=list(W=W_222))


# p=2, M=2, d=2, AR paramars same, non-diagonals zero, intercept
theta_222c <- c(0.33782, 0.183512, 0.472168, 0.095311, 0.201199, 0.600596, 0.237819,
                0.23529, 1.077816, -0.016343, 0.112771, 0.22199, 0.005582, 0.028126, 0.492844)
mat0 <- matrix(c(1, rep(0, 10), 1, rep(0, 8), 1, rep(0, 10), 1), nrow=2*2^2, byrow=FALSE)
C_222c <- rbind(mat0, mat0)
mod_222c <- GMVAR(gdpdef, p=2, M=2, d=2, params=theta_222c, conditional=TRUE, parametrization="intercept", constraints=C_222c)

# p=2, M=2, d=2, parametrization="mean", constraints=C_mat, same_means=list(1:2)
C_mat <- rbind(diag(2*2^2), diag(2*2^2))
params_222cm <- c(0.811034, 0.578587, 0.212084, 0.020444, -0.193005, 0.624671,
                  0.235827, 0.013962, 0.053267, 0.262703, 1.06105, -0.013519,
                  0.114109, 0.229542, 0.003092, 0.027266, 0.424341)
mod_222cm <- GMVAR(gdpdef, p=2, M=2, params=params_222cm, parametrization="mean", constraints=C_mat, same_means=list(1:2))

set.seed(1); sim_112 <- simulateGMVAR(mod_112, nsimu=1)
set.seed(2); sim_222 <- simulateGMVAR(mod_222c, nsimu=3)
set.seed(3); sim_222s <- simulateGMVAR(mod_222c, nsimu=2, init_values=matrix(rep(1, times=4), nrow=2))
set.seed(4); sim_112_2 <- simulateGMVAR(mod_112, nsimu=3, ntimes=3)
set.seed(5); sim_222_2 <- simulateGMVAR(mod_222c, nsimu=1, ntimes=2)
set.seed(6); sim_222cm <- simulateGMVAR(mod_222cm, nsimu=2, ntimes=2)

test_that("simulateGMVAR works correctly", {
  expect_equal(sim_112$sample[1,], c(1.797594, 1.629599), tolerance=1e-5)
  expect_equal(sim_112$component, 1)
  expect_equal(sim_112$mixing_weights, as.matrix(1))
  expect_equal(sim_222$sample[3,], c(1.2680622, 0.7005461), tolerance=1e-5)
  expect_equal(sim_222$component, c(2, 2, 2))
  expect_equal(sim_222$mixing_weights[3,], c(0.07284036, 0.92715964), tolerance=1e-5)
  expect_equal(sim_222s$sample[2,], c(0.4338438, 0.8722736), tolerance=1e-5)

  expect_equal(sim_112_2$sample[3, , 3], c(1.6326367, 0.4610806), tolerance=1e-5)
  expect_equal(sim_112_2$component[,1], c(1, 1 ,1))
  expect_equal(sim_112_2$mixing_weights[, , 3], c(1, 1, 1))
  expect_equal(sim_222_2$sample[1, , 2], c(0.5796451, 0.7611220), tolerance=1e-5)
  expect_equal(sim_222_2$component[1,], c(2, 2))
  expect_equal(sim_222_2$mixing_weights[, , 2], c(0.04113062, 0.95886938), tolerance=1e-5)

  expect_equal(sim_222cm$sample[1:8], c(2.17141141, 1.13700895, 0.06482195, 0.45695849, 3.57077564,
                                        1.60079169, 0.09547996, 0.80902621), tolerance=1e-5)
  expect_equal(sim_222cm$component[1:4], c(1, 1, 1, 1))
  expect_equal(sim_222cm$mixing_weights[1:8], c(9.998214e-01, 9.999746e-01, 1.786218e-04, 2.543522e-05,
                                                9.998214e-01, 9.999999e-01, 1.786218e-04, 6.528877e-08), tolerance=1e-5)
})
