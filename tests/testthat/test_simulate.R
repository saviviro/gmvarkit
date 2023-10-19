context("simulate")
library(gmvarkit)

# NOTE that some elements of these tests use random elements obtained from simulation algorithms

## A(M)(p)_(p)(M)(d)

# p=1, M=1, d=2, parametrization="mean"
phi10_112 <- c(0.75, 0.8)
A11_112 <- matrix(c(0.29, 0.02, -0.14, 0.9), nrow=2, byrow=FALSE)
Omega1_112 <- matrix(c(0.60, 0.01, 0.01, 0.07), nrow=2, byrow=FALSE)
theta_112 <- c(phi10_112, vec(A11_112), vech(Omega1_112))
mod_112 <- GSMVAR(gdpdef, p=1, M=1, d=2, params=theta_112, conditional=TRUE, parametrization="mean", constraints=NULL)

mod_112t <- GSMVAR(gdpdef, p=1, M=1, d=2, params=c(theta_112, 3), model="StMVAR", parametrization="mean")


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
mod_222s <- GSMVAR(gdpdef, p=2, M=2, d=2, params=theta_222s, structural_pars=list(W=W_222))

mod_222gss <- GSMVAR(gdpdef, p=2, M=c(1, 1), d=2, params=c(theta_222s, 20), model="G-StMVAR", structural_pars=list(W=W_222))

# p=2, M=2, d=2, AR paramars same, non-diagonals zero, intercept
theta_222c <- c(0.33782, 0.183512, 0.472168, 0.095311, 0.201199, 0.600596, 0.237819,
                0.23529, 1.077816, -0.016343, 0.112771, 0.22199, 0.005582, 0.028126, 0.492844)
mat0 <- matrix(c(1, rep(0, 10), 1, rep(0, 8), 1, rep(0, 10), 1), nrow=2*2^2, byrow=FALSE)
C_222c <- rbind(mat0, mat0)
mod_222c <- GSMVAR(gdpdef, p=2, M=2, d=2, params=theta_222c, constraints=C_222c)

mod_222tc <- GSMVAR(gdpdef, p=2, M=2, d=2, params=c(theta_222c, 20, 25), model="StMVAR", constraints=C_222c)

# p=2, M=2, d=2, parametrization="mean", constraints=C_mat, same_means=list(1:2)
C_mat <- rbind(diag(2*2^2), diag(2*2^2))
params_222cm <- c(0.811034, 0.578587, 0.212084, 0.020444, -0.193005, 0.624671,
                  0.235827, 0.013962, 0.053267, 0.262703, 1.06105, -0.013519,
                  0.114109, 0.229542, 0.003092, 0.027266, 0.424341)
mod_222cm <- GSMVAR(gdpdef, p=2, M=2, params=params_222cm, parametrization="mean", constraints=C_mat, same_means=list(1:2))

# p=2, M=2, d=2, model="GMVAR", parametrization="mean", constraints=C_mat, same_means=list(1:2),
# weight_constraints=0.4, structural_pars=list(W=W_222, fixed_lambdas=c(7, 4))
params_222cmwsF <- c(0.811034, 0.578587, 0.212084, 0.020444, -0.193005, 0.624671,
                     0.235827, 0.013962, 0.053267, 0.262703, Wvec(W_222))
mod_222cmwsF <- GSMVAR(gdpdef, p=2, M=2, params=params_222cmwsF,  model="GMVAR", parametrization="mean",
                       constraints=C_mat, same_means=list(1:2), weight_constraints=0.4,
                       structural_pars=list(W=W_222, fixed_lambdas=c(7, 4)))

sim_112 <- simulate.gsmvar(mod_112, nsim=1, seed=1)
sim_112t <- simulate.gsmvar(mod_112t, nsim=2, seed=1)
sim_222 <- simulate.gsmvar(mod_222c, nsim=3, seed=2)
sim_222s <- simulate.gsmvar(mod_222s, nsim=2, seed=3, init_values=matrix(rep(1, times=4), nrow=2))
sim_222gss <- simulate.gsmvar(mod_222gss, nsim=2, seed=3, init_regimes=2)
sim_112_2 <- simulate.gsmvar(mod_112, nsim=3, seed=4, ntimes=3)
sim_222_2 <- simulate.gsmvar(mod_222c, nsim=1, seed=5, ntimes=2)
sim_222tc <- simulate.gsmvar(mod_222tc, nsim=2, seed=5, ntimes=2)
sim_222cm <- simulate.gsmvar(mod_222cm, nsim=2, seed=6, ntimes=2)
sim_222cmwsF <- simulate.gsmvar(mod_222cmwsF, nsim=2, seed=6, ntimes=2, init_values=gdpdef)

test_that("simulate.gsmvar works correctly", {
  expect_equal(sim_112$sample[1,], c(1.571209, 1.040196), tolerance=1e-5)
  expect_equal(sim_112$component, 1)
  expect_equal(sim_112$mixing_weights, as.matrix(1))
  expect_equal(sim_112t$sample[2,], c(0.6914065, 2.6814250), tolerance=1e-5)
  expect_equal(sim_112t$component, c(1, 1))
  expect_equal(sim_112t$mixing_weights, as.matrix(c(1, 1)))
  expect_equal(sim_222$sample[3,], c(1.2680622, 0.7005461), tolerance=1e-5)
  expect_equal(sim_222$component, c(2, 2, 2))
  expect_equal(sim_222$mixing_weights[3,], c(0.07284036, 0.92715964), tolerance=1e-5)
  expect_equal(sim_222s$sample[2,], c(1.3121072, 0.7512483), tolerance=1e-5)
  expect_equal(sim_222gss$sample[2,], c(-0.5964731, 1.9836836), tolerance=1e-5)
  expect_equal(sim_222gss$mixing_weights[2,], c(1.246476e-06, 9.999988e-01), tolerance=1e-5)
  expect_equal(sim_222gss$component, c(2, 2), tolerance=1e-5)
  expect_equal(sim_112_2$sample[3, , 3], c(0.1946473, 1.0668018), tolerance=1e-5)
  expect_equal(sim_112_2$component[,1], c(1, 1 ,1))
  expect_equal(sim_112_2$mixing_weights[, , 3], c(1, 1, 1))
  expect_equal(sim_222_2$sample[1, , 2], c(0.5796451, 0.7611220), tolerance=1e-5)
  expect_equal(sim_222_2$component[1,], c(2, 2))
  expect_equal(sim_222_2$mixing_weights[, , 2], c(0.04113062, 0.95886938), tolerance=1e-5)
  expect_equal(c(sim_222tc$sample[2, , ]), c(0.3599588, 0.7637256, 1.9319898, 1.2949551), tolerance=1e-5)
  expect_equal(sim_222tc$component[1,], c(1, 1))
  expect_equal(c(sim_222tc$mixing_weights[2, , ]), c(0.93108341, 0.06891659, 0.98868343, 0.01131657), tolerance=1e-5)
  expect_equal(sim_222cm$sample[1:8], c(2.17141141, 1.13700895, 0.06482195, 0.45695849, 3.57077564,
                                        1.60079169, 0.09547996, 0.80902621), tolerance=1e-5)
  expect_equal(sim_222cm$component[1:4], c(1, 1, 1, 1))
  expect_equal(sim_222cm$mixing_weights[1:8], c(9.998214e-01, 9.999746e-01, 1.786218e-04, 2.543522e-05,
                                                9.998214e-01, 9.999999e-01, 1.786218e-04, 6.528877e-08), tolerance=1e-5)

  expect_equal(sim_222cmwsF$sample[1:8], c(0.1033549, -1.4098509, 0.5003297, 0.6083894, 0.2050738, 0.7318845,
                                           0.5346906, 0.6693941), tolerance=1e-5)
  expect_equal(sim_222cmwsF$component[1:4], c(1, 2, 1, 1))
  expect_equal(sim_222cmwsF$mixing_weights[1:8], c(0.91561446, 0.84410853, 0.08438554, 0.15589147, 0.91561446,
                                                   0.86236485, 0.08438554, 0.13763515), tolerance=1e-5)
})
