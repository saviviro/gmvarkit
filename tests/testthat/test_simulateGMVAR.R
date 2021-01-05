context("simulateGMVAR")
library(gmvarkit)

# NOTE that some elements of these tests use random elements obtained from simulation algorithms

data <- cbind(10*eurusd[,1], 100*eurusd[,2])

## A(M)(p)_(p)(M)(d)

# p=1, M=1, d=2, parametrization="mean"
phi10_112 <- c(1.07, 127.71)
A11_112 <- matrix(c(0.99, 0.00, -0.01, 0.99), nrow=2)
Omega1_112 <- matrix(c(4.05, 2.22, 8.87, 2.22), nrow=2)
theta_112 <- c(phi10_112, vec(A11_112), vech(Omega1_112))
mod_112 <- GMVAR(data, p=1, M=1, d=2, params=theta_112, conditional=TRUE, parametrization="mean", constraints=NULL)

# p=2, M=2, d=2, AR paramars same, non-diagonals zero, intercept
theta_222c <- c(0.3552775, 3.1929675, -0.1143198, 2.8294743, 1.2633425, 1.3375150, -0.2919742, -0.3624010,
                5.5971764, 3.4559442, 9.6221422, 0.9820759, -0.3267521, 5.2358855, 0.6501600)
mat0 <- matrix(c(1, rep(0, 10), 1, rep(0, 8), 1, rep(0, 10), 1), nrow=2*2^2, byrow=FALSE) # Laske paperilla monimutkaisemmat rajoitteet
C_222c <- rbind(mat0, mat0)
mod_222c <- GMVAR(data, p=2, M=2, d=2, params=theta_222c, conditional=TRUE, parametrization="intercept", constraints=C_222c)

# p=2, M=2, d=2, parametrization="mean", constraints=C_mat, same_means=list(1:2)
C_mat <- rbind(diag(2*2^2), diag(2*2^2))
params_222cm <- c(-9.12048853805255, 123.142757183508, 1.2658425363326, 0.0675545389606989, 0.0331264235657607,
                  1.33370494344656, -0.285882557831441, -0.0769144929653558, -0.0382772162867802, -0.351635998882842,
                  5.8625623309659, 3.57488618757834, 9.70846346569286, 0.869261580580846, -0.248703862116217,
                  5.17613656742281, 0.439575388572472)
mod_222cm <- GMVAR(data, p=2, M=2, params=params_222cm, parametrization="mean", constraints=C_mat, same_means=list(1:2))

set.seed(13)
sim_112 <- simulateGMVAR(mod_112, nsimu=1)
sim_222 <- simulateGMVAR(mod_222c, nsimu=3)
sim_112_2 <- simulateGMVAR(mod_112, nsimu=3, ntimes=3)
sim_222_2 <- simulateGMVAR(mod_222c, nsimu=1, ntimes=2)

set.seed(1); sim_222cm <- simulateGMVAR(mod_222cm, nsimu=2, ntimes=2)

test_that("simulateGMVAR works correctly", {
  expect_equal(sim_112$sample[1,], c(-5.90633, 113.67962), tolerance=1e-5)
  expect_equal(sim_112$component, 1)
  expect_equal(sim_112$mixing_weights, as.matrix(1))
  expect_equal(sim_222$sample[3,], c(-0.8352874, 119.1928212), tolerance=1e-5)
  expect_equal(sim_222$component, c(2, 2, 2))
  expect_equal(sim_222$mixing_weights[3,], c(0.1462826, 0.8537174), tolerance=1e-5)

  expect_equal(sim_112_2$sample[3, , 3], c(-10.24494, 109.41560), tolerance=1e-5)
  expect_equal(sim_112_2$component[,1], c(1, 1 ,1))
  expect_equal(sim_112_2$mixing_weights[, , 3], c(1, 1, 1))
  expect_equal(sim_222_2$sample[1, , 2], c(4.265082, 115.896110), tolerance=1e-5)
  expect_equal(sim_222_2$component[1,], c(1, 1))
  expect_equal(sim_222_2$mixing_weights[, , 2], c(0.96265621, 0.03734379), tolerance=1e-5)

  expect_equal(sim_222cm$sample[1:8], c(-10.496086, -10.881399, 145.642259, 151.055338,
                                        -8.321608, -8.195829, 145.239271, 144.538774), tolerance=1e-5)
  expect_equal(sim_222cm$component[1:4], c(2, 2, 2, 2))
  expect_equal(sim_222cm$mixing_weights[1:8], c(0.2172295, 0.1220861, 0.7827705, 0.8779139,
                                                0.2172295, 0.3708393, 0.7827705, 0.6291607), tolerance=1e-5)
})
