context("GIRF")
library(gmvarkit)

# NOTE that some elements of these tests use random elements obtained from simulation algorithms

data <- cbind(10*eurusd[,1], 100*eurusd[,2])

## A(M)(p)_(p)(M)(d)

# Structural GMVAR(2, 2), d=2 model identified with sign-constraints:
params222s <- c(1.428, -0.808, 1.029, 5.84, 1.314, 0.145, 0.094, 1.292,
                -0.389, -0.07, -0.109, -0.281, 1.248, 0.077, -0.04, 1.266, -0.272,
                -0.074, 0.034, -0.313, 0.903, 0.718, -0.324, 2.079, 7, 1.44, 0.742)
W_222 <- matrix(c(1, NA, -1, 1), nrow=2, byrow=FALSE)
mod222s <- GMVAR(data, p=2, M=2, params=params222s, structural_pars=list(W=W_222))

girf1 <- GIRF(mod222s, N=2, R2=2, R1=2, seeds=1:2)
girf2 <- GIRF(mod222s, which_shocks=2, shock_size=1, N=1, R2=1, R1=1, seeds=1,
              include_mixweights=FALSE, init_values=mod222s$data, ci=0.1)
girf3 <- GIRF(mod222s, N=2, R2=1, R1=1, which_cumulative=1:2, seeds=1)

test_that("GIRF works correctly", {
  expect_equal(unname(girf1$girf_res[[1]]$point_est[3,]), c(2.7647446, 2.6555342, -0.2440685, 0.2440685), tolerance=1e-4)
  expect_equal(unname(girf1$girf_res[[2]]$point_est[1,]), c(-0.358230, 2.298642, 0.000000, 0.000000), tolerance=1e-4)
  expect_equal(unname(girf1$girf_res[[2]]$conf_ints[3, , 2]), c(2.728152, 2.788786, 3.435544, 3.496178), tolerance=1e-4)
  expect_equal(unname(girf1$girf_res[[1]]$conf_ints[1, 1, ]), c(1.409993, 1.121124, 0.000000, 0.000000), tolerance=1e-4)

  expect_equal(unname(girf2$girf_res[[1]]$point_est[1,]), c(0.1317841, -0.8456145), tolerance=1e-4)
  expect_equal(unname(girf2$girf_res[[1]]$point_est[2,]), c(-3.333572, -1.151782), tolerance=1e-4)
  expect_equal(unname(girf2$girf_res[[1]]$conf_ints[2, 2, ]), c(-3.333572, -1.151782), tolerance=1e-4)

  expect_equal(unname(girf3$girf_res[[1]]$point_est[,1]), c(1.643883, 3.926812, 6.325788), tolerance=1e-4)
  expect_equal(unname(girf3$girf_res[[1]]$point_est[,2]), c(1.307096, 3.234228, 5.572740), tolerance=1e-4)
  expect_equal(unname(girf3$girf_res[[1]]$conf_ints[1, 2, ]), c(1.643883, 1.307096, 0.000000, 0.000000), tolerance=1e-4)

})
