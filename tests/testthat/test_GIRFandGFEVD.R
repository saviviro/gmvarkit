context("GIRF and GFEVD")
library(gmvarkit)

# NOTE that these tests contain random elements and the results might change is simulateGMVAR is modified

## A(M)(p)_(p)(M)(d)

# Structural GMVAR(2, 2), d=2 model identified with sign-constraints:
params22s <- c(0.36, 0.121, 0.484, 0.072, 0.223, 0.059, -0.151, 0.395, 0.406,
               -0.005, 0.083, 0.299, 0.218, 0.02, -0.119, 0.722, 0.093, 0.032,
               0.044, 0.191, 0.057, 0.172, -0.46, 0.016, 3.518, 5.154, 0.58)
W22 <- matrix(c(1, 1, -1, 1), nrow=2, byrow=FALSE)
mod22s <- GMVAR(gdpdef, p=2, M=2, params=params22s, structural_pars=list(W=W22))

girf1 <- GIRF(mod22s, N=2, R2=2, R1=2, seeds=1:2, plot=FALSE)
girf2 <- GIRF(mod22s, which_shocks=2, shock_size=1, N=1, R2=1, R1=1, seeds=1,
              include_mixweights=FALSE, init_values=mod22s$data, ci=0.1, plot=FALSE)
girf3 <- GIRF(mod22s, N=2, R2=1, R1=1, which_cumulative=1:2, seeds=1, plot=FALSE)
girf4 <- GIRF(mod22s, N=2, R2=1, R1=1, scale=c(1, 1, 1), seeds=1, plot=FALSE)

test_that("GIRF works correctly", {
  expect_equal(unname(girf1$girf_res[[1]]$point_est[3,]), c(-0.1014679, 0.1171628, -0.1047082, 0.1047082), tolerance=1e-4)
  expect_equal(unname(girf1$girf_res[[2]]$point_est[1,]), c(-0.57599848, 0.02003473, 0.00000000, 0.00000000), tolerance=1e-4)
  expect_equal(unname(girf1$girf_res[[2]]$conf_ints[3, , 2]), c(-0.009430567, -0.008801493, -0.002091365, -0.001462290), tolerance=1e-4)
  expect_equal(unname(girf1$girf_res[[1]]$conf_ints[1, 1, ]), c(0.06830714, 0.20611978, 0.00000000, 0.00000000), tolerance=1e-4)

  expect_equal(unname(girf2$girf_res[[1]]$point_est[1,]), c(0.124624506, -0.004334765), tolerance=1e-4)
  expect_equal(unname(girf2$girf_res[[1]]$point_est[2,]), c(0.028445814, 0.005640614), tolerance=1e-4)
  expect_equal(unname(girf2$girf_res[[1]]$conf_ints[2, 2, ]), c(0.028445814, 0.005640614), tolerance=1e-4)

  expect_equal(unname(girf3$girf_res[[1]]$point_est[,1]), c(0.11195283, 0.08590711, -0.41032582), tolerance=1e-4)
  expect_equal(unname(girf3$girf_res[[1]]$point_est[,2]), c(0.3378226, 0.4778677, 0.6470493), tolerance=1e-4)
  expect_equal(unname(girf3$girf_res[[1]]$conf_ints[1, 2, ]), c(0.1119528, 0.3378226, 0.0000000, 0.0000000), tolerance=1e-4)

  expect_equal(unname(girf4$girf_res[[1]]$point_est[,1]), c(1.0000000, -0.2326491, -4.4325177), tolerance=1e-4)
  expect_equal(unname(girf4$girf_res[[1]]$point_est[,2]), c(3.017544, 1.250930, 1.511186), tolerance=1e-4)
  expect_equal(unname(girf4$girf_res[[1]]$conf_ints[1, 2, ]), c(1.000000, 3.017544, 0.000000, 0.000000), tolerance=1e-4)
})


# Estimating the GFEVD using all possible histories in the data as the
# initial values:
mod22s_2 <- add_data(mod22s$data[1:50,], mod22s)
gfevd1 <- GFEVD(mod22s_2, N=3, R1=1, initval_type="data", seeds=1:49)
gfevd2 <- GFEVD(mod22s_2, N=2, R1=5, R2=10, initval_type="random", init_regimes=1:2, seeds=1:10)
gfevd3 <- GFEVD(mod22s_2, N=4, R1=6, initval_type="fixed",
                init_values=matrix(c(1, 1, 2, 2), nrow=2),
                include_mixweights=TRUE, seeds=1)

test_that("GFEVD works correctly", {
  expect_equal(unname(gfevd1$gfevd_res[3, , 1]), c(0.01382024, 0.98617976), tolerance=1e-3)
  expect_equal(unname(gfevd1$gfevd_res[1, , 2]), c(0.98893981, 0.01106019), tolerance=1e-3)
  expect_equal(unname(gfevd2$gfevd_res[3, , 1]), c(0.00793223, 0.99206777), tolerance=1e-3)
  expect_equal(unname(gfevd2$gfevd_res[2, , 2]), c(0.98524421, 0.01475579), tolerance=1e-3)
  expect_equal(unname(gfevd3$gfevd_res[2, , 1]), c(0.01929081, 0.98070919), tolerance=1e-3)
  expect_equal(unname(gfevd3$gfevd_res[1, , 2]), c(0.993477615, 0.006522385), tolerance=1e-3)
  expect_equal(unname(gfevd3$gfevd_res[5, , 3]), c(0.2543561, 0.7456439), tolerance=1e-3)
  expect_equal(unname(gfevd3$gfevd_res[4, , 4]), c(0.7643379, 0.2356621), tolerance=1e-3)

})
