context("GIRF and GFEVD")
library(gmvarkit)

# NOTE that these tests contain random elements and the results might change is simulateGMVAR is modified

## A(M)(p)_(p)(M)(d)

# Structural GMVAR(2, 2), d=2 model identified with sign-constraints:
params22s <- c(0.36, 0.121, 0.484, 0.072, 0.223, 0.059, -0.151, 0.395, 0.406,
               -0.005, 0.083, 0.299, 0.218, 0.02, -0.119, 0.722, 0.093, 0.032,
               0.044, 0.191, 0.057, 0.172, -0.46, 0.016, 3.518, 5.154, 0.58)
W22 <- matrix(c(1, 1, -1, 1), nrow=2, byrow=FALSE)
mod22s <- GSMVAR(gdpdef, p=2, M=2, params=params22s, structural_pars=list(W=W22))
mod22ts <- GSMVAR(gdpdef, p=2, M=2, params=c(params22s, 20, 30), model="StMVAR", structural_pars=list(W=W22))
mod22gss <- GSMVAR(gdpdef, p=2, M=c(1, 1), params=c(params22s, 30), model="G-StMVAR", structural_pars=list(W=W22))

girf1 <- GIRF(mod22s, N=2, R2=2, R1=2, seeds=1:2, plot_res=FALSE)
girf2 <- GIRF(mod22s, which_shocks=2, shock_size=1, N=1, R2=1, R1=1, seeds=1,
              include_mixweights=FALSE, init_values=mod22s$data, ci=0.1, plot_res=FALSE)
girf3 <- GIRF(mod22s, N=2, R2=1, R1=1, which_cumulative=1:2, seeds=1, plot_res=FALSE)
girf4 <- GIRF(mod22s, N=2, R2=1, R1=1, scale=c(1, 1, 1), seeds=1, plot_res=FALSE)

girf5 <- GIRF(mod22ts, N=2, R2=2, R1=2, seeds=1:2, plot_res=FALSE)
girf6 <- GIRF(mod22ts, which_shocks=1, shock_size=1, N=1, R2=1, R1=1, seeds=1, which_cumulative=2,
              include_mixweights=FALSE, init_values=mod22s$data, ci=0.1, plot_res=FALSE)
girf7 <- GIRF(mod22gss, N=10, R2=6, R1=2, scale=c(2, 2, -1), seeds=1:6, plot_res=FALSE)
girf8 <- GIRF(mod22gss, N=10, R2=5, R1=2, scale=c(1, 2, 0.3), scale_type="peak", seeds=1:5, plot_res=FALSE)



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

  expect_equal(unname(girf5$girf_res[[1]]$point_est[3,]), c(1.669579e-01, -7.251490e-02, 2.065481e-05, -2.065481e-05), tolerance=1e-4)
  expect_equal(unname(girf5$girf_res[[2]]$point_est[1,]), c(-0.9514335, -0.1542835, 0.0000000, 0.0000000), tolerance=1e-4)
  expect_equal(unname(girf5$girf_res[[2]]$conf_ints[3, , 2]), c(-0.4840127, -0.4551930, -0.1477824, -0.1189627), tolerance=1e-4)
  expect_equal(unname(girf5$girf_res[[1]]$conf_ints[1, 1, ]), c(-0.7153089, -0.1987499, 0.0000000, 0.0000000), tolerance=1e-4)

  expect_equal(unname(girf6$girf_res[[1]]$point_est[1,]), c(0.4467698, 0.1979680), tolerance=1e-4)
  expect_equal(unname(girf6$girf_res[[1]]$point_est[2,]), c(0.2507398, 0.3915768), tolerance=1e-4)
  expect_equal(unname(girf6$girf_res[[1]]$conf_ints[2, 2, ]), c(0.2507398, 0.3915768), tolerance=1e-4)

  expect_equal(unname(girf7$girf_res[[1]]$point_est[11,]), c(0.032326891, 0.240915655, -0.001879026, 0.001879026), tolerance=1e-4)
  expect_equal(unname(girf7$girf_res[[2]]$point_est[11,]), c(-3.5349067, -0.4698592, -3.1878321, 3.1878321), tolerance=1e-4)
  expect_equal(unname(girf7$girf_res[[1]]$conf_ints[10, 2, ]), c(-0.5697602, -0.0763127, -0.2645382, -0.2644381), tolerance=1e-4)
  expect_equal(unname(girf7$girf_res[[2]]$conf_ints[10, 4, ]), c(4.016071, 3.616466, 4.504059, 20.080489), tolerance=1e-4)

  expect_equal(unname(girf8$girf_res[[1]]$point_est[11,]), c(-0.04691326, 0.07027594, 0.11717142, -0.11717142), tolerance=1e-4)
  expect_equal(unname(girf8$girf_res[[2]]$point_est[11,]), c(-0.05080895, 0.03214937, 0.09800529, -0.09800529), tolerance=1e-4)
  expect_equal(unname(girf8$girf_res[[1]]$conf_ints[10, 2, ]), c(-0.3768587, -0.1879175, -0.2536318, -0.3432418), tolerance=1e-4)
  expect_equal(unname(girf8$girf_res[[2]]$conf_ints[10, 4, ]), c(0.2626703, 0.3669083, 0.5024262, 0.1319560), tolerance=1e-4)
})


# Estimating the GFEVD using all possible histories in the data as the
# initial values:
mod22s_2 <- add_data(mod22s$data[1:50,], mod22s)
gfevd1 <- GFEVD(mod22s_2, N=3, R1=1, initval_type="data", seeds=1:49)
gfevd2 <- GFEVD(mod22s_2, N=2, R1=5, R2=10, initval_type="random", init_regimes=1:2, seeds=1:10)
gfevd3 <- GFEVD(mod22s_2, N=4, R1=6, initval_type="fixed",
                init_values=matrix(c(1, 1, 2, 2), nrow=2),
                include_mixweights=TRUE, seeds=1)

mod22ts_2 <- add_data(mod22ts$data[1:50,], mod22ts)
mod22gss_2 <- add_data(mod22gss$data[1:50,], mod22gss)

gfevd4 <- GFEVD(mod22ts_2, N=2, R1=2, initval_type="data", seeds=1:49)
gfevd5 <- GFEVD(mod22gss_2, N=3, R1=3, R2=10, initval_type="random", init_regimes=2, seeds=1:10)
gfevd6 <- GFEVD(mod22s_2, N=2, R1=2, initval_type="fixed",
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
  expect_equal(c(unname(gfevd4$gfevd_res[3, , ])), c(0.003385848, 0.996614152, 0.970176207, 0.029823793), tolerance=1e-3)
  expect_equal(c(unname(gfevd4$gfevd_res[2, , ])), c(0.003247969, 0.996752031, 0.978639508, 0.021360492), tolerance=1e-3)
  expect_equal(c(unname(gfevd5$gfevd_res[3, , ])), c(0.06959060, 0.93040940, 0.96216059, 0.03783941), tolerance=1e-3)
  expect_equal(c(unname(gfevd5$gfevd_res[1, , ])), c(0.02707136, 0.97292864, 0.97663429, 0.02336571), tolerance=1e-3)

  expect_equal(c(unname(gfevd6$gfevd_res[1, , 1:2])), c(0.092954461, 0.907045539, 0.998705168, 0.001294832), tolerance=1e-3)
  expect_equal(c(unname(gfevd6$gfevd_res[3, , ])), c(0.0894120535, 0.9105879465, 0.9990029740, 0.0009970260, 0.0007465000,
                                                     0.9992535000, 0.0007465005, 0.9992534995), tolerance=1e-3)
})
