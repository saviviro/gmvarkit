context("GIRF and GFEVD")
library(gmvarkit)

# NOTE that these tests contain random elements and the results might change is simulateGMVAR is modified

## A(M)(p)_(p)(M)(d)

# Reduced form linear p=1, M=1, d=2
phi10_112 <- c(0.65, 0.7)
A11_112 <- matrix(c(0.29, 0.02, -0.14, 0.9), nrow=2, byrow=FALSE)
Omega1_112 <- matrix(c(0.60, 0.01, 0.01, 0.07), nrow=2, byrow=FALSE)

theta_112t <- c(phi10_112, vec(A11_112), vech(Omega1_112), 3)
mod11t <- GSMVAR(gdpdef, p=1, M=1, params=theta_112t, model="StMVAR")


# Reduced form p=1, M=2, d=2 model identified by Cholesky
phi10_122 <- c(0.55, 0.11)
A11_122 <- matrix(c(0.34, 0.05, -0.01, 0.72), nrow=2, byrow=FALSE)
Omega1_122 <- matrix(c(0.58, 0.01, 0.01, 0.06), nrow=2, byrow=FALSE)

phi20_122 <- c(0.17, 0.25)
A21_122 <- A11_122
Omega2_122 <- matrix(c(0.50, -0.01, -0.01, 0.20), nrow=2, byrow=FALSE)

alpha1_122 <- 0.60
upsilon1_122 <- c(phi10_122, vec(A11_122), vech(Omega1_122))
upsilon2_122 <- c(phi20_122, vec(A21_122), vech(Omega2_122))
theta_122 <- c(upsilon1_122, upsilon2_122, alpha1_122)

mod12 <- GSMVAR(gdpdef, p=1, M=2, params=theta_122)
mod12t <- GSMVAR(gdpdef, p=1, M=2, params=c(theta_122, 3, 12), model="StMVAR")
mod12gs <- GSMVAR(gdpdef, p=1, M=c(1, 1), params=c(theta_122, 3), model="G-StMVAR")


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

# Reduced form mod
girf9 <- GIRF(mod12, N=2, R2=2, R1=2, seeds=1:2, plot_res=FALSE)
girf10 <- GIRF(mod12t, which_shocks=2, shock_size=1, N=1, R2=1, R1=1, seeds=1,
              include_mixweights=FALSE, init_values=mod12$data, ci=0.2, plot_res=FALSE)
girf11 <- GIRF(mod12gs, N=2, R2=3, R1=4, which_cumulative=1:2, seeds=11:13, plot_res=FALSE)
girf12 <- GIRF(mod11t, N=2, R2=1, R1=1, scale=c(1, 1, 1), seeds=5, plot_res=FALSE)


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

  expect_equal(unname(girf5$girf_res[[1]]$point_est[3,]), c(-0.0020855158, 0.0702707899, -0.0009616155, 0.0009616155), tolerance=1e-4)
  expect_equal(unname(girf5$girf_res[[2]]$point_est[1,]), c(-0.8900453, 0.0309581, 0.0000000, 0.0000000), tolerance=1e-4)
  expect_equal(unname(girf5$girf_res[[2]]$conf_ints[3, , 2]), c(-0.0516066644, -0.0472164259, -0.0003872155, 0.0040030230), tolerance=1e-4)
  expect_equal(unname(girf5$girf_res[[1]]$conf_ints[1, 1, ]), c(-0.05380848, -0.16236944, 0.00000000, 0.00000000), tolerance=1e-4)

  expect_equal(unname(girf6$girf_res[[1]]$point_est[1,]), c(0.06994922, 0.21107483), tolerance=1e-4)
  expect_equal(unname(girf6$girf_res[[1]]$point_est[2,]), c(-0.01780591, 0.29782252), tolerance=1e-4)
  expect_equal(unname(girf6$girf_res[[1]]$conf_ints[2, 2, ]), c(-0.01780591, 0.29782252), tolerance=1e-4)

  expect_equal(unname(girf7$girf_res[[1]]$point_est[11,]), c(0.021456821, 0.006198903, 0.057140310, -0.057140310), tolerance=1e-4)
  expect_equal(unname(girf7$girf_res[[2]]$point_est[11,]), c(-2.012376, 3.994086, -5.961740, 5.961740), tolerance=1e-4)
  expect_equal(unname(girf7$girf_res[[1]]$conf_ints[10, 2, ]), c(-0.07265945, -0.12302023, -0.07657784, -0.13035412), tolerance=1e-4)
  expect_equal(unname(girf7$girf_res[[2]]$conf_ints[10, 4, ]), c(1.26414556, 15.81447560, -0.08471768, 6.51993375), tolerance=1e-4)

  expect_equal(unname(girf8$girf_res[[1]]$point_est[11,]), c(0.03833223, -0.04198505, 0.04963786, -0.04963786), tolerance=1e-4)
  expect_equal(unname(girf8$girf_res[[2]]$point_est[11,]), c(0.03162255, -0.09872999, 0.09132550, -0.09132550), tolerance=1e-4)
  expect_equal(unname(girf8$girf_res[[1]]$conf_ints[10, 2, ]), c(-0.3827951, -0.2052065, -0.4112630, -0.1960736), tolerance=1e-4)
  expect_equal(unname(girf8$girf_res[[2]]$conf_ints[10, 4, ]), c(0.276011894, -0.016043945, 0.162048134, -0.007608629), tolerance=1e-4)

  expect_equal(unname(girf9$girf_res[[1]]$point_est[2,]), c(0.47129726, 0.03285889, 0.17361305, -0.17361305), tolerance=1e-4)
  expect_equal(unname(girf9$girf_res[[2]]$point_est[3,]), c(-0.01661809, 0.13433957, -0.20154580, 0.20154580), tolerance=1e-4)
  expect_equal(unname(girf9$girf_res[[1]]$conf_ints[1, 1, ]), c(0.6019734824, 0.0008577321, 0.0000000000, 0.0000000000), tolerance=1e-4)
  expect_equal(unname(girf9$girf_res[[2]]$conf_ints[3, 4, ]), c(0.0239029, 0.2047201, -0.1177271, 0.2853645), tolerance=1e-4)

  expect_equal(unname(girf10$girf_res[[1]]$point_est[1,]), c(0.00000000, 0.09606795), tolerance=1e-4)
  expect_equal(unname(girf10$girf_res[[1]]$conf_ints[2, 2, ]), c(-0.01311336, 0.07134762), tolerance=1e-4)

  expect_equal(unname(girf11$girf_res[[1]]$point_est[1,]), c(1.09778453, -0.02195569, 0.00000000, 0.00000000), tolerance=1e-4)
  expect_equal(unname(girf11$girf_res[[2]]$point_est[2,]), c(0.08784187, 0.50221556, -0.01577926, 0.01577926), tolerance=1e-4)
  expect_equal(unname(girf11$girf_res[[1]]$conf_ints[2, 3, ]), c(3.502104325, 0.006528326, 0.345689728, 0.001160202), tolerance=1e-4)
  expect_equal(unname(girf11$girf_res[[2]]$conf_ints[3, 2, ]), c(-0.3986564, 0.1166099, -0.1710915, -0.1577015), tolerance=1e-4)

  expect_equal(unname(girf12$girf_res[[1]]$point_est[3,]), c(0.01806770, 0.03962662), tolerance=1e-4)
  expect_equal(unname(girf12$girf_res[[1]]$conf_ints[2, 2, ]), c(0.24309474, 0.04184041), tolerance=1e-4)
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

mod11t_2 <- add_data(mod11t$data[1:50,], mod11t)
mod12_2 <- add_data(mod12$data[1:50,], mod12)
mod12t_2 <- add_data(mod12t$data[1:50,], mod12t)
mod12gs_2 <- add_data(mod12gs$data[1:50,], mod12gs)

gfevd7 <- GFEVD(mod11t_2, N=3, R1=1, initval_type="data", seeds=1:50)
gfevd8 <- GFEVD(mod12_2, N=2, R1=5, R2=8, initval_type="random", init_regimes=1:2, seeds=1:8)
gfevd9 <- GFEVD(mod12t_2, N=4, R1=6, initval_type="fixed",
                init_values=matrix(c(1, 1, 2, 2), nrow=2),
                include_mixweights=TRUE, seeds=1)
gfevd10 <- GFEVD(mod12gs_2, N=2, R1=2, initval_type="data", seeds=1:50)


test_that("GFEVD works correctly", {
  expect_equal(unname(gfevd1$gfevd_res[3, , 1]), c(0.01382024, 0.98617976), tolerance=1e-3)
  expect_equal(unname(gfevd1$gfevd_res[1, , 2]), c(0.98893981, 0.01106019), tolerance=1e-3)
  expect_equal(unname(gfevd2$gfevd_res[3, , 1]), c(0.00793223, 0.99206777), tolerance=1e-3)
  expect_equal(unname(gfevd2$gfevd_res[2, , 2]), c(0.98524421, 0.01475579), tolerance=1e-3)
  expect_equal(unname(gfevd3$gfevd_res[2, , 1]), c(0.01929081, 0.98070919), tolerance=1e-3)
  expect_equal(unname(gfevd3$gfevd_res[1, , 2]), c(0.993477615, 0.006522385), tolerance=1e-3)
  expect_equal(unname(gfevd3$gfevd_res[5, , 3]), c(0.2543561, 0.7456439), tolerance=1e-3)
  expect_equal(unname(gfevd3$gfevd_res[4, , 4]), c(0.7643379, 0.2356621), tolerance=1e-3)

  expect_equal(c(unname(gfevd4$gfevd_res[3, , ])), c(0.01213561, 0.98786439, 0.98451084, 0.01548916), tolerance=1e-3)
  expect_equal(c(unname(gfevd4$gfevd_res[2, , ])), c(0.01167646, 0.98832354, 0.98816800, 0.01183200), tolerance=1e-3)
  expect_equal(c(unname(gfevd5$gfevd_res[3, , ])), c(0.006614715, 0.993385285, 0.985963802, 0.014036198), tolerance=1e-3)
  expect_equal(c(unname(gfevd5$gfevd_res[1, , ])), c(0.006968765, 0.993031235, 0.981418598, 0.018581402), tolerance=1e-3)
  expect_equal(c(unname(gfevd6$gfevd_res[1, , 1:2])), c(0.092954461, 0.907045539, 0.998705168, 0.001294832), tolerance=1e-3)
  expect_equal(c(unname(gfevd6$gfevd_res[3, , ])), c(0.0894120535, 0.9105879465, 0.9990029740, 0.0009970260, 0.0007465000,
                                                     0.9992535000, 0.0007464997, 0.9992535003), tolerance=1e-3)

  expect_equal(unname(gfevd7$gfevd_res[3, , 1]), c(0.996256034, 0.003743966), tolerance=1e-3)
  expect_equal(unname(gfevd7$gfevd_res[1, , 2]), c(0.004839729, 0.995160271), tolerance=1e-3)
  expect_equal(unname(gfevd8$gfevd_res[3, , 1]), c(0.996935931, 0.003064069), tolerance=1e-3)
  expect_equal(unname(gfevd8$gfevd_res[2, , 2]), c(0.0005259957, 0.9994740043), tolerance=1e-3)
  expect_equal(unname(gfevd9$gfevd_res[2, , 1]), c(9.999490e-01, 5.104191e-05), tolerance=1e-3)
  expect_equal(unname(gfevd9$gfevd_res[1, , 2]), c(0.0007162443, 0.9992837557), tolerance=1e-3)
  expect_equal(unname(gfevd9$gfevd_res[5, , 3]), c(0.5799644, 0.4200356), tolerance=1e-3)
  expect_equal(unname(gfevd9$gfevd_res[4, , 4]), c(0.7555897, 0.2444103), tolerance=1e-3)
  expect_equal(c(unname(gfevd10$gfevd_res[3, , ])), c(0.990785908, 0.009214092, 0.001832945, 0.998167055), tolerance=1e-3)
  expect_equal(c(unname(gfevd10$gfevd_res[2, , ])), c(0.991911836, 0.008088164, 0.002087102, 0.997912898), tolerance=1e-3)
})
