context("gradient and Hessian")
library(gmvarkit)

foo1 <- function(x) x^2
foo2 <- function(x, a=1, b=1) a*x[1]^2 - b*x[2]^2
foo3 <- function(x) x[1]^2 + log(x[2]) - x[3]^3

test_that("calc_gradient works correctly", {
  expect_equal(calc_gradient(x=0, fn=foo1), 0, tolerance=1e-5)
  expect_equal(calc_gradient(x=1, fn=foo1), 2, tolerance=1e-5)
  expect_equal(calc_gradient(x=-2, fn=foo1), -4, tolerance=1e-5)

  expect_equal(calc_gradient(x=c(0, 0), fn=foo2), c(0, 0), tolerance=1e-5)
  expect_equal(calc_gradient(x=c(0, 0), fn=foo2, a=2, b=3), c(0, 0), tolerance=1e-5)
  expect_equal(calc_gradient(x=c(1, 2), fn=foo2), c(2, -4), tolerance=1e-5)
  expect_equal(calc_gradient(x=c(1, 2), fn=foo2, a=2, b=3), c(4, -12), tolerance=1e-5)

  expect_equal(calc_gradient(x=c(1, 2, 3), fn=foo3), c(2.0, 0.5, -27.0), tolerance=1e-4)
  expect_equal(calc_gradient(x=c(1, 2, 3), fn=foo3, varying_h=c(0.1, 0.2, 5)), c(2.0000000, 0.5016767, -52.0000000), tolerance=1e-4)
})


test_that("calc_hessian works correctly", {
  expect_equal(calc_hessian(x=0, fn=foo1), as.matrix(2), tolerance=1e-4)
  expect_equal(calc_hessian(x=1, fn=foo1), as.matrix(2), tolerance=1e-4)
  expect_equal(calc_hessian(x=-2, fn=foo1), as.matrix(2), tolerance=1e-4)

  expect_equal(calc_hessian(x=c(0, 0), fn=foo2), diag(c(2, -2)), tolerance=1e-4)
  expect_equal(calc_hessian(x=c(0, 0), fn=foo2, a=2, b=3), diag(c(4, -6)), tolerance=1e-4)
  expect_equal(calc_hessian(x=c(1, 2), fn=foo2), diag(c(2, -2)), tolerance=1e-4)
  expect_equal(calc_hessian(x=c(1, 2), fn=foo2, a=2, b=3), diag(c(4, -6)), tolerance=1e-4)

  expect_equal(calc_hessian(x=c(1, 2, 3), fn=foo3), diag(c(1.99998, -0.2500222, -18.00002)), tolerance=1e-4)
  expect_equal(calc_hessian(x=c(1, 2, 3), fn=foo3, varying_h=c(0.1, 0.02, 1)), diag(c(2.000000, -0.2500500, -18.00000)), tolerance=1e-4)
})

params112 <- c(0.1, 0.2, 0.1, 0.2, 0.3, 0.4, 0.1, 0, 0.2)
params112t <- c(params112, 1001)
params112t_2 <- c(params112, 10001)

test_that("get_varying_h works correctly", {
  expect_equal(get_varying_h(M=1, params=params112, model="GMVAR"), rep(6e-6, times=9), tol=1e-7)
  expect_equal(get_varying_h(M=1, params=params112t, model="StMVAR"), c(rep(6e-6, times=9), 1), tol=1e-7)
  expect_equal(get_varying_h(M=1, params=params112t_2, model="StMVAR"), c(rep(6e-6, times=9), 10), tol=1e-7)
})

mod112 <- GSMVAR(gdpdef, p=1, M=1, params=params112, model="GMVAR")
mod112t <- suppressWarnings(GSMVAR(gdpdef, p=1, M=1, params=params112t, model="StMVAR"))

test_that("get_foc works correctly", {
  expect_equal(get_foc(mod112), c(820.141040, 156.493420, 998.092351, -64.297637, 264.061712, 342.167590,
                                  8704.787988, -816.327571, -1.652966), tol=1e-1)
  expect_equal(get_foc(mod112t), c(802.620914, 141.015935, 967.384188, -65.006605, 258.535179, 321.811796,
                                   8284.248737, -761.540532, -45.888875, -0.021505), tol=1e-1)
})

test_that("get_soc works correctly", {
  expect_equal(get_soc(mod112), c(170.6500, -215.1800, -326.7867, -1193.5314, -1801.9150, -4532.0571,
                                  -5148.3655, -99218.4250, -187936.4309), tol=1)
})
