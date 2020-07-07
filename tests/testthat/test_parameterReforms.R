context("Parameter reforms")
library(gmvarkit)

data0 <- t(t(eurusd))
data2 <- unname(data0)
data3 <- cbind(1:252, data2)
n_obs <- nrow(data3) # 252
dat2_1 <- reform_data(data2, p=1)
dat2_2 <- reform_data(data2, p=2)
dat2_3 <- reform_data(data2, p=3)
dat3_1 <- reform_data(data3, p=1)
dat3_2 <- reform_data(data3, p=2)

test_that("reform_data works correctly", {
  expect_equal(dat2_1, data2)

  expect_equal(nrow(dat2_2), n_obs - 2 + 1)
  expect_equal(dat2_2[1,], c(data2[2,], data2[1,]))
  expect_equal(dat2_2[23,], c(data2[24,], data2[23,]))
  expect_equal(dat2_2[nrow(dat2_2),], c(data2[252,], data2[251,]))

  expect_equal(nrow(dat2_3), n_obs - 3 + 1)
  expect_equal(dat2_3[1,], c(data2[3,], data2[2,], data2[1,]))
  expect_equal(dat2_3[100,], c(data2[102,], data2[101,], data2[100,]))
  expect_equal(dat2_3[nrow(dat2_3),], c(data2[252,], data2[251,], data2[250,]))

  expect_equal(dat3_1, data3)

  expect_equal(nrow(dat3_2), n_obs - 2 + 1)
  expect_equal(dat3_2[1,], c(data3[2,], data3[1,]))
  expect_equal(dat3_2[13,], c(data3[14,], data3[13,]))
  expect_equal(dat3_2[nrow(dat3_2),], c(data3[252,], data3[251,]))
})


## A(M)(p)_(p)(M)(d)

# p=1, M=1, d=2
phi10_112 <- c(1.03, 2.36)
A11_112 <- matrix(c(1.25, 0.06, 0.04, 1.34), nrow=2, byrow=FALSE)
Omega1_112 <- matrix(c(0.93, -0.15, -0.15, 5.20), nrow=2, byrow=FALSE)
theta_112 <- c(phi10_112, vec(A11_112), vech(Omega1_112))

# p=1, M=2, d=2
phi10_122 <- c(1.03, 2.36)
A11_122 <- matrix(c(1, -0.06, -0.04, 1), nrow=2, byrow=FALSE)
Omega1_122 <- matrix(c(0.93, -0.15, -0.15, 5.20), nrow=2, byrow=FALSE)
phi20_122 <- c(1.79, 3.00)
A21_122 <- A11_122 + 0.02
Omega2_122 <- matrix(c(5.88, 3.56, 3.56, 9.80), nrow=2, byrow=FALSE)
alpha1_122 <- 0.37
upsilon1_122 <- c(phi10_122, vec(A11_122), vech(Omega1_122))
upsilon2_122 <- c(phi20_122, vec(A21_122), vech(Omega2_122))
theta_122 <- c(upsilon1_122, upsilon2_122, alpha1_122)

W_122 <- matrix(c(-0.89, -0.71, 0.36, -2.10), nrow=2, ncol=2, byrow=FALSE)
lambdas_122 <- c(7.20, 1.30)
theta_122s <- c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), vec(W_122), lambdas_122, alpha1_122) # SGMVAR

# p=2, M=2, d=2
phi10_222 <- c(1.03, 2.36)
A11_222 <- matrix(c(1.25, 0.06, 0.04, 1.34), nrow=2, byrow=FALSE)
A12_222 <- matrix(c(-0.29, -0.08, -0.05, -0.36), nrow=2, byrow=FALSE)
Omega1_222 <- matrix(c(0.93, -0.15, -0.15, 5.20), nrow=2, byrow=FALSE)
phi20_222 <- c(1.79, 3.00)
A21_222 <- matrix(c(1.20, 0.05, 0.05, 1.30), nrow=2, byrow=FALSE)
A22_222 <- matrix(c(-0.30, -0.10, -0.05, -0.40), nrow=2, byrow=FALSE)
Omega2_222 <- matrix(c(5.88, 3.56, 3.56, 9.80), nrow=2, byrow=FALSE)
alpha1_222 <- 0.37
upsilon1_222 <- c(phi10_222, vec(A11_222), vec(A12_222), vech(Omega1_222))
upsilon2_222 <- c(phi20_222, vec(A21_222), vec(A22_222), vech(Omega2_222))
theta_222 <- c(upsilon1_222, upsilon2_222, alpha1_222)

# p=3, M=3, d=2
phi10_332 <- c(1.03, 2.36)
A11_332 <- matrix(c(1.25, 0.06, 0.04, 1.34), nrow=2, byrow=FALSE)
A12_332 <- matrix(c(-0.29, -0.08, -0.05, -0.36), nrow=2, byrow=FALSE)
A13_332 <- matrix(c(1, 0.21, 0.12, 2), nrow=2, byrow=FALSE)
Omega1_332 <- matrix(c(0.93, -0.15, -0.15, 5.20), nrow=2, byrow=FALSE)
phi20_332 <- c(1.79, 3.00)
A21_332 <- matrix(c(1.20, 0.05, 0.05, 1.30), nrow=2, byrow=FALSE)
A22_332 <- matrix(c(-0.30, -0.10, -0.05, -0.40), nrow=2, byrow=FALSE)
A23_332 <- matrix(c(0.30, 0.10, 0.05, 0.40), nrow=2, byrow=FALSE)
Omega2_332 <- matrix(c(5.88, 3.56, 3.56, 9.80), nrow=2, byrow=FALSE)
phi30_332 <- c(1.79, 3.00)
A31_332 <- matrix(c(1.30, 0.03, 0.08, 1.33), nrow=2, byrow=FALSE)
A32_332 <- matrix(c(-0.50, -0.20, -0.01, -0.40), nrow=2, byrow=FALSE)
A33_332 <- matrix(c(0.50, 0.20, 0.01, 0.40), nrow=2, byrow=FALSE)
Omega3_332 <- matrix(c(5.00, 3.00, 3.00, 9.00), nrow=2, byrow=FALSE)
alpha1_332 <- 0.5
alpha2_332 <- 0.3
upsilon1_332 <- c(phi10_332, vec(A11_332), vec(A12_332), vec(A13_332), vech(Omega1_332))
upsilon2_332 <- c(phi20_332, vec(A21_332), vec(A22_332), vec(A23_332), vech(Omega2_332))
upsilon3_332 <- c(phi30_332, vec(A31_332), vec(A32_332), vec(A33_332), vech(Omega3_332))
theta_332 <- c(upsilon1_332, upsilon2_332, upsilon3_332, alpha1_332, alpha2_332)

# p=1, M=2, d=3
phi10_123 <- c(1.1, 2.2, 3.3)
A11_123 <- matrix(c(1, 0.21, 0.31, 0.12, 2, 0.32, 0.13, 0.23, 3), nrow=3, byrow=FALSE)
Omega1_123 <- matrix(c(1, 0.22, 0.33, 0.22, 2, 0.44, 0.33, 0.44, 3), nrow=3, byrow=FALSE)
phi20_123 <- c(1.11, 2.22, 3.33)
A21_123 <- matrix(c(-1, -0.21, -0.31, -0.12, -2, -0.32, -0.13, -0.23, -3), nrow=3, byrow=FALSE)
Omega2_123 <- matrix(c(1.1, 0.222, 0.333, 0.222, 2.2, 0.444, 0.333, 0.444, 3.3), nrow=3, byrow=FALSE)
alpha1_123 <- 0.6
upsilon1_123 <- c(phi10_123, vec(A11_123), vech(Omega1_123))
upsilon2_123 <- c(phi20_123, vec(A21_123), vech(Omega2_123))
theta_123 <- c(upsilon1_123, upsilon2_123, alpha1_123)

# p=2, M=1, d=3
phi10_213 <- c(1.1, 2.2, 3.3)
A11_213 <- matrix(c(1, 0.21, 0.31, 0.12, 2, 0.32, 0.13, 0.23, 3), nrow=3, byrow=FALSE)
A12_213 <- matrix(c(-1, -0.21, -0.31, -0.12, -2, -0.32, -0.13, -0.23, -3), nrow=3, byrow=FALSE)
Omega1_213 <- matrix(c(1, 0.22, 0.33, 0.22, 2, 0.44, 0.33, 0.44, 3), nrow=3, byrow=FALSE)
upsilon1_213 <- c(phi10_213, vec(A11_213), vec(A12_213), vech(Omega1_213))
theta_213 <- upsilon1_213

## A(M)(p)_(p)(M)(d)
rbind_diags <- function(p, M, d) {
  I <- diag(p*d^2)
  Reduce(rbind, replicate(M, I, simplify=FALSE))
}

# Constraining AR-parameters to be the same for all regimes
C_112 <- rbind_diags(p=1, M=1, d=2)
theta_112c <- c(phi10_112, vec(A11_112), vech(Omega1_112))
C_122 <- rbind_diags(p=1, M=2, d=2)
theta_122c <- c(phi10_122, phi20_122, vec(A11_122), vech(Omega1_122), vech(Omega2_122), alpha1_122)
theta_122c_expanded <- c(phi10_122, vec(A11_122), vech(Omega1_122), phi20_122, vec(A11_122), vech(Omega2_122), alpha1_122)
C_222 <- rbind_diags(p=2, M=2, d=2)
theta_222c <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vech(Omega1_222), vech(Omega2_222), alpha1_222)
theta_222c_expanded <- c(phi10_222, vec(A11_222), vec(A12_222), vech(Omega1_222), phi20_222, vec(A11_222), vec(A12_222),
                         vech(Omega2_222), alpha1_222)
C_332 <- rbind_diags(p=3, M=3, d=2)
theta_332c <- c(phi10_332, phi20_332, phi30_332, vec(A11_332), vec(A12_332), vec(A13_332), vech(Omega1_332), vech(Omega2_332),
                vech(Omega3_332), alpha1_332, alpha2_332)
theta_332c_expanded <- c(phi10_332, vec(A11_332), vec(A12_332), vec(A13_332), vech(Omega1_332), phi20_332, vec(A11_332),
                         vec(A12_332), vec(A13_332), vech(Omega2_332), phi30_332, vec(A11_332), vec(A12_332), vec(A13_332),
                         vech(Omega3_332), alpha1_332, alpha2_332)
C_123 <- rbind_diags(p=1, M=2, d=3)
theta_123c <- c(phi10_123, phi20_123, vec(A11_123), vech(Omega1_123), vech(Omega2_123), alpha1_123)
theta_123c_expanded <- c(phi10_123, vec(A11_123), vech(Omega1_123), phi20_123, vec(A11_123), vech(Omega2_123), alpha1_123)
C_213 <- rbind_diags(p=2, M=1, d=3)
theta_213c <- c(phi10_213, vec(A11_213), vec(A12_213), vech(Omega1_213))

# p=2, M=2, d=2, constraint AR-parameters to be the same for all regimes
# and constraint the of-diagonal elements of AR-matrices to be zero.
mat0 <- matrix(c(1, rep(0, 10), 1, rep(0, 8), 1, rep(0, 10), 1), nrow=2*2^2, byrow=FALSE)
C_222_2 <- rbind(mat0, mat0)
A21_222_c2 <- A11_222_c2 <- matrix(c(1.26, 0, 0, 1.34), nrow=2, byrow=FALSE)
A22_222_c2 <- A12_222_c2 <- matrix(c(-0.29, 0, 0, -0.36), nrow=2, byrow=FALSE)
phi10_222_c2 <- c(-0.11, 2.83)
phi20_222_c2 <- c(0.36, 3.19)
Omega1_222_c2 <- matrix(c(0.98, -0.33, -0.33, 5.24), nrow=2, byrow=FALSE)
Omega2_222_c2 <- matrix(c(5.60, 3.46, 3.46, 9.62), nrow=2, byrow=FALSE)
alpha1_222_c2 <- 0.35
theta_222_c2 <- c(phi10_222_c2, phi20_222_c2, 1.26, 1.34, -0.29, -0.36, vech(Omega1_222_c2),
                  vech(Omega2_222_c2), alpha1_222_c2)
theta_222_c2_expanded <- c(phi10_222_c2, vec(A11_222_c2), vec(A12_222_c2), vech(Omega1_222_c2),
                           phi20_222_c2, vec(A21_222_c2), vec(A22_222_c2), vech(Omega2_222_c2),
                           alpha1_222_c2)

test_that("reform_constrained_pars works correctly", {
  expect_equal(reform_constrained_pars(p=1, M=1, d=2, params=theta_112c, constraints=C_112), theta_112)
  expect_equal(reform_constrained_pars(p=1, M=2, d=2, params=theta_122c, constraints=C_122), theta_122c_expanded)
  expect_equal(reform_constrained_pars(p=2, M=2, d=2, params=theta_222, constraints=NULL), theta_222)
  expect_equal(reform_constrained_pars(p=2, M=2, d=2, params=theta_222c, constraints=C_222), theta_222c_expanded)
  expect_equal(reform_constrained_pars(p=3, M=3, d=2, params=theta_332c, constraints=C_332), theta_332c_expanded)
  expect_equal(reform_constrained_pars(p=1, M=2, d=3, params=theta_123c, constraints=C_123), theta_123c_expanded)
  expect_equal(reform_constrained_pars(p=2, M=1, d=3, params=theta_213c, constraints=C_213), theta_213)

  expect_equal(reform_constrained_pars(p=2, M=2, d=2, params=theta_222_c2, constraints=C_222_2), theta_222_c2_expanded)
})


allA_112 <- pick_allA(p=1, M=1, d=2, params=theta_112)
allA_122 <- pick_allA(p=1, M=2, d=2, params=theta_122)
allA_222 <- pick_allA(p=2, M=2, d=2, params=theta_222)
allA_332 <- pick_allA(p=3, M=3, d=2, params=theta_332)
allA_123 <- pick_allA(p=1, M=2, d=3, params=theta_123)
allA_213 <- pick_allA(p=2, M=1, d=3, params=theta_213)


lower_part <- function(p, d) {
  cbind(diag(nrow=d*(p-1)), matrix(0, nrow=d*(p-1), ncol=d))
}

test_that("form_boldA works correctly", {
  expect_equal(form_boldA(p=1, M=1, d=2, all_A=allA_112)[, , 1], A11_112)

  expect_equal(form_boldA(p=1, M=2, d=2, all_A=allA_122)[, , 1], A11_122)
  expect_equal(form_boldA(p=1, M=2, d=2, all_A=allA_122)[, , 2], A21_122)

  expect_equal(form_boldA(p=2, M=2, d=2, all_A=allA_222)[, , 1], rbind(cbind(A11_222, A12_222), lower_part(p=2, d=2)))
  expect_equal(form_boldA(p=2, M=2, d=2, all_A=allA_222)[, , 2], rbind(cbind(A21_222, A22_222), lower_part(p=2, d=2)))

  expect_equal(form_boldA(p=3, M=3, d=2, all_A=allA_332)[, , 1], rbind(cbind(A11_332, A12_332, A13_332), lower_part(p=3, d=2)))
  expect_equal(form_boldA(p=3, M=3, d=2, all_A=allA_332)[, , 2], rbind(cbind(A21_332, A22_332, A23_332), lower_part(p=3, d=2)))
  expect_equal(form_boldA(p=3, M=3, d=2, all_A=allA_332)[, , 3], rbind(cbind(A31_332, A32_332, A33_332), lower_part(p=3, d=2)))

  expect_equal(form_boldA(p=1, M=2, d=3, all_A=allA_123)[, , 1], A11_123)
  expect_equal(form_boldA(p=1, M=2, d=3, all_A=allA_123)[, , 2], A21_123)

  expect_equal(form_boldA(p=2, M=1, d=3, all_A=allA_213)[, , 1], rbind(cbind(A11_213, A12_213), lower_part(p=2, d=3)))
})


alpha1_332_2 <- 0.1; alpha2_332_2 <- 0.4
theta_332_2 <- c(upsilon1_332, upsilon2_332, upsilon3_332, alpha1_332_2, alpha2_332_2)
alpha1_332_3 <- 0.5; alpha2_332_3 <- 0.2
theta_332_3 <- c(upsilon1_332, upsilon2_332, upsilon3_332, alpha1_332_3, alpha2_332_3)

upsilon1_342 <- upsilon1_332; upsilon2_342 <- upsilon2_332; upsilon3_342 <- upsilon3_332; upsilon4_342 <- upsilon3_342+0.07
alpha1_342 <- 0.1; alpha2_342 <- 0.2; alpha3_342 <- 0.3
theta_342 <- c(upsilon1_342, upsilon2_342, upsilon3_342, upsilon4_342, alpha1_342, alpha2_342, alpha3_342)
alpha1_342_2 <- 0.2; alpha2_342_2 <- 0.3; alpha3_342_2 <- 0.4
theta_342_2 <- c(upsilon1_342, upsilon2_342, upsilon3_342, upsilon4_342, alpha1_342_2, alpha2_342_2, alpha3_342_2)
alpha1_342_3 <- 0.3; alpha2_342_3 <- 0.4; alpha3_342_3 <- 0.1
theta_342_3 <- c(upsilon1_342, upsilon2_342, upsilon3_342, upsilon4_342, alpha1_342_3, alpha2_342_3, alpha3_342_3)

alpha1_123_2 <- 0.3
theta_123_2 <- c(upsilon1_123, upsilon2_123, alpha1_123_2)

test_that("sort_components works correctly", {
  expect_equal(sort_components(p=1, M=1, d=2, params=theta_112), theta_112)
  expect_equal(sort_components(p=1, M=2, d=2, params=theta_122), c(upsilon2_122, upsilon1_122, 1-alpha1_122))
  expect_equal(sort_components(p=2, M=2, d=2, params=theta_222), c(upsilon2_222, upsilon1_222, 1-alpha1_222))

  expect_equal(sort_components(p=3, M=3, d=2, params=theta_332), theta_332)
  expect_equal(sort_components(p=3, M=3, d=2, params=theta_332_2), c(upsilon3_332, upsilon2_332, upsilon1_332, 1-alpha1_332_2-alpha2_332_2, alpha2_332_2))
  expect_equal(sort_components(p=3, M=3, d=2, params=theta_332_3), c(upsilon1_332, upsilon3_332, upsilon2_332, alpha1_332_3, 1-alpha1_332_3-alpha2_332_3))

  expect_equal(sort_components(p=3, M=4, d=2, params=theta_342), c(upsilon4_342, upsilon3_342, upsilon2_342, upsilon1_342, 1-alpha1_342-alpha2_342-alpha3_342, alpha3_342, alpha2_342))
  expect_equal(sort_components(p=3, M=4, d=2, params=theta_342_2), c(upsilon3_342, upsilon2_342, upsilon1_342, upsilon4_342, alpha3_342_2, alpha2_342_2, alpha1_342_2))
  expect_equal(sort_components(p=3, M=4, d=2, params=theta_342_3), c(upsilon2_342, upsilon1_342, upsilon4_342, upsilon3_342, alpha2_342_3, alpha1_342_3, 1-alpha1_342_3-alpha2_342_3-alpha3_342_3))

  expect_equal(sort_components(p=1, M=2, d=3, params=theta_123), theta_123)
  expect_equal(sort_components(p=1, M=2, d=3, params=theta_123_2),  c(upsilon2_123, upsilon1_123, 1-alpha1_123_2))
  expect_equal(sort_components(p=2, M=1, d=3, params=theta_213), theta_213)
})


calc_mu <- function(p, M, d, params, constraints=NULL) {
  params <- reform_constrained_pars(p, M, d, params, constraints)
  all_A <- pick_allA(p=p, M=M, d=d, params=params)
  all_phi0 <- pick_phi0(p=p, M=M, d=d, params=params)
  vapply(1:M, function(m) solve(diag(d) - rowSums(all_A[, , , m, drop=FALSE], dims=2), all_phi0[,m]), numeric(d))
}

theta_112_mu <- change_parametrization(p=1, M=1, d=2, params=theta_112, change_to="mean")
theta_122_mu <- change_parametrization(p=1, M=2, d=2, params=theta_122, change_to="mean")
theta_222_mu <- change_parametrization(p=2, M=2, d=2, params=theta_222, change_to="mean")
theta_332_mu <- change_parametrization(p=3, M=3, d=2, params=theta_332, change_to="mean")
theta_342_mu <- change_parametrization(p=3, M=4, d=2, params=theta_342, change_to="mean")
theta_123_mu <- change_parametrization(p=1, M=2, d=3, params=theta_123, change_to="mean")
theta_213_mu <- change_parametrization(p=2, M=1, d=3, params=theta_213, change_to="mean")

theta_112c_mu <- change_parametrization(p=1, M=1, d=2, params=theta_112c, constraints=C_112, change_to="mean")
theta_222c_mu <- change_parametrization(p=2, M=2, d=2, params=theta_222c, constraints=C_222, change_to="mean")
theta_222c_mu2 <- change_parametrization(p=2, M=2, d=2, params=theta_222_c2, constraints=C_222_2, change_to="mean")
theta_123c_mu <- change_parametrization(p=1, M=2, d=3, params=theta_123c, constraints=C_123, change_to="mean")


test_that("change_parametrization works correctly", {
  expect_equal(pick_phi0(p=1, M=1, d=2, params=theta_112_mu), calc_mu(p=1, M=1, d=2, params=theta_112))
  expect_equal(change_parametrization(p=1, M=1, d=2, params=theta_112_mu, change_to="intercept"), theta_112)

  expect_equal(pick_phi0(p=1, M=2, d=2, params=theta_122_mu), calc_mu(p=1, M=2, d=2, params=theta_122))
  expect_equal(change_parametrization(p=1, M=2, d=2, params=theta_122_mu, change_to="intercept"), theta_122)

  expect_equal(pick_phi0(p=2, M=2, d=2, params=theta_222_mu), calc_mu(p=2, M=2, d=2, params=theta_222))
  expect_equal(change_parametrization(p=2, M=2, d=2, params=theta_222_mu, change_to="intercept"), theta_222)

  expect_equal(pick_phi0(p=3, M=3, d=2, params=theta_332_mu), calc_mu(p=3, M=3, d=2, params=theta_332))
  expect_equal(change_parametrization(p=3, M=3, d=2, params=theta_332_mu, change_to="intercept"), theta_332)

  expect_equal(pick_phi0(p=3, M=4, d=2, params=theta_342_mu), calc_mu(p=3, M=4, d=2, params=theta_342))
  expect_equal(change_parametrization(p=3, M=4, d=2, params=theta_342_mu, change_to="intercept"), theta_342)

  expect_equal(pick_phi0(p=1, M=2, d=3, params=theta_123_mu), calc_mu(p=1, M=2, d=3, params=theta_123))
  expect_equal(change_parametrization(p=1, M=2, d=3, params=theta_123_mu, change_to="intercept"), theta_123)

  expect_equal(pick_phi0(p=2, M=1, d=3, params=theta_213_mu), calc_mu(p=2, M=1, d=3, params=theta_213))
  expect_equal(change_parametrization(p=2, M=1, d=3, params=theta_213_mu, change_to="intercept"), theta_213)

  expect_equal(matrix(theta_112c_mu[1:(1*2)], nrow=2, byrow=FALSE), calc_mu(p=1, M=1, d=2, constraints=C_112, params=theta_112c))
  expect_equal(change_parametrization(p=1, M=1, d=2, params=theta_112c_mu, constraints=C_112, change_to="intercept"), theta_112c)

  expect_equal(matrix(theta_222c_mu[1:(2*2)], nrow=2, byrow=FALSE), calc_mu(p=2, M=2, d=2, constraints=C_222, params=theta_222c))
  expect_equal(change_parametrization(p=2, M=2, d=2, params=theta_222c_mu, constraints=C_222, change_to="intercept"), theta_222c)

  expect_equal(matrix(theta_222c_mu2[1:(2*2)], nrow=2, byrow=FALSE), calc_mu(p=2, M=2, d=2, constraints=C_222_2, params=theta_222_c2))
  expect_equal(change_parametrization(p=2, M=2, d=2, params=theta_222c_mu2, constraints=C_222_2, change_to="intercept"), theta_222_c2)

  expect_equal(matrix(theta_123c_mu[1:(2*3)], nrow=3, byrow=FALSE), calc_mu(p=1, M=2, d=3, constraints=C_123, params=theta_123c))
  expect_equal(change_parametrization(p=1, M=2, d=3, params=theta_123c_mu, constraints=C_123, change_to="intercept"), theta_123c)
})

theta_122_cr <- c(upsilon2_122, upsilon2_122, alpha1_122)
theta_222_cr <- c(upsilon1_222, upsilon1_222, alpha1_222)
theta_332_cr <- c(upsilon1_332, upsilon2_332, upsilon2_332, alpha1_332, alpha2_332)
theta_123_cr1 <- c(upsilon1_123, upsilon1_123, alpha1_123)
theta_123_cr2 <- c(upsilon2_123, upsilon2_123, alpha1_123)


test_that("change_regime works correctly", {
  expect_equal(change_regime(p=1, M=2, d=2, params=theta_122, m=1, regime_pars=upsilon2_122), theta_122_cr)
  expect_equal(change_regime(p=2, M=2, d=2, params=theta_222, m=2, regime_pars=upsilon1_222), theta_222_cr)
  expect_equal(change_regime(p=3, M=3, d=2, params=theta_332, m=3, regime_pars=upsilon2_332), theta_332_cr)
  expect_equal(change_regime(p=1, M=2, d=3, params=theta_123, m=2, regime_pars=upsilon1_123), theta_123_cr1)
  expect_equal(change_regime(p=1, M=2, d=3, params=theta_123, m=1, regime_pars=upsilon2_123), theta_123_cr2)
})

