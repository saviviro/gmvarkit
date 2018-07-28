context("Picking parameters from parameter vector")
library(gmvarkit)


## A(M)(p)_(p)(M)(d)

# p=1, M=1, d=2
phi10_112 <- c(1.03, 2.36)
A11_112 <- matrix(c(1.25, 0.06, 0.04, 1.34), nrow=2, byrow=FALSE)
Omega1_112 <- matrix(c(0.93, -0.15, -0.15, 5.20), nrow=2, byrow=FALSE)

theta_112 <- upsilon1_112 <- c(phi10_112, vec(A11_112), vech(Omega1_112))

# p=1, M=2, d=2
phi10_122 <- c(1.03, 2.36)
A11_122 <- matrix(c(1, -0.06, -0.04, 1), nrow=2, byrow=FALSE)
Omega1_122 <- matrix(c(0.93, -0.15, -0.15, 5.20), nrow=2, byrow=FALSE)

phi20_122 <- c(1.79, 3.00)
A21_122 <- A11_122
Omega2_122 <- matrix(c(5.88, 3.56, 3.56, 9.80), nrow=2, byrow=FALSE)

alpha1_122 <- 0.37
upsilon1_122 <- c(phi10_122, vec(A11_122), vech(Omega1_122))
upsilon2_122 <- c(phi20_122, vec(A21_122), vech(Omega2_122))
theta_122 <- c(upsilon1_122, upsilon2_122, alpha1_122)

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


test_that("pick_Ami works correctly", {
  expect_equal(pick_Ami(p=1, M=1, d=2, params=theta_112, m=1, i=1, unvec=TRUE), A11_112)

  expect_equal(pick_Ami(p=1, M=2, d=2, params=theta_122, m=1, i=1, unvec=TRUE), A11_122)
  expect_equal(pick_Ami(p=1, M=2, d=2, params=theta_122, m=2, i=1, unvec=TRUE), A21_122)

  expect_equal(pick_Ami(p=2, M=2, d=2, params=theta_222, m=1, i=1, unvec=TRUE), A11_222)
  expect_equal(pick_Ami(p=2, M=2, d=2, params=theta_222, m=2, i=2), A22_222)

  expect_equal(pick_Ami(p=3, M=3, d=2, params=theta_332, m=1, i=2, unvec=TRUE), A12_332)
  expect_equal(pick_Ami(p=3, M=3, d=2, params=theta_332, m=2, i=3, unvec=TRUE), A23_332)
  expect_equal(pick_Ami(p=3, M=3, d=2, params=theta_332, m=3, i=1, unvec=TRUE), A31_332)

  expect_equal(pick_Ami(p=1, M=2, d=3, params=theta_123, m=1, i=1, unvec=TRUE), A11_123)
  expect_equal(pick_Ami(p=1, M=2, d=3, params=theta_123, m=2, i=1, unvec=TRUE), A21_123)

  expect_equal(pick_Ami(p=2, M=1, d=3, params=theta_213, m=1, i=1, unvec=TRUE), A11_213)
  expect_equal(pick_Ami(p=2, M=1, d=3, params=theta_213, m=1, i=2, unvec=TRUE), A12_213)
})


test_that("pick_Am works correctly", {
  expect_equal(pick_Am(p=1, M=1, d=2, params=theta_112, m=1)[, , 1], A11_112)

  expect_equal(pick_Am(p=1, M=2, d=2, params=theta_122, m=1)[, , 1], A11_122)
  expect_equal(pick_Am(p=1, M=2, d=2, params=theta_122, m=2)[, , 1], A21_122)

  expect_equal(pick_Am(p=2, M=2, d=2, params=theta_222, m=1)[, , 2], A12_222)
  expect_equal(pick_Am(p=2, M=2, d=2, params=theta_222, m=2)[, , 1], A21_222)

  expect_equal(pick_Am(p=3, M=3, d=2, params=theta_332, m=1)[, , 3], A13_332)
  expect_equal(pick_Am(p=3, M=3, d=2, params=theta_332, m=2)[, , 2], A22_332)
  expect_equal(pick_Am(p=3, M=3, d=2, params=theta_332, m=3)[, , 1], A31_332)

  expect_equal(pick_Am(p=1, M=2, d=3, params=theta_123, m=1)[, , 1], A11_123)
  expect_equal(pick_Am(p=1, M=2, d=3, params=theta_123, m=2)[, , 1], A21_123)

  expect_equal(pick_Am(p=2, M=1, d=3, params=theta_213, m=1)[, , 1], A11_213)
  expect_equal(pick_Am(p=2, M=1, d=3, params=theta_213, m=1)[, , 2], A12_213)
})


test_that("pick_allA works correctly", {
  expect_equal(pick_allA(p=1, M=1, d=2, params=theta_112)[, , 1, 1], A11_112)

  expect_equal(pick_allA(p=1, M=2, d=2, params=theta_122)[, , 1, 1], A11_122)
  expect_equal(pick_allA(p=1, M=2, d=2, params=theta_122)[, , 1, 2], A21_122)

  expect_equal(pick_allA(p=2, M=2, d=2, params=theta_222)[, , 1, 1], A11_222)
  expect_equal(pick_allA(p=2, M=2, d=2, params=theta_222)[, , 2, 1], A12_222)
  expect_equal(pick_allA(p=2, M=2, d=2, params=theta_222)[, , 2, 2], A22_222)

  expect_equal(pick_allA(p=3, M=3, d=2, params=theta_332)[, , 3, 1], A13_332)
  expect_equal(pick_allA(p=3, M=3, d=2, params=theta_332)[, , 1, 2], A21_332)
  expect_equal(pick_allA(p=3, M=3, d=2, params=theta_332)[, , 2, 3], A32_332)
  expect_equal(pick_allA(p=3, M=3, d=2, params=theta_332)[, , 3, 3], A33_332)

  expect_equal(pick_allA(p=1, M=2, d=3, params=theta_123)[, , 1, 1], A11_123)
  expect_equal(pick_allA(p=1, M=2, d=3, params=theta_123)[, , 1, 2], A21_123)

  expect_equal(pick_allA(p=2, M=1, d=3, params=theta_213)[, , 1, 1], A11_213)
  expect_equal(pick_allA(p=2, M=1, d=3, params=theta_213)[, , 2, 1], A12_213)
})


test_that("pick_phi0 works correctly", {
  expect_equal(pick_phi0(p=1, M=1, d=2, params=theta_112)[, 1], phi10_112)

  expect_equal(pick_phi0(p=1, M=2, d=2, params=theta_122)[, 1], phi10_122)
  expect_equal(pick_phi0(p=1, M=2, d=2, params=theta_122)[, 2], phi20_122)

  expect_equal(pick_phi0(p=2, M=2, d=2, params=theta_222)[, 1], phi10_222)
  expect_equal(pick_phi0(p=2, M=2, d=2, params=theta_222)[, 2], phi20_222)

  expect_equal(pick_phi0(p=3, M=3, d=2, params=theta_332)[, 1], phi10_332)
  expect_equal(pick_phi0(p=3, M=3, d=2, params=theta_332)[, 2], phi20_332)
  expect_equal(pick_phi0(p=3, M=3, d=2, params=theta_332)[, 3], phi30_332)

  expect_equal(pick_phi0(p=1, M=2, d=3, params=theta_123)[, 1], phi10_123)
  expect_equal(pick_phi0(p=1, M=2, d=3, params=theta_123)[, 2], phi20_123)

  expect_equal(pick_phi0(p=2, M=1, d=3, params=theta_213)[, 1], phi10_213)
})


test_that("pick_all_phi0_A works correctly", {
  expect_equal(pick_all_phi0_A(p=1, M=1, d=2, params=theta_112), as.matrix(c(phi10_112, vec(A11_112))))
  expect_equal(pick_all_phi0_A(p=1, M=2, d=2, params=theta_122), cbind(c(phi10_122, vec(A11_122)),
                                                                       c(phi20_122, vec(A21_122))))
  expect_equal(pick_all_phi0_A(p=2, M=2, d=2, params=theta_222), cbind(c(phi10_222, vec(A11_222), vec(A12_222)),
                                                                       c(phi20_222, vec(A21_222), vec(A22_222))))
  expect_equal(pick_all_phi0_A(p=3, M=3, d=2, params=theta_332), cbind(c(phi10_332, vec(A11_332), vec(A12_332), vec(A13_332)),
                                                                       c(phi20_332, vec(A21_332), vec(A22_332), vec(A23_332)),
                                                                       c(phi30_332, vec(A31_332), vec(A32_332), vec(A33_332))))
  expect_equal(pick_all_phi0_A(p=1, M=2, d=3, params=theta_123), cbind(c(phi10_123, vec(A11_123)),
                                                                       c(phi20_123, vec(A21_123))))
  expect_equal(pick_all_phi0_A(p=2, M=1, d=3, params=theta_213), as.matrix(c(phi10_213, vec(A11_213), vec(A12_213))))
})


test_that("pick_Omegas works correctly", {
  expect_equal(pick_Omegas(p=1, M=1, d=2, params=theta_112)[, , 1], Omega1_112)

  expect_equal(pick_Omegas(p=1, M=2, d=2, params=theta_122)[, , 1], Omega1_122)
  expect_equal(pick_Omegas(p=1, M=2, d=2, params=theta_122)[, , 2], Omega2_122)

  expect_equal(pick_Omegas(p=2, M=2, d=2, params=theta_222)[, , 1], Omega1_222)
  expect_equal(pick_Omegas(p=2, M=2, d=2, params=theta_222)[, , 2], Omega2_222)

  expect_equal(pick_Omegas(p=3, M=3, d=2, params=theta_332)[, , 1], Omega1_332)
  expect_equal(pick_Omegas(p=3, M=3, d=2, params=theta_332)[, , 2], Omega2_332)
  expect_equal(pick_Omegas(p=3, M=3, d=2, params=theta_332)[, , 3], Omega3_332)

  expect_equal(pick_Omegas(p=1, M=2, d=3, params=theta_123)[, , 1], Omega1_123)
  expect_equal(pick_Omegas(p=1, M=2, d=3, params=theta_123)[, , 2], Omega2_123)

  expect_equal(pick_Omegas(p=2, M=1, d=3, params=theta_213)[, , 1], Omega1_213)
})


test_that("pick_alphas works correctly", {
  expect_equal(pick_alphas(p=1, M=1, d=2, params=theta_112), 1)
  expect_equal(pick_alphas(p=1, M=2, d=2, params=theta_122), c(alpha1_122, 1-alpha1_122))
  expect_equal(pick_alphas(p=2, M=2, d=2, params=theta_222), c(alpha1_222, 1-alpha1_222))
  expect_equal(pick_alphas(p=3, M=3, d=2, params=theta_332), c(alpha1_332, alpha2_332, 1-alpha1_332-alpha2_332))
  expect_equal(pick_alphas(p=1, M=2, d=3, params=theta_123), c(alpha1_123, 1-alpha1_123))
  expect_equal(pick_alphas(p=2, M=1, d=3, params=theta_213), 1)
})


test_that("pick_regime works correctly", {
  expect_equal(pick_regime(p=1, M=1, d=2, params=theta_112, m=1), upsilon1_112)

  expect_equal(pick_regime(p=1, M=2, d=2, params=theta_122, m=1), upsilon1_122)
  expect_equal(pick_regime(p=1, M=2, d=2, params=theta_122, m=2), upsilon2_122)

  expect_equal(pick_regime(p=2, M=2, d=2, params=theta_222, m=1), upsilon1_222)
  expect_equal(pick_regime(p=2, M=2, d=2, params=theta_222, m=2), upsilon2_222)


  expect_equal(pick_regime(p=3, M=3, d=2, params=theta_332, m=1), upsilon1_332)
  expect_equal(pick_regime(p=3, M=3, d=2, params=theta_332, m=2), upsilon2_332)
  expect_equal(pick_regime(p=3, M=3, d=2, params=theta_332, m=3), upsilon3_332)

  expect_equal(pick_regime(p=1, M=2, d=3, params=theta_123, m=1), upsilon1_123)
  expect_equal(pick_regime(p=1, M=2, d=3, params=theta_123, m=2), upsilon2_123)

  expect_equal(pick_regime(p=2, M=1, d=3, params=theta_213, m=1), upsilon1_213)
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
theta_123_mu <- change_parametrization(p=1, M=2, d=3, params=theta_123, change_to="mean")
theta_213_mu <- change_parametrization(p=2, M=1, d=3, params=theta_213, change_to="mean")


## A(M)(p)_(p)(M)(d)
rbind_diags <- function(p, M, d) {
  I <- diag(p*d^2)
  Reduce(rbind, replicate(M, I, simplify=FALSE))
}

# Constraining AR-parameters to be the same for all regimes
C_112 <- rbind_diags(p=1, M=1, d=2)
theta_112c <- c(phi10_112, vec(A11_112), vech(Omega1_112))
C_222 <- rbind_diags(p=2, M=2, d=2)
theta_222c <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vech(Omega1_222), vech(Omega2_222), alpha1_222)
C_123 <- rbind_diags(p=1, M=2, d=3)
theta_123c <- c(phi10_123, phi20_123, vec(A11_123), vech(Omega1_123), vech(Omega2_123), alpha1_123)

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

theta_112c_mu <- change_parametrization(p=1, M=1, d=2, params=theta_112c, constraints=C_112, change_to="mean")
theta_222c_mu <- change_parametrization(p=2, M=2, d=2, params=theta_222c, constraints=C_222, change_to="mean")
theta_222c_mu2 <- change_parametrization(p=2, M=2, d=2, params=theta_222_c2, constraints=C_222_2, change_to="mean")
theta_123c_mu <- change_parametrization(p=1, M=2, d=3, params=theta_123c, constraints=C_123, change_to="mean")

theta_112c_mu_exp <- reform_constrained_pars(p=1, M=1, d=2, params=theta_112c_mu, constraints=C_112)
theta_222c_mu_exp <- reform_constrained_pars(p=2, M=2, d=2, params=theta_222c_mu, constraints=C_222)
theta_222c_mu2_exp <- reform_constrained_pars(p=2, M=2, d=2, params=theta_222c_mu2, constraints=C_222_2)
theta_123c_mu_exp <- reform_constrained_pars(p=1, M=2, d=3, params=theta_123c_mu, constraints=C_123)


test_that("get_regime_means_int works correctly", {
  expect_equal(get_regime_means_int(p=1, M=1, d=2, params=theta_112, parametrization="intercept", constraints=NULL),
               pick_phi0(p=1, M=1, d=2, params=theta_112_mu))
  expect_equal(get_regime_means_int(p=1, M=1, d=2, params=theta_112_mu, parametrization="mean", constraints=NULL),
               calc_mu(p=1, M=1, d=2, params=theta_112))
  expect_equal(get_regime_means_int(p=1, M=2, d=2, params=theta_122, parametrization="intercept", constraints=NULL),
               pick_phi0(p=1, M=2, d=2, params=theta_122_mu))
  expect_equal(get_regime_means_int(p=1, M=2, d=2, params=theta_122_mu, parametrization="mean", constraints=NULL),
              calc_mu(p=1, M=2, d=2, params=theta_122))
  expect_equal(get_regime_means_int(p=2, M=2, d=2, params=theta_222, parametrization="intercept", constraints=NULL),
               pick_phi0(p=2, M=2, d=2, params=theta_222_mu))
  expect_equal(get_regime_means_int(p=2, M=2, d=2, params=theta_222_mu, parametrization="mean", constraints=NULL),
               calc_mu(p=2, M=2, d=2, params=theta_222))
  expect_equal(get_regime_means_int(p=3, M=3, d=2, params=theta_332, parametrization="intercept", constraints=NULL),
               pick_phi0(p=3, M=3, d=2, params=theta_332_mu))
  expect_equal(get_regime_means_int(p=3, M=3, d=2, params=theta_332_mu, parametrization="mean", constraints=NULL),
               calc_mu(p=3, M=3, d=2, params=theta_332))
  expect_equal(get_regime_means_int(p=1, M=2, d=3, params=theta_123, parametrization="intercept", constraints=NULL),
               pick_phi0(p=1, M=2, d=3, params=theta_123_mu))
  expect_equal(get_regime_means_int(p=1, M=2, d=3, params=theta_123_mu, parametrization="mean", constraints=NULL),
               calc_mu(p=1, M=2, d=3, params=theta_123))
  expect_equal(get_regime_means_int(p=2, M=1, d=3, params=theta_213, parametrization="intercept", constraints=NULL),
               pick_phi0(p=2, M=1, d=3, params=theta_213_mu))
  expect_equal(get_regime_means_int(p=2, M=1, d=3, params=theta_213_mu, parametrization="mean", constraints=NULL),
               calc_mu(p=2, M=1, d=3, params=theta_213))

  expect_equal(get_regime_means_int(p=1, M=1, d=2, params=theta_112c, parametrization="intercept", constraints=C_112),
               pick_phi0(p=1, M=1, d=2, params=theta_112c_mu_exp))
  expect_equal(get_regime_means_int(p=1, M=1, d=2, params=theta_112c_mu, parametrization="mean", constraints=C_112),
               calc_mu(p=1, M=1, d=2, params=theta_112c, constraints=C_112))
  expect_equal(get_regime_means_int(p=2, M=2, d=2, params=theta_222c, parametrization="intercept", constraints=C_222),
               pick_phi0(p=2, M=2, d=2, params=theta_222c_mu_exp))
  expect_equal(get_regime_means_int(p=2, M=2, d=2, params=theta_222c_mu, parametrization="mean", constraints=C_222),
               calc_mu(p=2, M=2, d=2, params=theta_222c, constraints=C_222))
  expect_equal(get_regime_means_int(p=2, M=2, d=2, params=theta_222_c2, parametrization="intercept", constraints=C_222_2),
               pick_phi0(p=2, M=2, d=2, params=theta_222c_mu2_exp))
  expect_equal(get_regime_means_int(p=2, M=2, d=2, params=theta_222c_mu2, parametrization="mean", constraints=C_222_2),
               calc_mu(p=2, M=2, d=2, params=theta_222_c2, constraints=C_222_2))
  expect_equal(get_regime_means_int(p=1, M=2, d=3, params=theta_123c, parametrization="intercept", constraints=C_123),
               pick_phi0(p=1, M=2, d=3, params=theta_123c_mu_exp))
  expect_equal(get_regime_means_int(p=1, M=2, d=3, params=theta_123c_mu, parametrization="mean", constraints=C_123),
               calc_mu(p=1, M=2, d=3, params=theta_123c, constraints=C_123))
})
