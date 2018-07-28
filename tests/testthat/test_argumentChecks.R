context("argument checking functions")
library(gmvarkit)

# The tests are very brief !!!

## A(M)(p)_(p)(M)(d)

# p=1, M=1, d=2
phi10_112 <- c(1.03, 2.36)
A11_112 <- matrix(c(1.25, 0.06, 0.04, 1.34), nrow=2, byrow=FALSE)
Omega1_112 <- matrix(c(0.93, -0.15, -0.15, 5.20), nrow=2, byrow=FALSE)

theta_112 <- c(phi10_112, vec(A11_112), vech(Omega1_112))

# p=1, M=2, d=2
phi10_122 <- c(1.03, 2.36)
A11_122 <- matrix(c(0.9, 0.3, -0.3, 0.9), nrow=2, byrow=FALSE)
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
A11_123 <- matrix(c(0.1, 0.21, 0.31, 0.12, 0.2, 0.32, 0.13, 0.23, 0.1), nrow=3, byrow=FALSE)
Omega1_123 <- matrix(c(1, 0.22, 0.33, 0.22, 2, 0.44, 0.33, 0.44, 3), nrow=3, byrow=FALSE)

phi20_123 <- c(1.11, 2.22, 3.33)
A21_123 <- matrix(c(-0.1, -0.21, -0.31, -0.12, -0.2, -0.32, -0.13, -0.23, -2.1), nrow=3, byrow=FALSE)
Omega2_123 <- matrix(c(1.1, 0.222, 0.333, 0.222, 2.2, 0.444, 0.333, 0.444, 3.3), nrow=3, byrow=FALSE)

alpha1_123 <- 0.6
upsilon1_123 <- c(phi10_123, vec(A11_123), vech(Omega1_123))
upsilon2_123 <- c(phi20_123, vec(A21_123), vech(Omega2_123))
theta_123 <- c(upsilon1_123, upsilon2_123, alpha1_123)

# p=2, M=1, d=3
phi10_213 <- c(1.1, 2.2, 3.3)
A11_213 <- matrix(c(1, 0.21, 0.31, 0.12, 1.3, 0.32, 0.13, 0.23, 1), nrow=3, byrow=FALSE)
A12_213 <- matrix(c(-0.1, -0.21, -0.31, -0.12, -0.2, -0.32, -0.13, -0.23, -0.3), nrow=3, byrow=FALSE)
Omega1_213 <- matrix(c(1, 0.22, 0.33, 0.22, 2, 0.44, 0.33, 0.44, 3), nrow=3, byrow=FALSE)

upsilon1_213 <- c(phi10_213, vec(A11_213), vec(A12_213), vech(Omega1_213))
theta_213 <- upsilon1_213


# Constraining AR-parameters to be the same for all regimes
rbind_diags <- function(p, M, d) {
  I <- diag(p*d^2)
  Reduce(rbind, replicate(M, I, simplify=FALSE))
}


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


test_that("is_stationary works correctly", {
  expect_false(is_stationary(p=1, M=1, d=2, params=theta_112))
  expect_true(is_stationary(p=2, M=2, d=2, params=theta_222))
  expect_false(is_stationary(p=3, M=3, d=2, params=theta_332))
  expect_false(is_stationary(p=1, M=2, d=3, params=theta_123))
  expect_false(is_stationary(p=2, M=1, d=3, params=theta_213))
})

test_that("in_paramspace works correctly", {
  expect_false(in_paramspace(p=1, M=1, d=2, params=theta_112))
  expect_true(in_paramspace(p=2, M=2, d=2, params=theta_222))
  expect_false(in_paramspace(p=3, M=3, d=2, params=theta_332))
  expect_false(in_paramspace(p=1, M=2, d=3, params=theta_123))
  expect_false(in_paramspace(p=2, M=1, d=3, params=theta_213))

  expect_false(in_paramspace(p=1, M=1, d=2, params=theta_112c, constraints=C_112))
  expect_true(in_paramspace(p=2, M=2, d=2, params=theta_222c, constraints=C_222))
  expect_true(in_paramspace(p=2, M=2, d=2, params=theta_222_c2, constraints=C_222_2))
})

test_that("check_parameters works correctly", {
  expect_error(check_parameters(p=1, M=1, d=2, params=theta_112))
  expect_error(check_parameters(p=2, M=3, d=2, params=theta_222))
  expect_error(check_parameters(p=3, M=3, d=2, params=theta_332))
  expect_error(check_parameters(p=1, M=2, d=3, params=theta_123))
  expect_error(check_parameters(p=2, M=1, d=3, params=theta_213))

  expect_error(check_parameters(p=1, M=1, d=2, params=theta_112c, constraints=C_112))
  expect_error(check_parameters(p=2, M=2, d=2, params=theta_222c, constraints=t(C_222)))
  expect_error(check_parameters(p=2, M=2, d=2, params=theta_222_c2, constraints=as.vector(C_222_2)))
})


test_that("check_constraints works correctly", {
  expect_error(check_constraints(p=1, M=1, d=2, constraints=cbind(C_112, C_112)))
  expect_error(check_constraints(p=1, M=1, d=2, constraints=rbind(C_112, C_112)))
  expect_error(check_constraints(p=1, M=2, d=2, constraints=cbind(C_122[,1:3], rep(0, 8))))
  expect_error(check_constraints(p=1, M=2, d=2, constraints=as.data.frame(C_122)))
  expect_error(check_constraints(p=2, M=2, d=2, constraints=C_222[-1,]))
  expect_error(check_constraints(p=2, M=2, d=2, constraints=C_222_2[-16,]))
  expect_error(check_constraints(p=3, M=3, d=2, constraints=t(C_332)))
  expect_error(check_constraints(p=1, M=2, d=3, constraints=cbind(C_123[,1:8], C_123[,8])))
  expect_error(check_constraints(p=2, M=1, d=3, constraints=as.vector(C_213)))
})



data_na <- eurusd; data_na[252, 2] <- NA
eurusd_mat <- t(t(eurusd))

test_that("check_data works correctly", {
  expect_equal(check_data(data=eurusd, p=10), eurusd)
  expect_equal(check_data(data=ts(eurusd_mat, start=c(1989, 1), deltat=1/12), p=1), eurusd)
  expect_equal(check_data(data=as.data.frame(eurusd_mat), p=2), eurusd_mat)

  expect_error(check_data(data=eurusd[,1, drop=FALSE], p=4))
  expect_error(check_data(data=as.data.frame(eurusd[,2]), p=4))
  expect_error(check_data(data=ts(eurusd[,2], start=c(1989, 1), deltat=1/12), p=3))
  expect_error(check_data(data=data_na, p=4))
  expect_error(check_data(data=eurusd, p=nrow(eurusd)))
})

test_that("check_pMd works correctly", {
  expect_error(check_pMd(p=1, M=1, d=1))
  expect_error(check_pMd(p=1, M=1.2, d=2))
  expect_error(check_pMd(p=0, M=1, d=2))
  expect_error(check_pMd(p=1, M=-1, d=2))
  expect_error(check_pMd(p=1.1, M=1, d=2))
  expect_error(check_pMd(p=1, M=1, d=2.2))
})

test_that("all_pos_ints works correctly", {
  expect_true(all_pos_ints(c(1, 2, 3)))
  expect_true(all_pos_ints(1))
  expect_true(all_pos_ints(list(1, 3, 100)))

  expect_false(all_pos_ints(c(1, 2, 0)))
  expect_false(all_pos_ints(-1))
  expect_false(all_pos_ints(0.1))
  expect_false(all_pos_ints(1.1))
  expect_false(all_pos_ints(list(1, 2, 3, 0.1)))





})
