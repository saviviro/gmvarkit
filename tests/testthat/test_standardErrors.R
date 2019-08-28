context("calculating approximate standard errors")
library(gmvarkit)

data <- cbind(10*eurusd[,1], 100*eurusd[,2])

## NOTE: For some reason these tests fail at win-builder with tolerance smaller than 1e-2 (while they pass locally).
# Apparently, when testing diagonal of the observed information matrix, the tests pass, but when testing
# diagonal of its inverse (or the standard errors), the tests fail (with average difference about 1e-3).
# Possibly a numerical error caused by (the default) imprecission of float-point presentation?

## A(M)(p)_(p)(M)(d)

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

# Constraint AR-parameters to be the same for all regimes
rbind_diags <- function(p, M, d) {
  I <- diag(p*d^2)
  Reduce(rbind, replicate(M, I, simplify=FALSE))
}

C_222 <- rbind_diags(p=2, M=2, d=2)
theta_222c <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vech(Omega1_222), vech(Omega2_222), alpha1_222)

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

test_that("standard_errors works correctly", {
   expect_equal(standard_errors(data, p=2, M=2, params=theta_222c, conditional=TRUE, parametrization="intercept",
                                constraints=C_222, minval=-99999),
                c(0.57323112, 1.42708542, 0.68407209, 1.79850113, 0.06503308, 0.09518451, 0.03464635, 0.06118081, 0.06339426,
                  0.09592758, 0.03487118, 0.06453526, 0.15615508, 0.27195919, 0.72882509, 0.72951956, 0.73219670, 1.18963949,
                  0.28865116), tolerance=1e-2)
   expect_equal(standard_errors(data, p=2, M=2, params=theta_222_c2, conditional=TRUE, parametrization="intercept",
                                constraints=C_222_2, minval=-99999),
                c(NA, 0.7887294, NA, 0.92849143, 0.05692548, 0.05744037, 0.05888322, 0.05761173, NA, NA, 0.85622084,
                  0.68689890, 0.72052650, 1.16661586, 0.40045562), tolerance=1e-2)
})
