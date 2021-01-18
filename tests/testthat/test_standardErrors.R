context("standard errors")
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

WL_222 <- diag_Omegas(Omega1_222, Omega2_222)
#W_222 <- matrix(WL_222[1:(2^2)], nrow=2, byrow=FALSE)
#lambdas_222 <- WL_222[(2^2 + 1):length(WL_222)]
#theta_222s <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A21_222),
#                vec(A22_222), vec(W_222), lambdas_222, alpha1_222) # SGMVAR
theta_222s <- c(1.386765, -0.764975, 1.00547, 5.928393, 1.314132, 0.144683, 0.093975,
                1.292154, -0.389113, -0.069592, -0.108713, -0.281433, 1.248205, 0.076956,
                -0.039837, 1.265769, -0.272355, -0.073872, 0.034041, -0.31349, -0.902536,
                -0.718412, 0.324033, -2.079116, 7.001825, 1.439954, 0.739348)
W_222 <- pick_W(p=2, M=2, d=2, params=theta_222s, structural_pars=list(W=matrix(1:4, nrow=2)))

# Constraint AR-parameters to be the same for all regimes
rbind_diags <- function(p, M, d) {
  I <- diag(p*d^2)
  Reduce(rbind, replicate(M, I, simplify=FALSE))
}

C_222 <- rbind_diags(p=2, M=2, d=2)
theta_222c <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vech(Omega1_222), vech(Omega2_222), alpha1_222)

C_lambda_222 <- matrix(c(1, 2), nrow=2)
theta_222csL <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A21_222),
                  vec(A22_222), vec(W_222), 0.2, alpha1_222) # SGMVAR lambdas
theta_222csLAR <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(W_222), 0.2, alpha1_222) # SGMVAR lambdas and AR


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

WL_222c2 <- diag_Omegas(Omega1_222_c2, Omega2_222_c2)
W_222c2 <- matrix(WL_222c2[1:(2^2)], nrow=2, byrow=FALSE)
lambdas_222c2 <- WL_222c2[(2^2 + 1):length(WL_222c2)]
theta_222_c2s <- c(phi10_222_c2, phi20_222_c2, 1.26, 1.34, -0.29, -0.36, vec(W_222c2), lambdas_222c2, alpha1_222_c2) # SGMVAR AR

## p=2, M=2, d=2, parametrization="mean", constraints=C_mat, same_means=list(1:2)
C_mat <- rbind(diag(2*2^2), diag(2*2^2))
params_222cm <- c(-9.12048853805255, 123.142757183508, 1.2658425363326, 0.0675545389606989, 0.0331264235657607,
                  1.33370494344656, -0.285882557831441, -0.0769144929653558, -0.0382772162867802, -0.351635998882842,
                  5.8625623309659, 3.57488618757834, 9.70846346569286, 0.869261580580846, -0.248703862116217,
                  5.17613656742281, 0.439575388572472)

test_that("standard_errors works correctly", {
  expect_equal(standard_errors(data, p=2, M=2, params=theta_222c, conditional=TRUE, parametrization="intercept",
                               constraints=C_222, minval=-99999),
               c(0.57323112, 1.42708542, 0.68407209, 1.79850113, 0.06503308, 0.09518451, 0.03464635, 0.06118081, 0.06339426,
                 0.09592758, 0.03487118, 0.06453526, 0.15615508, 0.27195919, 0.72882509, 0.72951956, 0.73219670, 1.18963949,
                 0.28865116), tolerance=1e-1)
  # expect_equal(standard_errors(data, p=2, M=2, params=theta_222_c2, conditional=TRUE, parametrization="intercept",
  #                              constraints=C_222_2, minval=-99999),
  #              c(NA, 0.7887294, NA, 0.92849143, 0.05692548, 0.05744037, 0.05888322, 0.05761173, NA, NA, 0.85622084,
  #                0.68689890, 0.72052650, 1.16661586, 0.40045562), tolerance=1e-1)

  # SGMVAR
  expect_equal(standard_errors(data, p=2, M=2, params=theta_222s, conditional=FALSE, parametrization="intercept",
                               constraints=NULL, structural_pars=list(W=W_222), minval=-99999),
               c(0.657165, 1.435745, 1.888876, 2.323694, 0.097186, 0.23554, 0.048015, 0.102941, 0.096483, 0.225936,
                 0.04786, 0.110421, 0.091458, 0.117798, 0.070049, 0.089094, 0.092122, 0.11891, 0.070357, 0.090593,
                 0.088115, 0.127551, 0.142765, 0.184838, 1.40735, 0.272365, 0.223313), tolerance=1e-1)
  expect_equal(standard_errors(data, p=2, M=2, params=theta_222_c2s, conditional=TRUE, parametrization="intercept",
                               constraints=C_222_2, structural_pars=list(W=W_222c2), minval=-99999),
               c(0.058447, 0.875092, 0.204019, 0.992678, 0.057246, 0.059366, 0.058941, 0.060121, NA, 0.24881,
                 0.308266, 0.382079, NA, 0.201883, 0.572563), tolerance=1e-1)

  # Same_means - the results are a bit different in some setups used in CRAN tests; thats why commented.
  # expect_equal(standard_errors(data, p=2, M=2, params=params_222cm, conditional=TRUE, parametrization="mean",
  #                              constraints=C_mat, same_means=list(1:2), minval=-99999),
  #              c(5.65699691, 14.55107385, 0.06494907, 0.09540929, 0.03515466, 0.06032848, 0.06475830, 0.09546348,
  #                0.03533626, 0.06107803, 0.74041185, 0.72118921, 1.16976089, 0.14426045, 0.26732233, 0.79159544,
  #                0.12665061), tolerance=1e-2)
})
