context("S3 Methods")
library(gmvarkit)

# NOTE that some elements of these tests use random elements obtained from simulation algorithms

## A(M)(p)_(p)(M)(d)

# p=1, M=1, d=2, parametrization="mean"
phi10_112 <- c(0.75, 0.8)
A11_112 <- matrix(c(0.29, 0.02, -0.14, 0.9), nrow=2, byrow=FALSE)
Omega1_112 <- matrix(c(0.60, 0.01, 0.01, 0.07), nrow=2, byrow=FALSE)
theta_112 <- c(phi10_112, vec(A11_112), vech(Omega1_112))
mod_112 <- GMVAR(gdpdef, p=1, M=1, d=2, params=theta_112, conditional=TRUE, parametrization="mean", constraints=NULL)


# p=2, M=2, d=2, no constraints
phi10_222 <- c(0.36, 0.12)
A11_222 <- matrix(c(0.22, 0.06, -0.15, 0.39), nrow=2, byrow=FALSE)
A12_222 <- matrix(c(0.41, -0.01, 0.08, 0.3), nrow=2, byrow=FALSE)
Omega1_222 <- matrix(c(0.21, 0.01, 0.01, 0.03), nrow=2, byrow=FALSE)

phi20_222 <- c(0.48, 0.07)
A21_222 <- matrix(c(0.22, 0.02, -0.12, 0.72), nrow=2, byrow=FALSE)
A22_222 <- matrix(c(0.09, 0.03, 0.04, 0.19), nrow=2, byrow=FALSE)
Omega2_222 <- matrix(c(1.10, 0.01, 0.01, 0.11), nrow=2, byrow=FALSE)

alpha1_222 <- 0.37
upsilon1_222 <- c(phi10_222, vec(A11_222), vec(A12_222), vech(Omega1_222))
upsilon2_222 <- c(phi20_222, vec(A21_222), vec(A22_222), vech(Omega2_222))
theta_222 <- c(upsilon1_222, upsilon2_222, alpha1_222)
mod_222 <- GMVAR(gdpdef, p=2, M=2, d=2, params=theta_222, conditional=TRUE, parametrization="intercept", constraints=NULL)

WL_222 <- diag_Omegas(Omega1_222, Omega2_222)
W_222 <- matrix(WL_222[1:(2^2)], nrow=2, byrow=FALSE)
lambdas_222 <- WL_222[(2^2 + 1):length(WL_222)]
theta_222s <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A21_222),
                vec(A22_222), vec(W_222), lambdas_222, alpha1_222) # SGMVAR
mod_222s <- GMVAR(gdpdef, p=2, M=2, d=2, params=theta_222s, conditional=TRUE, parametrization="intercept", constraints=NULL,
                  structural_pars=list(W=W_222))

# p=2, M=2, d=2, AR paramars same, non-diagonals zero, intercept
theta_222c <- c(0.33782, 0.183512, 0.472168, 0.095311, 0.201199, 0.600596, 0.237819,
                0.23529, 1.077816, -0.016343, 0.112771, 0.22199, 0.005582, 0.028126, 0.492844)
mat0 <- matrix(c(1, rep(0, 10), 1, rep(0, 8), 1, rep(0, 10), 1), nrow=2*2^2, byrow=FALSE)
C_222c <- rbind(mat0, mat0)
mod_222c <- GMVAR(gdpdef, p=2, M=2, d=2, params=theta_222c, conditional=TRUE, parametrization="intercept", constraints=C_222c)


# p=1, M=2, d=3, no constraints, rand_ind and simulated data
set.seed(13)
theta_123 <- random_ind2(p=1, M=2, d=3, mu_scale=c(-10, 0, 5), mu_scale2=1:3, omega_scale=1:3, ar_scale=1)
mod_123 <- GMVAR(p=1, M=2, d=3, params=theta_123, conditional=FALSE, parametrization="mean", constraints=NULL)
sim_123 <- simulateGMVAR(mod_123, nsimu=300)
data_123 <- sim_123$sample
mod_123 <- GMVAR(data_123, p=1, M=2, d=3, params=theta_123, conditional=FALSE, parametrization="mean", constraints=NULL)

set.seed(1)
pred222 <- predict.gmvar(mod_222, n_ahead=2, n_simu=10, pi=c(0.95, 0.80), plot_res=FALSE, pred_type="mean")
pred222s <- predict.gmvar(mod_222s, n_ahead=2, n_simu=10, pi=c(0.95, 0.80), plot_res=FALSE, pred_type="mean")
pred123 <- predict.gmvar(mod_123, n_ahead=1, n_simu=10, pi=0.99, pi_type="upper", pred_type="median", plot_res=FALSE)
tmp222 <- unname(pred222$pred[2,])

# p=2, M=2, d=2, parametrization="mean", constraints=C_mat, same_means=list(1:2)
C_mat <- rbind(diag(2*2^2), diag(2*2^2))
params_222cm <- c(0.811034, 0.578587, 0.212084, 0.020444, -0.193005, 0.624671,
                  0.235827, 0.013962, 0.053267, 0.262703, 1.06105, -0.013519,
                  0.114109, 0.229542, 0.003092, 0.027266, 0.424341)
mod_222cm <- GMVAR(gdpdef, p=2, M=2, params=params_222cm, parametrization="mean", constraints=C_mat, same_means=list(1:2))
set.seed(1); pred222cm <- predict.gmvar(mod_222cm, n_ahead=2, nsimu=1, pi=0.9, pi_type="two-sided", pred_type="mean", plot_res=FALSE)

test_that("predict works correctly", {
   expect_equal(predict.gmvar(mod_112, n_ahead=1, pred_type="cond_mean", plot_res=FALSE)$pred, c(0.7231782, 0.4431300), tolerance=1e-5)
   expect_equal(predict.gmvar(mod_222c, n_ahead=1, pred_type="cond_mean", plot_res=FALSE)$pred, c(0.7250053, 0.4209626), tolerance=1e-5)

   expect_equal(tmp222, c(0.6709308, 0.4618839), tolerance=1e-5)
   expect_equal(pred222$pred_ints[, 1, 1], c(0.07127095, -0.32711717), tolerance=1e-3)
   expect_equal(pred222$pred_ints[, 3, 2], c(0.5783407, 0.5995812), tolerance=1e-3)
   expect_equal(pred222$mix_pred_ints[, 1, 1], c(0.9352294, 0.8441136), tolerance=1e-3)

   expect_equal(unname(pred123$pred[1,]), c(-9.309146, 1.768539, 1.376839), tolerance=1e-5)
   expect_equal(pred123$pred_ints[ , 1, ], c(-8.154824, 3.236734, 3.383006), tolerance=1e-5)
   expect_equal(unname(pred123$mix_pred[1 ,]), c(1.947047e-10, 1.000000e+00), tolerance=1e-5)
   expect_equal(unname(pred123$mix_pred_ints[1 , 1, ]), c(1.947047e-10, 1.000000e+00), tolerance=1e-5)

   expect_equal(unname(pred222cm$pred[2,]), c(0.7434035, 0.4107316), tolerance=1e-5)
   expect_equal(unname(pred222cm$pred_ints[2, 2, ]), c(1.6821331, 0.7858234), tolerance=1e-5)
})

# p=2, M=2, d=2, parametrization="mean", constraints=C_mat, same_means=list(1:2)
C_mat <- rbind(diag(2*2^2), diag(2*2^2))
params_222cm <- c(0.811034, 0.578587, 0.212084, 0.020444, -0.193005, 0.624671,
                  0.235827, 0.013962, 0.053267, 0.262703, 1.06105, -0.013519,
                  0.114109, 0.229542, 0.003092, 0.027266, 0.424341)
mod_222cm <- GMVAR(gdpdef, p=2, M=2, params=params_222cm, parametrization="mean", constraints=C_mat, same_means=list(1:2))

test_that("summary method works correctly", {
  sum222cm <- summary(mod_222cm)
  expect_equal(sum222cm$abs_boldA_eigens[3,], c(0.3963773, 0.3963773), tolerance=1e-5)
  expect_equal(sum222cm$omega_eigens[1:4], c(1.06124296, 0.11391604, 0.22958925, 0.02721875), tolerance=1e-5)
})
