context("S3 Methods")
library(gmvarkit)

# NOTE that some elements of these tests use random elements obtained from simulation algorithms

## A(M)(p)_(p)(M)(d)

# p=1, M=1, d=2, parametrization="mean"
phi10_112 <- c(0.75, 0.8)
A11_112 <- matrix(c(0.29, 0.02, -0.14, 0.9), nrow=2, byrow=FALSE)
Omega1_112 <- matrix(c(0.60, 0.01, 0.01, 0.07), nrow=2, byrow=FALSE)
theta_112 <- c(phi10_112, vec(A11_112), vech(Omega1_112))
mod_112 <- GSMVAR(gdpdef, p=1, M=1, d=2, params=theta_112, conditional=TRUE, parametrization="mean")

mod_112t <- GSMVAR(gdpdef, p=1, M=1, d=2, params=c(theta_112, 3), model="StMVAR", parametrization="mean")

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
mod_222 <- GSMVAR(gdpdef, p=2, M=2, d=2, params=theta_222, conditional=TRUE, parametrization="intercept", constraints=NULL)

mod_222gs <- GSMVAR(gdpdef, p=2, M=c(1, 1), d=2, params=c(theta_222, 20), model="G-StMVAR",
                    conditional=TRUE, parametrization="intercept")


WL_222 <- diag_Omegas(Omega1_222, Omega2_222)
W_222 <- matrix(WL_222[1:(2^2)], nrow=2, byrow=FALSE)
lambdas_222 <- WL_222[(2^2 + 1):length(WL_222)]
theta_222s <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A21_222),
                vec(A22_222), vec(W_222), lambdas_222, alpha1_222) # SGMVAR
mod_222s <- GSMVAR(gdpdef, p=2, M=2, d=2, params=theta_222s, conditional=TRUE, parametrization="intercept", constraints=NULL,
                  structural_pars=list(W=W_222))

# p=2, M=2, d=2, AR paramars same, non-diagonals zero, intercept
theta_222c <- c(0.33782, 0.183512, 0.472168, 0.095311, 0.201199, 0.600596, 0.237819,
                0.23529, 1.077816, -0.016343, 0.112771, 0.22199, 0.005582, 0.028126, 0.492844)
mat0 <- matrix(c(1, rep(0, 10), 1, rep(0, 8), 1, rep(0, 10), 1), nrow=2*2^2, byrow=FALSE)
C_222c <- rbind(mat0, mat0)
mod_222c <- GSMVAR(gdpdef, p=2, M=2, d=2, params=theta_222c, conditional=TRUE, parametrization="intercept", constraints=C_222c)


# p=1, M=2, d=3, no constraints, rand_ind and simulated data
theta_123 <- c(-9.44567, -0.56054, 10.32549, 0.0965, 0.63617, 0.35771, 0.63339,
               0.2519, -0.32399, 0.56932, -0.47935, 0.32332, 1.04371, 0.08397,
               0.71741, 0.46644, 0.23572, 1.14101, -8.16384, 0.7148, 1.86377,
               0.2646, -0.07309, -0.78756, -0.86484, -0.16795, -0.26713,
               -0.0035, 0.6088, -0.19626, 0.36186, -0.16349, 0.06036, 0.58441,
               1.10884, 2.64874, 0.54711)
mod_123 <- GSMVAR(p=1, M=2, d=3, params=theta_123, conditional=FALSE, parametrization="mean", constraints=NULL)
sim_123 <- simulate.gsmvar(mod_123, nsim=300, seed=2)
data_123 <- sim_123$sample
mod_123 <- GSMVAR(data_123, p=1, M=2, d=3, params=theta_123, conditional=FALSE, parametrization="mean", constraints=NULL)

mod_123t <- GSMVAR(data_123, p=1, M=2, d=3, params=c(theta_123, 20, 30), model="StMVAR",
                   conditional=FALSE, parametrization="mean")

set.seed(1); pred112t <- predict.gsmvar(mod_112t, n_ahead=2, nsim=10, pi=c(0.80), plot_res=FALSE, pred_type="mean")
set.seed(1); pred222 <- predict.gsmvar(mod_222, n_ahead=2, nsim=10, pi=c(0.95, 0.80), plot_res=FALSE, pred_type="mean")
set.seed(1); pred222gs <- predict.gsmvar(mod_222gs, n_ahead=2, nsim=10, pi=c(0.95, 0.70), plot_res=FALSE, pred_type="mean")
set.seed(2); pred222s <- predict.gsmvar(mod_222s, n_ahead=2, nsim=10, pi=c(0.95, 0.80), plot_res=FALSE, pred_type="mean")
set.seed(3); pred123 <- predict.gsmvar(mod_123, n_ahead=1, nsim=10, pi=0.99, pi_type="upper", pred_type="median", plot_res=FALSE)
set.seed(3); pred123t <- predict.gsmvar(mod_123t, n_ahead=2, nsim=5, pi=0.99, pi_type="lower", pred_type="median", plot_res=FALSE)
tmp222 <- unname(pred222$pred[2,])

# p=2, M=2, d=2, parametrization="mean", constraints=C_mat, same_means=list(1:2)
C_mat <- rbind(diag(2*2^2), diag(2*2^2))
params_222cm <- c(0.811034, 0.578587, 0.212084, 0.020444, -0.193005, 0.624671,
                  0.235827, 0.013962, 0.053267, 0.262703, 1.06105, -0.013519,
                  0.114109, 0.229542, 0.003092, 0.027266, 0.424341)
mod_222cm <- GSMVAR(gdpdef, p=2, M=2, params=params_222cm, parametrization="mean", constraints=C_mat, same_means=list(1:2))
set.seed(1); pred222cm <- predict.gsmvar(mod_222cm, n_ahead=2, nsimu=1, pi=0.9, pi_type="two-sided", pred_type="mean", plot_res=FALSE)

test_that("predict works correctly", {
   expect_equal(predict.gsmvar(mod_112, n_ahead=1, pred_type="cond_mean", plot_res=FALSE)$pred, c(0.7231782, 0.4431300), tolerance=1e-5)
   expect_equal(predict.gsmvar(mod_222c, n_ahead=1, pred_type="cond_mean", plot_res=FALSE)$pred, c(0.7250053, 0.4209626), tolerance=1e-5)

   expect_equal(unname(pred112t$pred[2, ]), c(0.8970004, 0.4217662), tolerance=1e-3)
   expect_equal(pred112t$pred_ints[, 2, 2], c(0.5729150, 0.5105862), tolerance=1e-3)
   expect_equal(pred112t$pred_ints[, 1, 1], c(0.4116595, 0.3756889), tolerance=1e-3)
   expect_equal(pred112t$mix_pred_ints[, 1, 1], c(1, 1), tolerance=1e-3)

   expect_equal(tmp222, c(0.6709308, 0.4618839), tolerance=1e-5)
   expect_equal(pred222$pred_ints[, 1, 1], c(0.07127095, -0.32711717), tolerance=1e-3)
   expect_equal(pred222$pred_ints[, 3, 2], c(0.5783407, 0.5995812), tolerance=1e-3)
   expect_equal(pred222$mix_pred_ints[, 1, 1], c(0.9352294, 0.8441136), tolerance=1e-3)

   expect_equal(pred222s$pred_ints[, 2, 1], c(0.5230737, -0.3738420), tolerance=1e-3)
   expect_equal(pred222s$pred_ints[, 1, 2], c(0.09878082, 0.18761453), tolerance=1e-3)
   expect_equal(pred222s$mix_pred_ints[, 2, 2], c(0.06477058, 0.05950152), tolerance=1e-3)

   expect_equal(unname(pred222gs$pred[2, ]), c(0.7176830, 0.4060806), tolerance=1e-3)
   expect_equal(pred222gs$pred_ints[, 4, 2], c(0.7705175, 0.5784606), tolerance=1e-3)
   expect_equal(pred222gs$pred_ints[, 2, 2], c(0.1779556, 0.3322382), tolerance=1e-3)
   expect_equal(pred222gs$mix_pred_ints[, 2, 1], c(0.9260951, 0.7668288), tolerance=1e-3)

   expect_equal(unname(pred123$pred[1,]), c(-8.4121641, -0.3787007, 2.3372331), tolerance=1e-5)
   expect_equal(pred123$pred_ints[ , 1, ], c(-7.987103, 1.036073, 4.507460), tolerance=1e-5)
   expect_equal(unname(pred123$mix_pred[1 ,]), c(7.987841e-21, 1.000000e+00), tolerance=1e-5)
   expect_equal(unname(pred123$mix_pred_ints[1 , 1, ]), c(7.987841e-21, 1.000000e+00), tolerance=1e-5)

   expect_equal(unname(pred123t$pred[2,]), c(-7.782932, 1.377770, 2.024874), tolerance=1e-5)
   expect_equal(pred123t$pred_ints[2 , 1, ], c(-7.9877016, 0.7730701, 0.2280571), tolerance=1e-5)
   expect_equal(unname(pred123t$mix_pred[2,]), c(3.535557e-07, 9.999996e-01), tolerance=1e-5)
   expect_equal(unname(pred123t$mix_pred_ints[2 , 1, ]), c(1.487353e-07, 9.999948e-01), tolerance=1e-5)

   expect_equal(unname(pred222cm$pred[2,]), c(0.7434035, 0.4107316), tolerance=1e-5)
   expect_equal(unname(pred222cm$pred_ints[2, 2, ]), c(1.6821331, 0.7858234), tolerance=1e-5)
})

# p=2, M=2, d=2, parametrization="mean", constraints=C_mat, same_means=list(1:2)
C_mat <- rbind(diag(2*2^2), diag(2*2^2))
params_222cm <- c(0.811034, 0.578587, 0.212084, 0.020444, -0.193005, 0.624671,
                  0.235827, 0.013962, 0.053267, 0.262703, 1.06105, -0.013519,
                  0.114109, 0.229542, 0.003092, 0.027266, 0.424341)
mod_222cm <- GSMVAR(gdpdef, p=2, M=2, params=params_222cm, parametrization="mean", constraints=C_mat, same_means=list(1:2))

test_that("summary method works correctly", {
  sum222cm <- summary(mod_222cm)
  expect_equal(sum222cm$abs_boldA_eigens[3,], c(0.3963773, 0.3963773), tolerance=1e-5)
  expect_equal(sum222cm$omega_eigens[1:4], c(1.06124296, 0.11391604, 0.22958925, 0.02721875), tolerance=1e-5)
})
