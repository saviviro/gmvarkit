context("S3 Methods")
library(gmvarkit)

# NOTE that some elements of these tests use random elements obtained from simulation algorithms

data <- cbind(10*eurusd[,1], 100*eurusd[,2])

## A(M)(p)_(p)(M)(d)

# p=1, M=1, d=2, parametrization="mean"
phi10_112 <- c(1.07, 127.71)
A11_112 <- matrix(c(0.99, 0.00, -0.01, 0.99), nrow=2)
Omega1_112 <- matrix(c(4.05, 2.22, 8.87, 2.22), nrow=2)
theta_112 <- c(phi10_112, vec(A11_112), vech(Omega1_112))
mod_112 <- GMVAR(data, p=1, M=1, d=2, params=theta_112, conditional=TRUE, parametrization="mean", constraints=NULL)

# p=2, M=2, d=2, no constraints, GMVAR-paper
phi10_222 <- c(1.03, 2.36)
A11_222 <- matrix(c(1.25, 0.06, 0.04, 1.34), nrow=2, byrow=FALSE)
A12_222 <- matrix(c(-0.29, -0.08, -0.05, -0.36), nrow=2, byrow=FALSE)
Omega1_222 <- matrix(c(0.93, -0.15, -0.15, 5.20), nrow=2, byrow=FALSE)

phi20_222 <- c(1.79, 3.00)
A21_222 <- A11_222; A22_222 <- A12_222
Omega2_222 <- matrix(c(5.88, 3.56, 3.56, 9.80), nrow=2, byrow=FALSE)

alpha1_222 <- 0.37

upsilon1_222 <- c(phi10_222, vec(A11_222), vec(A12_222), vech(Omega1_222))
upsilon2_222 <- c(phi20_222, vec(A21_222), vec(A22_222), vech(Omega2_222))
theta_222 <- c(upsilon1_222, upsilon2_222, alpha1_222)
mod_222 <- GMVAR(data, p=2, M=2, d=2, params=theta_222, conditional=TRUE, parametrization="intercept", constraints=NULL)

# p=2, M=2, d=2, AR paramars same, non-diagonals zero, intercept
theta_222c <- c(0.3552775, 3.1929675, -0.1143198, 2.8294743, 1.2633425, 1.3375150, -0.2919742, -0.3624010,
                5.5971764, 3.4559442, 9.6221422, 0.9820759, -0.3267521, 5.2358855, 0.6501600)
mat0 <- matrix(c(1, rep(0, 10), 1, rep(0, 8), 1, rep(0, 10), 1), nrow=2*2^2, byrow=FALSE) # Laske paperilla monimutkaisemmat rajoitteet
C_222c <- rbind(mat0, mat0)
mod_222c <- GMVAR(data, p=2, M=2, d=2, params=theta_222c, conditional=TRUE, parametrization="intercept", constraints=C_222c)

# p=1, M=2, d=3, no constraints, rand_ind and simulated data
set.seed(13)
theta_123 <- random_ind2(p=1, M=2, d=3, mu_scale=c(-10, 0, 5), mu_scale2=1:3, omega_scale=1:3, ar_scale=1)
mod_123 <- GMVAR(p=1, M=2, d=3, params=theta_123, conditional=FALSE, parametrization="mean", constraints=NULL)
sim_123 <- simulateGMVAR(mod_123, nsimu=300)
data_123 <- sim_123$sample
# theta_123 <- c(-9.432281193, -0.505061517, 10.439237825, -0.007153524, 0.611600845, 0.316006743, 0.618395547, 0.242141051, -0.397496657,  0.607582014,
#                -0.588808209, 0.249527481, 1.000119135, 0.066231606, 0.772660170, 0.427032675, 0.173906741, 1.161863865, -8.261608811, 0.836296696,
#                2.010706334, 0.311276277, -0.049947762, -0.650465167, -0.857957924, -0.158176495, -0.218184205, -0.016438743, 0.606639158, -0.207387484,
#                0.393643727, -0.205999160, 0.008228229, 0.681252658, 1.261425042, 2.916105585, 0.685775441)
mod_123 <- GMVAR(data_123, p=1, M=2, d=3, params=theta_123, conditional=FALSE, parametrization="mean", constraints=NULL)

set.seed(1)
pred222 <- predict.gmvar(mod_222, n_ahead=2, n_simu=10, pi=c(0.95, 0.80), plot_res=FALSE, pred_type="mean")
pred123 <- predict.gmvar(mod_123, n_ahead=1, n_simu=10, pi=0.99, pi_type="upper", pred_type="median", plot_res=FALSE)
tmp222 <- unname(pred222$pred[2,])

# p=1 M=2, d=2, parametrization="mean", constraints=C_mat, same_means=list(1:2)
C_mat <- rbind(diag(2*2^2), diag(2*2^2))
params_222cm <- c(-9.12048853805255, 123.142757183508, 1.2658425363326, 0.0675545389606989, 0.0331264235657607,
                  1.33370494344656, -0.285882557831441, -0.0769144929653558, -0.0382772162867802, -0.351635998882842,
                  5.8625623309659, 3.57488618757834, 9.70846346569286, 0.869261580580846, -0.248703862116217,
                  5.17613656742281, 0.439575388572472)
mod_222cm <- GMVAR(data, p=2, M=2, params=params_222cm, parametrization="mean", constraints=C_mat, same_means=list(1:2))
set.seed(1); pred222cm <- predict.gmvar(mod_222cm, n_ahead=2, nsimu=1, pi=0.9, pi_type="two-sided", pred_type="mean", plot_res=FALSE)

test_that("predict works correctly", {
   expect_equal(predict.gmvar(mod_112, n_ahead=1, pred_type="cond_mean", plot_res=FALSE)$pred, c(2.625171, 145.951650), tolerance=1e-5)
   expect_equal(predict.gmvar(mod_222c, n_ahead=1, pred_type="cond_mean", plot_res=FALSE)$pred, c(2.643673, 144.594839), tolerance=1e-5)

   expect_equal(tmp222, c(2.202927, 144.450148), tolerance=1e-5)
   expect_equal(pred222$pred_ints[, 1, 1], c(-0.3797107, -4.2922563), tolerance=1e-3)
   expect_equal(pred222$pred_ints[, 3, 2], c(147.5922, 150.4961), tolerance=1e-3)
   expect_equal(pred222$mix_pred_ints[, 1, 1], c(0.005081765, 0.002005220), tolerance=1e-3)

   expect_equal(unname(pred123$pred[1,]), c(-9.203200, 2.321089, 1.795961), tolerance=1e-5)
   expect_equal(pred123$pred_ints[ , 1, ], c(-8.012776, 3.149987, 3.413862), tolerance=1e-5)
   expect_equal(unname(pred123$mix_pred[1 ,]), c(1.947047e-10, 1.000000e+00), tolerance=1e-5)
   expect_equal(unname(pred123$mix_pred_ints[1 , 1, ]), c(1.947047e-10, 1.000000e+00), tolerance=1e-5)

   expect_equal(unname(pred222cm$pred[2,]), c(1.303948, 143.158982), tolerance=1e-5)
   expect_equal(unname(pred222cm$pred_ints[2, 2, ]), c(7.594892, 151.985248), tolerance=1e-5)
})

# p=2, M=2, d=2, parametrization="mean", constraints=C_mat, same_means=list(1:2)
C_mat <- rbind(diag(2*2^2), diag(2*2^2))
params_222cm <- c(-9.12048853805255, 123.142757183508, 1.2658425363326, 0.0675545389606989, 0.0331264235657607,
                  1.33370494344656, -0.285882557831441, -0.0769144929653558, -0.0382772162867802, -0.351635998882842,
                  5.8625623309659, 3.57488618757834, 9.70846346569286, 0.869261580580846, -0.248703862116217,
                  5.17613656742281, 0.439575388572472)
mod_222cm <- GMVAR(data, p=2, M=2, params=params_222cm, parametrization="mean", constraints=C_mat, same_means=list(1:2))

test_that("summary method works correctly", {
  sum222cm <- summary(mod_222cm)
  expect_equal(sum222cm$abs_boldA_eigens[3,], c(0.3984175, 0.3984175), tolerance=1e-5)
  expect_equal(sum222cm$omega_eigens[1:4], c(11.8447678, 3.7262580, 5.1904506, 0.8549476), tolerance=1e-5)
})
