context("quantile residuals")
library(gmvarkit)

# NOTE that some test use random elements that might change if random_ind2 or simulateGMVAR is modified

## A(M)(p)_(p)(M)(d)

# p=1, M=1, d=2, parametrization="mean"
phi10_112 <- c(0.75, 0.8)
A11_112 <- matrix(c(0.29, 0.02, -0.14, 0.9), nrow=2, byrow=FALSE)
Omega1_112 <- matrix(c(0.60, 0.01, 0.01, 0.07), nrow=2, byrow=FALSE)
theta_112 <- c(phi10_112, vec(A11_112), vech(Omega1_112))
mod_112 <- GMVAR(gdpdef, p=1, M=1, d=2, params=theta_112, conditional=TRUE, parametrization="mean", constraints=NULL)

W_112 <- t(chol(Omega1_112))
theta_112sWC <- c(phi10_112, vec(A11_112), Wvec(W_112)) # SGMVAR, W constrained by Cholesky
mod_112s <- GMVAR(gdpdef, p=1, M=1, d=2, params=theta_112sWC, conditional=TRUE, parametrization="mean", constraints=NULL,
                  structural_pars=list(W=W_112))

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

# p=1, M=4, d=2
theta_142 <- c(0.680381, 0.099736, 0.502918, 0.080781, -0.627348, 0.674579,
               0.37666, 0.020433, 0.016747, 0.285857, 0.205766, 0.392568,
               0.047474, 0.317407, 0.683117, 0.415324, -0.059849, 0.079795,
               1.927008, 0.687905, 0.036475, -0.014841, -0.793236, 0.638711,
               1.281068, 0.017391, 0.135752, 1.716725, 0.668851, -0.184389,
               -0.155109, -1.746672, -0.348001, 0.43192, -0.028927, 0.002053,
               0.428343, 0.347302, 0.152716)
mod_142_int <- GMVAR(gdpdef, p=1, M=4, d=2, params=theta_142, conditional=TRUE, parametrization="intercept")
mod_142_mean <- GMVAR(gdpdef, p=1, M=4, d=2, params=theta_142, conditional=TRUE, parametrization="mean")


# p=2, M=2, d=2, SGMVAR AR params constrained to be the same in both regimes
rbind_diags <- function(p, M, d) {
  I <- diag(p*d^2)
  Reduce(rbind, replicate(M, I, simplify=FALSE))
}
C_222 <- rbind_diags(p=2, M=2, d=2)
C_lambda_222 <- matrix(c(7, 1), nrow=2)
theta_222cs <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(W_222), 1, alpha1_222)
mod_222csLAR <- GMVAR(gdpdef, p=2, M=2, d=2, params=theta_222cs, conditional=TRUE, parametrization="intercept",
                      constraints=C_222, structural_pars=list(W=W_222, C_lambda=C_lambda_222))

# p=2, M=2, d=2, AR parameters same for both regimes, non-diagonals zero, intercept
theta_222c <- c(0.33782, 0.183512, 0.472168, 0.095311, 0.201199, 0.600596, 0.237819,
                0.23529, 1.077816, -0.016343, 0.112771, 0.22199, 0.005582, 0.028126, 0.492844)
mat0 <- matrix(c(1, rep(0, 10), 1, rep(0, 8), 1, rep(0, 10), 1), nrow=2*2^2, byrow=FALSE)
C_222c <- rbind(mat0, mat0)
mod_222c <- GMVAR(gdpdef, p=2, M=2, d=2, params=theta_222c, conditional=TRUE, parametrization="intercept", constraints=C_222c)


# p=3, M=2, d=3, no constraints, rand_ind and simulated data
set.seed(42)
theta_323 <- random_ind2(p=3, M=2, d=3, mu_scale=c(-10, 0, 5), mu_scale2=1:3, omega_scale=1:3, ar_scale=1)
mod_323 <- GMVAR(p=3, M=2, d=3, params=theta_323, conditional=FALSE, parametrization="mean", constraints=NULL)
sim_323 <- simulateGMVAR(mod_323, nsimu=500)$sample
mod_323 <- add_data(data=sim_323, gmvar=mod_323)

# p=2, M=2, d=2, constraints=C_222, structural_pars=list(W=W_222, C_lambda=C_lambda_222), same_means=list(1:2)
theta_222csLAR_int <- c(phi10_222, vec(A11_222), vec(A12_222), vec(W_222), 0.2, alpha1_222)
mod_222csLAR_int <- GMVAR(gdpdef, p=2, M=2, params=theta_222csLAR_int, conditional=TRUE, parametrization="mean",
                          constraints=C_222, structural_pars=list(W=W_222, C_lambda=C_lambda_222), same_means=list(1:2))

# Get quantile residuals
res_112 <- quantile_residuals(mod_112)
res_222 <- quantile_residuals(mod_222)
res_222c <- quantile_residuals(mod_222c)
res_323 <- quantile_residuals(mod_323)

res_142_int <- quantile_residuals(mod_142_int)
res_142_mean <- quantile_residuals(mod_142_mean)

res_112s <- quantile_residuals(mod_112s)
res_222s <- quantile_residuals(mod_222s)
res_222csLAR <- quantile_residuals(mod_222csLAR)

res222csLAR_int <- quantile_residuals(mod_222csLAR_int)


test_that("quantile_residuals works correctly", {
  expect_equal(res_112[1,], c(1.402850, -0.674736), tolerance=1e-6)
  expect_equal(res_112[2,], c(-1.5576025, 0.5953575), tolerance=1e-6)
  expect_equal(res_112[100,], c(1.0476053, 0.7918417), tolerance=1e-6)
  expect_equal(res_112[243,], c(-0.4242621, 0.1649217), tolerance=1e-6)

  expect_equal(res_222[1:2, 1], c(-1.1163356, -0.3784424), tolerance=1e-6)
  expect_equal(res_222[100:102, 1], c(0.5437907, -0.5027021, -0.5797931), tolerance=1e-6)
  expect_equal(res_222[212:215, 1], c(-0.7499078, -1.2303636, -0.8677514, 0.8693726), tolerance=1e-6)
  expect_equal(res_222[1:2, 2], c(0.2045036, -0.1337020), tolerance=1e-6)
  expect_equal(res_222[100:102, 2], c(-0.14768871, 0.11550758, -0.09258657), tolerance=1e-6)
  expect_equal(res_222[212:215, 2], c(-0.1349727, 0.3764042, 0.4159688, -0.5631334), tolerance=1e-6)

  expect_equal(res_222c[1,], c(-1.3793624, 0.1456127), tolerance=1e-6)
  expect_equal(res_222c[100,], c(0.5482568, -0.2830299), tolerance=1e-6)

  expect_equal(res_323[1,], c(-0.04254252, 0.37966523, -1.43671034), tolerance=1e-6)
  expect_equal(res_323[13,], c(-1.1160317, -0.1883617, -0.4064072), tolerance=1e-6)
  expect_equal(res_323[150,], c(0.09401363, 0.44800577, 1.11962898), tolerance=1e-6)
  expect_equal(res_323[497,], c(0.3658767, 0.4336302, -0.6696751), tolerance=1e-6)

  expect_equal(res_142_int[122:125, 1], c(-0.2421819, -0.7722385, 1.0008730, -0.9905129), tolerance=1e-6)
  expect_equal(res_142_int[12:15, 2], c(0.4100877, -2.2630865, -0.6945604, -0.9637667), tolerance=1e-6)
  expect_equal(res_142_mean[132:135, 1], c(0.7112074, 0.3783539, -0.3482674, 0.2798411), tolerance=1e-6)
  expect_equal(res_142_int[2:5, 2], c(0.03564512, 0.07320782, -0.24077470, -0.44959987), tolerance=1e-6)

  # SGMVAR
  expect_equal(res_112s[3,], c(-0.42107591, -0.04453619), tol=1e-6)
  expect_equal(res_112s, res_112, tol=1e-6)
  expect_equal(res_222s, res_222, tol=1e-6)
  expect_equal(res_222csLAR[1,], c(-1.4608297, -0.1445763), tol=1e-6)
  expect_equal(res_222csLAR[242,], c(-0.7107413, -0.1398955), tol=1e-6)

  # Same_means
  expect_equal(res222csLAR_int[c(1, 5, 100, 200)], c(-2.6082444, -0.9789409, 0.8552812, 1.5705630), tolerance=1e-6)
})

test_that("quantile_residuals_int works correctly", {
  qr112 <- quantile_residuals_int(gdpdef, p=1, M=1, params=theta_112, conditional=TRUE, parametrization="mean", constraints=NULL)
  qr222csLAR <- quantile_residuals_int(gdpdef, p=2, M=2, params=theta_222cs, conditional=TRUE, parametrization="intercept",
                                       constraints=C_222, structural_pars=list(W=W_222, C_lambda=C_lambda_222))
  qr222csLAR_int <- quantile_residuals_int(gdpdef, p=2, M=2, params=theta_222csLAR_int, conditional=TRUE, parametrization="mean",
                                           constraints=C_222, structural_pars=list(W=W_222, C_lambda=C_lambda_222), same_means=list(1:2))
  expect_equal(qr112, res_112, tol=1e-6)
  expect_equal(qr222csLAR, res_222csLAR, tol=1e-6)
  expect_equal(qr222csLAR_int[1:3], c(-2.608244, -1.453003, 4.221412), tolerance=1e-6)
})








