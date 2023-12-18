context("quantile residuals")
library(gmvarkit)

# NOTE that some test use random elements that might change if simulate.gsmvar is modified
# ALSO: These tests do not really cover the numerical integration as is!
# Set: which_def <- numeric(0) in the appropriate place in the function quantile_residuals to test numerical integration.

## A(M)(p)_(p)(M)(d)

# p=1, M=1, d=2, parametrization="mean"
phi10_112 <- c(0.75, 0.8)
A11_112 <- matrix(c(0.29, 0.02, -0.14, 0.9), nrow=2, byrow=FALSE)
Omega1_112 <- matrix(c(0.60, 0.01, 0.01, 0.07), nrow=2, byrow=FALSE)
theta_112 <- c(phi10_112, vec(A11_112), vech(Omega1_112))
mod_112 <- GSMVAR(gdpdef, p=1, M=1, d=2, params=theta_112, conditional=TRUE, parametrization="mean", calc_std_errors=FALSE)

W_112 <- t(chol(Omega1_112))
theta_112sWC <- c(phi10_112, vec(A11_112), Wvec(W_112)) # SGMVAR, W constrained by Cholesky
mod_112s <- GSMVAR(gdpdef, p=1, M=1, d=2, params=theta_112sWC, conditional=TRUE, parametrization="mean",
                  structural_pars=list(W=W_112), calc_std_errors=FALSE)

mod_112t <- GSMVAR(gdpdef, p=1, M=1, d=2, params=c(theta_112, 10), model="StMVAR", conditional=TRUE, parametrization="mean", calc_std_errors=FALSE)
mod_112ts <- GSMVAR(gdpdef, p=1, M=1, d=2, params=c(theta_112sWC, 10), model="StMVAR", conditional=TRUE, parametrization="mean",
                   structural_pars=list(W=W_112), calc_std_errors=FALSE)

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
mod_222 <- GSMVAR(gdpdef, p=2, M=2, d=2, params=theta_222, conditional=TRUE, parametrization="intercept", calc_std_errors=FALSE)

WL_222 <- diag_Omegas(Omega1_222, Omega2_222)
W_222 <- matrix(WL_222[1:(2^2)], nrow=2, byrow=FALSE)
lambdas_222 <- WL_222[(2^2 + 1):length(WL_222)]
theta_222s <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A21_222),
                vec(A22_222), vec(W_222), lambdas_222, alpha1_222) # SGMVAR
mod_222s <- GSMVAR(gdpdef, p=2, M=2, d=2, params=theta_222s, conditional=TRUE, parametrization="intercept", calc_std_errors=FALSE,
                 structural_pars=list(W=W_222))

mod_222t <- GSMVAR(gdpdef, p=2, M=2, d=2, params=c(theta_222, 10, 20), model="StMVAR", conditional=TRUE, parametrization="intercept",
                  calc_std_errors=FALSE)
mod_222ts <- GSMVAR(gdpdef, p=2, M=2, d=2, params=c(theta_222s, 10, 20), model="StMVAR", conditional=TRUE, parametrization="intercept",
                   calc_std_errors=FALSE, structural_pars=list(W=W_222))
mod_222gs <- GSMVAR(gdpdef, p=2, M=c(1, 1), d=2, params=c(theta_222, 20), model="G-StMVAR", conditional=TRUE, parametrization="intercept",
                   calc_std_errors=FALSE)

# p=1, M=4, d=2
theta_142 <- c(0.680381, 0.099736, 0.502918, 0.080781, -0.627348, 0.674579,
               0.37666, 0.020433, 0.016747, 0.285857, 0.205766, 0.392568,
               0.047474, 0.317407, 0.683117, 0.415324, -0.059849, 0.079795,
               1.927008, 0.687905, 0.036475, -0.014841, -0.793236, 0.638711,
               1.281068, 0.017391, 0.135752, 1.716725, 0.668851, -0.184389,
               -0.155109, -1.746672, -0.348001, 0.43192, -0.028927, 0.002053,
               0.428343, 0.347302, 0.152716)
mod_142_int <- GSMVAR(gdpdef, p=1, M=4, d=2, params=theta_142, conditional=TRUE, parametrization="intercept", calc_std_errors=FALSE)
mod_142_mean <- GSMVAR(gdpdef, p=1, M=4, d=2, params=theta_142, conditional=TRUE, parametrization="mean", calc_std_errors=FALSE)

mod_142gs_1 <- GSMVAR(gdpdef, p=1, M=c(1, 3), d=2, params=c(theta_142, 20, 30, 40), model="G-StMVAR", conditional=TRUE,
                     parametrization="intercept", calc_std_errors=FALSE)
mod_142gs_1_mean <- GSMVAR(gdpdef, p=1, M=c(1, 3), d=2, params=c(theta_142, 20, 30, 40), model="G-StMVAR", conditional=TRUE,
                          parametrization="mean", calc_std_errors=FALSE)
mod_142gs_2 <- GSMVAR(gdpdef, p=1, M=c(3, 1), d=2, params=c(theta_142, 3), model="G-StMVAR", conditional=TRUE,
                     parametrization="intercept", calc_std_errors=FALSE)
mod_142gs_3 <- GSMVAR(gdpdef, p=1, M=c(2, 2), d=2, params=c(theta_142, 10, 12), model="G-StMVAR", conditional=TRUE,
                     parametrization="intercept", calc_std_errors=FALSE)
mod_142t <- GSMVAR(gdpdef, p=1, M=4, d=2, params=c(theta_142, 10, 20, 30, 40), model="StMVAR", conditional=TRUE,
                  parametrization="intercept", calc_std_errors=FALSE)


# p=2, M=2, d=2, SGMVAR AR params constrained to be the same in both regimes
rbind_diags <- function(p, M, d) {
  I <- diag(p*d^2)
  Reduce(rbind, replicate(M, I, simplify=FALSE))
}
C_222 <- rbind_diags(p=2, M=2, d=2)
C_lambda_222 <- matrix(c(7, 1), nrow=2)
theta_222cs <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(W_222), 1, alpha1_222)
mod_222csLAR <- GSMVAR(gdpdef, p=2, M=2, d=2, params=theta_222cs, conditional=TRUE, parametrization="intercept",
                      constraints=C_222, structural_pars=list(W=W_222, C_lambda=C_lambda_222), calc_std_errors=FALSE)

mod_222gscsLAR <- GSMVAR(gdpdef, p=2, M=c(1, 1), d=2, params=c(theta_222cs, 2.001), model="G-StMVAR", conditional=TRUE, parametrization="intercept",
                        constraints=C_222, structural_pars=list(W=W_222, C_lambda=C_lambda_222), calc_std_errors=FALSE)

# p=2, M=2, d=2, AR parameters same for both regimes, non-diagonals zero, intercept
theta_222c <- c(0.33782, 0.183512, 0.472168, 0.095311, 0.201199, 0.600596, 0.237819,
                0.23529, 1.077816, -0.016343, 0.112771, 0.22199, 0.005582, 0.028126, 0.492844)
mat0 <- matrix(c(1, rep(0, 10), 1, rep(0, 8), 1, rep(0, 10), 1), nrow=2*2^2, byrow=FALSE)
C_222c <- rbind(mat0, mat0)
mod_222c <- GSMVAR(gdpdef, p=2, M=2, d=2, params=theta_222c, conditional=TRUE, parametrization="intercept", constraints=C_222c,
                  calc_std_errors=FALSE)

mod_222tc <- GSMVAR(gdpdef, p=2, M=2, d=2, params=c(theta_222c, 5, 11), model="StMVAR", conditional=TRUE, parametrization="intercept", constraints=C_222c,
                   calc_std_errors=FALSE)

# p=3, M=2, d=3, no constraints, rand_ind and simulated data
set.seed(42)
theta_323 <- c(-8.62904, -1.1294, 6.08939, 0.30181, -0.01918, 0.20814, 2.09997, -0.76256, 0.34977, 1.08522, 0.1862, 0.56666,
               -0.20946, -0.19745, 0.41976, 0.37882, 0.36577, -0.64832, -1.56472, 0.62412, 0.67048, -0.21841, -0.34758,
               0.38573, 0.28527, 0.09343, -0.55505, -0.76295, 0.16713, -0.30743, 1.00487, 0.84722, 0.50618, 2.21109, -1.67388,
               3.40942, -10.78446, -1.70182, -2.24262, 0.72048, -0.12182, 0.60088, 0.32822, -0.04847, -0.83748, -0.89225,
               -1.14077, -0.29568, -0.61784, -0.36704, -0.79401, -0.92552, 0.51123, -0.52166, -0.45825, -0.24572, 0.3612,
               -0.10096, -0.58217, 0.17453, -0.19763, -0.2274, 0.77741, 1.09703, 1.24502, 0.74678, 0.90936, 0.71689, -0.086,
               2.57044, 0.81515, 0.45869, 0.84737)
mod_323 <- GSMVAR(p=3, M=2, d=3, params=theta_323, conditional=FALSE, parametrization="mean")
sim_323 <- simulate.gsmvar(mod_323, nsim=500)$sample
mod_323 <- add_data(data=sim_323, gsmvar=mod_323)

mod_323t <- GSMVAR(sim_323, p=3, M=2, params=c(theta_323, 10, 3), model="StMVAR", conditional=FALSE,
                  parametrization="mean", calc_std_errors=FALSE)
mod_323gs <- GSMVAR(sim_323, p=3, M=c(1, 1), params=c(theta_323, 3), model="G-StMVAR", conditional=FALSE,
                  parametrization="mean", calc_std_errors=FALSE)

# p=2, M=2, d=2, constraints=C_222, structural_pars=list(W=W_222, C_lambda=C_lambda_222), same_means=list(1:2)
theta_222csLAR_int <- c(phi10_222, vec(A11_222), vec(A12_222), vec(W_222), 0.2, alpha1_222)
mod_222csLAR_int <- GSMVAR(gdpdef, p=2, M=2, params=theta_222csLAR_int, conditional=TRUE, parametrization="mean",
                          constraints=C_222, structural_pars=list(W=W_222, C_lambda=C_lambda_222), same_means=list(1:2),
                          calc_std_errors=FALSE)

mod_222tcsLAR_int <- GSMVAR(gdpdef, p=2, M=2, params=c(theta_222csLAR_int, 20, 2.1), model="StMVAR", conditional=TRUE, parametrization="mean",
                          constraints=C_222, structural_pars=list(W=W_222, C_lambda=C_lambda_222), same_means=list(1:2),
                          calc_std_errors=FALSE)

# p=2, M=2, d=2, model="GMVAR", constraints=C_222, same_means=list(1:2), weight_constraints=0.6,
# structural_pars=list(W=W_222, fixed_lambdas=c(7, 3)), parametrizarion="mean"
theta_222cmwsF <- c(phi10_222, vec(A11_222), vec(A12_222), vec(W_222))
mod_222cmwsF <- GSMVAR(gdpdef,  p=2, M=2, params=theta_222cmwsF, model="GMVAR", constraints=C_222, same_means=list(1:2),
                       weight_constraints=0.6, parametrization="mean", structural_pars=list(W=W_222, fixed_lambdas=c(7, 3)))

# p=2, M=2, d=2, model="StMVAR", constraints=C_222, weight_constraints=0.6,
# structural_pars=list(W=W_222, fixed_lambdas=c(7, 3))
theta_222tcwsF <- c(phi10_222, c(1.79, 3.00), vec(A11_222), vec(A12_222), vec(W_222), 11, 12)
mod_222tcwsF <- GSMVAR(gdpdef,  p=2, M=2, params=theta_222tcwsF, model="StMVAR", constraints=C_222, weight_constraints=0.6,
                       structural_pars=list(W=W_222, fixed_lambdas=c(7, 3)))

# p=2, M=c(1, 1), d=2, model="G-StMVAR", constraints=C_222, structural_pars=list(W=W_222, fixed_lambdas=c(7, 3))
theta_222gscsF <- c(phi10_222, c(1.79, 3.00), vec(A11_222), vec(A12_222), vec(W_222), 0.6, 11)
mod_222gscsF <- GSMVAR(gdpdef,  p=2, M=c(1, 1), params=theta_222gscsF, model="G-StMVAR", constraints=C_222,
                       structural_pars=list(W=W_222, fixed_lambdas=c(7, 3)))



# Calculate quantile residuals
res_112 <- quantile_residuals(mod_112)
res_112t <- quantile_residuals(mod_112t) # StMVAR
res_222 <- quantile_residuals(mod_222)
res_222t <- quantile_residuals(mod_222t) # StMVAR
res_222gs <- quantile_residuals(mod_222gs) # G-StMVAR

res_142_int <- quantile_residuals(mod_142_int)
res_142_mean <- quantile_residuals(mod_142_mean)
res_142gs_1 <- quantile_residuals(mod_142gs_1) # G-StMVAR
res_142gs_1_mean <- quantile_residuals(mod_142gs_1_mean) # G-StMVAR
res_142gs_2 <- quantile_residuals(mod_142gs_2) # G-StMVAR
res_142gs_3 <- quantile_residuals(mod_142gs_3) # G-StMVAR
res_142t <- quantile_residuals(mod_142t) # StMVAR

res_222c <- quantile_residuals(mod_222c)
res_222tc <- quantile_residuals(mod_222tc) # StMVAR

res_323 <- quantile_residuals(mod_323)
res_323t <- quantile_residuals(mod_323t) # StMVAR
res_323gs <- quantile_residuals(mod_323gs) # G-StMVAR

res_112s <- quantile_residuals(mod_112s)
res_112ts <- quantile_residuals(mod_112ts) # SStMVAR
res_222s <- quantile_residuals(mod_222s)
res_222ts <- quantile_residuals(mod_222ts) # SStMVAR
res_222csLAR <- quantile_residuals(mod_222csLAR)
res_222gscsLAR <- quantile_residuals(mod_222gscsLAR) # G-StMVAR

res222csLAR_int <- quantile_residuals(mod_222csLAR_int)
res222tcsLAR_int <- quantile_residuals(mod_222tcsLAR_int) # StMVAR

res_222cmwsF <- quantile_residuals(mod_222cmwsF)
res_222tcwsF <- quantile_residuals(mod_222tcwsF)
res_222gscsF <- quantile_residuals(mod_222gscsF)


test_that("quantile_residuals works correctly", {
  expect_equal(res_112[1,], c(1.402850, -0.674736), tolerance=1e-6)
  expect_equal(res_112[2,], c(-1.5576025, 0.5953575), tolerance=1e-6)
  expect_equal(res_112[100,], c(1.0476053, 0.7918417), tolerance=1e-6)
  expect_equal(res_112[243,], c(-0.4242621, 0.1649217), tolerance=1e-6)

  expect_equal(res_112t[1,], c(1.3929479, -0.6148843), tolerance=1e-6)
  expect_equal(res_112t[12,], c(0.8126274, 0.4278665), tolerance=1e-6)
  expect_equal(res_112t[243,], c(-0.4870697, 0.1794448), tolerance=1e-6)

  expect_equal(res_222[1:2, 1], c(-1.1163356, -0.3784424), tolerance=1e-6)
  expect_equal(res_222[100:102, 1], c(0.5437907, -0.5027021, -0.5797931), tolerance=1e-6)
  expect_equal(res_222[212:215, 1], c(-0.7499078, -1.2303636, -0.8677514, 0.8693726), tolerance=1e-6)
  expect_equal(res_222[1:2, 2], c(0.2045036, -0.1337020), tolerance=1e-6)
  expect_equal(res_222[100:102, 2], c(-0.14768871, 0.11550758, -0.09258657), tolerance=1e-6)
  expect_equal(res_222[212:215, 2], c(-0.1349727, 0.3764042, 0.4159688, -0.5631334), tolerance=1e-6)

  expect_equal(res_222t[1, ], c(-1.1850393, 0.2037219), tolerance=1e-6)
  expect_equal(res_222t[13, ], c(0.0704471, -0.7364274), tolerance=1e-6)

  expect_equal(res_222gs[11, ], c(0.6028979, 0.4206634), tolerance=1e-6)
  expect_equal(res_222gs[200, ], c(-0.1850893, -0.4409537), tolerance=1e-6)

  expect_equal(res_222c[1,], c(-1.3793624, 0.1456127), tolerance=1e-6)
  expect_equal(res_222c[100,], c(0.5482568, -0.2830299), tolerance=1e-6)

  expect_equal(res_222tc[1,], c(-1.4847954, 0.2748462), tolerance=1e-6)
  expect_equal(res_222tc[11,], c(0.7904500, 0.6851541), tolerance=1e-6)

  expect_equal(res_323[1,], c(1.2789657, -0.4454311, 2.1568430), tolerance=1e-6)
  expect_equal(res_323[13,], c(0.9090539, 0.7672014, 1.4780458), tolerance=1e-6)
  expect_equal(res_323[150,], c(1.4215695, 0.6164883, 0.8540793), tolerance=1e-6)
  expect_equal(res_323[497,], c(0.4797046, -0.7006484, 0.7415371), tolerance=1e-6)

  expect_equal(res_323t[1,], c(0.7219346, -0.7109976, 1.4089146), tolerance=1e-6)
  expect_equal(res_323t[111,], c(1.2590007, -0.1740127, -0.4653121), tolerance=1e-6)

  expect_equal(res_323gs[1,], c(0.520961, -1.976524, 2.140768), tolerance=1e-6)
  expect_equal(res_323t[222,], c(-1.1211720, 0.8204242, -1.6024847), tolerance=1e-6) ##### huom jäi t

  expect_equal(res_142_int[122:125, 1], c(-0.2421819, -0.7722385, 1.0008730, -0.9905129), tolerance=1e-6)
  expect_equal(res_142_int[12:15, 2], c(0.4100877, -2.2630865, -0.6945604, -0.9637667), tolerance=1e-6)
  expect_equal(res_142_mean[132:135, 1], c(0.7112074, 0.3783539, -0.3482674, 0.2798411), tolerance=1e-6)
  expect_equal(res_142_int[2:5, 2], c(0.03564512, 0.07320782, -0.24077470, -0.44959987), tolerance=1e-6)

  expect_equal(res_142gs_1[2,], c(-1.69926459, 0.08292301), tolerance=1e-6)
  expect_equal(res_142gs_1_mean[3,], c(0.04384774, 0.54117696), tolerance=1e-6)
  expect_equal(res_142gs_2[4,], c(2.5858242, -0.3039629), tolerance=1e-6)
  expect_equal(res_142gs_3[5,], c(-2.2394528, -0.4245589), tolerance=1e-6)
  expect_equal(res_142t[50,], c(-0.1257199, -0.5621444), tolerance=1e-6)

  # Structural
  expect_equal(res_112s[3,], c(-0.42107591, -0.04453619), tol=1e-6)
  expect_equal(res_112ts[203,], c(0.4678634, 0.6003873), tol=1e-6)
  expect_equal(res_112s, res_112, tol=1e-6)
  expect_equal(res_222s, res_222, tol=1e-6)
  expect_equal(res_222ts[211,], c(0.1871982, 0.8858166), tol=1e-6)
  expect_equal(res_222gs[111,], c(-0.1466454, 0.9015567), tol=1e-6)
  expect_equal(res_222csLAR[1,], c(-1.4608297, -0.1445763), tol=1e-6)
  expect_equal(res_222gscsLAR[200,], c(0.3219058, -0.8329321), tol=1e-6)
  expect_equal(res_222csLAR[242,], c(-0.7107413, -0.1398955), tol=1e-6)

  # Same_means
  expect_equal(res222csLAR_int[c(1, 5, 100, 200)], c(-2.6082444, -0.9789409, 0.8552812, 1.5705630), tolerance=1e-6)
  expect_equal(res222tcsLAR_int[c(3, 7, 14, 201)], c(3.538059, 1.003902, -1.036117, 1.318223), tolerance=1e-6)

  # Fixed alphas and lambdas
  expect_equal(c(res_222cmwsF[c(1, 13),]), c(-1.1627924, 0.1856571, 0.3064827, -0.3428420), tolerance=1e-6)
  expect_equal(c(res_222tcwsF[c(1, 242),]), c(-2.4532966, -0.9148971, 0.3486117, -0.1578582), tolerance=1e-6)
  expect_equal(c(res_222gscsF[c(1, 2, 242),]), c(-3.4222887, -2.1161398, -0.7937896, 0.5470228, 0.8390555, -0.1438871), tolerance=1e-6)
})

test_that("quantile_residuals_int works correctly", {
  qr112 <- quantile_residuals_int(gdpdef, p=1, M=1, params=theta_112, conditional=TRUE, parametrization="mean")
  qr222csLAR <- quantile_residuals_int(gdpdef, p=2, M=2, params=theta_222cs, conditional=TRUE, parametrization="intercept",
                                       constraints=C_222, structural_pars=list(W=W_222, C_lambda=C_lambda_222))
  qr222csLAR_int <- quantile_residuals_int(gdpdef, p=2, M=2, params=theta_222csLAR_int, conditional=TRUE, parametrization="mean",
                                           constraints=C_222, structural_pars=list(W=W_222, C_lambda=C_lambda_222), same_means=list(1:2))

  qr222t <- quantile_residuals_int(gdpdef, p=2, M=2, params=c(theta_222, 3, 4), model="StMVAR", conditional=TRUE,
                                   parametrization="mean")
  qr222gs <- quantile_residuals_int(gdpdef, p=2, M=c(1, 1), params=c(theta_222, 4), model="G-StMVAR", conditional=TRUE,
                                    parametrization="mean")

  qr222tcwsF <- quantile_residuals_int(gdpdef, p=2, M=2, params=theta_222tcwsF, model="StMVAR", constraints=C_222, weight_constraints=0.6,
                                       structural_pars=list(W=W_222, fixed_lambdas=c(7, 3)), parametrization="intercept", conditional=TRUE)

  expect_equal(qr112, res_112, tolerance=1e-6)
  expect_equal(qr222csLAR, res_222csLAR, tolerance=1e-6)
  expect_equal(qr222csLAR_int[1:3], c(-2.608244, -1.453003, 4.221412), tolerance=1e-6)
  expect_equal(qr222t[c(1, 13, 6, 206)], c(-1.0250704, 0.5404009, -2.1362256, -0.2425385), tolerance=1e-6)
  expect_equal(qr222gs[c(1, 14, 7, 208)], c(-1.0049791, -0.7887779, 0.7620531, 0.6624912), tolerance=1e-6)
  expect_equal(qr222tcwsF, res_222tcwsF, tolerance=1e-6)
  expect_equal(qr222tcwsF[c(1, 2, 300)], c(-2.453297, -1.499268, 1.705310), tolerance=1e-6)

})








