context("quantile residual tests")
library(gmvarkit)

# These tests use random elements and some of the numerical values may change is simulatGSMAR is modified

## A(M)(p)_(p)(M)(d)

# p=1, M=1, d=2, parametrization="mean"
phi10_112 <- c(0.75, 0.8)
A11_112 <- matrix(c(0.29, 0.02, -0.14, 0.9), nrow=2, byrow=FALSE)
Omega1_112 <- matrix(c(0.60, 0.01, 0.01, 0.07), nrow=2, byrow=FALSE)
theta_112 <- c(phi10_112, vec(A11_112), vech(Omega1_112))
mod_112 <- GSMVAR(gdpdef, p=1, M=1, d=2, params=theta_112, conditional=TRUE, parametrization="mean", constraints=NULL)

W_112 <- t(chol(Omega1_112))
theta_112sWC <- c(phi10_112, vec(A11_112), Wvec(W_112)) # SGMVAR, W constrained by Cholesky
mod_112s <- GSMVAR(gdpdef, p=1, M=1, d=2, params=theta_112sWC, conditional=TRUE, parametrization="mean", constraints=NULL,
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
mod_222 <- GSMVAR(gdpdef, p=2, M=2, d=2, params=theta_222, conditional=TRUE, parametrization="intercept", constraints=NULL)

mod_222t <- GSMVAR(gdpdef, p=2, M=2, d=2, params=c(theta_222, 20, 30),  model="StMVAR",
                   conditional=TRUE, parametrization="intercept", constraints=NULL)
mod_222gs <- GSMVAR(gdpdef, p=2, M=c(1, 1), d=2, params=c(theta_222, 30),  model="G-StMVAR",
                   conditional=TRUE, parametrization="intercept", constraints=NULL)


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
set.seed(11)
theta_123 <- random_ind2(p=1, M=2, d=3, mu_scale=c(-10, 0, 5), mu_scale2=1:3, omega_scale=1:3, ar_scale=1)
mod_123 <- GSMVAR(p=1, M=2, d=3, params=theta_123, conditional=FALSE, parametrization="mean", constraints=NULL)
sim_123 <- simulate.gsmvar(mod_123, nsim=500)
data_123 <- sim_123$sample
mod_123 <- GSMVAR(data_123, p=1, M=2, d=3, params=theta_123, conditional=FALSE, parametrization="mean", constraints=NULL)

# p=2, M=2, d=2, parametrization="mean", constraints=C_mat, same_means=list(1:2)
C_mat <- rbind(diag(2*2^2), diag(2*2^2))
params_222cm <- c(0.811034, 0.578587, 0.212084, 0.020444, -0.193005, 0.624671,
                  0.235827, 0.013962, 0.053267, 0.262703, 1.06105, -0.013519,
                  0.114109, 0.229542, 0.003092, 0.027266, 0.424341)
mod_222cm <- GSMVAR(gdpdef, p=2, M=2, params=params_222cm, parametrization="mean", constraints=C_mat, same_means=list(1:2))

set.seed(1)
res_112 <- quantile_residual_tests(mod_112, lags_ac=1:2, lags_ch=1:2, nsim=300, print_res=FALSE)
res_222 <- quantile_residual_tests(mod_222, lags_ac=3, lags_ch=1, nsim=1, print_res=FALSE)
res_123 <- quantile_residual_tests(mod_123, lags_ac=1, lags_ch=2, nsim=1, print_res=FALSE)

set.seed(1); res_222t <- quantile_residual_tests(mod_222t, lags_ac=1, lags_ch=2, nsim=1, print_res=FALSE)
set.seed(1); res_222gs <- quantile_residual_tests(mod_222gs, lags_ac=2, lags_ch=1, nsim=1, print_res=FALSE)

set.seed(1); res_112s <- quantile_residual_tests(mod_112s, lags_ac=1:2, lags_ch=1:2, nsim=1, print_res=FALSE)
set.seed(1); res_222s <- quantile_residual_tests(mod_222s, lags_ac=1, lags_ch=2, nsim=300, print_res=FALSE)

set.seed(1); res_222cm <- quantile_residual_tests(mod_222cm, lags_ac=2, lags_ch=1, nsim=1, print_res=FALSE)

test_that("quantile_residual_tests - test_results - works correctly", {
  expect_equal(res_112$norm_res$test_stat, 154.3486, tolerance=1)
  expect_equal(res_112$ac_res$test_results$test_stat, c(21.13147, 28.86697), tolerance=1e-4)
  expect_equal(res_112$ch_res$test_results$test_stat, c(61.29593, 204.98674), tolerance=1e-4)

  expect_equal(res_222$norm_res$p_val, 0.01583773, tolerance=1e-4)
  expect_equal(res_222$ac_res$test_results$p_val, 0.1168246, tolerance=1e-4)
  expect_equal(res_222$ch_res$test_results$test_stat, 5.045936, tolerance=1e-4)

  expect_equal(res_222t$norm_res$p_val, 0.1554327, tolerance=1e-4)
  expect_equal(res_222t$ac_res$test_results$p_val, 0.2866306, tolerance=1e-4)
  expect_equal(res_222t$ch_res$test_results$p_val, 0.424535, tolerance=1e-4)

  expect_equal(res_222gs$norm_res$p_val, 0.05590034, tolerance=1e-4)
  expect_equal(res_222gs$ac_res$test_results$p_val, 0.07231058, tolerance=1e-4)
  expect_equal(res_222gs$ch_res$test_results$p_val, 0.2983806, tolerance=1e-4)

  expect_equal(res_123$norm_res$test_stat, 22.50893, tolerance=1e-4)
  expect_equal(res_123$ac_res$test_results$test_stat, 13.70742, tolerance=1e-4)
  expect_equal(res_123$ch_res$test_results$test_stat, 17.73436, tolerance=1e-4)

  # SGSMVAR
  expect_equal(res_112s$norm_res$test_stat, 17.45124, tolerance=1e-3)
  expect_equal(res_112s$ac_res$test_results$test_stat, c(15.53686, 29.04951), tolerance=1e-3)
  expect_equal(res_112s$ch_res$test_results$test_stat, c(10.21175, 24.40845), tolerance=1e-3)

  expect_equal(res_222s$norm_res$p_val, 0.008607749, tolerance=1e-4)
  expect_equal(res_222s$ac_res$test_results$p_val, 0.4742369, tolerance=1e-4)
  expect_equal(res_222s$ch_res$test_results$p_val, 0.4440993, tolerance=1e-4)

  # Same means
  expect_equal(res_222cm$norm_res$p_val, 0.2152119, tolerance=1e-4)
  expect_equal(res_222cm$ac_res$test_results$p_val, 0.05002536, tolerance=1e-4)
  expect_equal(res_222cm$ch_res$test_results$p_val, 0.270333, tolerance=1e-4)
})

test_that("quantile_residual_tests - ind_stats - works correctly", {
  expect_equal(res_112$ac_res$ind_stats$lag2, c(2.5963913, 0.1874684, 0.7390084, -0.9374877), tolerance=1e-4)
  expect_equal(res_112$ch_res$ind_stats$lag1, c(2.3595939, -0.1535945, 3.6163184, 5.7548554), tolerance=1e-4)

  expect_equal(res_222$ac_res$ind_stats$lag3, c(-1.022108, 1.337845, -1.045152, 1.079755), tolerance=1e-4)
  expect_equal(res_222$ch_res$ind_stats$lag1, c(-0.4215779, -1.0604109, -0.6596756, 1.5901816), tolerance=1e-4)

  expect_equal(res_222t$ac_res$ind_stats$lag1, c(-0.0817922, -0.7530136, 0.1631821, -1.6861624), tolerance=1e-4)
  expect_equal(res_222t$ch_res$ind_stats$lag2, c(-0.6548325, -0.8686823, 0.6983760, 1.3222143), tolerance=1e-4)

  expect_equal(res_222gs$ac_res$ind_stats$lag2, c(0.4558164, -0.3540209, 0.6928812, -3.3972828), tolerance=1e-4)
  expect_equal(res_222gs$ch_res$ind_stats$lag1, c(-0.3705498, -1.1441065, -0.8061067, 1.5904593), tolerance=1e-4)

  expect_equal(res_123$ac_res$ind_stats$lag1,
               c(0.35493268, 1.75181713, -0.31676151, 0.07667559, -2.97144849, -0.47202034, -0.26072843, -0.05403170, 0.47639484), tolerance=1e-4)
  expect_equal(res_123$ch_res$ind_stats$lag2,
               c(0.68604929, -2.36952048, -0.87931641, -1.33700310, -0.03940421, -0.68251561, -0.48329670, 1.29673759, -0.22832612), tolerance=1e-4)

  # SGSMVAR
  expect_equal(res_112s$ac_res$ind_stats$lag1, c(-2.20806502, 0.64582669, -0.04698889, -3.15105747), tolerance=1e-4)
  expect_equal(res_112s$ch_res$ind_stats$lag2, c(2.436755, 1.844220, 1.259395, 1.634222), tolerance=1e-4)

  expect_equal(res_222s$ac_res$ind_stats$lag1, c(-0.10349474, -1.21019700, -0.02516271, -1.59079039), tolerance=1e-4)
  expect_equal(res_222s$ch_res$ind_stats$lag2, c(-0.2116545, -0.2798438, 1.3557275, 1.5344440), tolerance=1e-4)
})


dim_g_norm <- 3*2 # 3*d, here d=2
g_norm <- function(r) { # "r" should be (T x d) quantile residual matrix
  d <- 2 ####
  T0 <- nrow(r)
  matrix((vapply(1:d, function(j) c(r[,j]^2 - 1, r[,j]^3, r[,j]^4 - 3), numeric(T0*3))), nrow=T0, ncol=dim_g_norm, byrow=FALSE)
}

# Function factory to produce function g for different lags
get_g <- function(lag) {
  d <- 2 ####
  function(r) {
    t(vapply((lag+1):nrow(r), function(t) vapply(1:lag, function(i1) tcrossprod(r[t,], r[t-i1,]), numeric(d^2)), numeric(lag*d^2)))
  }
} # Returns (T - lag x dim_g) matrix with values of g_t at each row, starting from t=lag+1 at the first row
g_ac1 <- get_g(1); dim_g_ac1 <- 1*2^2 #lag*d^2
g_ac2 <- get_g(2); dim_g_ac2 <- 2*2^2 #lag*d^2


test_that("get_test_Omega works correctly", {
  expect_equal(get_test_Omega(data=gdpdef, p=1, M=1, params=theta_112, model="GMVAR", conditional=TRUE, parametrization="mean",
                              constraints=NULL, same_means=NULL, g=g_norm, dim_g=dim_g_norm)[1,],
               c(1.1340811, 2.2652951, 15.5201108, 0.3014837, 0.6945552, 1.8175825), tolerance=1e-4)
  expect_equal(get_test_Omega(data=gdpdef, p=1, M=1, params=c(theta_112, 20), model="StMVAR", conditional=TRUE, parametrization="mean",
                              constraints=NULL, same_means=NULL, g=g_norm, dim_g=dim_g_norm)[6,],
               c(-0.6390766, -2.5941262, -20.5173315, 1.9402860, 5.7716302, 27.3629371), tolerance=1e-4)
  expect_equal(get_test_Omega(data=gdpdef, p=2, M=2, params=theta_222, model="GMVAR", conditional=TRUE, parametrization="intercept",
                              constraints=NULL, same_means=NULL, g=g_norm, dim_g=dim_g_norm)[,6],
               c(-0.8334590, -0.1859936, -5.9797771, 2.5482817, -0.9632214, 17.9082196), tolerance=1e-4)
  expect_equal(get_test_Omega(data=gdpdef, p=2, M=2, params=c(theta_222, 20, 30), model="StMVAR", conditional=TRUE, parametrization="intercept",
                              constraints=NULL, same_means=NULL, g=g_norm, dim_g=dim_g_norm)[3,],
               c(2.7171367, 1.9955770, 23.6659063, -0.7818911, 0.7877658, -2.7470121), tolerance=1e-4)
  expect_equal(get_test_Omega(data=gdpdef, p=2, M=c(1, 1), params=c(theta_222, 30), model="G-StMVAR", conditional=TRUE, parametrization="intercept",
                              constraints=NULL, same_means=NULL, g=g_norm, dim_g=dim_g_norm)[4,],
               c(-0.1684149, -0.4871592, -1.0006479, 1.0587743, -0.4385124, 5.5969645), tolerance=1e-4)
  expect_equal(get_test_Omega(data=gdpdef, p=2, M=2, params=theta_222c, model="GMVAR", conditional=TRUE, parametrization="intercept",
                              constraints=C_222c, same_means=NULL, g=g_norm, dim_g=dim_g_norm)[3,],
               c(3.9829341, 9.0409142, 47.4316793, -1.6580935, -0.6015382, -7.9351451), tolerance=1e-4)

  expect_equal(get_test_Omega(data=gdpdef, p=1, M=1, params=theta_112, model="GMVAR", conditional=TRUE, parametrization="mean",
                              constraints=NULL, same_means=NULL, g=g_ac1, dim_g=dim_g_ac1)[,1],
               c(0.193075103, 0.007521766, -0.175034038, -0.023736421), tolerance=1e-4)
  expect_equal(get_test_Omega(data=gdpdef, p=1, M=1, params=theta_112, model="GMVAR", conditional=TRUE, parametrization="mean",
                              constraints=NULL, same_means=NULL, g=g_ac2, dim_g=dim_g_ac2)[5,],
               c(-0.36497486, -0.03684782, 0.31026054, 0.15347275, 1.48304035, 0.26748799, -0.44855409, -0.26216707), tolerance=1e-4)

  expect_equal(get_test_Omega(data=gdpdef, p=2, M=2, params=theta_222, model="GMVAR", conditional=TRUE, parametrization="intercept",
                              constraints=NULL, same_means=NULL, g=g_ac1, dim_g=dim_g_ac1)[4,],
               c(0.056251016, -0.106900180, 0.006265508, 0.231235779), tolerance=1e-3)
  expect_equal(get_test_Omega(data=gdpdef, p=2, M=2, params=theta_222, model="GMVAR", conditional=TRUE, parametrization="intercept",
                              constraints=NULL, same_means=NULL, g=g_ac2, dim_g=dim_g_ac2)[8,],
               c(-0.01511927, 0.05053473, -0.02218766, 0.16288831, -0.01263294, -0.02870023, 0.10997098, 0.61942729), tolerance=1e-3)

  expect_equal(get_test_Omega(data=gdpdef, p=2, M=2, params=theta_222c, model="GMVAR", conditional=TRUE, parametrization="intercept",
                              constraints=C_222c, same_means=NULL, g=g_ac1, dim_g=dim_g_ac1)[2,],
               c(0.004225044, 0.829299155, -0.018035216, 0.010810051), tolerance=1e-4)
  expect_equal(get_test_Omega(data=gdpdef, p=2, M=2, params=theta_222c, model="GMVAR", conditional=TRUE, parametrization="intercept",
                              constraints=C_222c, same_means=NULL, g=g_ac2, dim_g=dim_g_ac2)[,2],
               c(0.007519279, 0.833487335, -0.015455803, 0.010043412, 0.079612914, 0.117168925, 0.095616855, 0.038757206), tolerance=1e-4)

  # SGSMVAR
  expect_equal(get_test_Omega(data=gdpdef, p=1, M=1, params=theta_112sWC, model="GMVAR", conditional=TRUE, parametrization="mean",
                              constraints=NULL, same_means=NULL, structural_pars=list(W=W_112), g=g_norm, dim_g=dim_g_norm)[2,],
               c(2.2652950, 30.1582430, 68.7565902, -0.1180411, 0.2864592, 0.4527952), tolerance=1e-4)
  expect_equal(get_test_Omega(data=gdpdef, p=1, M=1, params=theta_112sWC, model="GMVAR", conditional=TRUE, parametrization="mean",
                              constraints=NULL, same_means=NULL, structural_pars=list(W=W_112), g=g_ac1, dim_g=dim_g_ac1)[4,],
               c(-0.02373642, -0.07959005, -0.15221360, 1.29820691), tolerance=1e-4)
  expect_equal(get_test_Omega(data=gdpdef, p=2, M=2, params=theta_222s, model="GMVAR", conditional=TRUE, parametrization="intercept",
                              constraints=NULL, same_means=NULL, structural_pars=list(W=W_222), g=g_ac2, dim_g=dim_g_ac2)[,8],
               c(-0.01511927, 0.05053473, -0.02218765, 0.16288831, -0.01263294, -0.02870023, 0.10997098, 0.61942729), tolerance=1e-4)
})

