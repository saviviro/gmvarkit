context("quantile residual tests")
library(gmvarkit)

# NOTE these tests use random elements obtained from simulation algorithms

data <- cbind(10*eurusd[,1], 100*eurusd[,2])

## A(M)(p)_(p)(M)(d)

# p=1, M=1, d=2, parametrization="mean"
phi10_112 <- c(1.07, 127.71)
A11_112 <- matrix(c(0.99, 0.00, -0.01, 0.99), nrow=2)
Omega1_112 <- matrix(c(4.05, 2.22, 2.22, 2.22), nrow=2)
theta_112 <- c(phi10_112, vec(A11_112), vech(Omega1_112))
mod_112 <- GMVAR(data, p=1, M=1, d=2, params=theta_112, conditional=TRUE, parametrization="mean", constraints=NULL)

W_112 <- t(chol(Omega1_112))
theta_112sWC <- c(phi10_112, vec(A11_112), Wvec(W_112)) # SGMVAR, W constrained by Cholesky
mod_112s <- GMVAR(data, p=1, M=1, d=2, params=theta_112sWC, conditional=TRUE, parametrization="mean", constraints=NULL,
                  structural_pars=list(W=W_112))


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

WL_222 <- diag_Omegas(Omega1_222, Omega2_222)
W_222 <- matrix(WL_222[1:(2^2)], nrow=2, byrow=FALSE)
lambdas_222 <- WL_222[(2^2 + 1):length(WL_222)]
theta_222s <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A21_222),
                vec(A22_222), vec(W_222), lambdas_222, alpha1_222) # SGMVAR
mod_222s <- GMVAR(data, p=2, M=2, d=2, params=theta_222s, conditional=TRUE, parametrization="intercept", constraints=NULL,
                  structural_pars=list(W=W_222))


# p=2, M=2, d=2, AR paramars same, non-diagonals zero, intercept
theta_222c <- c(0.3552775, 3.1929675, -0.1143198, 2.8294743, 1.2633425, 1.3375150, -0.2919742, -0.3624010,
                5.5971764, 3.4559442, 9.6221422, 0.9820759, -0.3267521, 5.2358855, 0.6501600)
mat0 <- matrix(c(1, rep(0, 10), 1, rep(0, 8), 1, rep(0, 10), 1), nrow=2*2^2, byrow=FALSE)
C_222c <- rbind(mat0, mat0)
mod_222c <- GMVAR(data, p=2, M=2, d=2, params=theta_222c, conditional=TRUE, parametrization="intercept", constraints=C_222c)

# p=1, M=2, d=3, no constraints, rand_ind and simulated data
set.seed(11)
theta_123 <- random_ind2(p=1, M=2, d=3, mu_scale=c(-10, 0, 5), mu_scale2=1:3, omega_scale=1:3, ar_scale=1)
mod_123 <- GMVAR(p=1, M=2, d=3, params=theta_123, conditional=FALSE, parametrization="mean", constraints=NULL)
sim_123 <- simulateGMVAR(mod_123, nsimu=500)
data_123 <- sim_123$sample
mod_123 <- GMVAR(data_123, p=1, M=2, d=3, params=theta_123, conditional=FALSE, parametrization="mean", constraints=NULL)

# p=2, M=2, d=2, parametrization="mean", constraints=C_mat, same_means=list(1:2)
C_mat <- rbind(diag(2*2^2), diag(2*2^2))
params_222cm <- c(-9.12048853805255, 123.142757183508, 1.2658425363326, 0.0675545389606989, 0.0331264235657607,
                  1.33370494344656, -0.285882557831441, -0.0769144929653558, -0.0382772162867802, -0.351635998882842,
                  5.8625623309659, 3.57488618757834, 9.70846346569286, 0.869261580580846, -0.248703862116217,
                  5.17613656742281, 0.439575388572472)
mod_222cm <- GMVAR(data, p=2, M=2, params=params_222cm, parametrization="mean", constraints=C_mat, same_means=list(1:2))

set.seed(1)
res_112 <- quantile_residual_tests(mod_112, lags_ac=1:2, lags_ch=1:2, nsimu=300, print_res=FALSE)
res_222 <- quantile_residual_tests(mod_222, lags_ac=3, lags_ch=1, nsimu=1, print_res=FALSE)
res_123 <- quantile_residual_tests(mod_123, lags_ac=1, lags_ch=2, nsimu=1, print_res=FALSE)

set.seed(1); res_112s <- quantile_residual_tests(mod_112s, lags_ac=1:2, lags_ch=1:2, nsimu=1, print_res=FALSE)
set.seed(1); res_222s <- quantile_residual_tests(mod_222s, lags_ac=1, lags_ch=2, nsimu=300, print_res=FALSE)

set.seed(1); res_222cm <- quantile_residual_tests(mod_222cm, lags_ac=2, lags_ch=1, nsimu=1, print_res=FALSE)

test_that("quantile_residual_tests - test_results - works correctly", {
  expect_equal(res_112$norm_res$test_stat, 403037.1, tolerance=1)
  expect_equal(res_112$ac_res$test_results$test_stat, c(1009.7560, 990.2116), tolerance=1e-3)
  expect_equal(res_112$ch_res$test_results$test_stat, c(281816.3, 695870.9), tolerance=0.1)

  expect_equal(res_222$norm_res$p_val, 0.9026294, tolerance=1e-4)
  expect_equal(res_222$ac_res$test_results$p_val, 0.3135982, tolerance=1e-4)
  expect_equal(res_222$ch_res$test_results$test_stat, 1.767297, tolerance=1e-4)

  expect_equal(res_123$norm_res$test_stat, 22.50893, tolerance=1e-4)
  expect_equal(res_123$ac_res$test_results$test_stat, 13.70742, tolerance=1e-4)
  expect_equal(res_123$ch_res$test_results$test_stat, 17.73436, tolerance=1e-4)

  # SGMVAR
  expect_equal(res_112s$norm_res$test_stat, 115.5537, tolerance=1e-3)
  expect_equal(res_112s$ac_res$test_results$test_stat, c(32.74336, 35.58885), tolerance=1e-3)
  expect_equal(res_112s$ch_res$test_results$test_stat, c(16.27407, 22.02809), tolerance=1e-3)

  expect_equal(res_222s$norm_res$p_val, 0.6104193, tolerance=1e-4)
  expect_equal(res_222s$ac_res$test_results$p_val, 0.8085959, tolerance=1e-4)
  expect_equal(res_222s$ch_res$test_results$p_val, 0.6376873, tolerance=1e-4)

  # Same means
  expect_equal(res_222cm$norm_res$p_val, 0.02890172, tolerance=1e-4)
  expect_equal(res_222cm$ac_res$test_results$p_val, 0.2861438, tolerance=1e-4)
  expect_equal(res_222cm$ch_res$test_results$p_val, 0.9103965, tolerance=1e-4)
})

test_that("quantile_residual_tests - ind_stats - works correctly", {
  expect_equal(res_112$ac_res$ind_stats$lag2, c(-0.04276715, -0.80618069, -0.17001551, -0.05597274), tolerance=1e-4)
  expect_equal(res_112$ch_res$ind_stats$lag1, c(2.744658, 3.532080, 27.329656, 525.308165), tolerance=1e-4)

  expect_equal(res_222$ac_res$ind_stats$lag3, c(2.1920450, 0.3863022, -0.1567286, 1.2905688), tolerance=1e-4)
  expect_equal(res_222$ch_res$ind_stats$lag1, c(0.09774581, 0.02637839, -0.93709782, 0.75117796), tolerance=1e-4)

  expect_equal(res_123$ac_res$ind_stats$lag1,
               c(0.35493268, 1.75181713, -0.31676151, 0.07667559, -2.97144849, -0.47202034, -0.26072843, -0.05403170, 0.47639484), tolerance=1e-4)
  expect_equal(res_123$ch_res$ind_stats$lag2,
               c(0.68604929, -2.36952048, -0.87931641, -1.33700310, -0.03940421, -0.68251561, -0.48329670, 1.29673759, -0.22832612), tolerance=1e-4)

  # SGMVAR
  expect_equal(res_112s$ac_res$ind_stats$lag1, c(4.10238061, 0.95255874, 0.03488986, 3.79132960), tolerance=1e-4)
  expect_equal(res_112s$ch_res$ind_stats$lag2, c(1.3283931, -0.8333597, 0.5940387, 3.4563746), tolerance=1e-4)

  expect_equal(res_222s$ac_res$ind_stats$lag1, c(-0.1788567, 0.4186823, -0.6197361, 0.9011829), tolerance=1e-4)
  expect_equal(res_222s$ch_res$ind_stats$lag2, c(0.2758057, -0.1911319, 1.7897447, 1.0156843), tolerance=1e-4)
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
  expect_equal(get_test_Omega(data=data, p=1, M=1, params=theta_112, conditional=TRUE, parametrization="mean",
                              constraints=NULL, same_means=NULL, g=g_norm, dim_g=dim_g_norm)[1,],
               c(0.7057628, 2.2848667, 8.6067786, 0.9472146, -7.4855476, 42.0832749), tolerance=1e-4)
  expect_equal(get_test_Omega(data=data, p=2, M=2, params=theta_222, conditional=TRUE, parametrization="intercept",
                              constraints=NULL, same_means=NULL, g=g_norm, dim_g=dim_g_norm)[,6],
               c(-0.2246738, 0.2367933, -1.3629537, 2.5325211, -7.5225937, 38.2535773), tolerance=1e-4)
  expect_equal(get_test_Omega(data=data, p=2, M=2, params=theta_222c, conditional=TRUE, parametrization="intercept",
                              constraints=C_222c, same_means=NULL, g=g_norm, dim_g=dim_g_norm)[3,],
               c(1.5549257, 4.5538484, 28.3233066, 0.4093669, -0.6562638, 0.4822259), tolerance=1e-4)

  expect_equal(get_test_Omega(data=data, p=1, M=1, params=theta_112, conditional=TRUE, parametrization="mean",
                              constraints=NULL, same_means=NULL, g=g_ac1, dim_g=dim_g_ac1)[,1],
               c(1.0143162, -0.2785037, -0.1520407, 0.4920706), tolerance=1e-4)
  expect_equal(get_test_Omega(data=data, p=1, M=1, params=theta_112, conditional=TRUE, parametrization="mean",
                              constraints=NULL, same_means=NULL, g=g_ac2, dim_g=dim_g_ac2)[5,],
               c(0.26554361, -0.04612983, 0.41599389, 0.84506916, 1.18335341, 0.02087973, 0.23315441, 1.08043135), tolerance=1e-4)

  expect_equal(get_test_Omega(data=data, p=2, M=2, params=theta_222, conditional=TRUE, parametrization="intercept",
                              constraints=NULL, same_means=NULL, g=g_ac1, dim_g=dim_g_ac1)[4,],
               c(0.06105453, 0.03487609, 0.03182064, 0.30107295), tolerance=1e-3)
  expect_equal(get_test_Omega(data=data, p=2, M=2, params=theta_222, conditional=TRUE, parametrization="intercept",
                              constraints=NULL, same_means=NULL, g=g_ac2, dim_g=dim_g_ac2)[8,],
               c(-0.03558556, -0.07470666, -0.07362958, -0.27822363, 0.17240094, -0.01494992, 0.11873980, 1.05340983), tolerance=1e-3)

  expect_equal(get_test_Omega(data=data, p=2, M=2, params=theta_222c, conditional=TRUE, parametrization="intercept",
                              constraints=C_222c, same_means=NULL, g=g_ac1, dim_g=dim_g_ac1)[2,],
               c(0.19084136, 0.85543945, 0.09623485, -0.18203186), tolerance=1e-4)
  expect_equal(get_test_Omega(data=data, p=2, M=2, params=theta_222c, conditional=TRUE, parametrization="intercept",
                              constraints=C_222c, same_means=NULL, g=g_ac2, dim_g=dim_g_ac2)[,2],
               c(0.19172737, 0.85825300, 0.09773145, -0.18278803, 0.03803171, -0.07227835, -0.07058063, -0.08615313), tolerance=1e-4)

  # SGMVAR
  expect_equal(get_test_Omega(data=data, p=1, M=1, params=theta_112sWC, conditional=TRUE, parametrization="mean",
                              constraints=NULL, same_means=NULL, structural_pars=list(W=W_112), g=g_norm, dim_g=dim_g_norm)[2,],
               c(2.284867, 21.283539, 51.544341, 1.723343, 15.281678, -16.989623), tolerance=1e-4)
  expect_equal(get_test_Omega(data=data, p=1, M=1, params=theta_112sWC, conditional=TRUE, parametrization="mean",
                              constraints=NULL, same_means=NULL, structural_pars=list(W=W_112), g=g_ac1, dim_g=dim_g_ac1)[4,],
               c(0.4927337, -3.5265900, 0.9979045, 94.4919253), tolerance=1e-4)
  expect_equal(get_test_Omega(data=data, p=2, M=2, params=theta_222s, conditional=TRUE, parametrization="intercept",
                              constraints=NULL, same_means=NULL, structural_pars=list(W=W_222), g=g_ac2, dim_g=dim_g_ac2)[,8],
               c(-0.03558556, -0.07470666, -0.07362959, -0.27822363, 0.17240094, -0.01494992, 0.11873980, 1.05340983), tolerance=1e-4)
})

