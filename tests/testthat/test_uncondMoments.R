context("unconditional moments")
library(gmvarkit)


## A(M)(p)_(p)(M)(d)

# p=1, M=1, d=2
phi10_112 <- c(1.03, 2.36)
A11_112 <- matrix(c(0.25, 0.06, 0.04, 0.34), nrow=2, byrow=FALSE)
Omega1_112 <- matrix(c(0.93, -0.15, -0.15, 5.20), nrow=2, byrow=FALSE)

theta_112 <- upsilon1_112 <- c(phi10_112, vec(A11_112), vech(Omega1_112))

W_112 <- t(chol(Omega1_112))
theta_112sWC <- c(phi10_112, vec(A11_112), Wvec(W_112)) # SGMVAR, W constrained by Cholesky

theta_112t <- c(theta_112, 10)
theta_112tsWC <- c(theta_112sWC, 10)

# p=1, M=2, d=2
phi10_122 <- c(1.03, 2.36)
A11_122 <- matrix(c(0.1, -0.06, -0.04, 0.1), nrow=2, byrow=FALSE)
Omega1_122 <- matrix(c(0.93, -0.15, -0.15, 5.20), nrow=2, byrow=FALSE)

phi20_122 <- c(1.79, 3.00)
A21_122 <- A11_122
Omega2_122 <- matrix(c(5.88, 3.56, 3.56, 9.80), nrow=2, byrow=FALSE)

alpha1_122 <- 0.37
upsilon1_122 <- c(phi10_122, vec(A11_122), vech(Omega1_122))
upsilon2_122 <- c(phi20_122, vec(A21_122), vech(Omega2_122))
theta_122 <- c(upsilon1_122, upsilon2_122, alpha1_122)

WL_122 <- diag_Omegas(Omega1_122, Omega2_122)
W_122 <- matrix(WL_122[1:(2^2)], nrow=2, byrow=FALSE)
lambdas_122 <- WL_122[(2^2 + 1):length(WL_122)]
theta_122s <- c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), vec(W_122), lambdas_122, alpha1_122) # SGMVAR

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
W_222 <- matrix(WL_222[1:(2^2)], nrow=2, byrow=FALSE)
lambdas_222 <- WL_222[(2^2 + 1):length(WL_222)]
theta_222s <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A21_222),
                vec(A22_222), vec(W_222), lambdas_222, alpha1_222) # SGMVAR

theta_222t <- c(theta_222, 10, 20) # StMVAR
theta_222gs <- c(theta_222, 20) # G-StMVAR, M1=1, M2=1
theta_222ts <- c(theta_222s, 10, 20) # SStMVAR


# p=3, M=3, d=2
phi10_332 <- c(1.03, 2.36)
A11_332 <- matrix(c(0.25, 0.06, 0.04, 0.34), nrow=2, byrow=FALSE)
A12_332 <- matrix(c(-0.29, -0.08, -0.05, -0.36), nrow=2, byrow=FALSE)
A13_332 <- matrix(c(0.1, 0.21, 0.12, 0.2), nrow=2, byrow=FALSE)
Omega1_332 <- matrix(c(0.93, -0.15, -0.15, 5.20), nrow=2, byrow=FALSE)

phi20_332 <- c(1.79, 3.00)
A21_332 <- matrix(c(0.20, 0.05, 0.05, 0.30), nrow=2, byrow=FALSE)
A22_332 <- matrix(c(-0.30, -0.10, -0.05, -0.40), nrow=2, byrow=FALSE)
A23_332 <- matrix(c(0.30, 0.10, 0.05, 0.40), nrow=2, byrow=FALSE)
Omega2_332 <- matrix(c(5.88, 3.56, 3.56, 9.80), nrow=2, byrow=FALSE)

phi30_332 <- c(1.79, 3.00)
A31_332 <- matrix(c(0.30, 0.03, 0.08, 0.33), nrow=2, byrow=FALSE)
A32_332 <- matrix(c(-0.50, -0.20, -0.01, -0.40), nrow=2, byrow=FALSE)
A33_332 <- matrix(c(0.50, 0.20, 0.01, 0.40), nrow=2, byrow=FALSE)
Omega3_332 <- matrix(c(5.00, 3.00, 3.00, 9.00), nrow=2, byrow=FALSE)

alpha1_332 <- 0.5
alpha2_332 <- 0.3
upsilon1_332 <- c(phi10_332, vec(A11_332), vec(A12_332), vec(A13_332), vech(Omega1_332))
upsilon2_332 <- c(phi20_332, vec(A21_332), vec(A22_332), vec(A23_332), vech(Omega2_332))
upsilon3_332 <- c(phi30_332, vec(A31_332), vec(A32_332), vec(A33_332), vech(Omega3_332))
theta_332 <- c(upsilon1_332, upsilon2_332, upsilon3_332, alpha1_332, alpha2_332)

W_332 <- matrix(c(-0.8924620, 0, 0.3653923, -2.1643472), nrow=2, byrow=FALSE)
lambdas2_332 <- c(7.16, 1.30)
lambdas3_332 <- c(6.10, 1.22)
theta_332s <- c(phi10_332, phi20_332, phi30_332, vec(A11_332), vec(A12_332), vec(A13_332),
                vec(A21_332), vec(A22_332), vec(A23_332), vec(A31_332), vec(A32_332), vec(A33_332),
                vec(W_332), lambdas2_332, lambdas3_332, alpha1_332, alpha2_332) # SGMVAR
Omega1_332s <- tcrossprod(W_332)
Omega2_332s <- W_332%*%tcrossprod(diag(lambdas2_332), W_332)
Omega3_332s <- W_332%*%tcrossprod(diag(lambdas3_332), W_332)
theta_332_froms <-  c(phi10_332, vec(A11_332), vec(A12_332), vec(A13_332), vech(Omega1_332s), # based on SGMVAR
                      phi20_332, vec(A21_332), vec(A22_332), vec(A23_332), vech(Omega2_332s),
                      phi30_332, vec(A31_332), vec(A32_332), vec(A33_332), vech(Omega3_332s),
                      alpha1_332, alpha2_332)

theta_332sWC <- c(phi10_332, phi20_332, phi30_332, vec(A11_332), vec(A12_332), vec(A13_332),
                  vec(A21_332), vec(A22_332), vec(A23_332), vec(A31_332), vec(A32_332), vec(A33_332),
                  Wvec(W_332), lambdas2_332, lambdas3_332, alpha1_332, alpha2_332) # SGMVAR W constrained

theta_332gs <- c(theta_332, 20, 30) # G-StMVAR, M1=1, M2=2
theta_332gssWC <- c(theta_332sWC, 30) # SG-StMVAR, M1=2, M2=1

# p=1, M=2, d=3
phi10_123 <- c(1.1, 2.2, 3.3)
A11_123 <- matrix(c(0.1, 0.21, 0.31, 0.12, 0.2, 0.32, 0.13, 0.23, 0.3), nrow=3, byrow=FALSE)
Omega1_123 <- matrix(c(1, 0.22, 0.33, 0.22, 2, 0.44, 0.33, 0.44, 3), nrow=3, byrow=FALSE)

phi20_123 <- c(1.11, 2.22, 3.33)
A21_123 <- matrix(c(-0.1, -0.21, -0.31, -0.12, -0.2, -0.32, -0.13, -0.23, -0.3), nrow=3, byrow=FALSE)
Omega2_123 <- matrix(c(1.1, 0.222, 0.333, 0.222, 2.2, 0.444, 0.333, 0.444, 3.3), nrow=3, byrow=FALSE)

alpha1_123 <- 0.6
upsilon1_123 <- c(phi10_123, vec(A11_123), vech(Omega1_123))
upsilon2_123 <- c(phi20_123, vec(A21_123), vech(Omega2_123))
theta_123 <- c(upsilon1_123, upsilon2_123, alpha1_123)

WL_123 <- diag_Omegas(Omega1_123, Omega2_123)
W_123 <- matrix(WL_123[1:(3^2)], nrow=3, byrow=FALSE)
lambdas_123 <- WL_123[(3^2 + 1):length(WL_123)]
theta_123s <- c(phi10_123, phi20_123, vec(A11_123), vec(A21_123), vec(W_123), lambdas_123, alpha1_123) # SGMVAR

theta_123t <- c(theta_123, 10, 20) # StMVAR
theta_123gss <- c(theta_123s, 20) # SG-StMVAR

# p=2, M=1, d=3
phi10_213 <- c(1.1, 2.2, 3.3)
A11_213 <- matrix(c(0.1, 0.21, 0.31, 0.12, 0.2, 0.32, 0.13, 0.23, 0.3), nrow=3, byrow=FALSE)
A12_213 <- matrix(c(-0.1, -0.21, -0.31, -0.12, -0.2, -0.32, -0.13, -0.23, -0.3), nrow=3, byrow=FALSE)
Omega1_213 <- matrix(c(1, 0.22, 0.33, 0.22, 2, 0.44, 0.33, 0.44, 3), nrow=3, byrow=FALSE)

upsilon1_213 <- c(phi10_213, vec(A11_213), vec(A12_213), vech(Omega1_213))
theta_213 <- upsilon1_213

W_213 <- t(chol(Omega1_213))
theta_213s <- c(phi10_213, vec(A11_213), vec(A12_213), vec(W_213)) # SGMVAR

theta_213sWC <- c(phi10_213, vec(A11_213), vec(A12_213), Wvec(W_213)) # SGMVAR W constrained

theta_213t <- c(theta_213, 10) # StMVAR
theta_213tsWC <- c(theta_213sWC, 10) # SStMVAR

calc_mu <- function(p, M, d, params, model=c("GMVAR", "StMVAR", "G-StMVAR"), constraints=NULL, structural_pars=NULL) {
  model <- match.arg(model)
  params <- reform_constrained_pars(p=p, M=M, d=d, params=params, model=model, constraints=constraints, structural_pars=structural_pars)
  all_A <- pick_allA(p=p, M=M, d=d, params=params, structural_pars=structural_pars)
  all_phi0 <- pick_phi0(p=p, M=M, d=d, params=params, structural_pars=structural_pars)
  vapply(1:sum(M), function(m) solve(diag(d) - rowSums(all_A[, , , m, drop=FALSE], dims=2), all_phi0[,m]), numeric(d))
}

theta_112_mu <- change_parametrization(p=1, M=1, d=2, params=theta_112, change_to="mean")
theta_112t_mu <- change_parametrization(p=1, M=1, d=2, params=theta_112t, model="StMVAR", change_to="mean")
theta_122_mu <- change_parametrization(p=1, M=2, d=2, params=theta_122, change_to="mean")
theta_222_mu <- change_parametrization(p=2, M=2, d=2, params=theta_222, change_to="mean")
theta_222t_mu <- change_parametrization(p=2, M=2, d=2, params=theta_222t, model="StMVAR", change_to="mean")
theta_222gs_mu <- change_parametrization(p=2, M=c(1, 1), d=2, params=theta_222t, model="G-StMVAR", change_to="mean")
theta_332_mu <- change_parametrization(p=3, M=3, d=2, params=theta_332, change_to="mean")
theta_332gs_mu <- change_parametrization(p=3, M=c(1, 2), d=2, params=theta_332gs, model="G-StMVAR", change_to="mean")
theta_123_mu <- change_parametrization(p=1, M=2, d=3, params=theta_123, change_to="mean")
theta_123t_mu <- change_parametrization(p=1, M=2, d=3, params=theta_123t, model="StMVAR", change_to="mean")
theta_213_mu <- change_parametrization(p=2, M=1, d=3, params=theta_213, change_to="mean")
theta_213t_mu <- change_parametrization(p=2, M=1, d=3, params=theta_213t, model="StMVAR", change_to="mean")

theta_112sWC_mu <- change_parametrization(p=1, M=1, d=2, params=theta_112sWC, structural_pars=list(W=W_112), change_to="mean")
theta_112tsWC_mu <- change_parametrization(p=1, M=1, d=2, params=theta_112tsWC, model="StMVAR", structural_pars=list(W=W_112), change_to="mean")
theta_222ts_mu <- change_parametrization(p=2, M=2, d=2, params=theta_222ts, model="StMVAR", structural_pars=list(W=W_222), change_to="mean")



## A(M)(p)_(p)(M)(d)
rbind_diags <- function(p, M, d) {
  I <- diag(p*d^2)
  Reduce(rbind, replicate(M, I, simplify=FALSE))
}

# Constraining AR-parameters to be the same for all regimes

C_122 <- rbind_diags(p=1, M=2, d=2)

# p=1, M=1, d=2
C_112 <- rbind_diags(p=1, M=1, d=2)
theta_112c <- c(phi10_112, vec(A11_112), vech(Omega1_112))

theta_112csWAR <- c(phi10_112, vec(A11_112), Wvec(W_112)) # SGMVAR W and AR

theta_112tc <- c(theta_112c, 10) # StMVAR

theta_112tcsWAR <- c(theta_112csWAR, 10) # StMVAR

# p=2, M=2, d=2
C_222 <- rbind_diags(p=2, M=2, d=2)
theta_222c <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vech(Omega1_222), vech(Omega2_222), alpha1_222)

C_lambda_222 <- matrix(c(1, 2), nrow=2)
theta_222csL <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A21_222),
                  vec(A22_222), vec(W_222), 0.2, alpha1_222) # SGMVAR lambdas
theta_222csLAR <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(W_222), 0.2, alpha1_222) # SGMVAR lambdas and AR

theta_222gscsL <- c(theta_222csL, 20) # SG-StMVAR, M1=1, M2=1
#theta_222gscsL_expanded <- c(theta_222csL_expanded, 20) # SG-StMVAR, M1=1, M2=1

theta_222tcsLAR <- c(theta_222csLAR, 10, 20) # SStMVAR
#theta_222tcsLAR_expanded <- c(theta_222csLAR_expanded, 10, 20) # SStMVAR

# p=1, M=2, d=3
C_123 <- rbind_diags(p=1, M=2, d=3)
theta_123c <- c(phi10_123, phi20_123, vec(A11_123), vech(Omega1_123), vech(Omega2_123), alpha1_123)

C_lambda_123 <- matrix(c(1, 1, 0, 0, 0, 1), nrow=3, byrow=FALSE)
theta_123csL <- c(phi10_123, phi20_123, vec(A11_123), vec(A21_123), vec(W_123), 1, 2, alpha1_123) # SGMVAR lambdas
theta_123csLAR <- c(phi10_123, phi20_123, vec(A11_123), vec(W_123), 1, 2, alpha1_123) # SGMVAR lambdas and AR

theta_123tc <- c(theta_123c, 10, 20) # StMVAR
#theta_123tc_expanded <- c(theta_123c_expanded, 10, 20) # StMVAR

theta_123tcsL <- c(theta_123csL, 10, 20) # StMVAR
#theta_123tcL_expanded <- c(theta_123csL_expanded, 10, 20) # StMVAR


# p=2, M=2, d=2, constraint AR-parameters to be the same for all regimes
# and constraint the of-diagonal elements of AR-matrices to be zero.
mat0 <- matrix(c(1, rep(0, 10), 1, rep(0, 8), 1, rep(0, 10), 1), nrow=2*2^2, byrow=FALSE)
C_222_2 <- rbind(mat0, mat0)
A21_222_c2 <- A11_222_c2 <- matrix(c(0.26, 0, 0, 0.34), nrow=2, byrow=FALSE)
A22_222_c2 <- A12_222_c2 <- matrix(c(-0.29, 0, 0, -0.36), nrow=2, byrow=FALSE)
phi10_222_c2 <- c(-0.11, 2.83)
phi20_222_c2 <- c(0.36, 3.19)
Omega1_222_c2 <- matrix(c(0.98, -0.33, -0.33, 5.24), nrow=2, byrow=FALSE)
Omega2_222_c2 <- matrix(c(5.60, 3.46, 3.46, 9.62), nrow=2, byrow=FALSE)
alpha1_222_c2 <- 0.35
theta_222_c2 <- c(phi10_222_c2, phi20_222_c2, 0.26, 0.34, -0.29, -0.36, vech(Omega1_222_c2),
                  vech(Omega2_222_c2), alpha1_222_c2)

WL_222c2 <- diag_Omegas(Omega1_222_c2, Omega2_222_c2)
W_222c2 <- matrix(WL_222c2[1:(2^2)], nrow=2, byrow=FALSE)
lambdas_222c2 <- WL_222c2[(2^2 + 1):length(WL_222c2)]
theta_222_c2s <- c(phi10_222_c2, phi20_222_c2, 1.26, 1.34, -0.29, -0.36, vec(W_222c2), lambdas_222c2, alpha1_222_c2) # SGMVAR AR

theta_112c_mu <- change_parametrization(p=1, M=1, d=2, params=theta_112c, constraints=C_112, change_to="mean")
theta_112tc_mu <- change_parametrization(p=1, M=1, d=2, params=theta_112tc, model="StMVAR", constraints=C_112, change_to="mean")

theta_222c_mu <- change_parametrization(p=2, M=2, d=2, params=theta_222c, constraints=C_222, change_to="mean")
theta_222c_mu2 <- change_parametrization(p=2, M=2, d=2, params=theta_222_c2, constraints=C_222_2, change_to="mean")
theta_123c_mu <- change_parametrization(p=1, M=2, d=3, params=theta_123c, constraints=C_123, change_to="mean")
theta_123tc_mu <- change_parametrization(p=1, M=2, d=3, params=theta_123tc, model="StMVAR", constraints=C_123, change_to="mean")

theta_112c_mu_exp <- reform_constrained_pars(p=1, M=1, d=2, params=theta_112c_mu, constraints=C_112)
theta_112tc_mu_exp <- reform_constrained_pars(p=1, M=1, d=2, params=theta_112tc_mu, model="StMVAR", constraints=C_112)

theta_222c_mu_exp <- reform_constrained_pars(p=2, M=2, d=2, params=theta_222c_mu, constraints=C_222)
theta_222c_mu2_exp <- reform_constrained_pars(p=2, M=2, d=2, params=theta_222c_mu2, constraints=C_222_2)
theta_123c_mu_exp <- reform_constrained_pars(p=1, M=2, d=3, params=theta_123c_mu, constraints=C_123)
theta_123tc_mu_exp <- reform_constrained_pars(p=1, M=2, d=3, params=theta_123tc_mu, model="StMVAR", constraints=C_123)

theta_222s_mu <- change_parametrization(p=2, M=2, d=2, params=theta_222s, structural_pars=list(W=W_222), change_to="mean")

theta_112cs_mu <- change_parametrization(p=1, M=1, d=2, params=theta_112csWAR, constraints=C_112,
                                         structural_pars=list(W=W_112), change_to="mean")
theta_222cs_mu <- change_parametrization(p=2, M=2, d=2, params=theta_222csLAR, constraints=C_222,
                                         structural_pars=list(W=W_222, C_lambda=C_lambda_222), change_to="mean")
theta_112cs_mu_exp <- reform_constrained_pars(p=1, M=1, d=2, params=theta_112cs_mu, constraints=C_112, structural_pars=list(W=W_112))
theta_222cs_mu_exp <- reform_constrained_pars(p=2, M=2, d=2, params=theta_222cs_mu, constraints=C_222,
                                              structural_pars=list(W=W_222, C_lambda_222))

# p=2, M=2, d=2, constraints=C_222, same_means=list(1:2)
theta_222c_int <- c(phi10_222, vec(A11_222), vec(A12_222), vech(Omega1_222), vech(Omega2_222), alpha1_222)
theta_222c_int_expanded <- c(phi10_222, vec(A11_222), vec(A12_222), vech(Omega1_222), phi10_222, vec(A11_222), vec(A12_222),
                             vech(Omega2_222), alpha1_222)

theta_222gsc_int <- c(theta_222c_int, 20) # G-StMVAR, M1=1, M2=1
theta_222gsc_int_expanded <- c(theta_222c_int_expanded, 20) # G-StMVAR, M1=1, M2=1

# p=2, M=2, d=2, constraints=C_222, structural_pars=list(W=W_222, C_lambda=C_lambda_222), same_means=list(1:2)
theta_222csLAR_int <- c(phi10_222, vec(A11_222), vec(A12_222), vec(W_222), 0.2, alpha1_222)
theta_222csLAR_int_expanded <-  c(phi10_222, phi10_222, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222),
                                  vec(W_222), 0.2, 2*0.2, alpha1_222)

theta_222tcsLAR_int <- c(theta_222csLAR_int, 10, 20) # SStMVAR
theta_222tcsLAR_int_expanded <- c(theta_222csLAR_int_expanded, 10, 20) # SStMVAR

theta_222gscsLAR_int <- c(theta_222csLAR_int, 20) # SG-StMVAR, M1=1, M2=1
theta_222gscsLAR_int_expanded <- c(theta_222csLAR_int_expanded, 20) # SG-StMVAR, M1=1, M2=1

## Fixed alphas and lambdas

# p=1, M=2, d=2, model="GMVAR", constraints=C_122, weight_constraints=0.6, structural_pars=list(W=W_122, fixed_lambdas=c(4, 3))
theta_122cwsF <- c(phi10_122, phi20_122, vec(A11_122), vec(W_122))

# p=2, M=2, d=2, model="StMVAR", constraints=C_222, weight_constraints=0.6, structural_pars=list(W=W_222, fixed_lambdas=c(4, 3))
theta_222tcwsF <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(W_222), 11, 12)

# p=2, M=c(1, 1), d=2, model="G-StMVAR", constraints=C_222, same_means=list(1:2), weight_constraints=0.6,
# structural_pars=list(W=W_222, fixed_lambdas=c(4, 3))
theta_222gscmwsF <- c(phi10_222, vec(A11_222), vec(A12_222), vec(W_222), 11)

test_that("get_regime_means_int works correctly", {
  expect_equal(get_regime_means_int(p=1, M=1, d=2, params=theta_112, parametrization="intercept", constraints=NULL),
               pick_phi0(p=1, M=1, d=2, params=theta_112_mu))
  expect_equal(get_regime_means_int(p=1, M=1, d=2, params=theta_112_mu, parametrization="mean", constraints=NULL),
               calc_mu(p=1, M=1, d=2, params=theta_112))
  expect_equal(get_regime_means_int(p=1, M=1, d=2, params=theta_112t, model="StMVAR", parametrization="intercept", constraints=NULL),
               pick_phi0(p=1, M=1, d=2, params=theta_112t_mu))
  expect_equal(get_regime_means_int(p=1, M=1, d=2, params=theta_112t_mu, model="StMVAR", parametrization="mean", constraints=NULL),
               calc_mu(p=1, M=1, d=2, params=theta_112t, model="StMVAR"))
  expect_equal(get_regime_means_int(p=1, M=2, d=2, params=theta_122, parametrization="intercept", constraints=NULL),
               pick_phi0(p=1, M=2, d=2, params=theta_122_mu))
  expect_equal(get_regime_means_int(p=1, M=2, d=2, params=theta_122_mu, parametrization="mean", constraints=NULL),
              calc_mu(p=1, M=2, d=2, params=theta_122))
  expect_equal(get_regime_means_int(p=2, M=2, d=2, params=theta_222, parametrization="intercept", constraints=NULL),
               pick_phi0(p=2, M=2, d=2, params=theta_222_mu))
  expect_equal(get_regime_means_int(p=2, M=2, d=2, params=theta_222_mu, parametrization="mean", constraints=NULL),
               calc_mu(p=2, M=2, d=2, params=theta_222))
  expect_equal(get_regime_means_int(p=2, M=2, d=2, params=theta_222t, model="StMVAR", parametrization="intercept", constraints=NULL),
               pick_phi0(p=2, M=2, d=2, params=theta_222t_mu))
  expect_equal(get_regime_means_int(p=2, M=2, d=2, params=theta_222t_mu, model="StMVAR", parametrization="mean", constraints=NULL),
               calc_mu(p=2, M=2, d=2, params=theta_222t, model="StMVAR"))
  expect_equal(get_regime_means_int(p=2, M=c(1, 1), d=2, params=theta_222gs, model="G-StMVAR", parametrization="intercept", constraints=NULL),
               pick_phi0(p=2, M=c(1, 1), d=2, params=theta_222gs_mu))
  expect_equal(get_regime_means_int(p=2, M=c(1, 1), d=2, params=theta_222gs_mu, model="StMVAR", parametrization="mean", constraints=NULL),
               calc_mu(p=2, M=c(1, 1), d=2, params=theta_222t, model="G-StMVAR"))
  expect_equal(get_regime_means_int(p=3, M=3, d=2, params=theta_332, parametrization="intercept", constraints=NULL),
               pick_phi0(p=3, M=3, d=2, params=theta_332_mu))
  expect_equal(get_regime_means_int(p=3, M=3, d=2, params=theta_332_mu, parametrization="mean", constraints=NULL),
               calc_mu(p=3, M=3, d=2, params=theta_332))
  expect_equal(get_regime_means_int(p=3, M=c(1, 2), d=2, params=theta_332gs, model="G-StMVAR", parametrization="intercept", constraints=NULL),
               pick_phi0(p=3, M=c(1, 2), d=2, params=theta_332gs_mu))
  expect_equal(get_regime_means_int(p=3, M=c(1, 2), d=2, params=theta_332gs_mu, model="G-StMVAR", parametrization="mean", constraints=NULL),
               calc_mu(p=3, M=c(1, 2), d=2, params=theta_332, model="G-StMVAR"))

  expect_equal(get_regime_means_int(p=1, M=2, d=3, params=theta_123, parametrization="intercept", constraints=NULL),
               pick_phi0(p=1, M=2, d=3, params=theta_123_mu))
  expect_equal(get_regime_means_int(p=1, M=2, d=3, params=theta_123_mu, parametrization="mean", constraints=NULL),
               calc_mu(p=1, M=2, d=3, params=theta_123))
  expect_equal(get_regime_means_int(p=1, M=2, d=3, params=theta_123t, model="StMVAR", parametrization="intercept", constraints=NULL),
               pick_phi0(p=1, M=2, d=3, params=theta_123t_mu))
  expect_equal(get_regime_means_int(p=1, M=2, d=3, params=theta_123t_mu, model="StMVAR", parametrization="mean", constraints=NULL),
               calc_mu(p=1, M=2, d=3, params=theta_123, model="StMVAR"))
  expect_equal(get_regime_means_int(p=2, M=1, d=3, params=theta_213, parametrization="intercept", constraints=NULL),
               pick_phi0(p=2, M=1, d=3, params=theta_213_mu))
  expect_equal(get_regime_means_int(p=2, M=1, d=3, params=theta_213_mu, parametrization="mean", constraints=NULL),
               calc_mu(p=2, M=1, d=3, params=theta_213))
  expect_equal(get_regime_means_int(p=2, M=1, d=3, params=theta_213t, model="StMVAR", parametrization="intercept", constraints=NULL),
               pick_phi0(p=2, M=1, d=3, params=theta_213t_mu))
  expect_equal(get_regime_means_int(p=2, M=1, d=3, params=theta_213t_mu, model="StMVAR", parametrization="mean", constraints=NULL),
               calc_mu(p=2, M=1, d=3, params=theta_213t, model="StMVAR"))

  expect_equal(get_regime_means_int(p=1, M=1, d=2, params=theta_112c, parametrization="intercept", constraints=C_112),
               pick_phi0(p=1, M=1, d=2, params=theta_112c_mu_exp))
  expect_equal(get_regime_means_int(p=1, M=1, d=2, params=theta_112c_mu, parametrization="mean", constraints=C_112),
               calc_mu(p=1, M=1, d=2, params=theta_112c, constraints=C_112))
  expect_equal(get_regime_means_int(p=1, M=1, d=2, params=theta_112tc, model="StMVAR", parametrization="intercept", constraints=C_112),
               pick_phi0(p=1, M=1, d=2, params=theta_112tc_mu_exp))
  expect_equal(get_regime_means_int(p=1, M=1, d=2, params=theta_112tc_mu, model="StMVAR", parametrization="mean", constraints=C_112),
               calc_mu(p=1, M=1, d=2, params=theta_112tc, model="StMVAR", constraints=C_112))
  expect_equal(get_regime_means_int(p=2, M=2, d=2, params=theta_222c, parametrization="intercept", constraints=C_222),
               pick_phi0(p=2, M=2, d=2, params=theta_222c_mu_exp))
  expect_equal(get_regime_means_int(p=2, M=2, d=2, params=theta_222c_mu, parametrization="mean", constraints=C_222),
               calc_mu(p=2, M=2, d=2, params=theta_222c, constraints=C_222))
  expect_equal(get_regime_means_int(p=2, M=2, d=2, params=theta_222_c2, parametrization="intercept", constraints=C_222_2),
               pick_phi0(p=2, M=2, d=2, params=theta_222c_mu2_exp))
  expect_equal(get_regime_means_int(p=2, M=2, d=2, params=theta_222c_mu2, parametrization="mean", constraints=C_222_2),
               calc_mu(p=2, M=2, d=2, params=theta_222_c2, constraints=C_222_2))
  expect_equal(get_regime_means_int(p=1, M=2, d=3, params=theta_123c, parametrization="intercept", constraints=C_123),
               pick_phi0(p=1, M=2, d=3, params=theta_123c_mu_exp))
  expect_equal(get_regime_means_int(p=1, M=2, d=3, params=theta_123c_mu, parametrization="mean", constraints=C_123),
               calc_mu(p=1, M=2, d=3, params=theta_123c, constraints=C_123))
  expect_equal(get_regime_means_int(p=1, M=2, d=3, params=theta_123tc, model="StMVAR", parametrization="intercept", constraints=C_123),
               pick_phi0(p=1, M=2, d=3, params=theta_123tc_mu_exp))
  expect_equal(get_regime_means_int(p=1, M=2, d=3, params=theta_123tc_mu, model="StMVAR", parametrization="mean", constraints=C_123),
               calc_mu(p=1, M=2, d=3, params=theta_123tc, model="StMVAR", constraints=C_123))

  # SGSMVAR
  expect_equal(get_regime_means_int(p=1, M=1, d=2, params=theta_112sWC, parametrization="intercept", constraints=NULL,
                                    structural_pars=list(W=W_112)),
               pick_phi0(p=1, M=1, d=2, params=theta_112sWC_mu, structural_pars=list(W=W_112)))
  expect_equal(get_regime_means_int(p=1, M=1, d=2, params=theta_112sWC_mu, parametrization="mean", constraints=NULL,
                                    structural_pars=list(W=W_112)),
               calc_mu(p=1, M=1, d=2, params=theta_112sWC, structural_pars=list(W=W_112)))
  expect_equal(get_regime_means_int(p=2, M=2, d=2, params=theta_222s, parametrization="intercept", constraints=NULL,
                                    structural_pars=list(W=W_222)),
               pick_phi0(p=2, M=2, d=2, params=theta_222s_mu, structural_pars=list(W=W_222)))
  expect_equal(get_regime_means_int(p=2, M=2, d=2, params=theta_222s_mu, parametrization="mean", constraints=NULL,
                                    structural_pars=list(W=W_222)),
               calc_mu(p=2, M=2, d=2, params=theta_222s, structural_pars=list(W=W_222)))

  expect_equal(get_regime_means_int(p=1, M=1, d=2, params=theta_112csWAR, parametrization="intercept", constraints=C_112,
                                    structural_pars=list(W=W_112)),
               pick_phi0(p=1, M=1, d=2, params=theta_112cs_mu_exp, structural_pars=list(W=W_112)))
  expect_equal(get_regime_means_int(p=1, M=1, d=2, params=theta_112cs_mu, parametrization="mean", constraints=C_112,
                                    structural_pars=list(W=W_112)),
               calc_mu(p=1, M=1, d=2, params=theta_112csWAR, constraints=C_112, structural_pars=list(W=W_112)))
  expect_equal(get_regime_means_int(p=2, M=2, d=2, params=theta_222csLAR, parametrization="intercept", constraints=C_222,
                                    structural_pars=list(W=W_222, C_lambda=C_lambda_222)),
               pick_phi0(p=2, M=2, d=2, params=theta_222cs_mu_exp, structural_pars=list(W=W_222, C_lambda=C_lambda_222)))
  expect_equal(get_regime_means_int(p=2, M=2, d=2, params=theta_222cs_mu, parametrization="mean", constraints=C_222,
                                    structural_pars=list(W=W_222, C_lambda=C_lambda_222)),
               calc_mu(p=2, M=2, d=2, params=theta_222csLAR, constraints=C_222, structural_pars=list(W=W_222, C_lambda=C_lambda_222)))

  # Same_means
  expect_equal(get_regime_means_int(p=2, M=2, d=2, params=theta_222c_int, parametrization="mean", constraints=C_222, same_means=list(1:2)),
               unname(cbind(phi10_222, phi10_222)))
  expect_equal(get_regime_means_int(p=2, M=2, d=2, params=theta_222csLAR_int, parametrization="mean", constraints=C_222,
                                    structural_pars=list(W=W_222, C_lambda=C_lambda_222), same_means=list(1:2)),
               unname(cbind(phi10_222, phi10_222)))
  expect_equal(get_regime_means_int(p=2, M=2, d=2, params=theta_222tcsLAR_int, model="StMVAR", parametrization="mean", constraints=C_222,
                                    structural_pars=list(W=W_222, C_lambda=C_lambda_222), same_means=list(1:2)),
               unname(cbind(phi10_222, phi10_222)))
  expect_equal(get_regime_means_int(p=2, M=c(1, 1), d=2, params=theta_222tcsLAR_int, model="G-StMVAR", parametrization="mean", constraints=C_222,
                                    structural_pars=list(W=W_222, C_lambda=C_lambda_222), same_means=list(1:2)),
               unname(cbind(phi10_222, phi10_222)))

  # Fixed alphas and lambdas
  expect_equal(c(get_regime_means_int(p=1, M=2, d=2, params=theta_122cwsF, model="GMVAR", constraints=C_122, weight_constraints=0.6,
                                    structural_pars=list(W=W_122, fixed_lambdas=c(4, 3)))),
               c(1.030956, 2.553492, 1.846211, 3.210253), tolerance=1e-5)
  expect_equal(c(get_regime_means_int(p=2, M=2, d=2, params=theta_222tcwsF, model="StMVAR", constraints=C_222, weight_constraints=0.6,
                                      structural_pars=list(W=W_222, fixed_lambdas=c(4, 3)))),
               c(-5.000000, 123.000000, 9.666667, 140.333333), tolerance=1e-5)
  expect_equal(c(get_regime_means_int(p=2, M=c(1, 1), d=2, params=theta_222gscmwsF, model="G-StMVAR", constraints=C_222, same_means=list(1:2),
                                      weight_constraints=0.6, structural_pars=list(W=W_222, fixed_lambdas=c(4, 3)))),
               c(-5, 123, -5, 123), tolerance=1e-5)
})


test_that("get_regime_autocovs_int works correctly", {
  expect_equal(get_regime_autocovs_int(p=1, M=1, d=2, params=theta_112, constraints=NULL)[, 2, 2, 1],
               c(0.2201706, 1.9959185), tolerance=1e-6)
  expect_equal(get_regime_autocovs_int(p=1, M=1, d=2, params=theta_112t, model="StMVAR", constraints=NULL)[, 2, 2, 1],
               c(0.2201706, 1.9959185), tolerance=1e-6)
  expect_equal(get_regime_autocovs_int(p=1, M=2, d=2, params=theta_122, constraints=NULL)[, 1, 1, 2],
               c(5.926843, 3.528684), tolerance=1e-5)
  expect_equal(get_regime_autocovs_int(p=2, M=2, d=2, params=theta_222, constraints=NULL)[, 1, 3, 2],
               c(36.50580, 13.77247), tolerance=1e-5)
  expect_equal(get_regime_autocovs_int(p=2, M=2, d=2, params=theta_222t, model="StMVAR", constraints=NULL)[, 1, 3, 2],
               c(36.50580, 13.77247), tolerance=1e-5)
  expect_equal(get_regime_autocovs_int(p=2, M=c(1, 1), d=2, params=theta_222gs, model="G-StMVAR", constraints=NULL)[, 1, 3, 2],
               c(36.50580, 13.77247), tolerance=1e-5)
  expect_equal(get_regime_autocovs_int(p=3, M=3, d=2, params=theta_332, constraints=NULL)[2, , 4, 3],
               c(2.89655, 4.28663), tolerance=1e-6)
  expect_equal(get_regime_autocovs_int(p=3, M=c(1, 2), d=2, params=theta_332gs, model="G-StMVAR", constraints=NULL)[2, , 4, 3],
               c(2.89655, 4.28663), tolerance=1e-6)
  expect_equal(get_regime_autocovs_int(p=1, M=2, d=3, params=theta_123, constraints=NULL)[1, , 2, 2],
               c(-0.3194244, -0.5931677, -0.8712061), tolerance=1e-6)
  expect_equal(get_regime_autocovs_int(p=1, M=2, d=3, params=theta_123t, model="StMVAR", constraints=NULL)[1, , 2, 2],
               c(-0.3194244, -0.5931677, -0.8712061), tolerance=1e-6)
  expect_equal(get_regime_autocovs_int(p=2, M=1, d=3, params=theta_213, constraints=NULL)[, 3, 3, 1],
               c(-0.582898, -1.038395, -1.413783), tolerance=1e-6)
  expect_equal(get_regime_autocovs_int(p=2, M=1, d=3, params=theta_213t, model="StMVAR", constraints=NULL)[, 3, 3, 1],
               c(-0.582898, -1.038395, -1.413783), tolerance=1e-6)
  expect_equal(get_regime_autocovs_int(p=2, M=2, d=2, params=theta_222c, constraints=C_222)[, 1, 3, 2],
               c(102.509805, 3.535469), tolerance=1e-6)


  # SGSMVAR
  expect_equal(get_regime_autocovs_int(p=1, M=1, d=2, params=theta_112sWC, constraints=NULL, structural_pars=list(W=W_112))[, 2, 2, 1],
               c(0.2201706, 1.9959185), tolerance=1e-6)
  expect_equal(get_regime_autocovs_int(p=1, M=1, d=2, params=theta_112tsWC, model="StMVAR", constraints=NULL, structural_pars=list(W=W_112))[, 2, 2, 1],
               c(0.2201706, 1.9959185), tolerance=1e-6)
  expect_equal(get_regime_autocovs_int(p=1, M=2, d=2, params=theta_122s, constraints=NULL, structural_pars=list(W=W_122))[, 1, 1, 2],
               c(5.926843, 3.528684), tolerance=1e-5)
  expect_equal(get_regime_autocovs_int(p=2, M=2, d=2, params=theta_222s, constraints=NULL, structural_pars=list(W=W_222))[, 1, 3, 2],
               c(36.50580, 13.77247), tolerance=1e-5)
  expect_equal(get_regime_autocovs_int(p=2, M=2, d=2, params=theta_222ts, model="StMVAR", constraints=NULL, structural_pars=list(W=W_222))[, 1, 3, 2],
               c(36.50580, 13.77247), tolerance=1e-5)
  expect_equal(get_regime_autocovs_int(p=3, M=3, d=2, params=theta_332sWC, constraints=NULL, structural_pars=list(W=W_332))[2, , 4, 3],
               c(1.396341, 2.174643), tolerance=1e-5)
  expect_equal(get_regime_autocovs_int(p=3, M=c(2, 1), d=2, params=theta_332gssWC, model="G-StMVAR", constraints=NULL, structural_pars=list(W=W_332))[2, , 4, 3],
               c(1.396341, 2.174643), tolerance=1e-5)

  expect_equal(get_regime_autocovs_int(p=1, M=2, d=3, params=theta_123s, constraints=NULL, structural_pars=list(W=W_123))[1, , 2, 2],
               c(-0.3194244, -0.5931677, -0.8712061), tolerance=1e-6)
  expect_equal(get_regime_autocovs_int(p=1, M=2, d=3, params=theta_123gss, model="G-StMVAR", constraints=NULL, structural_pars=list(W=W_123))[1, , 2, 2],
               c(-0.3194244, -0.5931677, -0.8712061), tolerance=1e-6)
  expect_equal(get_regime_autocovs_int(p=2, M=1, d=3, params=theta_213sWC, constraints=NULL, structural_pars=list(W=W_213))[, 3, 3, 1],
               c(-0.582898, -1.038395, -1.413783), tolerance=1e-6)
  expect_equal(get_regime_autocovs_int(p=1, M=2, d=3, params=theta_123csLAR, constraints=C_123, structural_pars=list(W=W_123, C_lambda=C_lambda_123))[, 1, 2, 1],
               c(0.2987136, 0.5610925, 0.8093760), tolerance=1e-6)

  # Same_means
  expect_equal(get_regime_autocovs_int(p=2, M=2, d=2, params=theta_222c_int, constraints=C_222, same_means=list(1:2))[, 2, 2, 1],
               c(-48.59455, 246.67447), tolerance=1e-5)
  expect_equal(get_regime_autocovs_int(p=2, M=2, d=2, params=theta_222c_int, constraints=C_222, same_means=list(1:2))[, 1, 1, 2],
               c(111.265016, 7.133962), tolerance=1e-5)
  expect_equal(get_regime_autocovs_int(p=2, M=c(1, 1), d=2, params=theta_222c_int, model="G-StMVAR", constraints=C_222, same_means=list(1:2))[, 1, 1, 2],
               c(111.265016, 7.133962), tolerance=1e-5)
  expect_equal(get_regime_autocovs_int(p=2, M=2, d=2, params=theta_222csLAR_int, constraints=C_222,
                                       structural_pars=list(W=W_222, C_lambda=C_lambda_222), same_means=list(1:2))[, 2, 1, 2],
               c(-21.57249, 97.57641), tolerance=1e-5)
  expect_equal(get_regime_autocovs_int(p=2, M=2, d=2, params=theta_222csLAR_int, constraints=C_222,
                                       structural_pars=list(W=W_222, C_lambda=C_lambda_222), same_means=list(1:2))[, 1, 2, 1],
               c(29.21410, -47.35956), tolerance=1e-5)
  expect_equal(get_regime_autocovs_int(p=2, M=c(1, 1), d=2, params=theta_222csLAR_int, model="G-StMVAR", constraints=C_222,
                                       structural_pars=list(W=W_222, C_lambda=C_lambda_222), same_means=list(1:2))[, 1, 2, 1],
               c(29.21410, -47.35956), tolerance=1e-5)

  # Fixed lambdas and alphas
  expect_equal(get_regime_autocovs_int(p=1, M=2, d=2, params=theta_122cwsF, model="GMVAR", constraints=C_122, weight_constraints=0.6,
                                       structural_pars=list(W=W_122, fixed_lambdas=c(4, 3)))[, 2, 2, 1],
               c(-0.2282206, 0.5365515), tolerance=1e-5)
  expect_equal(get_regime_autocovs_int(p=2, M=2, d=2, params=theta_222tcwsF, model="StMVAR", constraints=C_222, weight_constraints=0.6,
                                       structural_pars=list(W=W_222, fixed_lambdas=c(4, 3)))[, 2, 2, 1],
               c(-48.59455, 246.67447), tolerance=1e-5)
  expect_equal(get_regime_autocovs_int(p=2, M=c(1, 1), d=2, params=theta_222gscmwsF, model="G-StMVAR", constraints=C_222,
                                       same_means=list(1:2), weight_constraints=0.6,
                                       structural_pars=list(W=W_222, fixed_lambdas=c(4, 3)))[, 2, 2, 1],
               c(-48.59455, 246.67447), tolerance=1e-5)
})


test_that("uncond_moments_int works correctly", {
  expect_equal(uncond_moments_int(p=1, M=1, d=2, params=theta_112, parametrization="intercept", constraints=NULL)$uncond_mean,
               c(1.571661, 3.718636), tolerance=1e-6)
  expect_equal(uncond_moments_int(p=1, M=1, d=2, params=theta_112, parametrization="intercept", constraints=NULL)$autocors[, 1, 1],
               c(1.00000000, -0.02484575), tolerance=1e-6)
  expect_equal(uncond_moments_int(p=1, M=1, d=2, params=theta_112t, model="StMVAR", parametrization="intercept", constraints=NULL)$uncond_mean,
               c(1.571661, 3.718636), tolerance=1e-6)
  expect_equal(uncond_moments_int(p=1, M=1, d=2, params=theta_112t, model="StMVAR", parametrization="intercept", constraints=NULL)$autocors[, 1, 1],
               c(1.00000000, -0.02484575), tolerance=1e-6)
  expect_equal(uncond_moments_int(p=1, M=2, d=2, params=theta_122, parametrization="intercept", constraints=NULL)$uncond_mean,
               c(1.544567, 2.967251), tolerance=1e-6)
  expect_equal(uncond_moments_int(p=1, M=2, d=2, params=theta_122, parametrization="intercept", constraints=NULL)$autocors[, 2, 2],
               c(0.002322773, 0.095293066), tolerance=1e-6)
  expect_equal(uncond_moments_int(p=2, M=2, d=2, params=theta_222, parametrization="intercept", constraints=NULL)$uncond_mean,
               c(9.4270, 58.7715), tolerance=1e-4)
  expect_equal(uncond_moments_int(p=2, M=2, d=2, params=theta_222, parametrization="intercept", constraints=NULL)$autocors[1, , 3],
               c(0.9623361, -0.8532178), tolerance=1e-6)
  expect_equal(uncond_moments_int(p=2, M=2, d=2, params=theta_222, parametrization="intercept", constraints=NULL)$autocors[1, , 3],
               c(0.9623361, -0.8532178), tolerance=1e-6)
  expect_equal(uncond_moments_int(p=2, M=c(1, 1), d=2, params=theta_222gs, model="G-StMVAR", parametrization="intercept", constraints=NULL)$uncond_mean,
               c(9.4270, 58.7715), tolerance=1e-4)
  expect_equal(uncond_moments_int(p=2, M=c(1, 1), d=2, params=theta_222gs, model="G-StMVAR", parametrization="intercept", constraints=NULL)$autocors[1, , 3],
               c(0.9623361, -0.8532178), tolerance=1e-6)
  expect_equal(uncond_moments_int(p=3, M=3, d=2, params=theta_332, parametrization="intercept", constraints=NULL)$uncond_mean,
               c(2.108140, 3.872402), tolerance=1e-6)
  expect_equal(uncond_moments_int(p=3, M=3, d=2, params=theta_332, parametrization="intercept", constraints=NULL)$autocors[2, , 2],
               c(0.08213133, 0.17093255), tolerance=1e-6)
  expect_equal(uncond_moments_int(p=3, M=c(1, 2), d=2, params=theta_332gs, model="G-StMVAR", parametrization="intercept", constraints=NULL)$uncond_mean,
               c(2.108140, 3.872402), tolerance=1e-6)
  expect_equal(uncond_moments_int(p=3, M=c(1, 2), d=2, params=theta_332gs, model="G-StMVAR", parametrization="intercept", constraints=NULL)$autocors[2, , 2],
               c(0.08213133, 0.17093255), tolerance=1e-6)
  expect_equal(uncond_moments_int(p=1, M=2, d=3, params=theta_123, parametrization="intercept", constraints=NULL)$uncond_mean,
               c(2.263278, 4.278151, 6.265532), tolerance=1e-6)
  expect_equal(uncond_moments_int(p=1, M=2, d=3, params=theta_123, parametrization="intercept", constraints=NULL)$autocors[, 3, 1],
               c(0.7763858, 0.8178659, 1.0000000), tolerance=1e-6)
  expect_equal(uncond_moments_int(p=1, M=2, d=3, params=theta_123t, model="StMVAR", parametrization="intercept", constraints=NULL)$uncond_mean,
               c(2.263278, 4.278151, 6.265532), tolerance=1e-6)
  expect_equal(uncond_moments_int(p=1, M=2, d=3, params=theta_123t, model="StMVAR", parametrization="intercept", constraints=NULL)$autocors[, 3, 1],
               c(0.7763858, 0.8178659, 1.0000000), tolerance=1e-6)
  expect_equal(uncond_moments_int(p=2, M=1, d=3, params=theta_213, parametrization="intercept", constraints=NULL)$uncond_mean,
               c(1.1, 2.2, 3.3), tolerance=1e-1)
  expect_equal(uncond_moments_int(p=2, M=1, d=3, params=theta_213, parametrization="intercept", constraints=NULL)$autocors[3, , 3],
               c(-0.2324331, -0.2750250, -0.2901025), tolerance=1e-6)
  expect_equal(uncond_moments_int(p=2, M=2, d=2, params=theta_222c, parametrization="intercept", constraints=C_222)$uncond_mean,
               c(4.24, 133.92), tolerance=1e-2)
  expect_equal(uncond_moments_int(p=2, M=2, d=2, params=theta_222c, parametrization="intercept", constraints=C_222)$autocors[, 2, 1],
               c(0.1983422, 1.0000000), tolerance=1e-6)

  # SGSMVAR
  expect_equal(uncond_moments_int(p=1, M=1, d=2, params=theta_112sWC, parametrization="intercept", constraints=NULL, structural_pars=list(W=W_112))$uncond_mean,
               c(1.571661, 3.718636), tolerance=1e-6)
  expect_equal(uncond_moments_int(p=1, M=1, d=2, params=theta_112sWC, parametrization="intercept", constraints=NULL, structural_pars=list(W=W_112))$autocors[, 1, 1],
               c(1.00000000, -0.02484575), tolerance=1e-6)
  expect_equal(uncond_moments_int(p=1, M=1, d=2, params=theta_112tsWC, model="StMVAR", parametrization="intercept", constraints=NULL, structural_pars=list(W=W_112))$uncond_mean,
               c(1.571661, 3.718636), tolerance=1e-6)
  expect_equal(uncond_moments_int(p=1, M=1, d=2, params=theta_112tsWC, model="StMVAR", parametrization="intercept", constraints=NULL, structural_pars=list(W=W_112))$autocors[, 1, 1],
               c(1.00000000, -0.02484575), tolerance=1e-6)
  expect_equal(uncond_moments_int(p=1, M=2, d=2, params=theta_122s, parametrization="intercept", constraints=NULL, structural_pars=list(W=W_122))$uncond_mean,
               c(1.544567, 2.967251), tolerance=1e-6)
  expect_equal(uncond_moments_int(p=1, M=2, d=2, params=theta_122s, parametrization="intercept", constraints=NULL, structural_pars=list(W=W_122))$autocors[, 2, 2],
               c(0.002322773, 0.095293066), tolerance=1e-6)
  expect_equal(uncond_moments_int(p=2, M=2, d=2, params=theta_222s, parametrization="intercept", constraints=NULL, structural_pars=list(W=W_222))$uncond_mean,
               c(9.4270, 58.7715), tolerance=1e-4)
  expect_equal(uncond_moments_int(p=2, M=2, d=2, params=theta_222s, parametrization="intercept", constraints=NULL, structural_pars=list(W=W_222))$autocors[1, , 3],
               c(0.9623361, -0.8532178), tolerance=1e-6)
  expect_equal(uncond_moments_int(p=3, M=3, d=2, params=theta_332s, parametrization="intercept", constraints=NULL, structural_pars=list(W=W_332))$uncond_mean,
               c(2.108140, 3.872402), tolerance=1e-6)
  expect_equal(uncond_moments_int(p=3, M=3, d=2, params=theta_332s, parametrization="intercept", constraints=NULL, structural_pars=list(W=W_332))$autocors[2, , 2],
               c(0.4336740, 0.4434898), tolerance=1e-6)
  expect_equal(uncond_moments_int(p=1, M=2, d=3, params=theta_123s, parametrization="intercept", constraints=NULL, structural_pars=list(W=W_123))$uncond_mean,
               c(2.263278, 4.278151, 6.265532), tolerance=1e-6)
  expect_equal(uncond_moments_int(p=1, M=2, d=3, params=theta_123s, parametrization="intercept", constraints=NULL, structural_pars=list(W=W_123))$autocors[, 3, 1],
               c(0.7763858, 0.8178659, 1.0000000), tolerance=1e-6)
  expect_equal(uncond_moments_int(p=1, M=c(1, 1), d=3, params=theta_123gss, model="G-StMVAR", parametrization="intercept", constraints=NULL, structural_pars=list(W=W_123))$uncond_mean,
               c(2.263278, 4.278151, 6.265532), tolerance=1e-6)
  expect_equal(uncond_moments_int(p=1, M=c(1, 1), d=3, params=theta_123gss, model="G-StMVAR", parametrization="intercept", constraints=NULL, structural_pars=list(W=W_123))$autocors[, 3, 1],
               c(0.7763858, 0.8178659, 1.0000000), tolerance=1e-6)
  expect_equal(uncond_moments_int(p=2, M=1, d=3, params=theta_213sWC, parametrization="intercept", constraints=NULL, structural_pars=list(W=W_213))$uncond_mean,
               c(1.1, 2.2, 3.3), tolerance=1e-1)
  expect_equal(uncond_moments_int(p=2, M=1, d=3, params=theta_213sWC, parametrization="intercept", constraints=NULL, structural_pars=list(W=W_213))$autocors[3, , 3],
               c(-0.2324331, -0.2750250, -0.2901025), tolerance=1e-6)
  expect_equal(uncond_moments_int(p=2, M=1, d=3, params=theta_213tsWC, model="StMVAR", parametrization="intercept", constraints=NULL, structural_pars=list(W=W_213))$uncond_mean,
               c(1.1, 2.2, 3.3), tolerance=1e-1)
  expect_equal(uncond_moments_int(p=2, M=1, d=3, params=theta_213tsWC, model="StMVAR", parametrization="intercept", constraints=NULL, structural_pars=list(W=W_213))$autocors[3, , 3],
               c(-0.2324331, -0.2750250, -0.2901025), tolerance=1e-6)
  expect_equal(uncond_moments_int(p=2, M=2, d=2, params=theta_222csLAR, parametrization="intercept", constraints=C_222, structural_pars=list(W=W_222, C_lambda=C_lambda_222))$uncond_mean,
               c(4.24, 133.92), tolerance=1e-2)
  expect_equal(uncond_moments_int(p=2, M=2, d=2, params=theta_222csLAR, parametrization="intercept", constraints=C_222, structural_pars=list(W=W_222, C_lambda=C_lambda_222))$autocors[, 2, 1],
               c(0.2277718, 1.0000000), tolerance=1e-6)

  # Same_means
  expect_equal(uncond_moments_int(p=2, M=2, d=2, params=theta_222c_int, parametrization="mean", constraints=C_222, same_means=list(1:2))$uncond_mean,
               c(1.03, 2.36), tolerance=1e-6)
  expect_equal(uncond_moments_int(p=2, M=2, d=2, params=theta_222c_int, parametrization="mean", constraints=C_222, same_means=list(1:2))$autocors[, 1, 2],
               c(0.97148078, -0.08139881), tolerance=1e-6)
  expect_equal(uncond_moments_int(p=2, M=c(1, 1), d=2, params=theta_222gsc_int, model="G-StMVAR", parametrization="mean", constraints=C_222, same_means=list(1:2))$uncond_mean,
               c(1.03, 2.36), tolerance=1e-6)
  expect_equal(uncond_moments_int(p=2, M=c(1, 1), d=2, params=theta_222gsc_int, model="G-StMVAR", parametrization="mean", constraints=C_222, same_means=list(1:2))$autocors[, 1, 2],
               c(0.97148078, -0.08139881), tolerance=1e-6)
  expect_equal(uncond_moments_int(p=2, M=2, d=2, params=theta_222csLAR_int, parametrization="mean", constraints=C_222,
                                  structural_pars=list(W=W_222, C_lambda=C_lambda_222), same_means=list(1:2))$uncond_mean,
               c(1.03, 2.36), tolerance=1e-6)
  expect_equal(uncond_moments_int(p=2, M=2, d=2, params=theta_222csLAR_int, parametrization="mean", constraints=C_222,
                                  structural_pars=list(W=W_222, C_lambda=C_lambda_222), same_means=list(1:2))$autocors[, 2, 1],
               c(-0.6139921, 1.0000000), tolerance=1e-6)
  expect_equal(uncond_moments_int(p=2, M=2, d=2, params=theta_222tcsLAR_int, model="StMVAR", parametrization="mean", constraints=C_222,
                                  structural_pars=list(W=W_222, C_lambda=C_lambda_222), same_means=list(1:2))$uncond_mean,
               c(1.03, 2.36), tolerance=1e-6)
  expect_equal(uncond_moments_int(p=2, M=2, d=2, params=theta_222tcsLAR_int, model="StMVAR", parametrization="mean", constraints=C_222,
                                  structural_pars=list(W=W_222, C_lambda=C_lambda_222), same_means=list(1:2))$autocors[, 2, 1],
               c(-0.6139921, 1.0000000), tolerance=1e-6)

  # Fixed lambdas and alphas
  expect_equal(uncond_moments_int(p=1, M=2, d=2, params=theta_122cwsF, model="GMVAR", constraints=C_122, weight_constraints=0.6,
                                  structural_pars=list(W=W_122, fixed_lambdas=c(4, 3)))$autocors[, 2, 1],
               c(0.01365825, 1.00000000), tolerance=1e-6)
  expect_equal(uncond_moments_int(p=2, M=2, d=2, params=theta_222tcwsF, model="StMVAR", constraints=C_222, weight_constraints=0.6,
                                  structural_pars=list(W=W_222, fixed_lambdas=c(4, 3)))$uncond_mean,
               c(0.8666667, 129.9333333), tolerance=1e-6)
  expect_equal(uncond_moments_int(p=2, M=c(1, 1), d=2, params=theta_222gscmwsF, model="G-StMVAR", constraints=C_222, same_means=list(1:2),
                                  weight_constraints=0.6, structural_pars=list(W=W_222, fixed_lambdas=c(4, 3)))$autocors[, 2, 1],
               c(-0.5008924, 1.0000000), tolerance=1e-6)
})


test_that("non_int uncond moment functions work", {
  mod122cwsF <- GSMVAR(p=1, M=2, d=2, params=theta_122cwsF, model="GMVAR", constraints=C_122, weight_constraints=0.6,
                       structural_pars=list(W=W_122, fixed_lambdas=c(4, 3)))
  mod222tcwsF <- GSMVAR(p=2, M=2, d=2, params=theta_222tcwsF, model="StMVAR", constraints=C_222, weight_constraints=0.6,
                        structural_pars=list(W=W_222, fixed_lambdas=c(4, 3)))
  mod222gscmwsF <- GSMVAR(p=2, M=c(1, 1), d=2, params=theta_222gscmwsF, model="G-StMVAR", constraints=C_222, same_means=list(1:2),
                          parametrization="mean", weight_constraints=0.6, structural_pars=list(W=W_222, fixed_lambdas=c(4, 3)))


  mod122 <- GSMVAR(p=1, M=2, d=2, params=theta_122)
  mod112csWAR <- GSMVAR(p=1, M=1, d=2, params=theta_112csWAR, structural_pars=list(W=W_112))
  mod112tcsWAR <- GSMVAR(p=1, M=1, d=2, params=c(theta_112csWAR, 10), model="StMVAR", structural_pars=list(W=W_112))
  mod222csLAR <- GSMVAR(p=2, M=2, d=2, params=theta_222csLAR, constraints=C_222,
                       structural_pars=list(W=W_222, C_lambda=C_lambda_222))
  mod222gscsLAR <- GSMVAR(p=2, M=c(1, 1), d=2, params=c(theta_222csLAR, 10), model="G-StMVAR", constraints=C_222,
                        structural_pars=list(W=W_222, C_lambda=C_lambda_222))
  mod222c_int <- GSMVAR(p=2, M=2, d=2, params=theta_222c_int, parametrization="mean", constraints=C_222, same_means=list(1:2))
  mod222tc_int <- GSMVAR(p=2, M=2, d=2, params=c(theta_222c_int, 10, 20), model="StMVAR",
                         parametrization="mean", constraints=C_222, same_means=list(1:2))


  unc122 <- uncond_moments(mod122)
  unc112csWAR <- uncond_moments(mod112csWAR)
  unc112tcsWAR <- uncond_moments(mod112tcsWAR)
  unc222csLAR <- uncond_moments(mod222csLAR)
  unc222gscsLAR <- uncond_moments(mod222gscsLAR)
  unc222c_int <- uncond_moments(mod222c_int)
  unc222tc_int <- uncond_moments(mod222tc_int)

  unc122cwsF <- uncond_moments(mod122cwsF)
  unc222tcwsF <- uncond_moments(mod222tcwsF)
  unc222gscmwsF <- uncond_moments(mod222gscmwsF)

  expect_equal(unc122cwsF, uncond_moments_int(p=1, M=2, d=2, params=theta_122cwsF, model="GMVAR", constraints=C_122,
                                              weight_constraints=0.6, structural_pars=list(W=W_122, fixed_lambdas=c(4, 3))), tolerance=1e-6)
  expect_equal(unc222tcwsF, uncond_moments_int(p=2, M=2, d=2, params=theta_222tcwsF, model="StMVAR", constraints=C_222, weight_constraints=0.6,
                                               structural_pars=list(W=W_222, fixed_lambdas=c(4, 3))), tolerance=1e-6)
  expect_equal(unc222gscmwsF, uncond_moments_int(p=2, M=c(1, 1), d=2, params=theta_222gscmwsF, model="G-StMVAR", constraints=C_222,
                                                 same_means=list(1:2), weight_constraints=0.6, parametrization="mean",
                                                 structural_pars=list(W=W_222, fixed_lambdas=c(4, 3))), tolerance=1e-6)



  expect_equal(unc122, uncond_moments_int(p=1, M=2, d=2, params=theta_122), tolerance=1e-6)
  expect_equal(unc112csWAR, uncond_moments_int(p=1, M=1, d=2, params=theta_112csWAR,
                                               structural_pars=list(W=W_112)), tolerance=1e-6)
  expect_equal(unc112tcsWAR, uncond_moments_int(p=1, M=1, d=2, params=c(theta_112tcsWAR, 10), model="StMVAR",
                                                structural_pars=list(W=W_112)), tolerance=1e-6)
  expect_equal(unc222csLAR$autocovs[1, 2, ], c(27.91380, 27.59758, 27.16074), tolerance=1e-4)
  expect_equal(unc222csLAR$autocors[2, , 1], c(0.2277718, 1.0000000), tolerance=1e-4)
  expect_equal(unc222csLAR$uncond_mean, c(4.24, 133.92), tolerance=1e-4)

  expect_equal(unc222gscsLAR$autocovs[1, 2, ], c(27.91380, 27.59758, 27.16074), tolerance=1e-4)
  expect_equal(unc222gscsLAR$autocors[2, , 1], c(0.2277718, 1.0000000), tolerance=1e-4)
  expect_equal(unc222gscsLAR$uncond_mean, c(4.24, 133.92), tolerance=1e-4)
  expect_equal(unc222c_int,
               uncond_moments_int(p=2, M=2, d=2, params=theta_222c_int,
                                  parametrization="mean", constraints=C_222, same_means=list(1:2)), tolerance=1e-6)
  expect_equal(unc222tc_int,
               uncond_moments_int(p=2, M=2, d=2, params=c(theta_222c_int, 10, 20), model="StMVAR",
                                  parametrization="mean", constraints=C_222, same_means=list(1:2)), tolerance=1e-6)

  reg_means122 <- get_regime_means(mod122)
  reg_means112csWAR <- get_regime_means(mod112csWAR)
  reg_means112tcsWAR <- get_regime_means(mod112tcsWAR)
  reg_means222csLAR <- get_regime_means(mod222csLAR)
  reg_means222gscsLAR <- get_regime_means(mod222gscsLAR)
  reg_means222c_int <- get_regime_means(mod222c_int)
  reg_means222tc_int <- get_regime_means(mod222tc_int)

  expect_equal(reg_means122[2, ], c(2.553492, 3.210253), tolerance=1e-5)
  expect_equal(reg_means112csWAR[, 1], c(1.571661, 3.718636), tolerance=1e-5)
  expect_equal(reg_means112tcsWAR[, 1], c(1.571661, 3.718636), tolerance=1e-5)
  expect_equal(reg_means222csLAR[, 2], c(9.666667, 140.333333), tolerance=1e-5)
  expect_equal(reg_means222gscsLAR[, 2], c(9.666667, 140.333333), tolerance=1e-5)
  expect_equal(reg_means222c_int[, 2], c(1.03, 2.36), tolerance=1e-5)
  expect_equal(reg_means222tc_int[, 2], c(1.03, 2.36), tolerance=1e-5)

  reg_autocovs122 <- get_regime_autocovs(mod122)
  reg_autocovs112csWAR <- get_regime_autocovs(mod112csWAR)
  reg_autocovs112tcsWAR <- get_regime_autocovs(mod112tcsWAR)
  reg_autocovs222csLAR <- get_regime_autocovs(mod222csLAR)
  reg_autocovs222gscsLAR <- get_regime_autocovs(mod222gscsLAR)
  reg_autocovs222c_int <- get_regime_autocovs(mod222c_int)
  reg_autocovs222tc_int <- get_regime_autocovs(mod222tc_int)

  expect_equal(reg_autocovs122[2 ,2 ,1 , ], c(5.258146, 9.877770), tolerance=1e-5)
  expect_equal(reg_autocovs112csWAR[1, ,2 , 1], c(0.2477767, 0.2201706), tolerance=1e-5)
  expect_equal(reg_autocovs112tcsWAR[1, ,2 , 1], c(0.2477767, 0.2201706), tolerance=1e-5)
  expect_equal(reg_autocovs222csLAR[1, , 2, 2], c(9.302809, -21.716886), tolerance=1e-5)
  expect_equal(reg_autocovs222gscsLAR[1, , 2, 2], c(9.302809, -21.716886), tolerance=1e-5)
  expect_equal(reg_autocovs222c_int[, 2, 2, 1], c(-48.59455, 246.67447), tolerance=1e-5)
  expect_equal(reg_autocovs222tc_int[, 2, 2, 1], c(-48.59455, 246.67447), tolerance=1e-5)
})


alt_pcovmat <- function(p, d, all_A, all_Omega) {
  # Calculate the (dp x dp) covariance matrix Sigma_{m,p} (Lutkepohl 2005, eq. (2.1.39))
  all_boldA <- form_boldA(p=p, M=1, d=d, all_A=all_A)
  I_dp2 <- diag(nrow=(d*p)^2)
  ZER_lower <- matrix(0, nrow=d*(p - 1), ncol=d*p)
  ZER_right <- matrix(0, nrow=d, ncol=d*(p - 1))
  kronmat <- I_dp2 - kronecker(all_boldA[, , 1], all_boldA[, , 1])
  sigma_epsm <- rbind(cbind(all_Omega[, , 1], ZER_right), ZER_lower)
  Sigma_m <- solve(kronmat, vec(sigma_epsm))
  matrix(Sigma_m, nrow=d*p, ncol=d*p)
}

# pMd
params112 <- c(-0.114, 4.704, 0.23, -0.522, 0.014, 0.217, 0.04, -0.047, 1.649)
all_A112 <- pick_allA(p=1, M=1, d=2, params=params112)
all_Omega112 <- pick_Omegas(p=1, M=1, d=2, params=params112)
all_boldA112 <- form_boldA(p=1, M=1, d=2, all_A=all_A112)

params115 <- c(0.936, 3.738, 4.135, 9.663, 3.955, -0.579, -0.009, -0.187, 0.019, -0.241, -0.02, -0.128, -0.036,
               -0.167, -0.218, 0.042, -0.259, -0.029, -0.03, -0.133, 0.377, -0.321, 0.017, -0.205, -0.29, -0.01,
               0.483, 0.286, -0.035, -0.216, 2.373, -0.613, 1.99, -2.255, 0.08, 1.14, 0.319, 0.342, -0.723, 2.832,
               -1.437, -0.769, 4.711, 2.736, 6.54)
all_A115 <- pick_allA(p=1, M=1, d=5, params=params115)
all_Omega115 <- pick_Omegas(p=1, M=1, d=5, params=params115)

params212 <- c(1.551, 4.498, 0.146, -0.079, -0.045, -0.487, 0.271, -0.087, -0.697, -0.221, 0.25, 0.512, 2.507)
all_A212 <- pick_allA(p=2, M=1, d=2, params=params212)
all_Omega212 <- pick_Omegas(p=2, M=1, d=2, params=params212)

params213 <- c(1.076, 1.777, -1.985, -0.141, 0.234, 0.088, -0.311, 0.112, 0.234, -0.421, -0.217, 0.017, -0.067,
               -0.104, 0.067, 0.075, 0.454, -0.249, 0.034, -0.247, -0.231, 0.499, -0.218, -1.044, 0.422, 0.986, 3.175)
all_A213 <- pick_allA(p=2, M=1, d=3, params=params213)
all_Omega213 <- pick_Omegas(p=2, M=1, d=3, params=params213)
all_boldA213 <- form_boldA(p=2, M=1, d=3, all_A=all_A213)

params312 <- c(0.712, 4.819, -0.147, 0.365, -0.159, -0.055, 0.52, -0.257, -0.057, -0.313, 0.373, 0.182, -0.109,
               -0.111, 1.326, -0.895, 0.66)
all_A312 <- pick_allA(p=3, M=1, d=2, params=params312)
all_Omega312 <- pick_Omegas(p=3, M=1, d=2, params=params312)

params714 <- c(2.211, 3.498, 5.362, 4.747, -0.114, -0.073, -0.138, 0.186, -0.062, -0.015, -0.102, -0.276, 0.203,
               0.084, 0.087, 0.049, 0.257, 0.262, 0.037, 0.087, 0.1, -0.243, -0.109, -0.275, 0.199, -0.125, -0.081,
               -0.129, 0.323, 0.426, -0.021, 0.323, 0.062, 0.057, 0.076, -0.045, -0.303, 0.093, 0.072, -0.209, 0.219,
               0.068, -0.131, -0.381, -0.03, -0.101, -0.083, 0.26, -0.116, 0.011, -0.043, 0.099, 0.114, 0.152, -0.209,
               0.003, 0.258, 0.047, 0.457, -0.034, -0.291, 0.325, 0.033, -0.004, -0.017, 0.244, 0.027, 0.054, 0.095,
               0.079, 0.258, -0.101, 0.09, 0.014, -0.123, -0.079, 0.178, 0.27, 0.173, 0.161, -0.075, -0.118, -0.245,
               0.08, -0.172, 0.261, -0.074, -0.046, 0.137, -0.081, -0.225, -0.055, 0.136, -0.366, 0.04, 0.162, -0.137,
               0.292, 0.087, 0.074, -0.026, -0.137, -0.011, 0.068, 0.083, 0.11, 0.139, 0.091, -0.297, -0.116, 0.166,
               -0.064, 0.132, -0.062, 0.251, -0.042, 0.681, -0.328, 0.944, 0.585, 2.154, 0.684, -2.075, 4.601, -1.727, 3.392)
all_A714 <- pick_allA(p=7, M=1, d=4, params=params714)
all_Omega714 <- pick_Omegas(p=7, M=1, d=4, params=params714)
all_boldA714 <- form_boldA(p=7, M=1, d=4, all_A=all_A714)

params516 <- c(0.225461, 0.2072, -0.143174, -0.089739, -0.018469, -0.050561, -0.170449, 0.056068, 0.008218, 0.104494,
               -0.054231, -0.065922, -0.101427, 0.303574, 0.09502, 0.173378, 0.154364, -0.050695, 0.125909, 0.423252,
               -0.090888, 0.008861, 0.033686, 0.229432, -0.204607, 0.12374, -0.246613, -0.272826, -0.079953, -0.155406,
               0.147952, -0.045087, 0.135298, 0.024628, 0.112773, -0.087763, -0.040376, -0.072829, 0.006114, -0.006132,
               -0.081802, -0.099939, -0.026151, -0.080008, -0.066702, 0.160827, 0.10961, 0.032163, 0.278354, -0.151524,
               -0.159157, 0.069496, 0.2282, 0.141799, -0.042404, -0.035857, 0.001668, 0.051633, 0.177808, -0.108081,
               -0.003776, 0.086418, -0.177295, -0.007459, 0.001096, -0.08415, 0.077412, -0.153839, -0.222861, -0.221952,
               0.029957, 0.081076, 0.187479, 0.236314, 0.069664, -0.034321, 0.138018, -0.090784, 0.078337, 0.11314,
               0.072923, -0.069141, 0.155786, 0.095607, -0.013958, 0.25297, 0.143652, 0.085351, -0.191761, 0.110038,
               0.122346, -0.068535, 0.090091, 0.069526, -0.131725, -0.173419, 0.080161, 0.013833, -0.163593, -0.001283,
               -0.011312, -0.020348, -0.004064, 0.089246, 0.115923, -0.001136, 0.108662, 0.044961, 0.150245, -0.355093,
               -0.065505, -0.08575, -0.084668, 0.024026, -0.114004, -0.00582, -0.014616, -0.11759, 0.217696, 0.088498,
               0.078428, -0.24494, 0.108445, -0.085994, 0.010226, -0.012086, 0.062526, 0.032096, 0.15084, 0.187635,
               0.125639, 0.05251, 0.111083, 0.091175, -0.042892, 0.016866, -0.035713, -0.028081, -0.018831, -0.174531,
               -0.084074, 0.047865, 0.122736, 0.11565, -0.131177, 0.187818, 0.123726, -0.065527, 0.019519, 0.037157,
               -0.368623, 0.195931, -0.015609, -0.162628, -0.024976, -0.138307, -0.112833, -0.01342, -0.152109, -0.127772,
               -0.246636, 0.073784, -0.110782, -0.007831, -0.003678, 0.048329, -0.294818, 0.267569, -0.020779, 0.251712,
               -0.087998, -0.010173, -0.012708, -0.014015, 0.206771, 0.188524, -0.057268, -0.050239, 0.098685, -0.18502,
               0.33091, -0.154566, 0.044506, -0.18042, -0.087658, 0.090132, 0.934282, -0.050769, 0.148077, 0.137113,
               -0.141831, 0.107824, 0.549019, -0.137383, 0.08692, -0.036844, 0.092005, 0.782197, 0.003268, 0.06118,
               0.088202, 0.701364, -0.142686, 0.224613, 0.735951, -0.054841, 0.979966)
all_A516 <- pick_allA(p=5, M=1, d=6, params=params516)
all_Omega516 <- pick_Omegas(p=5, M=1, d=6, params=params516)

params412 <- c(3.201, 1.754, -2.525, 1.181, -1.138, 0.693, -3.132, 5.712, -1.882, 2.174, -1.832, 6.152, -1.5, 3.525,
               -0.297, 1.859, -0.786, 3.065, 2.544, 2.446, 2.805)
all_A412 <- pick_allA(p=4, M=1, d=2, params=params412)
all_Omega412 <- pick_Omegas(p=4, M=1, d=2, params=params412)

params413 <- c(-0.756, 0.121, 2.396, -1.019, -3.577, -2.476, 1.372, 3.627, 1.099, -1.98, -3.385, 0.109, 0.723, 2.879,
               1.895, -2.923, -4.605, -1.746, 2.225, 3.212, 1.138, -2.113, -4.057, -1.349, 1.13, 1.123, -0.01, -0.357,
               -0.428, -0.677, 1.995, 2.943, 1.059, -0.368, 0.026, 0.268, -0.308, -0.642, -0.022, 1.339, -1.746, 0.607,
               2.511, -1.196, 1.17)
all_A413 <- pick_allA(p=4, M=1, d=3, params=params413)
all_Omega413 <- pick_Omegas(p=4, M=1, d=3, params=params413)

params512 <- c(-1.413, 1.365, 1.333, 0.491, 0.073, -0.233, -1.222, 0.302, 0.046, 1.803, 1.099, -0.423, -1.833, 0.412,
               -0.238, -0.545, 2.679, -1.216, 0.006, -0.221, -0.975, -0.543, 1.86, 0.459, 1.203)
all_A512 <- pick_allA(p=5, M=1, d=2, params=params512)
all_Omega512 <- pick_Omegas(p=5, M=1, d=2, params=params512)

params314 <- c(-0.888, 3.587, 1.846, 6.059, 0.54, 0.847, -0.409, -0.053, -1.328, 1.203, -2.687, 0.804, 0.005, -0.067,
               -1.156, 0.765, 0.662, 0.269, 0.584, 0.17, 0.377, -0.612, 1.724, -1.784, 1.445, -0.185, 2.275, -0.267,
               -1.149, 0.441, -1.468, 0.835, 0.361, -0.781, 1.687, 0.039, -0.113, 0.042, -0.81, -0.882, 0.102, -0.316,
               0.774, 0.458, 0.653, 0.406, -0.119, 0.939, -0.481, 0.169, -0.835, 0.549, 0.817, -0.224, -0.043, -0.638,
               3.028, 0.801, 0.083, 0.762, 1.942, 8.485)
all_A314 <- pick_allA(p=3, M=1, d=4, params=params314)
all_Omega314 <- pick_Omegas(p=3, M=1, d=4, params=params314)


test_that("VAR_pcovmat works correctly", {
  expect_equal(VAR_pcovmat(p=1, d=2, all_Am=array(all_A112[, , , 1], dim=c(2, 2, 1)), Omega_m=all_Omega112[, , 1]),
               alt_pcovmat(p=1, d=2, all_A=all_A112, all_Omega=all_Omega112), tolerance=1e-5)

  expect_equal(VAR_pcovmat(p=1, d=5, all_Am=array(all_A115[, , , 1], dim=c(5, 5, 1)), Omega_m=all_Omega115[, , 1]),
               alt_pcovmat(p=1, d=5, all_A=all_A115, all_Omega=all_Omega115), tolerance=1e-5)

  expect_equal(VAR_pcovmat(p=2, d=2, all_Am=all_A212[, , , 1], Omega_m=all_Omega212[, , 1]),
               alt_pcovmat(p=2, d=2, all_A=all_A212, all_Omega=all_Omega212), tolerance=1e-5)

  expect_equal(VAR_pcovmat(p=2, d=3, all_Am=all_A213[, , , 1], Omega_m=all_Omega213[, , 1]),
               alt_pcovmat(p=2, d=3, all_A=all_A213, all_Omega=all_Omega213), tolerance=1e-5)

  expect_equal(VAR_pcovmat(p=3, d=2, all_Am=all_A312[, , , 1], Omega_m=all_Omega312[, , 1]),
               alt_pcovmat(p=3, d=2, all_A=all_A312, all_Omega=all_Omega312), tolerance=1e-5)

  expect_equal(VAR_pcovmat(p=7, d=4, all_Am=all_A714[, , , 1], Omega_m=all_Omega714[, , 1]),
               alt_pcovmat(p=7, d=4, all_A=all_A714, all_Omega=all_Omega714), tolerance=1e-5)

  expect_equal(VAR_pcovmat(p=5, d=6, all_Am=all_A516[, , , 1], Omega_m=all_Omega516[, , 1]),
               alt_pcovmat(p=5, d=6, all_A=all_A516, all_Omega=all_Omega516), tolerance=1e-5)

  expect_equal(VAR_pcovmat(p=4, d=2, all_Am=all_A412[, , , 1], Omega_m=all_Omega412[, , 1]),
               alt_pcovmat(p=4, d=2, all_A=all_A412, all_Omega=all_Omega412), tolerance=1e-5)

  expect_equal(VAR_pcovmat(p=4, d=3, all_Am=all_A413[, , , 1], Omega_m=all_Omega413[, , 1]),
               alt_pcovmat(p=4, d=3, all_A=all_A413, all_Omega=all_Omega413), tolerance=1e-5)

  expect_equal(VAR_pcovmat(p=5, d=2, all_Am=all_A512[, , , 1], Omega_m=all_Omega512[, , 1]),
               alt_pcovmat(p=5, d=2, all_A=all_A512, all_Omega=all_Omega512), tolerance=1e-5)

  expect_equal(VAR_pcovmat(p=3, d=4, all_Am=all_A314[, , , 1], Omega_m=all_Omega314[, , 1]),
               alt_pcovmat(p=3, d=4, all_A=all_A314, all_Omega=all_Omega314), tolerance=1e-5)
})


test_that("get_Sigmas works correctly", {
  expect_equal(get_Sigmas(p=1, M=1, d=2, all_A=all_A112, all_Omega=all_Omega112, all_boldA=all_boldA112),
               array(alt_pcovmat(p=1, d=2, all_A=all_A112, all_Omega=all_Omega112), dim=c(2, 2, 1)), tolerance=1e-5)

  expect_equal(get_Sigmas(p=2, M=1, d=3, all_A=all_A213, all_Omega=all_Omega213, all_boldA=all_boldA213),
               array(alt_pcovmat(p=2, d=3, all_A=all_A213, all_Omega=all_Omega213), dim=c(6, 6, 1)), tolerance=1e-5)

  expect_equal(get_Sigmas(p=7, M=1, d=4, all_A=all_A714, all_Omega=all_Omega714, all_boldA=all_boldA714),
               array(alt_pcovmat(p=7, d=4, all_A=all_A714, all_Omega=all_Omega714), dim=c(28, 28, 1)), tolerance=1e-5)
})
