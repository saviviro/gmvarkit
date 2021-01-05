context("generate parameters")
library(gmvarkit)

data <- cbind(10*eurusd[,1], 100*eurusd[,2])

params122 <- c(0.623, -0.129, 0.959, 0.089, -0.006, 1.006, 1.746, 0.804, 5.804, 3.245, 7.913,
               0.952, -0.037, -0.019, 0.943, 6.926, 3.982, 12.135, 0.789) # p=1, M=2, d=2
set.seed(1); W_313 <- chol(tcrossprod(matrix(rnorm(3*3), nrow=3)))
params313s <- random_ind2(p=3, M=1, d=3, mu_scale=1:3, mu_scale2=1:3, W_scale=1:3, lambda_scale=1:3, structural_pars=list(W=W_313)) # W idenfitied with Cholesky
set.seed(1); W_122 <- tcrossprod(matrix(rnorm(2*2), nrow=2))
C_lambda_122 <- matrix(1, nrow=1)
params122s <- random_ind2(p=1, M=2, d=2, mu_scale=1:2, mu_scale2=1:2, W_scale=1:2, lambda_scale=1,
                          structural_pars=list(W=W_122, C_lambda=C_lambda_122))

params222s_sm <- random_ind2(p=2, M=2, d=2, mu_scale=1:2, mu_scale2=1:2, W_scale=1:2, lambda_scale=1,
                             same_means=list(1:2),
                             structural_pars=list(W=W_122, C_lambda=C_lambda_122))


test_that("generate parameters don't throw errors", {
  test_length0 <- function(x, length_x) expect_equal(length(x), length_x)
  test_length0(random_ind(p=2, M=3, d=2, mu_scale=1:2, mu_scale2=1:2, omega_scale=1:2), 41)
  test_length0(random_ind(p=3, M=1, d=3, mu_scale=1:3, mu_scale2=1:3, W_scale=1:3, lambda_scale=1:3,
                          structural_pars=list(W=chol(tcrossprod(matrix(rnorm(3*3), nrow=3))))), 36)
  test_length0(random_ind(p=1, M=2, d=2, mu_scale=1:2, mu_scale2=1:2, W_scale=1:2, lambda_scale=1,
                          structural_pars=list(W=tcrossprod(matrix(rnorm(2*2), nrow=2)), C_lambda=matrix(1, nrow=1))), 18)
  test_length0(random_ind(p=2, M=2, d=2, mu_scale=1:2, mu_scale2=1:2, W_scale=1:2, lambda_scale=1, same_means=list(1:2),
                          structural_pars=list(W=tcrossprod(matrix(rnorm(2*2), nrow=2)), C_lambda=matrix(1, nrow=1))), 24)

  test_length0(smart_ind(p=1, M=2, d=2, params=params122, which_random=1, mu_scale=1:2, mu_scale2=1:2, omega_scale=1:2), 19)
  test_length0(smart_ind(p=1, M=2, d=2, params=params122s, which_random=2, mu_scale=1:2, mu_scale2=1:2, W_scale=1:2,
                         lambda_scale=2, structural_pars=list(W=W_122, C_lambda=C_lambda_122)), 18)
  test_length0(smart_ind(p=3, M=1, d=3, params=params313s, which_random=1, mu_scale=1:2, mu_scale2=1:2, W_scale=1:3,
                         lambda_scale=1:2, structural_pars=list(W=W_313)), 36)
  test_length0(smart_ind(p=2, M=2, d=2, params=params222s_sm, which_random=2, mu_scale=1:2, mu_scale2=1:2, W_scale=1:2,
                         lambda_scale=2, same_means=list(1:2), structural_pars=list(W=W_122, C_lambda=C_lambda_122)), 24)

  test_length0(random_ind2(p=3, M=1, d=3, mu_scale=1:3, mu_scale2=1:3, omega_scale=1:3), 36)
  test_length0(random_ind2(p=3, M=1, d=3, mu_scale=1:3, mu_scale2=1:3, W_scale=1:3, lambda_scale=1:3,
                           structural_pars=list(W=chol(tcrossprod(matrix(rnorm(3*3), nrow=3))))), 36)
  test_length0(random_ind2(p=1, M=2, d=2, mu_scale=1:2, mu_scale2=1:2, W_scale=1:2, lambda_scale=1,
                           structural_pars=list(W=tcrossprod(matrix(rnorm(2*2), nrow=2)), C_lambda=matrix(1, nrow=1))), 18)
  test_length0(random_ind2(p=2, M=2, d=2, mu_scale=1:2, mu_scale2=1:2, W_scale=1:2, lambda_scale=1, same_means=list(1:2),
                           structural_pars=list(W=tcrossprod(matrix(rnorm(2*2), nrow=2)), C_lambda=matrix(1, nrow=1))), 24)
  test_length0(random_ind2(p=1, M=2, d=2, mu_scale=1:2, mu_scale2=1:2, W_scale=1:2, lambda_scale=1, omega_scale=1:2,
                           structural_pars=NULL, same_means=list(1:2)), 17)

  test_length0(random_coefmats(d=2, how_many=1, scale=1), 4)
  test_length0(random_coefmats2(p=2, d=3), 18)
  test_length0(random_covmat(d=3, omega_scale=1:3), 6)
  test_length0(random_covmat(d=3, M=3, omega_scale=1:3, W_scale=1:3, lambda_scale=1:3,
                             structural_pars=list(W=chol(tcrossprod(matrix(rnorm(3*3), nrow=3))))), 12)
  test_length0(random_covmat(d=2, M=2, omega_scale=1:2, W_scale=1:2, lambda_scale=1,
                             structural_pars=list(W=tcrossprod(matrix(rnorm(2*2), nrow=2)), C_lambda=matrix(1, nrow=1))), 5)
  test_length0(smart_covmat(d=2, Omega=diag(2), accuracy=1), 3)
  test_length0(smart_covmat(d=3, M=2, W_and_lambdas=c(rnorm(3*3), abs(rnorm(3*2))), accuracy=10,
                             structural_pars=list(1)), 15)
})
