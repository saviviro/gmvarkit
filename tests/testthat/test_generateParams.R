context("generate parameters")
library(gmvarkit)

data <- cbind(10*eurusd[,1], 100*eurusd[,2])

params122 <- c(0.623, -0.129, 0.959, 0.089, -0.006, 1.006, 1.746, 0.804, 5.804, 3.245, 7.913,
               0.952, -0.037, -0.019, 0.943, 6.926, 3.982, 12.135, 0.789) # p=1, M=2, d=2

test_that("generate parameters don't throw errors", {
  test_length0 <- function(x, length_x) expect_equal(length(x), length_x)
  test_length0(random_ind(p=2, M=3, d=2, mu_scale=1:2, mu_scale2=1:2, omega_scale=1:2), 41)
  test_length0(smart_ind(p=1, M=2, d=2, params=params122, which_random=1, mu_scale=1:2, mu_scale2=1:2, omega_scale=1:2), 19)
  test_length0(random_ind2(p=3, M=1, d=3, mu_scale=1:3, mu_scale2=1:3, omega_scale=1:3), 36)
  test_length0(random_coefmats(d=2, how_many=1, scale=1), 4)
  test_length0(random_coefmats2(p=2, d=3), 18)
  test_length0(random_covmat(d=3, omega_scale=1:3), 6)
  test_length0(smart_covmat(d=2, Omega=diag(2), accuracy=1), 3)
})
