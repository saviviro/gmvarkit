context("linear_IRF")
library(gmvarkit)

## GMVAR, p=1, M=2, d=2 model with linear AR dynamics.
# recursive identification, IRF based on the first regime:
params12cm <- c(0.805, 0.546, 0.304, 0.024, -0.132, 0.852, 0.287,
                0.005, 0.025, 1.12, -0.017, 0.143, 0.647)
mod12cm <- GSMVAR(data=gdpdef, p=1, M=2, params=params12cm,
                  constraints=rbind(diag(1*2^2), diag(1*2^2)),
                  same_means=list(1:2), parametrization="mean")

# Identification by heteroskedasticity, bootstrapped confidence intervals and
# and scaled instantaneous effects of the shocks. Note that in actual
# empirical application, a larger number of bootstrap reps should be used.
mod12cms <- gsmvar_to_sgsmvar(mod12cm)

test_that("linear_IRF works correctly", {
  irf1 <- linear_IRF(mod12cm, regime=1, N=3, scale=cbind(c(1, 1, 1), c(2, 2, 1)),
                     which_cumulative=2)
  irf2 <- linear_IRF(mod12cms, regime=1, N=3, ci=0.68, bootstrap_reps=2,
                     ncalls=1, seeds=1:2, ncores=1)

  expect_equal(c(irf1$point_est[1:2, 1:2, 4]), c(0.02099899, 0.13304379, -0.14178912, 3.18684486), tolerance=1e-3)
  expect_equal(c(irf2$point_est[1:2, 1:2, 1]), c(0.11472142, -0.15217808, 0.52329628, 0.04291658), tolerance=1e-3)
  expect_equal(dim(irf2$conf_int), c(2, 2, 4, 2), tolerance=1e-3)
  # expect_equal((irf2$conf_int[1:2, 1:2, 4, 1:2]), c(-0.003625922, -0.041964108, -0.002235934, 0.006063775,
  #                                                    0.002667661, 0.050257831, -0.001106741, 0.021349308), tolerance=1e-3)
  # Variation in small numerical errror variation across different machines may result in different estimation and
  # thus bootstrap results

})
