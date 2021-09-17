context("Wald and LR tests")
library(gmvarkit)

## A(M)(p)_(p)(M)(d)

# Structural GMVAR(2, 2), d=2 model identified with sign-constraints:
params222s <- c(0.36, 0.121, 0.484, 0.072, 0.223, 0.059, -0.151, 0.395, 0.406,
                -0.005, 0.083, 0.299, 0.218, 0.02, -0.119, 0.722, 0.093, 0.032,
                0.044, 0.191, 0.057, 0.172, -0.46, 0.016, 3.518, 5.154, 0.58)
W_222 <- matrix(c(1, 1, -1, 1), nrow=2, byrow=FALSE)
mod222s <- GMVAR(gdpdef, p=2, M=2, params=params222s, structural_pars=list(W=W_222))

A1 <- rbind(c(rep(0, times=5), 1, rep(0, times=21)),
           c(rep(0, times=6), 1, rep(0, times=20)))
c1 <- c(0.1, 0)
wald1 <- Wald_test(mod222s, A1, c1)

test_that("Wald_test works correctly", {
  expect_equal(unname(wald1$parameter), 2)
  expect_equal(unname(wald1$statistic), 1.818517, tolerance=1e-1)
  expect_equal(wald1$p.value, 0.4028229, tolerance=1e-1)
})


# The same model as above but with the AR parameters restricted to be the
# same in both regimes.
C_mat <- rbind(diag(2*2^2), diag(2*2^2))
params222sc <- c(0.418, 0.153, 0.513, 0.057, 0.204, 0.028, -0.169, 0.591, 0.241,
                 0.014, 0.091, 0.248, 0.345, 0.31, -0.974, 0.12, 0.256, 0.199, 0.501)
mod222sc <- GMVAR(gdpdef, p=2, M=2, params=params222sc, constraints=C_mat, structural_pars=list(W=W_222))

lr1 <- LR_test(mod222s, mod222sc)

test_that("LR_test works correctly", {
  expect_equal(unname(lr1$parameter), 8)
  expect_equal(unname(lr1$statistic), 14.18665, tolerance=1e-4)
  expect_equal(lr1$p.value, 0.07702821, tolerance=1e-4)
})
