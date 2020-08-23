context("Wald and LR tests")
library(gmvarkit)

data <- cbind(10*eurusd[,1], 100*eurusd[,2])

## A(M)(p)_(p)(M)(d)

# Structural GMVAR(2, 2), d=2 model identified with sign-constraints:
params222s <- c(1.428, -0.808, 1.029, 5.84, 1.314, 0.145, 0.094, 1.292,
                -0.389, -0.07, -0.109, -0.281, 1.248, 0.077, -0.04, 1.266, -0.272,
                -0.074, 0.034, -0.313, 0.903, 0.718, -0.324, 2.079, 7, 1.44, 0.742)
W_222 <- matrix(c(1, NA, -1, 1), nrow=2, byrow=FALSE)
mod222s <- GMVAR(data, p=2, M=2, params=params222s, structural_pars=list(W=W_222))

A1 <- rbind(c(rep(0, times=5), 1, rep(0, times=21)),
           c(rep(0, times=6), 1, rep(0, times=20)))
c1 <- c(0.2, 0)
wald1 <- Wald_test(mod222s, A1, c1)

test_that("Wald_test works correctly", {
  expect_equal(wald1$df, 2)
  expect_equal(wald1$test_stat, 4.147194, tolerance=1e-4)
  expect_equal(wald1$p_value, 0.1257327, tolerance=1e-4)
})


# The same model as above but with the AR parameters restricted to be the
# same in both regimes.
C_mat <- rbind(diag(2*2^2), diag(2*2^2))
params222sc <- c(1.786, 3.001, 1.031, 2.357, 1.25, 0.06, 0.036, 1.335,
                 -0.29, -0.083, -0.047, -0.356, 2.389, 1.924, -0.419, 2.469,
                 0.14, 0.768, 0.632)
mod222sc <- GMVAR(data, p=2, M=2, params=params222sc, constraints=C_mat, structural_pars=list(W=W_222))

lr1 <- LR_test(mod222s, mod222sc)

test_that("LR_test works correctly", {
  expect_equal(lr1$df, 8)
  expect_equal(lr1$test_stat, 18.91746, tolerance=1e-4)
  expect_equal(lr1$p_value, 0.01530736, tolerance=1e-4)
})
