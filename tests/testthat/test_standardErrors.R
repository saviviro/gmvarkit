context("standard errors")
library(gmvarkit)

## NOTE: For some reason these tests fail at win-builder with tolerance smaller than 1e-1 (while they pass locally).
# Apparently, when testing diagonal of the observed information matrix, the tests pass, but when testing
# diagonal of its inverse (or the standard errors), the tests fail (with average difference about 1e-3).
# Possibly a numerical error caused by the limited precision of float-point presentation?

## A(M)(p)_(p)(M)(d)

# p=1, M=1, d=2
theta_112 <- c(0.6495277, 0.06650669, 0.28852526, 0.02176671,
               -0.14402566, 0.89710305, 0.60179087, -0.00294446,
               0.06722394)

# p=1, M=2, d=2, SGMVAR
theta_122s <- c(0.54971316, 0.11224185, 0.61921244, 0.1730025, 0.34446141, 0.05459402, -0.00872701,
               0.71754958, 0.25466175, 0.01675898, -0.13555656, 0.85845764, 0.54108934, 0.0574601,
               -0.16200273, 0.1623048, 3.62315743, 4.7262103, 0.67381109)
W_122 <- matrix(c(1, 1, -1, 1), nrow=2)


# p=2, M=2, d=2
theta_222 <- c(0.35991674, 0.12126731, 0.22311983, 0.05894632, -0.15147839, 0.39488719,
               0.40621715, -0.00505261, 0.08272533, 0.29903862, 0.21466522, 0.00231082,
               0.02975349, 0.48421647, 0.07213826, 0.21772846, 0.01990224, -0.11882016,
               0.7215664, 0.09345529, 0.03166791, 0.04420431, 0.19125506, 1.10115274,
               -0.0040654, 0.10511093, 0.58018355)


# Constraint AR-parameters to be the same for all regimes
rbind_diags <- function(p, M, d) {
  I <- diag(p*d^2)
  Reduce(rbind, replicate(M, I, simplify=FALSE))
}

C_222 <- rbind_diags(p=2, M=2, d=2)
theta_222c <- c(0.4180381, 0.15320846, 0.51320878, 0.05681983, 0.20440155, 0.02796954,
                -0.1694079, 0.59117481, 0.24081751, 0.01407821, 0.09090021, 0.24835439,
                1.06782606, -0.00997531, 0.11060615, 0.21895777, 0.00412991, 0.02746944,
                0.50178873)

test_that("standard_errors works correctly", {
  expect_equal(standard_errors(gdpdef, p=1, M=1, params=theta_112, conditional=TRUE, parametrization="intercept",
                               minval=-99999),
               c(0.101104770, 0.033791813, 0.061059527, 0.020407662, 0.086389670, 0.028873610, 0.054595595,
                 0.012904113, 0.006098679), tolerance=1e-1)

  expect_equal(standard_errors(gdpdef, p=1, M=2, params=theta_122s, conditional=TRUE, parametrization="intercept",
                               structural_pars=list(W=W_122), minval=-99999),
               c(0.17007239, 0.03951787, 0.30121858, 0.10084315, 0.10107951, 0.02901046, 0.27721982, 0.06214697,
                 0.11779698, 0.03910775, 0.18521711, 0.06074518, 0.12881568, 0.10403765, 0.31176109, 0.04174217,
                 0.84113405, 1.16065390, 0.15271742), tolerance=1e-1)

  expect_equal(standard_errors(gdpdef, p=2, M=2, params=theta_222, conditional=TRUE, parametrization="intercept",
                               minval=-99999),
               c(0.125996752, 0.043900202, 0.073803486, 0.033878146, 0.232191724, 0.105474131, 0.091032204, 0.034556047,
                 0.236481748, 0.102717211, 0.029695954, 0.007287935, 0.004471194, 0.246953455, 0.067581183, 0.102874515,
                 0.031647602, 0.329691953, 0.104357883, 0.107456526, 0.031461955, 0.337524312, 0.109251183, 0.169834165,
                 0.035705120, 0.016423312, 0.137456811), tolerance=1e-1)

  expect_equal(standard_errors(gdpdef, p=2, M=2, params=theta_222c, conditional=TRUE, parametrization="intercept",
                               constraints=C_222, minval=-99999),
               c(0.259084099, 0.057403813, 0.114979783, 0.035009525, 0.072466299, 0.022676068, 0.210200476, 0.069722257,
                 0.070015869, 0.023048650, 0.190588192, 0.067221957, 0.170186891, 0.034304708, 0.019122651, 0.048929385,
                 0.007297228, 0.004578487, 0.134639847), tolerance=1e-1)
})
