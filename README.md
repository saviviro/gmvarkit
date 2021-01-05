
<!-- README.md is generated from README.Rmd. Please edit that file -->

# gmvarkit

<!-- badges: start -->

<!-- badges: end -->

The goal of gmvarkit is to provide tools to analyse structural and
reduced form Gaussian mixture vector autoregressive (GMVAR) model.
`gmvarkit` provides functions for unconstrained and constrained maximum
likelihood estimation of the model parameters, quantile residual based
model diagnostics, simulation from the processes, forecasting,
estimation of generalized impulse response function, and more.

## Installation

You can install the released version of gmvarkit from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("gmvarkit")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("saviviro/gmvarkit")
```

## Example

## Simple example

This is a basic example on how to use `gmvarkit` in time series
analysis. The example data is the same that is used by Kalliovirta et
al. (2016) in their paper introducing the GMVAR model. The estimation
process is computationally demanding and takes advantage of parallel
computing. After estimating the model, it is shown by simple examples
how to conduct some further analysis.

``` r
# These examples use the data 'eurusd' which comes with the package, but in a scaled form.
data(eurusd, package="gmvarkit")
data <- cbind(10*eurusd[,1], 100*eurusd[,2])
colnames(data) <- colnames(eurusd)

## Reduced form GMVAR model ##

# Estimate a GMVAR(2,2) model: 16 estimation rounds and seeds for reproducible results
fit <- fitGMVAR(data, p=2, M=2, ncalls=16, seeds=1:16, ncores=4)
fit

# Estimate a GMVAR(2,2) model with the unconditional means of the two regimes 
# restricted to be the same:
fitm <- fitGMVAR(data, p=2, M=2, parametrization="mean", same_means=list(1:2),
                 ncalls=16, seeds=1:16, ncores=4)
fitm

# Estimate a GMVAR(2,2) model with autoregressive parameters restricted to be the same for all regimes
C_mat <- rbind(diag(2*2^2), diag(2*2^2))
fitc <- fitGMVAR(data, p=2, M=2, constraints=C_mat, ncalls=16, seeds=1:16, ncores=4)
fitc

# Test the above constraints on the AR parameters with likelihood ratio test:
LR_test(fit, fitc)

# Further information on the estimated model:
plot(fitc)
summary(fitc)
print_std_errors(fitc)
get_foc(fitc) # The first order condition
get_soc(fitc) # The second order condition (eigenvalues of approximated Hessian)
profile_logliks(fitc) # Profile log-likelihood functions

# Quantile residual diagnostics
diagnostic_plot(fitc) # type=c("all", "series", "ac", "ch", "norm")
qrt <- quantile_residual_tests(fitc, nsimu=10000)

# Simulate a sample path from the estimated process
sim <- simulateGMVAR(fitc, nsimu=100)
plot.ts(sim$sample)

# Forecast future values of the process
predict(fitc, n_ahead=10)


## Structural GMVAR model ##

# Estimate structural GMVAR(2,2) model identified with sign constraints:
W_22 <- matrix(c(1, 1, -1, 1), nrow=2, byrow=FALSE)
fit22s <- fitGMVAR(data, p=2, M=2, structural_pars=list(W=W_22),
                   ncalls=20, seeds=1:20, ncores=4)
fit22s

# Alternatively, if there are two regimes (M=2), a stuctural model can 
# be build based on the reduced form model:
fit22s_2 <- gmvar_to_sgmvar(fit)
fit22s_2

# Columns of the matrix W can be permutated and all signs in any column
# can be swapped without affecting the implied reduced form model, as 
# long as the lambda parameters are also rearranged accordingly: 
fit22s_3 <- reorder_W_columns(fit22s_2, perm=c(2, 1))
fit22s_3

fit22s_4 <- swap_W_signs(fit22s_3, which_to_swap=1)
fit22s_4

# Estimate generalized impulse response function (GIRF) with starting values
# generated from the stationary distribution of the process:
girf1 <- GIRF(fit22s, N=60, ci=c(0.95, 0.8), R1=200, R2=200)

# Estimate GIRF with starting values generated from the stationary distribution
# of the first regime:
girf2 <- GIRF(fit22s, N=60, ci=c(0.95, 0.8), init_regimes=1, R1=200, R2=200)

# Estimate GIRF with starting values given by the last p observations of the
# data:
girf3 <- GIRF(fit22s, N=60, init_values=fit22s$data, R1=1000)

# Test with Wald test whether the diagonal elements of the first AR coefficient
# matrix of the second regime are identical:
# fit22s has parameter vector of length 27 with the diagonal elements  of the
# first A-matrix of the second regime are in elements 13 and 16.
A <- matrix(c(rep(0, times=12), 1, 0, 0, -1, rep(0, times=27-16)), nrow=1, ncol=27)
c <- 0
Wald_test(fit22s, A, c)

# The same functions used in the demonstration of the reduced form model also
# work with structural models.
```

## References

  - Kalliovirta L., Meitz M. and Saikkonen P. (2016) Gaussian mixture
    vector autoregression. *Journal of Econometrics*, **192**, 485-498.
  - Kalliovirta L. and Saikkonen P. (2010) Reliable Residuals for
    Multivariate Nonlinear Time Series Models. *Unpublished Revision of
    HECER Discussion Paper No. 247*.
  - Virolainen S. 2020. Structural Gaussian mixture vector
    autoregressive model. Unpublished working paper, available as
    arXiv:2007.04713.
