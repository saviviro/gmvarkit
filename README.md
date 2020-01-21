
<!-- README.md is generated from README.Rmd. Please edit that file -->

# gmvarkit

<!-- badges: start -->

<!-- badges: end -->

The goal of gmvarkit is to provide tools to analyse the Gaussian mixture
vector autoregressive (GMVAR) model. `gmvarkit` provides functions for
unconstrained and constrained maximum likelihood estimation of the model
parameters, quantile residual based model diagnostics, simulation from
the processes, and forecasting.

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

This is a basic example how to estimate a GMVAR model to data. The
example data is the same that is used by Kalliovirta et al. (2016) in
their paper introducing the GMVAR model. The estimation process is
computationally demanding and takes advantage of parallel computing.
After estimating the model, it’s shown by simple examples how to conduct
some further
analysis.

``` r
# These examples use the data 'eurusd' which comes with the package, but in a scaled form.
data(eurusd, package="gmvarkit")
data <- cbind(10*eurusd[,1], 100*eurusd[,2])
colnames(data) <- colnames(eurusd)

# Estimate a GMVAR(2,2) model: 20 estimation rounds and seeds for reproducible results
fit <- fitGMVAR(data, p=2, M=2, ncalls=20, seeds=1:20)
fit

# Estimate a GMVAR(2,2) model with autoregressive parameters restricted to be the same for all regimes
C_mat <- rbind(diag(2*2^2), diag(2*2^2))
fitc <- fitGMVAR(data, p=2, M=2, constraints=C_mat, ncalls=20, seeds=21:40)
fitc

# Further information on the estimated model:
plot(fitc)
summary(fitc)
print_std_errors(fitc)
get_foc(fitc) # The first order condition
get_soc(fitc) # The second order condition (eigenvalues of approximated Hessian)
profile_logliks(fitc) # Profile log-likelihood functions

# Quantile residual diagnostics
diagnostic_plot(fitc, type="series") # type=c("series", "ac", "ch", "norm")
qrt <- quantile_residual_tests(fitc)

# Simulate a sample path form the estimated process
sim <- simulateGMVAR(fitc, nsimu=10)

# Forecast future values of the process
predict(fitc, n_ahead=10)
```

## References

  - Kalliovirta L., Meitz M. and Saikkonen P. (2016) Gaussian mixture
    vector autoregression. *Journal of Econometrics*, **192**, 485-498.
  - Kalliovirta L. and Saikkonen P. (2010) Reliable Residuals for
    Multivariate Nonlinear Time Series Models. *Unpublished Revision of
    HECER Discussion Paper No. 247*.
