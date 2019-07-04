<!-- README.md is generated from README.Rmd. Please edit that file -->
gmvarkit
========

The goal of gmvarkit is to provide tools to work with Gaussian Mixture Vector Autoregressive (GMVAR) model. Most importantly gmvarkit provides function `fitGMVAR` for two phase maximum likelihood estimation, but it also constains functions for quantile residual tests, graphical diagostics, forecasting and simulations. Also applying general linear constraints to the autoregressive parameters is supported.

Simple example
--------------

This is a basic example how to estimate a GMVAR model to an example data. The estimation process is computationally heavy and uses parallel computing.

``` r
# These examples use the data 'eurusd' which comes with the package, but in a scaled form.
data <- cbind(10*eurusd[,1], 100*eurusd[,2])
colnames(data) <- colnames(eurusd)

# GMVAR(2,2) model
fit22 <- fitGMVAR(data, p=2, M=2)
fit22

# GMVAR(2,2) model with autoregressive parameters restricted to be the same for all regimes
C_mat <- rbind(diag(2*2^2), diag(2*2^2))
fit22c <- fitGMVAR(data, p=2, M=2, constraints=C_mat)
fit22c
```

References
----------

-   Kalliovirta L., Meitz M. and Saikkonen P. (2016) Gaussian mixture vector autoregression. *Journal of Econometrics*, **192**, 485-498.
-   Kalliovirta L. and Saikkonen P. (2010) Reliable Residuals for Multivariate Nonlinear Time Series Models. *Unpublished Revision of HECER Discussion Paper No. 247*.
