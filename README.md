
<!-- README.md is generated from README.Rmd. Please edit that file -->

# gmvarkit

<!-- badges: start -->
<!-- badges: end -->

The goal of gmvarkit is to provide tools to analyse structural and
reduced form Gaussian mixture vector autoregressive (GMVAR) model,
Student’s t mixtue vector autoregressive (StMVAR) model, and Gaussian
and Student’s t mixture vector autoregressive (G-StMVAR) model.
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
analysis. The estimation process is computationally demanding and takes
advantage of parallel computing. After estimating the model, it is shown
by simple examples how to conduct some further analysis. StMVAR and
G-StMVAR model are briefly covered after giving examples related to the
reduced form and structural GMVAR model.

``` r
# These examples use the data 'gdpdef' which comes with the package, and contains the quarterly percentage growth rate
# of real U.S. GDP and quarterly percentage growth rate of U.S. GDP implicit price deflator, covering the period 
# from 1959Q1 to 2019Q4.
data(gdpdef, package="gmvarkit")

## Reduced form GMVAR model ##

# Estimate a GMVAR(2, 2) model: 20 estimation rounds and seeds for reproducible results
# (note: many empirical applications require more estimation rounds, e.g., hundreds
# or thousands).
fit <- fitGSMVAR(gdpdef, p=2, M=2, ncalls=20, seeds=1:20, ncores=2)
fit

# Estimate a GSMVAR(2, 2) model with autoregressive parameters restricted to be the same for all regimes
C_mat <- rbind(diag(2*2^2), diag(2*2^2))
fitc <- fitGMVAR(gdpdef, p=2, M=2, constraints=C_mat, ncalls=20, seeds=1:20, ncores=2)
fitc

# Estimate a GMVAR(2, 2) model with autoregressive parameters and the unconditional means
# restricted to be the same in both regimes (only the covariance matrix varies)
fitcm <- fitGSMVAR(gdpdef, p=2, M=2, parametrization="mean", constraints=C_mat, same_means=list(1:2),
                  ncalls=20, seeds=1:20, ncores=2)
fitcm 

# Test the above constraints on the AR parameters with likelihood ratio test:
LR_test(fit, fitc)

# Further information on the estimated model:
plot(fit)
summary(fit)
print_std_errors(fit)
get_foc(fit) # The first order condition (gradient of the log-likelihood function)
get_soc(fit) # The second order condition (eigenvalues of approximated Hessian)
profile_logliks(fit) # Profile log-likelihood functions

# Note: models can be built based the results from any estimation round 
# conveniently with the function 'alt_gmvar'.

# Quantile residual diagnostics
diagnostic_plot(fit) # type=c("all", "series", "ac", "ch", "norm")
qrt <- quantile_residual_tests(fit, nsim=10000)

# Simulate a sample path from the estimated process
sim <- simulate(fit, nsim=100)
plot.ts(sim$sample)

# Forecast future values of the process
predict(fit, n_ahead=12)


## Structural GMVAR model ##

# gmvarkit implements two identification methods: recursive identification and
# identification by heteroskedasticity. Reduced form models can be directly used
# as recursively identified structural models. The below examples consider
# identification by heteroskedasticity. 

# Estimate structural GMVAR(2,2) model identified with sign constraints:
W22 <- matrix(c(1, 1, -1, 1), nrow=2, byrow=FALSE)
fit22s <- fitGSMVAR(gdpdef, p=2, M=2, structural_pars=list(W=W22),
                    ncalls=20, seeds=1:20, ncores=2)
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

fit22s_4 <- swap_W_signs(fit22s_3, which_to_swap=2)
fit22s_4 # The same model as fit22s

all.equal(fit22s$loglik, fit$loglik)
all.equal(fit22s$loglik, fit22s_2$loglik)
all.equal(fit22s$loglik, fit22s_3$loglik)
all.equal(fit22s$loglik, fit22s_4$loglik)

# Check out also the function: estimate_sgsmvar
# for estimating overidentified structural models.

# Estimate generalized impulse response function (GIRF) with starting values
# generated from the stationary distribution of the process:
girf1 <- GIRF(fit22s, N=20, ci=c(0.95, 0.8), R1=200, R2=200, ncores=2)

# Estimate GIRF with starting values generated from the stationary distribution
# of the first regime:
girf2 <- GIRF(fit22s, N=20, ci=c(0.95, 0.8), init_regimes=1, R1=200, R2=200, ncores=2)

# Estimate GIRF with starting values given by the last p observations of the
# data:
girf3 <- GIRF(fit22s, N=20, init_values=fit22s$data, R1=1000, ncores=2)

# Estimate generalized forecast error variance decmposition (GFEVD) with the
# initial values being all possible lenght p the histories in the data:
gfevd1 <- GFEVD(fit22s, N=20, R1=100, initval_type="data", ncores=2)
plot(gfevd1)

# Estimate GFEVD with the initial values generated from the stationary
# distribution of the second regime:
gfevd2 <- GFEVD(fit22s, N=20, R1=100, R2=100, initval_type="random",
                init_regimes=2, ncores=2)
plot(gfevd2)

# Estimate GFEVD with fixed starting values that are the unconditional mean
# of the process: 
myvals <- rbind(fit22s$uncond_moments$uncond_mean,
                fit22s$uncond_moments$uncond_mean)
gfevd3 <- GFEVD(fit22s, N=48, R1=250, initval_type="fixed",
                init_values=myvals, ncores=2)
plot(gfevd3)

# Check also the function linear_IRF for estimating linear impulse response 
# functions based on a specific regime.

# Test with Wald test whether the diagonal elements of the first AR coefficient
# matrix of the second regime are identical:
# fit22s has parameter vector of length 27 with the diagonal elements  of the
# first A-matrix of the second regime are in elements 13 and 16.
A <- matrix(c(rep(0, times=12), 1, 0, 0, -1, rep(0, times=27-16)), nrow=1, ncol=27)
c <- 0
Wald_test(fit22s, A, c)

# The same functions used in the demonstration of the reduced form model also
# work with structural models.


## StMVAR and G-StMVAR models ##

# Fit a StMVAR(2, 2) model
fit22t <- fitGSMVAR(gdpdef, p=2, M=2, model="StMVAR", ncalls=1, seeds=1)

# Printout shows that there is an overly large degrees of freedom estimate
# in the 2nd regime!
fit22t

# So we change it to GMVAR type by switching to the appropariate G-StMVAR model
fit22gs <- stmvar_to_gstmvar(fit22t) 

# Printout of the appropriate G-StMVAR model that is based on the above StMVAR model.
fit22gs 

# Switch to the stastitically identified structural model
fit22gss <- gsmvar_to_sgsmvar(fit22gs)

# Printout of the structural model that is implied by the reduced form model
fit22gss
```

## References

- Kalliovirta L., Meitz M. and Saikkonen P. (2016) Gaussian mixture
  vector autoregression. *Journal of Econometrics*, **192**, 485-498.
- Kalliovirta L. and Saikkonen P. (2010) Reliable Residuals for
  Multivariate Nonlinear Time Series Models. *Unpublished Revision of
  HECER Discussion Paper No. 247*.
- Virolainen S. (forthcoming). A statistically identified structural
  vector autoregression with endogenously switching volatility regime.
  *Journal of Business & Economic Statistics*,
  <doi:10.1080/07350015.2024.2322090>.
- Virolainen S. (2022). Gaussian and Student’s t mixture vector
  autoregressive model with application to the assymetric effects of
  monetary policy shocks in the Euro area. Unpublished working paper,
  available as arXiv:2109.13648.
