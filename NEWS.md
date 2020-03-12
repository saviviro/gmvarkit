# gmvarkit 1.1.0

* Added a `NEWS.md` file to track changes to the package.
* New exported functions: uncond_moments, cond_moments, get_regime_means, get_regime_autocovs, get_omega_eigens, get_soc.
* simulateGMVAR now provides better tools for forecasting. This update includes non-backward compatible changes for the return values if the argument ntimes is set to be larger than one. In additional to the samples, it now returns a list containing the mixing weights and component that was used to generate each observation.
* diagnostic_plot normality plot now plots histograms instead of kernel estimates.
* In the predict method arguments "ci" and "ci_type" were changed to "pi" and "pi_type" to signify "prediction interval"" as it's more correct expression than "confidence interval". Also the default prediction method is now median, and not mean.
* Generally more functionality for conditional and unconditional moments.
* Increased numerical tolerance for considering models to be stationary to attain better numerical stability.
* Updates on documentation, fixed typos.

# gmvarkit 1.1.1

* Changed the default number of CPU cores employed by the estimation function fitGMVAR to be at most two due to CRAN policy.
* Added argument "seeds" to fitGMVAR allowing one to set the random number generator seed for each call to the genetic algorithm.
* New exported function: alt_gmvar
* Updated the vignette to include the definition of the GMVAR model

# gmvarkit 1.1.2

* New functions: 'profile_logliks' (plot profile log-likelihood functions), 'get_foc' (gradient at the estimates, remember also 'get_soc')
* Now standard errors are printed correctly for models imposing all kinds of constraints (in earlier versions standard errors for such constrained AR parameters that involved sums or multiplication were incorrect).
* Minor update on the genetic algorithm.
* Minor update on the print and summary-print methods. 
* Improvements on the comments and documentation.

# gmvarkit 1.1.3

* Updated the predict method to include forecasts for the mixing weights. Also updated the return value.
* The headlines in 'profile_logliks' are now correct for mean-parametrized models (there was phi parameter instread of mu).
* Improved the vignette. 
* The default population size in the genetic algorithm is now 50*ceiling(sqrt(npars)), was 10*npars.
* In the function quantile_residual_tests the default argument 'nsimu' is now 1 so that the tests are based on the given data only (and not on simulation).
* Corrected an error-causing bug in the predict method when the argument 'nt' was not specified.
