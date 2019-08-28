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
