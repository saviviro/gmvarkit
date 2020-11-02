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

* New exported functions: 'profile_logliks' (plot profile log-likelihood functions), 'get_foc' (gradient at the estimates, recall also 'get_soc')
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

# gmvarkit 1.2.0

* A structural GMVAR model is introduced! The SGMVAR model can be estimated with the function 'fitGMVAR' or constructed by hand with the function 'GMVAR', while also other existing functions such as 'profile_logliks', 'quantile_residual_tests', 'diagnostic_plot', etc. work with the SGMVAR model as well. See the cited paper by Virolainen (2020) or the vignette for the model definition and identification conditions.
* New exported function: 'GIRF' for estimating generalized impulse function for the variables in a SGMVAR model.
* New exported function: 'Wald_test' for conducting a Wald test testing validity of parameter constraints.
* New exported function: 'LR_test' for conducting a likelihood ratio test testing validity of parameter constraints.
* New exported function: 'gmvar_to_sgmvar' for building a structural GMVAR model based on a two-regime reduced form GMVAR model.
* New exported function: 'reorderd_W_column' for reordering the columns of W matrix of a structural GMVAR model.
* New exported function: 'swap_W_signs' for swapping all signs the pointed columns of W matrix (and also B-matrix) of a structural GMVAR model.
* Bug fix: the prediction intervals for mixing weights were incorrect when calculating upper or lower prediction intervals with only one level of significance.
* Minor computation speed improvements.
* Non-backward-compatible change: the functions 'get_boldA_eigens' and 'get_omega_eigens' now return matrices and not lists.
* Included the structural GMVAR model in a vignette. 

# gmvarkit 1.2.1

* Fixed CRAN issues concerning certain unit tests (not visible to users).

# gmvarkit 1.2.2

* The plot method of class gmvar objects now also plots marginal stationary densities along with kernel density estimates of the individual time series.
* New exported function: 'cond_moment_plot' for plotting in-sample conditional means and variances
* Changed 'diagnostic_plot' to automatically show all four figures by default and improved the quantile residual time series plot. 
* Added the argument 'which_largest' to the function 'alt_gmvar'
* In GIRF the default shock size is now one, which is then amplified by the B-matrix according to the conditional standard deviation of the model.
* Fixed a bug that caused error when trying to estimate a structural model with more than two regimes and no zero restrictions.
* Fixed a bug that caused error when estimating GIRF with only one initial regime that is not the first regime.
* Fixed a bug that caused error when printing standard errors of a SGMVAR model with constraints on the lambda parameters.

# gmvarkit 1.2.3

* The function GIRF now enables to directly calculate cumulative impulse responses.
* Updated some of the examples, readme, and vignette.
* Fixed a bug that caused error in estimation of GIRF with very large shock size.
* Fixed a bug that caused error in estimation of AR constrained models when initial population is used in the genetic algorithm.
