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

# gmvarkit 1.3.0

* Major speed improvement!
* Changing the numerical tolerances for stationarity of the AR parameters and positive definiteness of the error term covariance matrices is now possible everywhere expect in the main estimation function and genetic algorithm (to avoid errors in the estimation). This might help in the analysis (std errors, quantile residual tests, further iterative estimation, etc) of models with parameter values that are very close to the border of the parameter space, but too small numerical tolerance may also cause errors. 
* In the function GIRF, the argument 'variables' was renamed as 'which_shocks'.
* New exported function: "update_numtols" which enables one to update the numerical tolerances for stationarity of the AR parameters and positive definiteness of the error term covariance matrices in existing class 'gmvar' model.
* Log-likelihood function should not return -Inf anymore but with extremely bad estimates (for some t) the value is not precise (before in such cases it was -Inf). This should lead to better estimation results in some cases.
* Fixed a problem in the estimation algorithm that occurred when estimating a structural model with zero constraints in the W matrix. 
* Increased the default maxit from 300 to 500 in fitGMVAR, the default maxlag from 10 to 12 in diagnostic_plot, and decreased the default initial ar_scale from 1 to 0.2 in GAfit.
* There is now randomly varying "ar_scale" in the genetic algorithm but initial population ar_scale can still be adjusted.
* There is now possibility to set random preliminary smart mutations in the genetic algorithm.
* The package Brobdingnag is not imported anymore.
* The default lags used in quantile residual tests are now 1, 3, 6, and 12.
* Updated the documentation of fitGMVAR: more instructions on what to do when the algorithm fails to create an initial population were added.
* Some updates in the vignette, examples, and other documentation.

# gmvarkit 1.3.1 

* Bug fix which required that Brobdingnag was added back to the imported packages. There was in some rare cases a problem with the exact log-likelihood function that was introduced in the version 1.3.0 (the fix was introduced in the same day as the bug, however).

# gmvarkit 1.4.0

* Major speed improvement (when d*p > 12)!
* New feature: possibility to constrain the (unconditional) mean parameters to the same among some regimes with the argument 'same_means'. This feature is available for mean-parametrized models only (due to technical reasons).
* The default number of estimation rounds was increased to 'floor(10 + 30*log(M))'.
* Added more adjustable parameters to the genetic algorithm.
* Some changes related to estimation with the genetic algorithm (the estimation results might be different with this version, if AR or lambda constraints are employed, or in some cases when there are no overidentifying constraints in the W matrix). 
* Fixed a bug in the predict method that appeared in the very rare cases when one uses the exact one-step conditional mean without confidence intervals as the the prediction method with mean-parametrization.
* In this version, the estimation results with a given seed are different to those in the previous one due to a hidden problem that was introduced. 

# gmvarkit 1.4.1

* Fixed CRAN check issues and examples regarding LR and Wald tests.
* Fixed a problem in the estimation procedure that was introduced in the previous version.
* Updated some of the examples so that the running time for all of them is now shorter.
* The genetic algorithm now sorts regimes of the structural models by mixing weight parameters to decreasing order by redecomposing the error term covariance matrices if the first regime changes (before, only regimes 2,..,M were sorted). As a result, the MLE is now found with higher probability in each estimation round.
* In this version, the estimation results with a given seed are different to those in the previous versions due to the updates (see above).

# gmvarkit 1.4.2

* New exported function: GFEVD for estimating the generalized forecast error variance decomposition.
* Print and summary methods for gmvar objects now display the number of parameters and observations.
* Now also iterate_more and and alt_gmvar return the results from all estimation rounds.
* There is now possibility to choose not to calculate approximate standard errors when using gmvar_to_sgmvar, because it is computationally demanding for large models.
* Fixed a bug in the GIRF print method: it sometimes referred to wrong shocks when GIRF was estimated for only a subset of the shocks or the shocks were not in an increasing order.
* Fixed a bug in that caused an error in some functions when the model was adjusted to have smaller stationarity/posdef tolerance than the default one and the model was outside the default tolerance.
* Now fitGMVAR, iterate_more, and alt_gmvar warn if some regime has near-unit-roots or near-singular error term covariance matrix.
* Now the argument shock_size in GIRF should be a scalar value.
* Fixed a bug in the estimation algorithm that occurred when estimating structural models. This might have caused an error or just the model not to estimate that well. A small adjustment was also done. Because there was a change in the algorithm, the structural model estimations done with the previous versions are generally not reproducible with this version.
* Fixed the documentation for the argument ar_scale2.
* Fixed the argument "precission" in profile_logliks to "precision".
* Now fitGMVAR, GFEVD, and GIRF don't call closeAllConnections: instead, they only close the connections they opened.
* Added new datasets: usamone_prec and usamone

# gmvarkit 1.5.0

* The function GIRF now allows to scale the GIRFs of some shocks to normalize the magnitude of the instantaneous movement of some variable.
* Changed the default ncalls in fitGMVAR to 100. 
* Updated the plot method for class girf objects.
* Re-used the class 'htest' for the objects returned by the functions Wald_test and LR_test.
* Removed confusing NAs from the standard error prints when the related statistics are not parametrized in the first place. 
* The function GMVAR now throws an error if there are more parameters in the model than d*nrow(data).
* Internal functions are now removed from the user manual. 
* Added the data gdpdef, updated the examples for this data, and removed the old data.

# gmvarkit 2.0.0

* This update should be mostly backward compatible, but some changes made in the argument names are not.
* gmvarkit now accommodates new models: the StMVAR model and the G-StMVAR model, as well their structural versions (see the vignette or the references in the package description). 
* Changed the model class from 'gmvar' to 'gsmvar' to accommodate also StMVAR and G-StMVAR models.
* Renamed functions: GMVAR -> GSMVAR; fitGMVAR -> fitGSMVAR; alt_gmvar -> alt_gsmvar; gmvar_to_sgmvar -> gsmvar_to_sgsmvars; also all the class 'gmvar' methods were changed to class 'gsmvar' methods.
* New improved vignette.
* The old simulation function 'simulateGMVAR' is now deprecated. Now, we use the class 'gsmvar' simulation method 'simulate.gsmvar' instead.
* simulate.gsmvar now allows to generate initial values from the stationary distribution of a specific regime or from a mixture distribution of any set of regimes. Some argument names needed to be changed to make the method CRAN compatible.
* New exported function 'stmvar_to_gstmvar': estimate a G-StMVAR model based on a StMVAR model with large degrees of freedom parameters.
* Now the function 'quantile_residual_tests' takes use of parallel computing to shorten the computation time.
* Changes to the default arguments 'M' and 'maxit' of the function  'fitGSMVAR'
* In the method predict.gsmvar, changed the argument name n_simu to nsim, since nsim is now used in the simulation method as well.
* Changed the argument name nsimu to nsim in quantile_residual_tests as well, so that the argument name is the same as in the predict and simulate methods (while this specific name was required by CRAN compatibility in the simulation method). 
* Fixed a bug in the argument scale of the function GIRF. Fixed also a bug that caused an error in some cases.
* Now the plot method for class 'gsmvar' objects enables to plot only a single figure: series plot or density plot.
* Note that this version might produce different results with the same random number generator seed than the previous versions.

# gmvarkit 2.0.1

* Now GIRF also returns all the girfs from each Monte Carlo repetition.
* Added the data "usamone" that is used in the cited paper by Virolainen introducing the structural GMVAR model.
* Updated the references.
* Updated the vignette.

# gmvarkit 2.0.2

* We do not recycle the class "htest" for Wald and LR tests anymore because the print method was buggy and beyond the control of the author. Now Wald and LR tests use the class "hypotest" with a dedicated print method.
* Changed the default ncalls in fitGSMVAR to (M + 1)^5.
* Fixed a bug that sometimes caused an error in the estimation algorithm because the random degrees of freedom parameters were too close to their strict lower bound 2. This change may have an effect on the estimation results of StMVAR and G-StMVAR models with specific seed.
* Adjusted the tolerance when the estimation functions warns about near-non-stationary estimates or near-singular error term covariance matrices.
* Updated the vignette.

# gmvarkit 2.0.3

* Fixed a bug that caused an error when trying to estimate model a with restrictions on the mean parameters or a G-StMVAR model by setting up initial population in the genetic algorithm.
* Fixed bug that made some of the functions unable to calculate first and second differences for structural models
* Added new functionality to the function gsmvar_to_sgsmvar: now it is possible to employ Cholesky identification if M=1.
* Added new function: estimate_sgsmvar which can be used to conveniently estimate overidentified structural models.
* Updated the data 'usamone' to contain the observations until the end of 2021.
* Added the data 'euromone' that is used in the (G-)StMVAR paper of Virolainen.
* Fixed CRAN issues.

# gvmarkit 2.0.4

* Updated the references

# gmvarkit 2.0.5

* Fixed a CRAN issue

# gmvarkit 2.0.6

* Fixed a CRAN issue

# gmvarkit 2.0.7

* Fixed a bug in simulations from Student's t regimes
* Updated the print method of gsmvarpred objects

# gmvarkit 2.0.8

* Fixed a bug in simulations from Student's t regimes

# gmvarkit 2.0.9

* Added the function 'Rao_test' implementing the Rao's score test.

# gmvarkit 2.0.10

* Fixed the package overview help file and removed redundant argument from the documentation of warn_eigens. 

# gmvarkit 2.1.0

* Added support for recursively identified structural GSMVAR models.
* Added the function linear_IRF to calculate linear impulse response functions based on a single regime. Bootstrapped confidence bounds are also available for models that impose linear autoregressive dynamics. 
* Added a possibility to constrain the mixing weight parameters alphas to fixed constants.
* Added a possibility to constrain the lambda parameters of structural models to fixed constant.
* Added the function Pearson_residuals to calculate the standardized Pearson residuals or the raw residuals. 
* Added the dataset 'usamon'.
* fitGSMVAR can now filter out (many of the) inappropriate estimates by setting filter_estimates=TRUE. 
* fitGSMVAR can now estimate without parallel computing and printout by specifying use_parallel=FALSE.
* Updated the vignette (e.g., added a description of the Monte Carlo algorithm implemented to estimate GIRFs). 
* Fixed the .Rd files for random_df and smart_df.
* Fixed a bug in the estimation of GIRF and GFEVD for StMVAR and G-StMVAR models 
* Fixed a bug that occurred using stmvar_to_gstmvar for structural models.
* Fixed a bug that occurred when using swap_W_signs with C_lambda constraints.
* Fixed a bug in Rao test when using G-StMVAR model
* Fixed a bug in get_foc, get_gradient, get_soc, and get_hessian that returned only NAs for models with same_means.

# gmvarkit 2.1.1

* Bug fix for estimate_sgmvar for the case when relaxing a zero constraints.
* Slightly adjusted the setting of estimate_sgmvar.
* Added a strict upper bound of one for ar_scale to avoid problems caused by numerical inaccuracies caused by imprecise machine accuracy.
* Added note to the GAfit's documentation for using ar_scale with large p or d. Also bounded upper_ar_scale by 1-p*d/150 when p*d>40 to avoid numerical issues (but it does not go below 0.05). 

# gmvarkit 2.1.2

* fitGSMVAR now automatically filters out inappropriate estimates by default. 
* Updated the reference structural models identified by heterokedasticity.
* Fixed a lot of typos etc in the vignette. 
* Updates to readme. 

# gmvarkit 2.1.3

* Added table of contents to the vignette.
* Fixed typos etc in the vignette.
