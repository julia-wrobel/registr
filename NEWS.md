# registr 0.1.0

* Added a `NEWS.md` file to track changes to the package.


# registr 1.0.0

* Preparing for first release to CRAN
* Added Erin McDonnell as an author


# registr 1.1.0

* Erin McDonnel added non-parametric updates
* Updated vignette as well
* Use of periodic b-splines
* get her to implement gradient for parametric warping stuff

# registr 2.0.0

* registr is now able to handle incomplete curves (see new vignette)
* Added 'two-step GFPCA' as alternative GFPCA approach
* New plot function for GFPCA results
* Parallelized the registration call
* Added user-specified template functions by argument 'Y_template'
* Added Gamma and Poisson family
* Added Alexander Bauer an an author

# registr 2.1.0

* Changed convergence threshold in register_fpca from 0.00001 to 0.0001 as a more reasonable threshold in many data situations.
* Improved robustness of 'cov_hall' by (i) first centering the curves before taking the cross product and smoothing it and by (ii) ensuring positive (semi-)definiteness of the covariance matrix with 'Matrix::nearPD'.
* Improved speed of 'cov_hall' for large data by using 'mgcv::bam' for smoothing.
* Added argument 'npc_varExplained' to functions 'fpca_gauss' and 'bfpca' to choose the number of FPCs based on the share of explained variation.
* All 'verbose' arguments now take values between 0-4 to better control the level of detail of info messages.
* A lot of minor fixes and refinements.
* Updated potentially invalid URLs


Preparing things below

* mgcv spline choices
* penalization for optimization (not sure how to choose best sigma)

# registr 2.1.1

* Bug fixes to address source code pointing to old GCC compiler
* Updated email address of maintainer
* Minor documentation fixes

# registr 2.2.0

* New submission after package was archived