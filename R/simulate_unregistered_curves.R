#' Simulate mean curve
#' 
#' This function generates mean for simulated accelerometer data.
#' 
#' @param grid Grid of x values over which to evaluate the function.
#' 
#' @return A numeric vector.
#' 
mean_curve = function(grid) {
	# define mean on evenly spaced grid and return interpolation for subject grid
	tstar = seq(0, 1, length.out = length(grid))
	BS = splines::bs(tstar, df = 9, intercept = TRUE, degree = 3)
	coefs = c(-2, -2, -1, 3, 1, -3, 2, 2 , -2)
	mean = (BS %*% coefs)
	
	approx(tstar, mean, xout = grid)$y
}

#' Simulate amplitude variance
#' 
#' This function generates amplitudes for simulated accelerometer data.
#' 
#' @importFrom stats approx
#' 
#' @param grid Grid of x values over which to evaluate the function.
#' 
#' @return A numeric vector.
amp_curve = function(grid) {
	tstar = seq(0, 1, length.out = length(grid))		
	BS = splines::bs(tstar, df = 9, intercept = TRUE, degree = 3)
	amp_coefs = c(.1, .1, .1, 0.5, 0.5, 0, -0.5, -0.5, .05)
	amp = (BS %*% amp_coefs)
	
	approx(tstar, amp, xout = grid)$y
}

#' Generate subject-specific grid (t_star)
#' 
#' This function creates subject-specific time grid
#' 
#' @importFrom stats approx
#' 
#' @param coefs Spline basis coefficients for reconstructing the subject-specific grid. 
#' @param D Number of grid points per subject.
#' @return A numeric vector.
grid_subj_create = function(coefs, D) {
	BS = splines::bs(1:D, df = 3, intercept = FALSE, degree = 3)
	coefs = cumsum(coefs) / sum(coefs) 
	(BS %*% coefs)
}

#' Simulate unregistered curves
#' 
#' This function simulates unregistered curves, providing the time values for both 
#' the unregistered curves (t_star) and the registered curves (t). Curves all have one peak, the location
#' of which is shifted on the unregistered domain, meant to mimic accelerometer data.  
#'  
#' @param I Number of subjects. Defaults is 50.
#' @param D Number of grid points per subject. Default is 100.
#' @param lambda Standard deviation for subject-specific amplitudes.
#' @param seed Seed for reproducibility. Default is 1988.
#' @param phase_variation If TRUE, creates phase variation 
#' (registered curves are observed on uneven grid). If FALSE, no phase variation.
#' 
#' @author Julia Wrobel \email{jw3134@@cumc.columbia.edu},
#' Jeff Goldsmith \email{ajg2202@@cumc.columbia.edu}
#' @importFrom magrittr %>%
#' @importFrom stats rbinom runif plogis
#' 
#' @return A simulated dataframe with variables id, value, index, latent_mean, and t. Index is the domain
#' on which curves are unregistered and t is the domain on which curves are registered.
#' @export
#'
simulate_unregistered_curves = function(I = 50, D = 100, lambda = 1, seed = 1988,
																				phase_variation = TRUE) {
	
	## NULLify global values called by tidyverse functions
	value = index = NULL
	
	grid = seq(0, 1, length = D)
	
	set.seed(seed)
	
	## generate features for curves
	c_true = rnorm(I, mean = 0, sd = sqrt(lambda)) 
	
	## generate curves
	Yi_obs = Yi_latent = pi_true = t_subj = matrix(NA, I, D)
	Yi_regis.true = matrix(NA, I, D)
	for (i in 1:I) {
		if(phase_variation){
			t_subj[i,] =  grid_subj_create(runif(3, 0, 1), D = D) %>% as.vector
		}else{
			t_subj[i,] =  grid
		}
		Yi_latent[i,] = mean_curve(grid = t_subj[i,]) + 
			c_true[i] * amp_curve(grid = t_subj[i,])
		pi_true[i,] = plogis(Yi_latent[i,])
		Yi_regis.true[i,] = mean_curve(grid = grid) + 
			c_true[i] * amp_curve(grid = grid)
		for (j in 1:D) {
			Yi_obs[i,j] = rbinom(1, 1, pi_true[i,j])
		}
	}
	
	t_vec = as.vector(t(t_subj))
	tstar_vec = rep(grid, I)
	#subject = rep(1:I, each = D * 4)
	#visit = rep(rep(1:4, each = D), I)
	subject = rep(1:I, each = D)
	
	simulated_data = data.frame(
		id = subject,
		index = tstar_vec,
		value = as.vector(t(Yi_obs)),
		latent_mean = as.vector(t(Yi_latent)),
		t = t_vec
	)
	

	
	simulated_data
}
