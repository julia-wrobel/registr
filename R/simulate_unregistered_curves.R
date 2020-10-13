#' Simulate mean curve
#' 
#' This function generates mean for simulated accelerometer data.
#' 
#' @param grid Grid of x values over which to evaluate the function.
#' @param period Controls the period of the mean curve
#' @param spline_based If FALSE curve is constructed using sine and cosine functions,
#' if TRUE, curve is constructed using B-spline basis.
#' 
#' @return A numeric vector.
#' 
mean_curve = function(grid, period = 2*pi, spline_based = FALSE) {
	if(spline_based){
		
		# define mean on evenly spaced grid and return interpolation for subject grid
		tstar = seq(0, 1, length.out = length(grid))
		BS = splines::bs(tstar, df = 9, intercept = TRUE, degree = 3)
		coefs = c(-2, -2, -1, 3, 1, -3, 2, 2 , -2)
		mean = (BS %*% coefs)
		
		approx(tstar, mean, xout = grid)$y
		
	}else{
		-1.5 * (sin(period*grid) + cos(2*pi*grid) )  
	}
}

#' Simulate amplitude variance
#' 
#' This function generates amplitudes for simulated accelerometer data.
#' 
#' @importFrom stats approx
#' 
#' @param grid Grid of x values over which to evaluate the function.
#' @param period Controls the period of the mean curve
#' @param spline_based If FALSE curve is constructed using sine and cosine functions,
#' if TRUE, curve is constructed using B-spline basis.
#' 
#' @return A numeric vector.
amp_curve = function(grid, period = 2*pi, spline_based = FALSE) {
	if(spline_based) {

		tstar = seq(0, 1, length.out = length(grid))		
		BS = splines::bs(tstar, df = 9, intercept = TRUE, degree = 3)
		amp_coefs = c(0, 0, 0, 0.5, 0.5, 0, -0.5, -0.5, 0)
		amp = (BS %*% amp_coefs)
		
		approx(tstar, amp, xout = grid)$y
		
	}else{
		-(sin(2*pi*grid) + cos(period*grid) )  / sqrt(322)
	}
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
#' @param period Controls the period of the mean curve
#' @param spline_based If FALSE curve is constructed using sine and cosine functions,
#' if TRUE, curve is constructed using B-spline basis.
#' @param phase_variation If TRUE, creates phase variation 
#' (registered curves are observed on uneven grid). If FALSE, no phase variation.
#' 
#' @author Julia Wrobel \email{julia.wrobel@@cuanschutz.edu},
#' Jeff Goldsmith \email{ajg2202@@cumc.columbia.edu}
#' @importFrom magrittr %>%
#' @importFrom stats rbinom runif plogis
#' 
#' @return A simulated dataframe with variables id, value, index, latent_mean, and t. Index is the domain
#' on which curves are unregistered and t is the domain on which curves are registered.
#' @export
#'
simulate_unregistered_curves = function(I = 50, D = 100, lambda = 15, seed = 1988,
																				period = 2 * pi, spline_based = FALSE, 
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
		Yi_latent[i,] = mean_curve(grid = t_subj[i,], 
															 period = period, spline_based = spline_based) + 
			c_true[i] * amp_curve(grid = t_subj[i,], 
														period = period, spline_based = spline_based)
		pi_true[i,] = plogis(Yi_latent[i,])
		Yi_regis.true[i,] = mean_curve(grid = grid, 
																	 period = period, spline_based = spline_based) + 
			c_true[i] * amp_curve(grid = grid, period = period, spline_based = spline_based)
		for (j in 1:D) {
			Yi_obs[i,j] = rbinom(1, 1, pi_true[i,j])
		}
	}
	
	t_vec = as.vector(t(t_subj))
	tstar_vec = rep(grid, I)
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
