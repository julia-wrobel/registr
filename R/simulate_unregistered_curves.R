#' Simulate mean curve
#' 
#' This function generates mean for simulated accelerometer data.
#' 
#' @param grid Grid of x values over which to evaluate the function.
mean_curve = function(grid) {
	1.5 * (0 - sin(2*grid*pi) - cos(2*grid*pi) )  
}

#' Simulate amplitude variance
#' 
#' This function generates amplitudes for simulated accelerometer data.
#' 
#' @param grid Grid of x values over which to evaluate the function.
amp_curve = function(grid) {
	(0 - sin(2*grid*pi) - cos(2*grid*pi) )  / sqrt(322)
}

#' Generate subject-specific grid (t_star)
#' 
#' This function creates subject-specific time grid
#' 
#' @param coefs Spline basis coefficients for reconstructing the subject-specific grid. 
#' @param D Number of grid points per subject.
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
#' @param seed Seed for reprodicibility. Default is 1988.
#' @param vary_D Indicates if grid length vary by subject. If FALSE all subjects have grid length D.
#' 
#' @author Jeff Goldsmith \email{ajg2202@@cumc.columbia.edu}
#' @importFrom magrittr %>%
#' @importFrom stats rbinom runif
#' @importFrom boot inv.logit
#' @importFrom mvtnorm rmvnorm
#' 
#' @export
#'
simulate_unregistered_curves = function(I = 50, D = 100, lambda = 15, seed = 1988) {
	
	grid = seq(0, 1, length = D)
	
	set.seed(seed)
	
	## generate features for curves
	c_true = rmvnorm(I, mean = rep(0, 1), sigma = diag(lambda, 1, 1))
	
	## generate curves
	Yi_obs = Yi_latent = pi_true = t_subj = matrix(NA, I, D)
	Yi_regis.true = matrix(NA, I, D)
	for (i in 1:I) {
		t_subj[i,] = runif(3, 0, 1) %>% grid_subj_create(., D = D) %>% as.vector
		Yi_latent[i,] = mean_curve(grid = t_subj[i,]) + c_true[i] * amp_curve(grid = t_subj[i,])
		pi_true[i,] = inv.logit(Yi_latent[i,])
		Yi_regis.true[i,] = mean_curve(grid = grid) + c_true[i] * amp_curve(grid = grid)
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
		tstar = tstar_vec,
		t = t_vec,
		value = as.vector(t(Yi_obs)),
		latent_prob = as.vector(t(Yi_latent))
	)
	
	simulated_data
}
