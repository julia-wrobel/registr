#' 
#' Simulate mean
#' 
#' This function generates mean for simulated functional data.
#' 
#' @param grid Grid of x values over which to evaluate the function.
#' 
mean_sim = function(grid){
	sin(grid * 3 * pi) 
}


#' 
#' Simulate PC1
#' 
#' This function generates the first principal component for simulated functional data.
#' 
#' @param grid Grid of x values over which to evaluate the function.
#'
psi1_sim = function(grid){
	sin(grid * 3 * pi)
}


#' 
#' Simulate PC2
#' 
#' This function generates the second principal component for simulated functional data.
#' 
#' @param grid Grid of x values over which to evaluate the function.
#'
psi2_sim = function(grid){
	cos(grid * 3 * pi)
}


#' 
#' Simulate functional data
#' 
#' This function simulates functional data. The data it outputs is generated from a mean function
#' and two orthogonal principal component basis functions. The mean and principal components are 
#' based on sine and cosine functions. Subject-specific scores for each PC are drawn from normal 
#' distributions with standard deviation lambda1 and lambda2.
#' 
#' @param lambda1 Standard deviation for PC1 scores.
#' @param lambda2 Standard deviation for PC2 scores.
#' @param I Number of subjects. Defaults is 50.
#' @param D Number of grid points per subject. Default is 100.
#' @param seed Seed for reprodicibility. Default is 1988.
#' 
#' @author Julia Wrobel \email{jw3134@@cumc.columbia.edu}
#' @importFrom magrittr %>%
#' @importFrom tidyr gather
#' @importFrom stats rnorm rbinom
#' @importFrom boot inv.logit
#' 
#' @import dplyr
#' @export
#'
simulate_functional_data = function(lambda1 = 2, lambda2 = 1, I = 50, D = 100, seed = 1988){
	###
	grid = seq(0, 1, length.out = D)
	set.seed(seed)
	
	Y = matrix(NA, I, D)
	colnames(Y) = paste0("time", 1:D)
	for(i in 1:I){
		Y[i,] = mean_sim(grid) + psi1_sim(grid) * rnorm(1, 0, lambda1) + psi2_sim(grid) * rnorm(1, 0, lambda2)
	}
	
	last_time = paste0("time", D)
	Y = as.data.frame(Y) %>%	
		mutate(id = row_number()) %>%
		gather(key, value, time1:last_time) %>%
		arrange(id) %>%
		mutate(index = rep(grid, I),
					 prob = inv.logit(value),
					 value = rbinom(I*D, 1, prob)) 
		
	
	return(list(
		Y = Y,
		psi1 = psi1_sim(grid),
		psi2 = psi2_sim(grid),
		alpha = mean_sim(grid)
	))
	###
}
