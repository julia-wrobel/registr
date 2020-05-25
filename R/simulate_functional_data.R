#' Simulate mean
#' 
#' This function generates mean for simulated functional data.
#' 
#' @param grid Grid of x values over which to evaluate the function.
#' 
mean_sim = function(grid){
	sin(grid * 3 * pi) 
}


#' Simulate PC1
#' 
#' This function generates the first principal component for simulated functional data.
#' 
#' @param grid Grid of x values over which to evaluate the function.
#'
psi1_sim = function(grid){
	sin(grid * 3 * pi)
}


#' Simulate PC2
#' 
#' This function generates the second principal component for simulated functional data.
#' 
#' @param grid Grid of x values over which to evaluate the function.
#'
psi2_sim = function(grid){
	cos(grid * 3 * pi)
}


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
#' @param seed Seed for reproducibility. Default is 1988.
#' @param vary_D Indicates if grid length vary by subject. If FALSE all subjects have grid length D.
#' 
#' @author Julia Wrobel \email{julia.wrobel@@cuanschutz.edu}
#' @importFrom magrittr %>%
#' @importFrom tidyr gather
#' @importFrom stats rnorm rbinom runif plogis
#' 
#' @import dplyr
#' 
#' @return A list containing:
#' \item{Y}{Simulated dataframe with variables id, value, index, and latent_mean.}
#' \item{psi1}{True values for first principal component.}
#' \item{psi2}{True values for second principal component.}
#' \item{alpha}{True values for population-level mean.}
#' @export
#'
#' @return A list containing:
#' \item{Y}{A dataframe of simulated data.}
#' \item{psi1}{The first simulated eigenfunction.}
#' \item{psi2}{The second simulated eigenfunction.}
#' \item{alpha}{The population mean.}
#'
simulate_functional_data = function(lambda1 = 2, lambda2 = 1, I = 50, D = 100, seed = 1988,
																		vary_D = FALSE){
	
	## NULLify global values called by tidyverse functions
	value = key = time1 = latent_mean =  NULL
	
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
					 latent_mean = value,
					 value = rbinom(I*D, 1, plogis(latent_mean))) %>%
		dplyr::select(-key)
	
	if(vary_D){
		D_vec = D + as.integer(runif(I, -(D/10), D/10))
		Y = data.frame(
			id = rep(1:I,  D_vec),
			index = rep(NA, sum(D_vec)),
			value = rep(NA, sum(D_vec))
		)
		
		for(i in 1:I){
			grid = seq(0, 1, length.out = D_vec[i])
			Y_i = mean_sim(grid) + psi1_sim(grid) * rnorm(1, 0, lambda1) + psi2_sim(grid) * rnorm(1, 0, lambda2)
			
			Y[which(Y$id == i), "index"] = grid
			Y[which(Y$id == i), "value"] = Y_i
		}
		
		Y = Y %>% 
			mutate(latent_mean = value,
						 value = rbinom( sum(D_vec), 1, plogis(latent_mean)) ) 
	}
	
	# orthogonalize PCs
	psi_svd = svd(cbind( psi1_sim(grid),  psi2_sim(grid)))
	efunctions = psi_svd$u
	evalues = ( psi_svd$d ) ^ 2
	
	return(list(
		Y = Y,
		psi1 = efunctions[, 1],
		psi2 = efunctions[, 2], 
		alpha = mean_sim(grid)
	))
	###
}