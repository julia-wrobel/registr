
#' Correct slightly improper parameter vectors
#' 
#' Internal function. In the joint iterations between registration and GFPCA,
#' the optimization with \code{constrOptim()} in the registration step sometimes
#' leads to slightly improper solutions, which cause the optimization to
#' throw an error in the following optimization step. This function corrects
#' the parameter vector if one of the following slight inconsistencies occurs
#' that can mess with the optimization of \code{constrOptim()}: \cr
#' - two neighboring values of the parameter vector are too similar \cr
#' - the initial values of the parameter vector are smaller than \code{t_min},
#' the minimum of the underlying time domain \cr
#' - the last values of the parameter vector are greater than \code{t_max},
#' the maximum of the underlying time domain \cr
#' - one parameter value is slightly greater than its following value, i.e.
#' the parameter vector is not monotone.
#' 
#' @param beta Parameter vector.
#' @param t_min,t_max Minimum and maximum of the underlying time domain in the
#' registration step.
#' 
#' @author Alexander Bauer \email{alexander.bauer@@stat.uni-muenchen.de}
#' 
#' @return A slightly changed parameter vector that ensures a proper solution
#' in the optimization of the registration step.
#' 
#' @examples
#' beta_improper = c(0.24, 1.000047, 1.000002)
#' registr:::ensure_proper_beta(beta_improper, t_min = 0, t_max = 1)
ensure_proper_beta = function(beta, t_min, t_max) {
	
	# problem 1: constrOptim sometimes runs into problems if two neighboring
	#            values are too similar
	# solution:  ensure a minimal difference between neighboring values
	while (any(diff(beta) < 1e-6)) {
		index = which(diff(beta) < 1e-6)[1]
		beta[index]     = min(beta[index], beta[index + 1])
		beta[index + 1] = beta[index] + 1e-5
	}
	
	# problem 2: the initial betas are sometimes minimally smaller than t_min
	# solution:  add some value to the innitial values
	if (any(beta < t_min)) {
		indices_tooSmall = which(beta <= t_min)
		min_beta         = min(beta)
		beta[indices_tooSmall] = beta[indices_tooSmall] + (t_min - min_beta)
	}

	# problem 3: the last betas are sometimes minimally greater than t_max
	# solution:  subtract some value from the last values
	if (any(beta > t_max)) {
		indices_tooBig = which(beta >= t_max)
		max_beta       = max(beta)
		beta[indices_tooBig] = beta[indices_tooBig] + (t_max - max_beta)
	}
	
	# problem 4: one value is minimally greater than its following value
	# solution:  swap the two values
	if (any(diff(beta) < 0)) {
		indices_diff_negative       = which(diff(beta) < 0)
		indices_beta_negative       = unique(unlist(lapply(indices_diff_negative, function(i) { c(i, i + 1) })))
		beta_negative_sorted        = sort(beta[indices_beta_negative])
		beta[indices_beta_negative] = beta_negative_sorted
	}
	
	return(beta)
}
