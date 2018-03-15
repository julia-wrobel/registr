#' Loss function for registration step optimization
#'
#' @param Y vector of observed points.
#' @param tstar observed time points
#' @param mean_coefs spline coefficient vector for mean curve.
#' @param knots knot locations for B-spline basis used to estimate mean and FPC basis function.
#' @param family \code{gaussian} or \code{binomial}.
#' @param t_min minimum value to be evaluated on the time domain. 
#' @param t_max maximum value to be evaluated on the time domain. 
#' @param alpha_beta alpha beta coefficients for warping function h.
#' 
#' @importFrom boot inv.logit
#' @export
#'
loss_h2 = function(Y, tstar, mean_coefs, knots, family, t_min, t_max, alpha_beta){
	
	# probably can simpliy the arguments here
	# do I actually want to be starting with the old warping functions each iteration?
	hinv_tstar = pbeta(tstar, alpha_beta[1], alpha_beta[2])
	Theta_phi = bs(hinv_tstar, knots = knots, intercept = TRUE)
	g_mu_t = Theta_phi %*% mean_coefs
	
	if (family == "gaussian") {
		return(sum((Y - g_mu_t) ^ 2))
	} else if (family == "binomial") {
		pi_h = inv.logit(g_mu_t)
		return(-1 * sum(Y * log(pi_h) + (1 - Y) * log(1 - pi_h) ))
	}
}