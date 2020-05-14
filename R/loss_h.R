#' Loss function for registration step optimization
#'
#' @param Y vector of observed points.
#' @param Theta_h B-spline basis for inverse warping functions.
#' @param mean_coefs spline coefficient vector for mean curve.
#' @param knots knot locations for B-spline basis used to estimate mean and FPC basis function.
#' @param beta.inner spline coefficient vector to be estimated for warping function h.
#' @param family \code{gaussian} or \code{binomial}.
#' @param t_min minimum value to be evaluated on the time domain. 
#' @param t_max maximum value to be evaluated on the time domain. 
#' @param parametric_warps If FALSE (default), inverse warping functions are 
#' estimated nonparametrically. If 'beta_cdf', they are assumed to have the form of a 
#' Beta(a,b) CDF. If 'piecewise' they follow a piecewise parameterized function.
#' 
#' @return The scalar value taken by the loss function.
#' 
#' @importFrom stats plogis
#' @export
#'

loss_h = function(Y, Theta_h, mean_coefs, knots, beta.inner, family, t_min, t_max, 
									parametric_warps = FALSE, prior_sd = NULL){
  
  beta = c(t_min, beta.inner, t_max)
  hinv_tstar = Theta_h %*% beta 
  
  Theta_phi = bs(hinv_tstar, knots = knots, intercept = TRUE)
  g_mu_t = Theta_phi %*% mean_coefs
  
  
  if (family == "gaussian") {
    
    loss = sum((Y - g_mu_t) ^ 2)
    #return(sum((Y - g_mu_t) ^ 2))
    
  } else if (family == "binomial") {
    pi_h = plogis(g_mu_t)
    loss = -1 * sum(Y * log(pi_h) + (1 - Y) * log(1 - pi_h) )
    #return(-1 * sum(Y * log(pi_h) + (1 - Y) * log(1 - pi_h) ))
  }
  
  if(parametric_warps == "monotone_prior"){
    print("hi julia")
    #loss = loss - (1- all(beta == cummax(beta))) * 100 # enforces monotonicity of spline coefficients
  
    Kh = length(beta)
    #beta_prior = coef(lm(seq(t_min, t_max, length.out = length(Y)) ~ 0 + Theta_h))[-c(1, Kh)]
    
    for(k in 2:(Kh-1)){
      is_mono = (beta[k-1] < beta[k] & beta[k] < beta[k + 1])
      loss = loss - log(is_mono) # should be -Inf if out of bounds..
    }
    
    #loss = loss - sum(dnorm(x = beta.inner, mean = beta_prior, sd = prior_sd, log = TRUE))  
  }

  if(!is.null(prior_sd)){
    Kh = length(beta)
    #beta_prior = seq(t_min, t_max, length.out = Kh)[-c(1, Kh)]
    
    beta_prior = coef(lm(seq(t_min, t_max, length.out = length(Y)) ~ 0 + Theta_h))[-c(1, Kh)]
    loss = loss - sum(dnorm(x = beta.inner, mean = beta_prior, sd = prior_sd, log = TRUE))  
  }

  return(loss)
}