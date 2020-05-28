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
#' @param prior_sd standard deviation on prior for warping functions.
#' 
#' @return The scalar value taken by the loss function.
#' 
#' @importFrom stats plogis
#' @export
#'

loss_h = function(Y, Theta_h, mean_coefs, knots, beta.inner, family, t_min, t_max, 
									#parametric_warps = FALSE, 
									prior_sd = 1){
  
  beta = c(t_min, beta.inner, t_max)
  hinv_tstar = Theta_h %*% beta 
  
  Theta_phi = bs(hinv_tstar, knots = knots, intercept = TRUE)
  g_mu_t = Theta_phi %*% mean_coefs
  
  if (family == "gaussian") {
    loss = -1 * sum((Y - g_mu_t) ^ 2)
    
  } else if (family == "binomial") {
    pi_h = plogis(g_mu_t)
    loss =  sum(Y * log(pi_h) + (1 - Y) * log(1 - pi_h) )
  }
  
  Kh = length(beta)
  beta_prior = coef(lm(seq(t_min, t_max, length.out = length(Y)) ~ 0 + Theta_h))[-c(1, Kh)]
  priorlik = sum(dnorm(x = beta.inner, mean = beta_prior, sd = prior_sd, log = TRUE))  
  loss = loss + priorlik
  
  # if(parametric_warps == "monotone_prior"){
  #   #is_monotonic = all(beta == cummax(beta))
  #   #loss = loss + log(is_monotonic) #+ priorlik 
  #   
  #   if(has_beta_prior){
  # 
  #   }
  #}

  return(-1 * loss)
}