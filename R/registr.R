#' Register curves from exponential family
#' 
#' Function used in the registration step of an FPCA-based approach for registering exponential-family functional data. 
#' This method uses constrained optimization to estimate spline coefficients for warping functions,
#' where the objective function for optimization comes from maximizing the EF likelihood 
#' subject to monotonicity constraints on the warping functions. You have to either specify \code{obj}, which is a fpca 
#' object from an earlier step, or \code{Y}, which can be a matrix in wide format or a dataframe in long format. If 
#' \code{Y} is a dataframe, the subjects IDs, times, and observations should have column names id, index, and value, 
#' respectively. If Y is input, then we automatically calculate a template using glm.
#'
#' @param obj current estimate of FPC objects (including spline coefs and scores). Can be NULL only if Y argument is selected.
#' @param Y data matrix, each row is a subject.
#' @param Kt number of B-spline basis functions used to estimate mean functions. Defaults to 10.
#' @param Kh number of B-spline basis functions used to estimate warping functions \emph{h}. Defaults to 5.
#' @param family \code{gaussian} or \code{binomial}.
#' @param gradient if TRUE, uses analytic gradient to calculate derivative. 
#' If FALSE, calculates gradient numerically.
#' @param beta Initial values for beta for each subject. If NULL, these are chosen using seq().
#' @param t_min minimum value to be evaluated on the time domain (useful if data are sparse and / or irregular). 
#' if `NULL`, taken to be minimum observed value.
#' @param t_max maximum value to be evaluated on the time domain (useful if data are sparse and / or irregular). 
#' if `NULL`, taken to be maximum observed value.
#' @param row_obj if NULL, the function cleans the data and calculates row indices. Keep this NULL if you are using 
#' standalone \code{registr} function.
#' @param ... additional arguments passed to or from other functions
#' 
#' @author Julia Wrobel \email{jw3134@@cumc.columbia.edu}
#' @export
#' 
#' @importFrom stats glm coef constrOptim quantile
#' 
#' @examples
#' 
#' \dontrun{
#'    register_step = registr(obj = NULL, Y = Y_sim, Kt = 10, Kh = 5, family = "binomial", 
#'    gradient = TRUE)
#' }
#'
registr = function(obj = NULL, Y = NULL, Kt = 10, Kh = 5, family = "binomial", gradient = TRUE,
									 beta = NULL, t_min = NULL, t_max = NULL, row_obj = NULL, ...){
  
  if(is.null(Y)) { Y = obj$Y}
  
  if(is.null(row_obj)){
    data = data_clean(Y)
    Y = data$Y
    rows = data$Y_rows
    I = data$I
  }else{
    rows = row_obj
    I = dim(rows)[1]
  }
  
  tstar = Y$index
  if (is.null(t_min)) {t_min = min(tstar)}
  if (is.null(t_max)) {t_max = max(tstar)}
  Theta_h = bs(tstar, df = Kh, intercept = FALSE)
  
  if(is.null(obj)){
    # define population mean
    global_knots = quantile(tstar, probs = seq(0, 1, length = Kt - 2))[-c(1, Kt - 2)]
    basis = bs(c(t_min, t_max, tstar), knots = global_knots, intercept = TRUE)[-(1:2),] 
    mean_coefs = coef(glm(Y$value ~ 0 + basis, family = family))
  }else{
    global_knots = obj$knots
    mean_coefs = obj$subject_coefs
  }

  #### Define optimization constraints
  constrs = constraints(Kh, t_min, t_max)
  ui = constrs$ui
  ci = constrs$ci

  ### Calculate warping functions  
  t_hat = rep(NA, dim(Y)[1])
  loss_subjects = rep(NA, I)
  beta_new = matrix(NA, Kh - 1, I)
	beta_0 = seq(t_min, t_max, length.out = Kh + 1)[-c(1, Kh + 1)]
  for (i in 1:I) {
    
    subject_rows = rows$first_row[i]:rows$last_row[i]
    
    Yi = Y$value[subject_rows]
    Di = length(Yi)
    
    tstar_i = tstar[subject_rows]
    Theta_h_i = Theta_h[subject_rows ,]
    
    if (is.null(beta)) {beta_i = beta_0} else {beta_i = beta[, i]}
    if (is.null(obj)) {mean_coefs_i = mean_coefs} else {mean_coefs_i = mean_coefs[, i]}
    if (gradient) {gradf = loss_h_gradient} else {gradf = NULL}
    
    beta_new[,i] = constrOptim(beta_i, loss_h, grad = gradf, ui = ui, ci = ci, Y = Yi, 
                              Theta_h = Theta_h_i, mean_coefs = mean_coefs_i, knots = global_knots,
                              family = family, t_min = t_min, t_max = t_max)$par
    
    beta_full_i = c(t_min, beta_new[,i], t_max)
    
    t_hat[subject_rows] = cbind(1, Theta_h_i) %*% beta_full_i
    loss_subjects[i] = loss_h(Yi, Theta_h_i, mean_coefs_i, global_knots, beta_new[,i], family = family, 
                              t_min = t_min, t_max = t_max)
  }
  Y$index = t_hat

  return(list(Y = Y, loss = sum(loss_subjects), beta = beta_new)) 
} # end registration function
