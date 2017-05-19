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
#' If FALSE, calculates gradient nuemrically.
#' @param t.min minimum value to be evaluated on the time domain (useful if data are sparse and / or irregular). 
#' if `NULL`, taken to be minimum observed value.
#' @param t.max maximum value to be evaluated on the time domain (useful if data are sparse and / or irregular). 
#' if `NULL`, taken to be maximum observed value.
#' @param row_obj if NULL, the function cleans the data and calculates row indices. Keep this NULL if you are using 
#' standalone \code{register} function.
#' 
#' @author Julia Wrobel \email{jw3134@@cumc.columbia.edu}
#' @export
#' 
#' @importFrom stats glm coef constrOptim quantile
#' 
#' @examples
#' 
#' \dontrun{
#'    register_step = register(obj = NULL, Y = Y_sim, Kt = 10, Kh = 5, family = "binomial", gradient = TRUE)
#' }
#'
register = function(obj = NULL, Y = NULL, Kt = 10, Kh = 5, family = "gaussian", gradient = TRUE,
                    t.min = NULL, t.max = NULL, row_obj = NULL){
  
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
  if (is.null(t.min)) {t.min = min(tstar)}
  if (is.null(t.max)) {t.max = max(tstar)}
  basis.tstar = bs(tstar, df = Kh, intercept = FALSE)
  
  if(is.null(obj)){
    # define population mean
    global_knots = quantile(tstar, probs = seq(0, 1, length = Kt - 2))[-c(1, Kt - 2)]
    basis = bs(c(t.min, t.max, tstar), knots = global_knots, intercept = TRUE)[-(1:2),] 
    mean.coefs = coef(glm(Y$value ~ 0 + basis, family = family))
  }else{
    global_knots = obj$knots
    mean.coefs = obj$subject_coefs
  }

  
  #### Define optimization constraints
  constrs = constraints(Kh)
  ui = constrs$ui
  ci = constrs$ci
  
  beta.0 = seq(t.min, t.max, length.out = Kh + 1)[-c(1, Kh + 1)]
  
  ### Calculate warping functions  
  t.hat = rep(NA, dim(Y)[1])
  loss.subjects = rep(NA, I)
  for (i in 1:I) {
    
    subject_rows = rows$first_row[i]:rows$last_row[i]
    
    Yi = Y$value[subject_rows]
    Di = length(Yi)
    
    tstar.i = tstar[subject_rows]
    basis.tstar.i = basis.tstar[subject_rows ,]
    
    if (is.null(obj)) {mean.coefs.i = mean.coefs} else {mean.coefs.i = mean.coefs[, i]}
    if (gradient) {gradf = loss_h_gradient} else {gradf = NULL}
    
    coefs.inner = constrOptim(beta.0, loss_h, grad = gradf, ui = ui, ci = ci, Y = Yi, 
                              basis.tstar = basis.tstar.i, mean.coefs = mean.coefs.i, knots = global_knots,
                              family = family, t.min = t.min, t.max = t.max)$par
    
    coefs.h = c(0, coefs.inner, 1)

    t.hat[subject_rows] = cbind(1, basis.tstar.i) %*% coefs.h
    loss.subjects[i] = loss_h(Yi, basis.tstar.i, mean.coefs.i, global_knots, coefs.inner, family = family, 
                              t.min = t.min, t.max = t.max)
  }
  Y$index = t.hat

  return(list(Y = Y, Kt = Kt, Kh = Kh, loss = sum(loss.subjects))) 
} # end registration function
