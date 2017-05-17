#' Register curves from exponential family
#' 
#' Function used in the registration step. Assumes you feed it an object from data_clean, mean coefficients.
#'
#' @param data output of data_clean function
#' @param mean.coefs If null, calculates the population mean from the data.
#' @param Kt number of B-spline basis functions used to estimate mean functions. Defaults to 10.
#' @param Kh number of B-spline basis functions used to estimate warping functions \emph{h}. Defaults to 5.
#' @param family \code{gaussian} or \code{binomial}.
#' @param gradient if TRUE, uses analytic gradient to calculate derivative. 
#' If FALSE, calculates gradient nuemrically.
#' @param t.min minimum value to be evaluated on the time domain (useful if data are sparse and / or irregular). 
#' if `NULL`, taken to be minimum observed value.
#' @param t.max maximum value to be evaluated on the time domain (useful if data are sparse and / or irregular). 
#' if `NULL`, taken to be maximum observed value.
#' 
#' 
#' @author Julia Wrobel \email{jw3134@@cumc.columbia.edu}
#' @export
#'
reg_basic = function(data = NULL, mean.coefs = NULL, Kt = 10, Kh = 5, family = "gaussian", gradient = TRUE,
                    t.min = NULL, t.max = NULL){
  if(is.null(obj)){
    # clean data
    data = data_clean(Y)
    Y = data$Y
    tstar = Y$index
    
    # define population mean
    global_knots = define_basis(tstar, t.min, t.max, Kt)$knots
    basis = define_basis(tstar, t.min, t.max, Kt)$basis
    mean.coefs = coef(glm(Y$value ~ 0 + basis, family = family))
  }else{
    data = data_clean(obj$Y)
    Y = data$Y
    tstar = Y$index
    mean.coefs = obj$subject_coefs
    global_knots = obj$knots
  }
  
  I = data$I
  Y_rows = data$Y_rows
  
  basis.tstar = bs(tstar, df = Kh, intercept = FALSE)
  
  #### Define optimization constraints
  constrs = constraints(Kh)
  ui = constrs$ui
  ci = constrs$ci
  
  beta.0 = seq(0, 1, length.out = Kh + 1)[-c(1, Kh + 1)]
  
  ### Calculate warping functions  
  t.hat = rep(NA, dim(Y)[1])
  for(i in 1:I){
    
    subject_rows = Y_rows$first_row[i]:Y_rows$last_row[i]
    
    Yi = Y$value[subject_rows]
    Di = length(Yi)
    
    tstar.i = tstar[subject_rows]
    basis.tstar.i = basis.tstar[subject_rows ,]
    
    if(is.null(obj)){mean.coefs.i = mean.coefs}else{mean.coefs.i = mean.coefs[, i]}
    if(gradient){gradf = loss_h_gradient}else{gradf = NULL}
    
    coefs.inner = constrOptim(beta.0, loss_h, grad = gradf, ui = ui, ci = ci, Y = Yi, 
                              basis.tstar = basis.tstar.i, mean.coefs = mean.coefs.i, knots = global_knots,
                              family = family)$par
    
    coefs.h = c(0, coefs.inner, 1)
    #t.hat[i,] = cbind(1, basis.tstar.i) %*% coefs.h
    
    t.hat[subject_rows] = cbind(1, basis.tstar.i) %*% coefs.h
      
  }
   Y$index = t.hat
  
  return(list(Y = Y, Kt = Kt, Kh = Kh)) ## what if you want to return all of the steps of the registration?
} # end registration function
