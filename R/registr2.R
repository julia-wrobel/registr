#' Register curves from exponential family
#' @author Julia Wrobel \email{jw3134@@cumc.columbia.edu}
#' @export
#' 
#' @importFrom stats glm coef constrOptim quantile
#' 
registr2 = function(obj = NULL, Y = NULL, Kt = 8, family = "binomial",
										gradient = TRUE, t_min = NULL, t_max = NULL, 
										row_obj = NULL, ...){
	
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
	
	
	if(is.null(obj)){
		# define population mean
		global_knots = quantile(tstar, probs = seq(0, 1, length = Kt - 2))[-c(1, Kt - 2)]
		basis = bs(c(t_min, t_max, tstar), knots = global_knots, intercept = TRUE)[-(1:2),] 
		mean_coefs = coef(glm(Y$value ~ 0 + basis, family = family))
	}else{
		global_knots = obj$knots
		mean_coefs = obj$subject_coefs
	}
	
	
	### Calculate warping functions  
	t_hat = rep(NA, dim(Y)[1])
	loss_subjects = rep(NA, I)
	alpha_beta_new = matrix(NA, I, 2)
	colnames(alpha_beta_new) = c("a", "b")
	alpha_beta_0 = c(1, 1)
	for (i in 1:I){
		subject_rows = rows$first_row[i]:rows$last_row[i]
		Yi = Y$value[subject_rows]
		Di = length(Yi)
		
		tstar_i = tstar[subject_rows]
		
		# change below
		mean_coefs_i = mean_coefs
		
		alpha_beta_new[i,] = optim(alpha_beta_0, loss_h2, Y = Yi, tstar = tstar_i,
															 mean_coefs = mean_coefs_i,knots = global_knots,
															 family = family, t_min = t_min, t_max = t_max)$par
		

		t_hat[subject_rows] = pbeta(tstar_i, alpha_beta_new[i, 1], alpha_beta_new[i, 2])

		
	} # end loop over subjects
	Y$index = t_hat
	
	alpha_beta = as.data.frame(alpha_beta_new)
	alpha_beta$id = unique(Y$id)
	
	return(list(Y = Y, alpha_beta = alpha_beta)) 
	
}# end function
	