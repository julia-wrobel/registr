#' Visualize results of registration
#' 
#' Function takes in registration object and plots unregistered curves, registered curves, 
#' and warping functions. Note that there are a lot of specifications right now that only make sense with 
#' a simulation.
#'
#' @param reg_obj registration object to be plotted.
#' @param fpca_obj fpca object that contains mean values (fitted values).
#' @param family \code{gaussian} or \code{binomial}.
#' @param alpha alpha value for ggplot.
#' @param tstar vector of original observed time values.
#' @param t.min minimum value to be evaluated on the time domain (useful if data are sparse and / or irregular). 
#' if `NULL`, taken to be minimum observed value.
#' @param t.max maximum value to be evaluated on the time domain (useful if data are sparse and / or irregular). 
#' if `NULL`, taken to be maximum observed value.
#' @param Kt number of B-spline basis functions used to estimate coefficient functions. Defaults to 10.
#' 
#' @author Julia Wrobel \email{jw3134@@cumc.columbia.edu}
#' 
#' @importFrom gridExtra grid.arrange
#' @importFrom boot inv.logit
#' @import ggplot2
#' @export
#'
#'

plot_registration = function(reg_obj, fpca_obj = NULL, family = "gaussian", alpha = 0.25, tstar, t.min = NULL, t.max = NULL){

  Y = reg_obj$Y$value
  Kt = reg_obj$Kt
  
  if(is.null(fpca_obj)){
    global_knots = quantile(tstar, probs = seq(0, 1, length = Kt - 2))[-c(1, Kt - 2)]
    basis = bs(c(t.min, t.max, tstar), knots = global_knots, intercept = TRUE)[-(1:2),] 
    
    mean.coefs = coef(glm(Y ~ 0 + basis, family = family))
    mu = basis %*% mean.coefs
    mean.df = NULL
  }else{
    mu = fpca_obj$Yhat$value
    t = reg_obj$Y$index
    mean.df = data.frame(mean = fpca_obj$mu, id = 1, time = tstar[1:length(fpca_obj$mu)], t = t[1:length(fpca_obj$mu)])
  }
  
  if(family == "binomial" ){
    mu = inv.logit(mu)
    #mean.df$mean = inv.logit(mean.df$mean)
  }
  

  Y.df = data.frame(
    id = reg_obj$Y$id,
    t.star = tstar,
    t.hat = reg_obj$Y$index,
    Y = Y,
    mu = mu
  )
    
    baseplot = ggplot(Y.df, aes(x = t.star, y = mu, group = id)) + theme_bw()  + ylab("")
    
    unreg = baseplot + geom_line(alpha = alpha) + xlab("t star values")  
      # + geom_point(data = mean.df, aes(x = time, y = mean), col = "red", size = 0.5)

    reg = baseplot + geom_line(aes(x = t.hat, y = mu), alpha = alpha) + xlab(" t hat values") +
      theme(axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank()) 
      #+ geom_point(data = mean.df, aes(x = time, y = mean), col = "red", size = 0.5)
     
    warp = baseplot + geom_line(aes(x = t.star, y = t.hat, group = id), alpha = alpha) +
      geom_point(aes(x = t.star, y = t.star), col = "red", size = 0.5) + 
      xlab("Observed time") + ylab("Warped time") +
      theme(axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank())
    
    grid.arrange(unreg, reg, warp, ncol = 3)
  
}