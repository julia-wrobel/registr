#' Visualize results of binary registration
#' 
#' Function takes in registration object and plots unregistered kernel smoothed curves, registered kernel smoothed curves, 
#' and warping functions. 
#' 
#' @param reg_obj registration object to be plotted.
#' @param alpha alpha value for ggplot.
#' @param tstar vector of original observed time values.
#' @param t.min minimum value to be evaluated on the time domain (useful if data are sparse and / or irregular). 
#' if `NULL`, taken to be minimum observed value.
#' @param t.max maximum value to be evaluated on the time domain (useful if data are sparse and / or irregular). 
#' if `NULL`, taken to be maximum observed value.
#' 
#' @author Julia Wrobel \email{jw3134@@cumc.columbia.edu}
#' 
#' @importFrom gridExtra grid.arrange
#' @importFrom stats ksmooth
#' @import ggplot2
#' @export
#'
#'
plot_binary = function(reg_obj, alpha = 0.25, tstar, t.min = NULL, t.max = NULL){
  
  Y = reg_obj$Y
  t = reg_obj$Y$index
  Kt = reg_obj$Kt
  
  data = data_clean(Y)
  Y_rows = data$Y_rows
  I = data$I
  
  Y.smooth = rep(NA, dim(Y)[1])
  for(i in 1:I){
    subject_rows = Y_rows$first_row[i]:Y_rows$last_row[i]
    Yi = Y$value[subject_rows]
    tstar.i = tstar[subject_rows]
    
    Y.smooth[subject_rows] = ksmooth(tstar.i, Yi, "normal", bandwidth = 1/Kt)$y
  }
  
  
  Y.df = data.frame(
    id = reg_obj$Y$id,
    t.star = tstar,
    t.hat = t,
    Y = Y$value,
    mu = Y.smooth
  )
  
  baseplot = ggplot(Y.df, aes(x = t.star, y = mu, group = id)) + theme_bw()  + ylab("")
  
  unreg = baseplot + geom_line(alpha = alpha) + xlab("t star values")  
  
  reg = baseplot + geom_line(aes(x = t.hat, y = mu), alpha = alpha) + xlab(" t hat values") +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) 
  
  warp = baseplot + geom_line(aes(x = t.star, y = t.hat, group = id), alpha = alpha) +
    geom_point(aes(x = t.star, y = t.star), col = "red", size = 0.5) + 
    xlab("Observed time") + ylab("Warped time") +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  
  grid.arrange(unreg, reg, warp, ncol = 3)
  
}