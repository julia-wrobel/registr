
#' Plot the results of a functional PCA
#' 
#' Plot FPCA results by visualizing the variation of the individual FPCs around
#' the global mean. based on an object created with function
#' \code{\link{fpca_gauss}}, \code{\link{bfpca}} or \code{\link{gfpca_twoStep}}. \cr \cr
#' The shares of explained variance are included in the plot titles if
#' \code{fpca_obj} contains an element \code{evalues_sum}. This is currently
#' only the case for FPCA objects created with \code{\link{gfpca_twoStep}}.
#' 
#' @param fpca_obj Object of class \code{"fpca"}.
#' @param plot_npc Optional number of the main FPCs to be plotted. Defaults to
#' all FPCs contained in \code{fpca_obj}.
#' @param var_factor Numeric factor with which the FPC's are multiplied
#' to display their variation in the plots. Defaults to 2, but can be set to
#' 2 times the standard deviation of the obtained FPC scores.
#' @param response_function Optional response function to be applied before
#' plotting the curves. Defaults to the identity function.
#' @param subtitle If TRUE (default) the parameter \code{var_factor}
#' is displayed in the plot subtitle.
#' @param xlim,ylim Optional numeric vectors with limits for the x and y axis.
#' @param xlab,ylab Optional titles for the x and y axis.
#' 
#' @import ggplot2
#' @importFrom cowplot plot_grid
#' @importFrom dplyr bind_rows
#' @export
#' 
#' @return Grid of \code{ggplot} plots, created with \code{cowplot::plot_grid}.
#' 
#' @author Alexander Bauer \email{alexander.bauer@@stat.uni-muenchen.de}
#' 
#' @examples
#' data(growth_incomplete)
#' 
#' fpca_obj = fpca_gauss(Y = growth_incomplete, npc = 2)
#' plot_fpca(fpca_obj)
#' 
plot_fpca = function(fpca_obj, plot_npc = fpca_obj$npc, var_factor = 2,
                     response_function = function(x) { x },
                     subtitle = TRUE, xlim = NULL, ylim = NULL,
                     xlab = "t [registered]", ylab = "y") {
  
  # data preparation
  mean_dat = data.frame(t    = fpca_obj$t_vec,
                        mean = as.vector(fpca_obj$mu))
  fpc_dat  = data.frame(id        = rep(1:fpca_obj$npc, each = length(fpca_obj$t_vec)),
                        t         = rep(fpca_obj$t_vec, times = fpca_obj$npc),
                        fpc_value = as.vector(fpca_obj$efunctions))
  
  # apply response function to mean
  mean_dat = mean_dat %>%
    mutate(mean_y = response_function(mean))
  
  # prepare eigenvalue information for the plot titles
  if (!is.null(fpca_obj$evalues_sum)) {
    ev_share = round(100 * fpca_obj$evalues / fpca_obj$evalues_sum, 1)
    ev_share[ev_share == 0] = "<1"
    ev_info  = paste0(" (",ev_share,"%)")
  }
  
  # create a list with one ggplot object per principal components
  plotDat_list = lapply(1:plot_npc, function(i) {
    # dataset with mean +/- <var_factor>*FPC
    plot_dat = fpc_dat %>%
      filter(id == i) %>%
      arrange(t) %>%
      dplyr::left_join(mean_dat, by = "t") %>%
      mutate(mean_plusFPC  = response_function(mean + var_factor*fpc_value),
             mean_minusFPC = response_function(mean - var_factor*fpc_value))
    plot_dat = data.frame(id    = unique(plot_dat$id),
                          t     = rep(plot_dat$t, times = 3),
                          value = c(plot_dat$mean_y,
                                    plot_dat$mean_plusFPC,
                                    plot_dat$mean_minusFPC),
                          curve = rep(c("mean curve",
                                        paste0("mean + ",var_factor,"*FPC"),
                                        paste0("mean - ",var_factor,"*FPC")), 
                                      each = nrow(plot_dat)),
                          stringsAsFactors = FALSE) %>%
      mutate(curve = factor(curve, levels = c("mean curve",
                                              paste0("mean + ",var_factor,"*FPC"),
                                              paste0("mean - ",var_factor,"*FPC"))))
  })
  
  if (is.null(xlim))
    xlim = dplyr::bind_rows(plotDat_list) %>% pull(t) %>% range()
  if (is.null(ylim))
    ylim = dplyr::bind_rows(plotDat_list) %>% pull(value) %>% range()
  
  # create list with one plot per FPC
  plot_list = lapply(1:plot_npc, function(i) {
    
    # plot '+' and '-' symbols along the curves
    n_symbols = 7 # only plot the symbols this many times along each line
    t_symbols = unique(plotDat_list[[i]]$t)[round(seq(1, length(unique(plotDat_list[[i]]$t)), 
                                                      length.out = n_symbols))]
    symbol_dat      = plotDat_list[[i]] %>% filter(t %in% t_symbols)
    symbolPlus_dat  = symbol_dat %>% filter(curve == paste0("mean + ",var_factor,"*FPC"))
    symbolMinus_dat = symbol_dat %>% filter(curve == paste0("mean - ",var_factor,"*FPC"))
    
    # plot
    ggplot(mapping = aes(x = t, y = value, col = curve)) +
      geom_line(data = plotDat_list[[i]]) +
      scale_color_manual(values = c("black", "gray80", "gray80")) +
      geom_text(data = symbolPlus_dat,  label = "+", size = 7, col = "gray40") +
      geom_text(data = symbolMinus_dat, label = "-", size = 7, col = "gray40") +
      xlim(xlim) + ylim(ylim) + xlab(xlab) + ylab(ylab) +
      ggtitle(paste0("FPC ", i, ifelse(!is.null(fpca_obj$evalues_sum), ev_info[i], ""),
                     ifelse(subtitle, paste0("\n(mean +/- ",var_factor,"*FPC)"), ""))) +
      theme(legend.position  = "none",
            panel.grid.minor = element_blank(),
            plot.title       = element_text(hjust = 0.5))
  })
  
  # plot grid
  cowplot::plot_grid(plotlist = plot_list, align = T)
}
