
#' Plot the results of a functional PCA
#' 
#' S3 plot method for class \code{fpca}.
#' Plot FPCA results by visualizing the variation of the individual FPCs around
#' the global mean. based on an object created with function
#' \code{\link{fpca_gauss}}, \code{\link{bfpca}} or \code{\link{gfpca_twoStep}}. \cr \cr
#' The shares of explained variance are included in the plot titles if
#' \code{x} contains an element \code{evalues_sum}.
#' 
#' @param x Object of class \code{"fpca"}.
#' @param plot_FPCs Optional index vector of the FPCs to be plotted.
#' Defaults to all FPCs contained in \code{x}.
#' @param sd_factor Numeric factor with which the standard deviations of each
#' FPC's scores are multiplied to display its variation in the plots.
#' Defaults to 2.
#' @param response_function Optional response function to be applied before
#' plotting the curves. Defaults to \code{NULL}, i.e. the identity function if
#' \code{x$family} is one of \code{c("gaussian","binomial")} or
#' \code{exp()} if \code{x$family} is one of \code{c("gamma","poisson")}.
#' @param add_symbols Indicator if '+' and '-' symbols should be added to the
#' plot to highlight the direction of the displayed FPCs. Defaults to TRUE.
#' @param subtitle If TRUE (default) the parameter \code{sd_factor}
#' is displayed in the plot subtitle.
#' @param xlim,ylim Optional numeric vectors with limits for the x and y axis.
#' @param xlab,ylab Optional titles for the x and y axis.
#' @param ... Additional arguments passed to \code{\link[ggplot2]{theme}}.
#' 
#' @importFrom dplyr bind_rows
#' @importFrom stats sd
#' @export
#' 
#' @return @return If multiple FPCs are plotted, returns a grid of \code{ggplot}
#' plots, created with \code{cowplot::plot_grid}. If only one FPC is plotted,
#' returns a single \code{ggplot} plot.
#' 
#' @author Alexander Bauer \email{alexander.bauer@@stat.uni-muenchen.de}
#' 
#' @examples
#' data(growth_incomplete)
#' 
#' fpca_obj = fpca_gauss(Y = growth_incomplete, npc = 2)
#' if (requireNamespace("ggplot2", quietly = TRUE) &&
#' requireNamespace("cowplot", quietly = TRUE)) {
#' library(ggplot2)
#' plot(fpca_obj)
#' }
#' 
plot.fpca = function(x, plot_FPCs = 1:x$npc, sd_factor = 2,
                     response_function = NULL, add_symbols = TRUE,
                     subtitle = TRUE, xlim = NULL, ylim = NULL,
                     xlab = "t [registered]", ylab = "y", ...) {
  id = NULL
  rm(list="id")
  if (!requireNamespace("ggplot2", quietly = TRUE) || 
      !requireNamespace("cowplot", quietly = TRUE)) {
    stop("ggplot2 and cowplot are required for plot.fpca function")
  }
  # some NULL variable definitions to appease CRAN package checks regarding the use of ggplot2
  fpc_value = type = value = NULL
  
  if (is.null(response_function)) {
    if (x$family %in% c("gamma","poisson")) { # log link
      response_function = exp
    } else {
      response_function = function(x) { x }
    }
  }
  
  # determine the standard deviation of each FPC's scores
  FPC_sdVec <- apply(x$scores, 2, stats::sd)
  
  # data preparation
  mean_dat = data.frame(t    = x$t_vec,
                        mean = as.vector(x$mu))
  fpc_dat  = data.frame(id        = rep(1:x$npc, each = length(x$t_vec)),
                        t         = rep(x$t_vec, times = x$npc),
                        fpc_value = as.vector(x$efunctions))
  
  # apply response function to mean
  mean_dat = mean_dat %>%
    mutate(mean_y = response_function(mean))
  
  # prepare eigenvalue information for the plot titles
  if (!is.null(x$evalues_sum)) {
    ev_share = round(100 * x$evalues / x$evalues_sum, 1)
    ev_share[ev_share == 0] = "<1"
    ev_info  = paste0(" (",ev_share,"%)")
  }
  
  # create a list with one ggplot object per principal components
  plotDat_list = lapply(plot_FPCs, function(i) {
    sd_factor_i <- sd_factor * FPC_sdVec[i]
    # dataset with mean +/- <sd_factor_i>*FPC
    plot_dat = fpc_dat %>%
      dplyr::filter(id == i) %>%
      dplyr::arrange(t) %>%
      dplyr::left_join(mean_dat, by = "t") %>%
      dplyr::mutate(mean_plusFPC  = response_function(mean + sd_factor_i*fpc_value),
                    mean_minusFPC = response_function(mean - sd_factor_i*fpc_value))
    plot_dat = data.frame(id    = unique(plot_dat$id),
                          t     = rep(plot_dat$t, times = 3),
                          value = c(plot_dat$mean_y,
                                    plot_dat$mean_plusFPC,
                                    plot_dat$mean_minusFPC),
                          type = rep(c("mean curve",
                                       paste0("mean + ",sd_factor,"*sd*FPC"),
                                       paste0("mean - ",sd_factor,"*sd*FPC")), 
                                     each = nrow(plot_dat)),
                          stringsAsFactors = FALSE) %>%
      dplyr::mutate(type = factor(type, levels = c("mean curve",
                                            paste0("mean + ",sd_factor,"*sd*FPC"),
                                            paste0("mean - ",sd_factor,"*sd*FPC"))))
  })
  
  if (is.null(xlim))
    xlim = dplyr::bind_rows(plotDat_list) %>% 
    dplyr::pull(t) %>% 
    range()
  if (is.null(ylim))
    ylim = dplyr::bind_rows(plotDat_list) %>% 
    dplyr::pull(value) %>% 
    range()
  
  # create list with one plot per FPC
  plot_list = lapply(plot_FPCs, function(i) {
    
    fpc_index <- match(i, plot_FPCs)
    
    # plot
    gg <- ggplot2::ggplot(mapping = ggplot2::aes(x = t, y = value, col = type, lty = type)) +
      ggplot2::geom_line(data = plotDat_list[[fpc_index]]) +
      ggplot2::scale_color_manual(values = c("black", "gray70", "gray60")) +
      ggplot2::scale_linetype_manual(values = c(1, 2, 3)) +
      ggplot2::xlim(xlim) + 
      ggplot2::ylim(ylim) + 
      ggplot2::xlab(xlab) + 
      ggplot2::ylab(ylab) +
      ggplot2::ggtitle(paste0("FPC ", i, ifelse(!is.null(x$evalues_sum), ev_info[i], ""),
                              ifelse(subtitle, paste0("\n(mean +/- ",sd_factor,"*sd*FPC)"), ""))) +
      ggplot2::theme(legend.position  = "none",
                     panel.grid.minor = ggplot2::element_blank(),
                     plot.title       = ggplot2::element_text(hjust = 0.5),
                     ...)
    
    if (add_symbols) { # plot '+' and '-' symbols along the curves
      n_symbols = 4 # only plot the symbols this many times along each line
      t_symbols = unique(plotDat_list[[fpc_index]]$t)[round(seq(1, length(unique(plotDat_list[[fpc_index]]$t)), 
                                                                length.out = n_symbols))]
      symbol_dat      = plotDat_list[[fpc_index]] %>% filter(t %in% t_symbols)
      symbolPlus_dat  = symbol_dat %>% filter(type == paste0("mean + ",sd_factor,"*sd*FPC"))
      symbolMinus_dat = symbol_dat %>% filter(type == paste0("mean - ",sd_factor,"*sd*FPC"))

      gg <- gg +
        ggplot2::geom_text(data = symbolPlus_dat,  label = "+", size = 7, col = "gray40") +
        ggplot2::geom_text(data = symbolMinus_dat, label = "-", size = 7, col = "gray40")
    }
    
    return(gg)
  })
  
  # plot grid
  if (length(plot_list) > 1) { # multiple FPCs are plotted
    cowplot::plot_grid(plotlist = plot_list, align = T)
  } else { # only one FPC is plotted
    plot_list[[1]]
  }
}
