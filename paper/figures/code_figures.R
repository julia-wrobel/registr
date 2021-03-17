
library(registr)
library(dplyr)
library(ggplot2)
theme_set(theme_minimal(base_size = 16))
library(cowplot)


# data preparation --------------------------------------------------------
dat <- registr::growth_incomplete

# observed curves
gg <- ggplot(dat, aes(index, value, group = id)) +
	geom_line(alpha = 0.22, col = "dodgerblue4") +
	xlab("t* [observed]") + ylab("First derivative") +
	xlim(c(0,18)) +
	ggtitle("Observed growth development") +
	theme(plot.title       = element_text(hjust = 0.5),
				panel.grid.minor = element_blank())
# (with added space to both sides, s.t. it renders nicely in markdown)
gg_empty <- ggplot()
cowplot::plot_grid(gg_empty, gg, gg_empty, nrow = 1, rel_widths = c(.335,.33,.335)) +
	ggsave("1_data.png", width = 12, height = 3)

# and again, slightly changed, for the joint plot
gg1 <- ggplot(dat, aes(index, value, group = id)) +
	geom_line(alpha = 0.22, col = "dodgerblue4") +
	xlab("t* [observed]") + ylab("First derivative") +
	# ggtitle("Observed curves") +
	xlim(c(0,18)) +
	theme(plot.title       = element_blank(),
				panel.grid.minor = element_blank())



# run the joint approach --------------------------------------------------
reg <- register_fpca(Y = dat, npc = 2, max_iterations = 100,
										 incompleteness = "full", lambda_inc = 0.01)

# registered curves
gg2 <- ggplot(reg$Y, aes(t_hat, value, group = id)) +
	geom_line(alpha = 0.2, col = "dodgerblue4") +
	xlab("t [registered]") + ylab("First derivative") +
	# ggtitle("Registered curves") +
	xlim(c(0,18)) +
	theme(plot.title       = element_blank(),
				axis.title.y     = element_blank(),
				axis.text.y      = element_blank(),
				axis.ticks.y     = element_blank(),
				panel.grid.minor = element_blank())

# estimated warping functions
gg3 <- ggplot(reg$Y, aes(tstar, t_hat, group = id)) +
	geom_line(alpha = 0.2, col = "dodgerblue4") +
	xlab("t* [observed]") + ylab("t [registered]") +
	# ggtitle("Estimated inverse warping functions") +
	xlim(c(0,18)) + ylim(c(0,18)) +
	theme(plot.title       = element_blank(),
				panel.grid.minor = element_blank())

cowplot::plot_grid(gg1, gg2, gg3, nrow = 1, rel_widths = c(0.37, 0.3, 0.33)) +
	ggsave("2_registration.png", width = 12, height = 3)

# estimated FPCs
# (with added space to both sides, s.t. it renders nicely in markdown)
ggf      <- plot(reg$fpca_obj, plot_FPCs = 2, ylab = "First derivative")
gg_empty <- ggplot()
cowplot::plot_grid(gg_empty, ggf, gg_empty, nrow = 1, rel_widths = c(.165,.67,.165)) +
	ggsave("3_FPCA.png", width = 12, height = 3)
	