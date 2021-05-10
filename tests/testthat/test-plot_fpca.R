context("plot.fpca")

test_that("plot.fpca returns a cowplot object with classes 'gg' and 'ggplot' for all FPCA functions",{
	
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("cowplot")
  # Gaussian FPCA
	data(growth_incomplete)
	expect_warning({
	  fpca_obj1a = fpca_gauss(Y = growth_incomplete, npc = 2)
	}, "fpca_gauss convergence not reached. Try increasing maxiter.")
	fpca_obj1b = gfpca_twoStep(Y = growth_incomplete, npc = 2, family = "gaussian")
	plot_obj1a = plot(fpca_obj1a)
	plot_obj1b = plot(fpca_obj1b)
	
	# Binomial FPCA
	Y = simulate_functional_data()$Y
	fpca_obj2a = bfpca(Y, npc = 2)
	fpca_obj2b = gfpca_twoStep(Y = Y, npc = 2, family = "binomial")
	plot_obj2a = plot(fpca_obj2a)
	plot_obj2b = plot(fpca_obj2b)
	
	# Gamma FPCA
	Y$value = Y$value + 1 # get strictly positive data
	fpca_obj3 = gfpca_twoStep(Y = Y, npc = 2, family = "gamma")
	plot_obj3 = plot(fpca_obj3)
	
	# Poisson FPCA
	fpca_obj4 = gfpca_twoStep(Y = Y, npc = 2, family = "poisson")
	plot_obj4 = plot(fpca_obj4)
	
	expect_s3_class(plot_obj1a, class = c("gg","ggplot"))
	expect_s3_class(plot_obj1b, class = c("gg","ggplot"))
	expect_s3_class(plot_obj2a, class = c("gg","ggplot"))
	expect_s3_class(plot_obj2b, class = c("gg","ggplot"))
	expect_s3_class(plot_obj3,  class = c("gg","ggplot"))
	expect_s3_class(plot_obj4,  class = c("gg","ggplot"))
})
