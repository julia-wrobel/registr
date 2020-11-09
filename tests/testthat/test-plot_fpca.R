context("plot_fpca")

test_that("plot_fpca returns a cowplot object with classes 'gg' and 'ggplot' for all FPCA functions",{
	
	# Gaussian FPCA
	data(growth_incomplete)
	fpca_obj  = fpca_gauss(Y = growth_incomplete)
	fpca_obj2 = gfpca_twoStep(Y = growth_incomplete, npc = 2, family = "gaussian")
	plot_obj  = plot_fpca(fpca_obj)
	plot_obj2 = plot_fpca(fpca_obj2)
	
	# Binomial FPCA
	Y = simulate_functional_data()$Y
	fpca_obj3 = bfpca(Y)
	plot_obj3 = plot_fpca(fpca_obj3)
	
	expect_s3_class(plot_obj, class = c("gg","ggplot"))
	expect_s3_class(plot_obj2, class = c("gg","ggplot"))
	expect_s3_class(plot_obj3, class = c("gg","ggplot"))
})
