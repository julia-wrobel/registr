context("constrOptim_helpers")

test_that("ensure_proper_beta corrects too similar neighboring values",{
	beta = c(0.1, 0.2, 0.36742637, 0.36742638, 0.8, 0.9)
	
	expect_equal(ensure_proper_beta(beta, t_min = 0, t_max = 1),
							 expected = c(0.1, 0.2, 0.36742637, 0.36743638, 0.8, 0.9))
})

test_that("ensure_proper_beta corrects values outside the domain",{
	beta1 = c(-0.0000017, -0.00002, 0.1, 0.8, 0.9)
	beta2 = c(0.1, 0.2, 0.9, 1.00002, 1.000017)
	
	expect_equal(ensure_proper_beta(beta1, t_min = 0, t_max = 1),
							 expected = c(0, 1e-5, 0.1, 0.8, 0.9))
	expect_equal(ensure_proper_beta(beta2, t_min = 0, t_max = 1),
							 expected = c(0.1, 0.2, 0.9, 0.99999, 1))
})

test_that("ensure_proper_beta corrects nonmonotone parameters",{
	beta = c(0.1, 0.35326874, 0.3532675, 0.8, 0.9)
	
	expect_equal(ensure_proper_beta(beta, t_min = 0, t_max = 1),
							 expected = c(0.1, 0.3532675, 0.3532775, 0.8, 0.9))
})