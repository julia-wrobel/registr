context("fpca_gauss")

test_that("bfpca output is a list with non-null items and class fpca",{
	Y = simulate_functional_data()$Y
	Y$value = Y$latent_mean
	
	fpca_object = fpca_gauss(Y, npc = 2, print.iter = TRUE)
	
	expect_equal(class(fpca_object), "fpca")
	expect_equal(fpca_object$family, "gaussian")
	
	expect_false(any(is.na(fpca_object$mu)))
	expect_false(any(is.na(fpca_object$efunctions)))
	expect_false(any(is.na(fpca_object$evalues)))
	expect_false(any(is.na(fpca_object$scores)))
})