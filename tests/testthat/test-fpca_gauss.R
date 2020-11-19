context("fpca_gauss")

test_that("fpca_gauss output is a list with non-null items and class fpca",{
	Y = simulate_functional_data(seed = 2020)$Y
	Y$value = Y$latent_mean
	
	fpca_object = fpca_gauss(Y, npc = 2, print.iter = TRUE)
	
	expect_equal(class(fpca_object), "fpca")
	expect_equal(fpca_object$family, "gaussian")
	
	expect_false(any(is.na(fpca_object$mu)))
	expect_false(any(is.na(fpca_object$efunctions)))
	expect_false(any(is.na(fpca_object$evalues)))
	expect_false(any(is.na(fpca_object$scores)))
})

test_that("fpca_gauss output has correct number of dimensions",{
	Y = simulate_functional_data(I = 100, D = 200, seed = 2020)$Y
	Kt = 8
	fpca_object = fpca_gauss(Y, npc = 2, Kt = Kt, print.iter = TRUE)
	
	expect_equal(dim(fpca_object$alpha), c(200, 1))
	expect_equal(dim(fpca_object$efunctions), c(200, 2))
	expect_equal(length(fpca_object$evalues), 2)
	expect_equal(dim(fpca_object$scores), c(100, 2))
	expect_equal(length(fpca_object$knots), Kt - 4)
	
	fpca_object = fpca_gauss(Y, npc = 1, Kt = Kt, print.iter = TRUE)
	expect_equal(dim(fpca_object$efunctions), c(200, 1))
})

test_that("fpca_gauss output has correct number of dimensions when periodic = TRUE",{
	Y = simulate_functional_data(I = 100, D = 200, seed = 2020)$Y
	Kt = 8
	fpca_object = fpca_gauss(Y, npc = 2, Kt = Kt, periodic = TRUE, print.iter = TRUE)
	
	expect_equal(dim(fpca_object$alpha), c(200, 1))
	expect_equal(dim(fpca_object$efunctions), c(200, 2))
	expect_equal(length(fpca_object$evalues), 2)
	expect_equal(dim(fpca_object$scores), c(100, 2))
	expect_equal(length(fpca_object$knots), Kt - 1)
	
	fpca_object = fpca_gauss(Y, npc = 1, Kt = Kt, periodic = TRUE, print.iter = TRUE)
	expect_equal(dim(fpca_object$efunctions), c(200, 1))
})

test_that("fpca_gauss has correct number of dimensions when applied on incomplete curves",{
	Y = registr::growth_incomplete
	Kt = 8
	fpca_object = fpca_gauss(Y, npc = 2, Kt = Kt, print.iter = TRUE)
	
	expect_equal(dim(fpca_object$alpha), c(length(unique(Y$index)), 1))
	expect_equal(dim(fpca_object$efunctions), c(length(unique(Y$index)), 2))
	expect_equal(length(fpca_object$evalues), 2)
	expect_equal(dim(fpca_object$scores), c(length(unique(Y$id)), 2))
	expect_equal(length(fpca_object$knots), Kt - 4)
	
	fpca_object = fpca_gauss(Y, npc = 1, Kt = Kt, print.iter = TRUE)
	expect_equal(dim(fpca_object$efunctions), c(length(unique(Y$index)), 1))
})
