context("fpca_gauss")

test_that("fpca_gauss output is a list with non-null items and class fpca",{
	Y = simulate_functional_data(seed = 2020)$Y
	Y$value = Y$latent_mean
	
	fpca_object  = fpca_gauss(Y, npc = 2)
	expect_warning({
	  fpca_object2 = fpca_gauss(Y, npc_varExplained = 0.9, maxiter = 2)
	}, "fpca_gauss convergence not reached. Try increasing maxiter.")
	
	expect_equal(class(fpca_object),  "fpca")
	expect_equal(class(fpca_object2), "fpca")
	expect_equal(fpca_object$family,  "gaussian")
	expect_equal(fpca_object2$family, "gaussian")
	
	expect_false(any(is.na(fpca_object$mu)))
	expect_false(any(is.na(fpca_object2$mu)))
	expect_false(any(is.na(fpca_object$efunctions)))
	expect_false(any(is.na(fpca_object2$efunctions)))
	expect_false(any(is.na(fpca_object$evalues)))
	expect_false(any(is.na(fpca_object2$evalues)))
	expect_false(any(is.na(fpca_object$scores)))
	expect_false(any(is.na(fpca_object2$scores)))
})

test_that("fpca_gauss output has correct number of dimensions",{
	Y = simulate_functional_data(I = 100, D = 200, seed = 2020)$Y
	Kt = 8
	expect_warning({
	  fpca_object  = fpca_gauss(Y, npc = 2, Kt = Kt)
	}, "fpca_gauss convergence not reached. Try increasing maxiter.")
	expect_warning({
	  fpca_object2 = fpca_gauss(Y, npc_varExplained = 0.9, maxiter = 2)
	}, "fpca_gauss convergence not reached. Try increasing maxiter.")
	
	expect_equal(dim(fpca_object$alpha),  c(100, 1))
	expect_equal(dim(fpca_object2$alpha), c(100, 1))
	expect_equal(dim(fpca_object$efunctions),  c(100, 2))
	expect_equal(dim(fpca_object2$efunctions), c(100, 2))
	expect_equal(length(fpca_object$evalues),  2)
	expect_equal(length(fpca_object2$evalues), 2)
	expect_equal(dim(fpca_object$scores),  c(100, 2))
	expect_equal(dim(fpca_object2$scores), c(100, 2))
	expect_equal(length(fpca_object$knots),  Kt - 4)
	expect_equal(length(fpca_object2$knots), 20 - 4)
	
	expect_warning({
	  fpca_object = fpca_gauss(Y, npc = 1, Kt = Kt)
	}, "fpca_gauss convergence not reached. Try increasing maxiter.")
	expect_equal(dim(fpca_object$efunctions), c(100, 1))
})

test_that("fpca_gauss output has correct number of dimensions when periodic = TRUE",{
	Y = simulate_functional_data(I = 100, D = 200, seed = 2020)$Y
	Kt = 8
	expect_warning({
	  fpca_object = fpca_gauss(Y, npc = 2, Kt = Kt, periodic = TRUE)
	}, "fpca_gauss convergence not reached. Try increasing maxiter.")
	
	expect_equal(dim(fpca_object$alpha), c(100, 1))
	expect_equal(dim(fpca_object$efunctions), c(100, 2))
	expect_equal(length(fpca_object$evalues), 2)
	expect_equal(dim(fpca_object$scores), c(100, 2))
	expect_equal(length(fpca_object$knots), Kt - 1)
	
	expect_warning({
	  fpca_object = fpca_gauss(Y, npc = 1, Kt = Kt, periodic = TRUE)
	}, "fpca_gauss convergence not reached. Try increasing maxiter.")
	expect_equal(dim(fpca_object$efunctions), c(100, 1))
})

test_that("fpca_gauss has correct number of dimensions when applied on incomplete curves",{
	Y = registr::growth_incomplete
	Kt = 8
	expect_warning({
	  fpca_object = fpca_gauss(Y, npc = 2, Kt = Kt)
	}, "fpca_gauss convergence not reached. Try increasing maxiter.")
	
	expect_equal(dim(fpca_object$alpha), c(100, 1))
	expect_equal(dim(fpca_object$efunctions), c(100, 2))
	expect_equal(length(fpca_object$evalues), 2)
	expect_equal(dim(fpca_object$scores), c(length(unique(Y$id)), 2))
	expect_equal(length(fpca_object$knots), Kt - 4)
	
	expect_warning({
	  fpca_object = fpca_gauss(Y, npc = 1, Kt = Kt)
	}, "fpca_gauss convergence not reached. Try increasing maxiter.")
	expect_equal(dim(fpca_object$efunctions), c(100, 1))
})

