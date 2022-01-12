context("bfpca")

test_that("bfpca only accepts binary input values",{
	Y = simulate_functional_data(seed = 2020, I = 10, D = 50)$Y
	Y$value = Y$value + 2
	
	expect_that(bfpca(Y, npc = 1), 
							throws_error("'binomial' family requires data with binary values of 0 or 1"))
	
	Y$value = Y$latent_mean
	
	expect_that(bfpca(Y, npc = 1), 
							throws_error("'binomial' family requires data with binary values of 0 or 1"))
})

test_that("bfpca output is a list with non-null items and class fpca",{
	Y = simulate_functional_data(seed = 2020, I = 10, D = 50)$Y
	expect_warning({
		bfpca_object = bfpca(Y, npc = 2)
	}, "BFPCA convergence not reached. Try increasing maxiter.")
	expect_warning({
	  bfpca_object2 = bfpca(Y, npc_varExplained = 0.8, maxiter = 2)
	}, "BFPCA convergence not reached. Try increasing maxiter.")
	
	expect_equal(class(bfpca_object),  "fpca")
	expect_equal(class(bfpca_object2), "fpca")
	expect_equal(bfpca_object$family,  "binomial")
	expect_equal(bfpca_object2$family, "binomial")
	
	expect_false(any(is.na(bfpca_object$mu)))
	expect_false(any(is.na(bfpca_object2$mu)))
	expect_false(any(is.na(bfpca_object$efunctions)))
	expect_false(any(is.na(bfpca_object2$efunctions)))
	expect_false(any(is.na(bfpca_object$evalues)))
	expect_false(any(is.na(bfpca_object2$evalues)))
	expect_false(any(is.na(bfpca_object$scores)))
	expect_false(any(is.na(bfpca_object2$scores)))
})

test_that("bfpca output has correct number of dimensions",{
	Y = simulate_functional_data(I = 100, D = 200, seed = 2020)$Y
	Kt = 8
	expect_warning({
		bfpca_object = bfpca(Y, npc = 4, Kt = Kt)
	}, "BFPCA convergence not reached. Try increasing maxiter.")
	expect_warning({
	  bfpca_object2 = bfpca(Y, npc_varExplained = 0.8, Kt = Kt, maxiter = 2)
	}, "BFPCA convergence not reached. Try increasing maxiter.")
	
	expect_equal(dim(bfpca_object$alpha),  c(100, 1))
	expect_equal(dim(bfpca_object2$alpha), c(100, 1))
	expect_equal(dim(bfpca_object$efunctions),  c(100, 4))
	expect_equal(dim(bfpca_object2$efunctions), c(100, 4))
	expect_equal(length(bfpca_object$evalues),  4)
	expect_equal(length(bfpca_object2$evalues), 4)
	expect_equal(dim(bfpca_object$scores),  c(100, 4))
	expect_equal(dim(bfpca_object2$scores), c(100, 4))
	expect_equal(length(bfpca_object$knots),  Kt - 4)
	expect_equal(length(bfpca_object2$knots), 20 - 4)
	
	bfpca_object = bfpca(Y, npc = 1, Kt = Kt)
	expect_equal(dim(bfpca_object$efunctions), c(100, 1))
})

test_that("bfpca works for time domains other than (0, 1)",{
	Y = simulate_functional_data(seed = 2020, I = 10, D = 50)$Y
	Y$index = Y$index + 1
	expect_warning({
		bfpca_object = bfpca(Y, npc = 2)
	}, "BFPCA convergence not reached. Try increasing maxiter.")
	t_min = min(Y$index)
	t_max = max(Y$index)
	
	expect_equal(range(Y$index), range(bfpca_object$Yhat$index))
	expect_equal(range(Y$index), range(bfpca_object$Y$index))
})


test_that("bfpca works for subjects with different grid lengths",{
	Y = simulate_functional_data(vary_D = TRUE, seed = 2020)$Y
	expect_warning({
	  bfpca_object = bfpca(Y, npc = 2)
	}, "BFPCA convergence not reached. Try increasing maxiter.")
	
	expect_equal(class(bfpca_object), "fpca")
	expect_true(length(unique(table(Y$id))) > 1)
})


test_that("bfpca works for different start seeds",{
	Y = simulate_functional_data(vary_D = TRUE, seed = 2020)$Y
	seeds = as.integer(runif(3, 10, 100))
	
	expect_warning({
	  expect_error(bfpca(Y, npc = 2, seed = seeds[1]), NA)
	}, "BFPCA convergence not reached. Try increasing maxiter.")
	expect_error(bfpca(Y, npc = 2, seed = seeds[2]), NA)
	bfpca(Y, npc = 2, seed = seeds[3])
})


test_that("bfpca output has correct number of dimensions when periodic = TRUE ",{
	Y = simulate_functional_data(I = 100, D = 200, seed = 2020)$Y
	Kt = 8
	expect_warning({
		bfpca_object = bfpca(Y, npc = 4, periodic = TRUE, Kt = Kt)
	}, "BFPCA convergence not reached. Try increasing maxiter.")
	
	expect_equal(dim(bfpca_object$alpha), c(100, 1))
	expect_equal(dim(bfpca_object$efunctions), c(100, 4))
	expect_equal(length(bfpca_object$evalues), 4)
	expect_equal(dim(bfpca_object$scores), c(100, 4))
	expect_equal(length(bfpca_object$knots), Kt - 1)
	
	expect_warning({
		bfpca_object = bfpca(Y, npc = 1, periodic = TRUE)
	}, "BFPCA convergence not reached. Try increasing maxiter.")
	expect_equal(dim(bfpca_object$efunctions), c(100, 1))
})

test_that("bfpca has correct number of dimensions when applied on incomplete curves",{
	# simulate incompleteness by cutting-off some ids at some point
	Y  = simulate_functional_data(I = 100, D = 200, seed = 2020)$Y
	Y  = Y[!(Y$id %in% unique(Y$id)[1:50]) | (Y$index < 0.5),]
	Kt = 8
	expect_warning({
	  fpca_object = bfpca(Y, npc = 2, Kt = Kt)
	}, "BFPCA convergence not reached. Try increasing maxiter.")
	
	expect_equal(dim(fpca_object$alpha), c(100, 1))
	expect_equal(dim(fpca_object$efunctions), c(100, 2))
	expect_equal(length(fpca_object$evalues), 2)
	expect_equal(dim(fpca_object$scores), c(length(unique(Y$id)), 2))
	expect_equal(length(fpca_object$knots), Kt - 4)
	
	expect_warning({
		fpca_object = bfpca(Y, npc = 1, Kt = Kt)
	}, "BFPCA convergence not reached. Try increasing maxiter.")
	expect_equal(dim(fpca_object$efunctions), c(100, 1))
})

