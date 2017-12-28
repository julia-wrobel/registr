context("bfpca")

test_that("bfpca only accepts binary input values",{
	#Y = simulate_functional_data()$Y
	# run with acceptable input
	# run with factor input
	# run with unacceptable input
})

test_that("bfpca output is a list with non-null items and class fpca",{
	Y = simulate_functional_data()$Y
	bfpca_object = bfpca(Y, npc = 2, print.iter = TRUE)
	
	#expect_output(str(bfpca_object), "List of 2") ## figure out final number of items
	
	expect_equal(class(bfpca_object), "fpca")
	expect_equal(bfpca_object$family, "binomial")
	
	expect_false(any(is.na(bfpca_object$mu)))
	expect_false(any(is.na(bfpca_object$efunctions)))
	expect_false(any(is.na(bfpca_object$evalues)))
	expect_false(any(is.na(bfpca_object$scores)))
})

test_that("bfpca output has correct number of dimensions",{
	Y = simulate_functional_data(I = 100, D = 200)$Y
	bfpca_object = bfpca(Y, npc = 4, print.iter = TRUE)
	
	expect_equal(dim(bfpca_object$alpha), c(200, 1))
	expect_equal(dim(bfpca_object$efunctions), c(200, 4))
	expect_equal(length(bfpca_object$evalues), 4)
	expect_equal(dim(bfpca_object$scores), c(100, 4))
	
	bfpca_object = bfpca(Y, npc = 1, print.iter = TRUE)
	expect_equal(dim(bfpca_object$efunctions), c(200, 1))
})


test_that("bfpca works for time values on different grids",{
	Y = simulate_functional_data()$Y
	Y$index = Y$index + 1
	bfpca_object = bfpca(Y, npc = 2, print.iter = TRUE)
	
	# works meaning doesn't break
	# test min and max so that t_min and t_max are necessary
	# this doesnt do what you want - mean alpha is evaluated on wrong grid
	
})


test_that("bfpca works for subjects with different grid lengths",{
	# works meaning doesn't break
	# test min and max so that t_min and t_max are necessary
})


test_that("bfpca function works as intended",{
	# set simulated dataset
	# check that true mean is close to simulated mean,
		# check that fitted values are close for some subject
	# ok to use function that generates values?
})