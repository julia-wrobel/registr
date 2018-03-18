context("bfpca")

test_that("bfpca only accepts binary input values",{
	Y = simulate_functional_data()$Y
	Y$value = Y$value + 2
	
	expect_that(bfpca(Y), 
							throws_error("'binomial' family requires data with binary values of 0 or 1"))
	
	Y$value = Y$latent_mean
	
	expect_that(bfpca(Y), 
							throws_error("'binomial' family requires data with binary values of 0 or 1"))
})

test_that("bfpca output is a list with non-null items and class fpca",{
	Y = simulate_functional_data()$Y
	bfpca_object = bfpca(Y, npc = 2, print.iter = TRUE)
	
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

test_that("bfpca works for time domains other than (0, 1)",{
	Y = simulate_functional_data()$Y
	Y$index = Y$index + 1
	bfpca_object = bfpca(Y, npc = 2, print.iter = TRUE)
	t_min = min(Y$index)
	t_max = max(Y$index)
	
	expect_equal(range(Y$index), range(bfpca_object$Yhat$index))
	expect_equal(range(Y$index), range(bfpca_object$Y$index))
})


test_that("bfpca works for subjects with different grid lengths",{
	Y = simulate_functional_data(vary_D = TRUE)$Y
	bfpca_object = bfpca(Y, npc = 2, print.iter = TRUE)
	
	expect_equal(class(bfpca_object), "fpca")
	expect_true(length(unique(table(Y$id))) > 1)
})


test_that("bfpca works for different start seeds",{
	Y = simulate_functional_data(vary_D = TRUE)$Y
	seeds = as.integer(runif(3, 10, 100))
	
	expect_error(bfpca(Y, npc = 2, print.iter = TRUE, seed = seeds[1]), NA)
	expect_error(bfpca(Y, npc = 2, print.iter = TRUE, seed = seeds[2]), NA)
	expect_error(bfpca(Y, npc = 2, print.iter = TRUE, seed = seeds[3]), NA)
})

test_that("bfpca iterations are strictly decreasing",{
	Y = simulate_functional_data()$Y
	bfpca_object = bfpca(Y, npc = 2, print.iter = TRUE)
	
	expect_true(all(diff(bfpca_object$error) < 0))
})