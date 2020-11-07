context("gfpca_twoStep")

test_that("coarsen_index only accepts positive integer 'relevant_digits'",{
	index_values = 1:3
	
	expect_that(coarsen_index(index_values, relevant_digits = 0),
							throws_error("'relevant_digits' must be a positive integer."))
})

test_that("coarsen_index correctly rounds positive values",{
	index_values1 = c(0.84729, 0.9379, 0.19328)
	index_values2 = c(307, 86938, 13)
	
	expect_equal(coarsen_index(index_values1, relevant_digits = 1),
							 expected = c(0.8, 0.9, 0.2))
	expect_equal(coarsen_index(index_values1, relevant_digits = 3),
							 expected = c(0.847, 0.938, 0.193))
	expect_equal(coarsen_index(index_values2, relevant_digits = 1),
							 expected = c(0, 90000, 0))
	expect_equal(coarsen_index(index_values2, relevant_digits = 4),
							 expected = c(310, 86940, 10))
})

test_that("coarsen_index correctly rounds negative values",{
	index_values1 = c(0.84729, -0.9379, 0.19328)
	index_values2 = c(-307, 86938, -13)
	
	expect_equal(coarsen_index(index_values1, relevant_digits = 1),
							 expected = c(0.8, -0.9, 0.2))
	expect_equal(coarsen_index(index_values1, relevant_digits = 3),
							 expected = c(0.847, -0.938, 0.193))
	expect_equal(coarsen_index(index_values2, relevant_digits = 1),
							 expected = c(0, 90000, 0))
	expect_equal(coarsen_index(index_values2, relevant_digits = 4),
							 expected = c(-310, 86940, -10))
})
