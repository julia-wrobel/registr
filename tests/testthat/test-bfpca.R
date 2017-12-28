context("bfpca")


test_that("bfpca only accepts binary input values",{
	# run with acceptable input
	# run with factor input
	# run with unacceptable input
})

test_that("bfpca object is class bfpca and not null",{
	
})


test_that("bfpca works for time values on different grids",{
	# works meaning doesn't break
	# test min and max so that t_min and t_max are necessary
})


test_that("bfpca works for subjects with different grid lengths",{
	# works meaning doesn't break
	# test min and max so that t_min and t_max are necessary
})

test_that("bfpca output has correct number of dimensions",{
	# set I, Kt, npc = 1 and npc = 3, test each output thing
})


test_that("bfpca function works as intended",{
	# set simulated dataset
	# check that true mean is close to simulated mean,
		# check that fitted values are close for some subject
	# ok to use function that generates values?
})