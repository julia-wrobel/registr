context("registr")


## binary tests
## gaussian tests

## do both binary and gaussian checks for the following

test_that("registr function returns non-null items", {
	## tests stuff it returns
})


test_that("registr function returns items of expected length", {
	## tests stuff it returns
	## test based on values of Kt, Kh that are input
})

test_that("registr function works for time domains other than (0, 1)", {
	
})


test_that("registr function works for subjects with different grid lengths",{
	Y = simulate_functional_data(vary_D = TRUE)$Y
	#bfpca_object = bfpca(Y, npc = 2, print.iter = TRUE)
	#expect_true(length(unique(table(Y$id))) > 1)
})
