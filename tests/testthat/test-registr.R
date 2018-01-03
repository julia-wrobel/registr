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
	Y = simulate_functional_data()$Y
	Y$index = Y$index + 5
	expect_error(registr(Y = Y, family = "binomial"), NA)
})


test_that("registr function works for subjects with different grid lengths",{
	Y = simulate_functional_data(seed = 40, vary_D = TRUE)$Y
	
	expect_true(length(unique(table(Y$id))) > 1)
	expect_error(registr(Y = Y, family = "binomial"), NA)
})
