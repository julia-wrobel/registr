context("registr")


## binary tests
## gaussian tests

## do both binary and gaussian checks for the following

test_that("registr function returns items of expected length", {
	I = 50
	Kt = 11
	Kh = 3
	Y = simulate_functional_data(I = I)$Y
	registr_object = registr(Y = Y, family = "binomial", Kt = Kt, Kh = Kh)
	
	expect_equal(dim(registr_object$beta), c(Kh-1, I))
	
	Y$value = Y$prob
	expect_error(registr(Y = Y, family = "gaussian", Kt = Kt, Kh = Kh), NA)
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
