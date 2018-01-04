context("registr")

## do both binary and gaussian checks for the following
test_that("registr function returns items of expected length", {
	I = 50
	Kt = 3
	Kh = 3
	Y = simulate_functional_data(I = I)$Y
	registr_object = registr(Y = Y, family = "binomial", Kt = Kt, Kh = Kh)
	
	expect_equal(dim(registr_object$beta), c(Kh-1, I))
	
	Y$value = Y$latent_mean
	expect_error(registr(Y = Y, family = "gaussian", Kt = Kt, Kh = Kh), NA)
})

test_that("registr function throws error when invalid parameters are input", {
	Y = simulate_functional_data()$Y

	expect_error(registr(Y = Y, family = "gaussian", Kt = 2),"Kt must be greater than or equal to 3.")
	expect_error(registr(Y = Y, family = "gaussian", Kh = 2),"Kh must be greater than or equal to 3.")
})

test_that("registr function works for time domains other than (0, 1)", {
	Y = simulate_functional_data()$Y
	Y$index = Y$index + 5
	expect_error(registr(Y = Y, family = "binomial"), NA)
	
	Y$value = Y$latent_mean
	expect_error(registr(Y = Y, family = "gaussian"), NA)
})


test_that("registr function works for subjects with different grid lengths",{
	Y = simulate_functional_data(seed = 40, vary_D = TRUE)$Y
	
	expect_true(length(unique(table(Y$id))) > 1)
	expect_error(registr(Y = Y, family = "binomial"), NA)
	
	Y$value = Y$latent_mean
	expect_error(registr(Y = Y, family = "gaussian"), NA)
})
