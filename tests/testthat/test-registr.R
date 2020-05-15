context("registr")

## do both binary and gaussian checks for the following
test_that("registr function returns items of expected length", {
	I = 50
	Kt = 4
	Kh = 4
	Y = simulate_functional_data(I = I)$Y
	registr_object = registr(Y = Y, family = "binomial", Kt = Kt, Kh = Kh)
	
	expect_equal(dim(registr_object$beta), c(Kh-2, I))
	
	Y$value = Y$latent_mean
	expect_error(registr(Y = Y, family = "gaussian", Kt = Kt, Kh = Kh), NA)
})

test_that("registr function throws error when invalid parameters are input", {
	Y = simulate_functional_data()$Y

	expect_error(registr(Y = Y, family = "gaussian", Kt = 2),"Kt must be greater than or equal to 4.")
	expect_error(registr(Y = Y, family = "gaussian", Kh = 2),"Kh must be greater than or equal to 4.")
	expect_error(registr(Y = Y, family = "gaussian", warping = "test"), "warping argument can only take values of nonparametric or piecewise_linear2")
})

test_that("registr function works for time domains other than (0, 1)", {
	Y = simulate_functional_data()$Y
	Y$index = Y$index + 5
	
	registr_object = registr(Y = Y, family = "binomial")
	expect_error(registr_object, NA)
	expect_equal(range(Y$index), range(registr_object$Y$index))
})

test_that("registr function works for subjects with different grid lengths",{
	Y = simulate_functional_data(seed = 40, vary_D = TRUE)$Y
	
	expect_true(length(unique(table(Y$id))) > 1)
	expect_error(registr(Y = Y, family = "binomial"), NA)
	
	Y$value = Y$latent_mean
	expect_error(registr(Y = Y, family = "gaussian"), NA)
})

test_that("registr function requires an index_scaled variable with parametric warping",{
	Y = simulate_functional_data()$Y

	expect_error(registr(Y = Y, family = "binomial", warping = "piecewise_linear2"), 
							 "For warping = piecewise_linear2, need an index_scaled variable that ranges from 0 to 1.")
})

test_that("registr function only uses gradient when periodic = FALSE and warping = nonparametric",{
	Y = simulate_functional_data()$Y
	data = data_clean(Y)
	Y = data$Y
	
	expect_warning(registr(Y = Y, family = "binomial", gradient = TRUE, periodic = TRUE), 
								 "gradient = TRUE is only available for periodic = FALSE. Setting gradient = FALSE.")
	expect_warning(registr(Y = Y, family = "binomial", gradient = TRUE, warping = "piecewise_linear2"), 
								 "gradient = TRUE is only available for warping = nonparametric. Setting gradient = FALSE.")
})

test_that("registr function only accepts warping types nonparametric and piecewise_linear2",{
	Y = simulate_functional_data()$Y
	
	expect_error(registr(Y = Y, family = "binomial", warping = "test"), 
							 "warping argument can only take values of nonparametric or piecewise_linear2")
})