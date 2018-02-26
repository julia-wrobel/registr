context("register_fpca")


# for all tests implement for both gaussian and binary data
test_that("code only accepts supported family of distributions", {
	Y = simulate_unregistered_curves()
	registr_object = register_fpca(Y, family = "binomial", max_iterations = 3)
	
	expect_error(register_fpca(Y, family = "gaussian", max_iterations = 3), NA)
	expect_error(register_fpca(Y, family = "poisson"),
							 "Package currently handles only 'binomial' or 'gaussian' families.")

	expect_error(register_fpca(Y, family = 25),
							 "Package currently handles only 'binomial' or 'gaussian' families.")
})

test_that("registering binary and gaussian data throws no errors",{
	Y = simulate_unregistered_curves(seed = 10001)
	expect_error(register_fpca(Y, family = "binomial"), NA)
	
	dat = Y$index
	expect_error(register_fpca(dat, family = "binomial"), 
							 "Input dataset must have variables 'id', 'index', and 'value'.")
})

test_that("register_fpca output is a list with non-null items and class registration",{
	Y = simulate_unregistered_curves(seed = 1001)
	registr_object = register_fpca(Y, family = "binomial", max_iterations = 10)
	
	expect_equal(class(registr_object), "registration")
	
	expect_false(any(is.na(registr_object$loss)))
	expect_false(any(is.na(registr_object$time_warps)))
	expect_false(any(is.na(registr_object$Y$t_hat)))
	
	# throws error in code:
	# Y = simulate_unregistered_curves(seed = 10001)
	# registr_object = register_fpca(Y, family = "binomial", max_iterations = 25)
	
	
})
