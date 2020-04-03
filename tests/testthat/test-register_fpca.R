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

test_that("register_fpca function priors must be vectors of length 2 and
					specified only when warping = piecewise_linear2 and family = binomial",{
	Y = simulate_unregistered_curves()
	data = data_clean(Y)
	Y = data$Y
	
	expect_error(register_fpca(Y = Y, family = "binomial", warping = "piecewise_linear2",
														 gradient = FALSE,
												     prior_1_x_mean_sd = c(0.5, 1), prior_1_y_mean_sd = c(0.5, 1),
												     prior_2_x_mean_sd = 1, prior_2_y_mean_sd = c(0.5, 1)),
							 "'prior_2_x_mean_sd' must be NULL or a vector of length 2")
	
	expect_error(register_fpca(Y = Y, family = "binomial", warping = "nonparametric",
														 prior_1_x_mean_sd = c(0.5, 1), prior_1_y_mean_sd = c(0.5, 1),
														 prior_2_x_mean_sd = 1, prior_2_y_mean_sd = c(0.5, 1)),
							 "'prior' arguments are only available for warping = piecewise_linear2")
	
	expect_error(register_fpca(Y = Y, family = "gaussian", warping = "piecewise_linear2",
														 gradient = FALSE,
														 prior_1_x_mean_sd = c(0.5, 1), prior_1_y_mean_sd = c(0.5, 1),
														 prior_2_x_mean_sd = 1, prior_2_y_mean_sd = c(0.5, 1)),
							 "'prior' arguments are only available for family = binomial")
})

test_that("register_fpca function forces gradient = FALSE when using warping = piecewise_linear2 or periodic = TRUE",{
	Y = simulate_unregistered_curves()
	data = data_clean(Y)
	Y = data$Y

	expect_warning(register_fpca(Y, family = "binomial",  warping = "piecewise_linear2", gradient = TRUE), 
								 "gradient = TRUE is only available for warping = nonparametric. Setting gradient = FALSE.")
	expect_warning(register_fpca(Y, family = "binomial",  periodic = TRUE, gradient = TRUE), 
								 "gradient = TRUE is only available for periodic = FALSE. Setting gradient = FALSE.")
})

test_that("register_fpca function with priors on the piecewise_linear2 warping functions works as expected",{
	Y = simulate_unregistered_curves()
	data = data_clean(Y)
	Y = data$Y
	
	prior_mean = 0.5
	
	test1 = register_fpca(Y = Y, family = "binomial", warping = "piecewise_linear2",
												gradient = FALSE,
												prior_1_x_mean_sd = c(prior_mean, 1), prior_1_y_mean_sd = c(prior_mean, 1),
												prior_2_x_mean_sd = c(prior_mean, 1), prior_2_y_mean_sd = c(prior_mean, 1))
	
	test2 = register_fpca(Y = Y, family = "binomial", warping = "piecewise_linear2",
												gradient = FALSE,
												prior_1_x_mean_sd = c(prior_mean, 0.1), prior_1_y_mean_sd = c(prior_mean, 0.1),
												prior_2_x_mean_sd = c(prior_mean, 0.1), prior_2_y_mean_sd = c(prior_mean, 0.1))
	
	expect_true(all(colMeans(abs(test2$beta[,c(1:4)] - prior_mean)) < colMeans(abs(test1$beta[,c(1:4)] - prior_mean))))
})

test_that("register_fpca for binary data with parametric warping functions and/or periodic basis functions throws no errors",{
	Y = simulate_unregistered_curves()
	data = data_clean(Y)
	Y = data$Y
	
	expect_error(register_fpca(Y, family = "binomial", periodic = TRUE,  warping = "nonparametric",     gradient = FALSE), NA)
	expect_error(register_fpca(Y, family = "binomial", periodic = FALSE, warping = "piecewise_linear2", gradient = FALSE), NA)
	expect_error(register_fpca(Y, family = "binomial", periodic = TRUE,  warping = "piecewise_linear2", gradient = FALSE), NA)
})
