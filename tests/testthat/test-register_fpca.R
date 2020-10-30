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

test_that("register_fpca function forces gradient = FALSE when using warping = piecewise_linear2 or periodic = TRUE",{
	Y = simulate_unregistered_curves()
	data = data_clean(Y)
	Y = data$Y

	expect_warning(register_fpca(Y, family = "binomial",  warping = "piecewise_linear2", gradient = TRUE), 
								 "gradient = TRUE is only available for warping = nonparametric. Setting gradient = FALSE.")
	expect_warning(register_fpca(Y, family = "binomial",  periodic = TRUE, gradient = TRUE), 
								 "gradient = TRUE is only available for periodic = FALSE. Setting gradient = FALSE.")
})

test_that("register_fpca function priors must be specified only when warping = piecewise_linear2 and family = binomial",{
	Y = simulate_unregistered_curves()
	data = data_clean(Y)
	Y = data$Y
	
	expect_error(register_fpca(Y = Y, family = "binomial", warping = "nonparametric", 
														 priors = TRUE, prior_sd = 1),
							 "priors are only available for warping = piecewise_linear2")
	
	expect_error(register_fpca(Y = Y, family = "gaussian", warping = "piecewise_linear2", gradient = FALSE, 
														 priors = TRUE, prior_sd = 1),
							 "priors are only available for family = binomial")
})

test_that("register_fpca function with priors on the piecewise_linear2 warping functions works as expected",{
	Y = simulate_unregistered_curves()
	data = data_clean(Y)
	Y = data$Y
	
	expect_error(register_fpca(Y = Y, family = "binomial", warping = "piecewise_linear2", 
														 gradient = FALSE, priors = TRUE, prior_sd = NULL),
							 "priors = TRUE but no prior_sd supplied.")

	expect_message(register_fpca(Y = Y, family = "binomial", warping = "piecewise_linear2", 
															 gradient = FALSE, priors = FALSE, prior_sd = 1),
								 "prior_sd supplied but priors = FALSE. No priors included.")
	
	test1 = register_fpca(Y = Y, family = "binomial", warping = "piecewise_linear2",
												gradient = FALSE, priors = TRUE, prior_sd = 1)
	
	test2 = register_fpca(Y = Y, family = "binomial", warping = "piecewise_linear2",
												gradient = FALSE, priors = TRUE, prior_sd = 0.1)

	expect_true(all(colMeans(abs(test2$beta[,c(1:2)] - 0.25)) < colMeans(abs(test1$beta[,c(1:2)] - 0.25))))
	expect_true(all(colMeans(abs(test2$beta[,c(3:4)] - 0.75)) < colMeans(abs(test1$beta[,c(3:4)] - 0.75))))
})

test_that("register_fpca for binary data with parametric warping functions and/or periodic basis functions throws no errors",{
	Y = simulate_unregistered_curves()
	data = data_clean(Y)
	Y = data$Y
	
	expect_error(register_fpca(Y, family = "binomial", periodic = TRUE,  warping = "nonparametric",     gradient = FALSE), NA)
	expect_error(register_fpca(Y, family = "binomial", periodic = FALSE, warping = "piecewise_linear2", gradient = FALSE), NA)
	expect_error(register_fpca(Y, family = "binomial", periodic = TRUE,  warping = "piecewise_linear2", gradient = FALSE), NA)
})

test_that("register_fpca with preserve_domain = TRUE leads to endpoints on the diagonal",{
	Y   = registr::growth_incomplete
	reg = register_fpca(Y = Y, family = "gaussian", preserve_domain = TRUE)
	
	t_max_observed   = tapply(X = reg$Y$tstar, INDEX = reg$Y$id, FUN = max)
	t_max_registered = tapply(X = reg$Y$index, INDEX = reg$Y$id, FUN = max)
	expect_equal(t_max_registered, expected = t_max_observed)
})

test_that("register_fpca with preserve_domain = FALSE and lambda_endpoint = 0: warping functions do not exceed the overall time domain",{
	Y   = registr::growth_incomplete
	reg = register_fpca(Y = Y, family = "gaussian", preserve_domain = TRUE)
	
	t_min = min(Y$index)
	t_max = max(Y$index)
	# round the registered times to prevent numerical / machine number errors
	expect_gte(round(min(reg$Y$t_hat), 10), expected = t_min)
	expect_lte(round(max(reg$Y$t_hat), 10), expected = t_max)
})

test_that("register_fpca with preserve_domain = FALSE: higher lambda_endpoint values lead to endpoints closer to the diagonal",{
	Y = registr::growth_incomplete
	reg1 = register_fpca(Y = Y, family = "gaussian", preserve_domain = FALSE,
											 lambda_endpoint = 1)
	reg2 = register_fpca(Y = Y, family = "gaussian", preserve_domain = FALSE,
											 lambda_endpoint = 10)
	
	# compare the MSE values of the warping function endpoints to the diagonal
	t_max_observed_1   = tapply(X = reg1$Y$tstar, INDEX = reg1$Y$id, FUN = max)
	t_max_registered_1 = tapply(X = reg1$Y$t_hat, INDEX = reg1$Y$id, FUN = max)
	t_max_observed_2   = tapply(X = reg2$Y$tstar, INDEX = reg2$Y$id, FUN = max)
	t_max_registered_2 = tapply(X = reg2$Y$t_hat, INDEX = reg2$Y$id, FUN = max)
	MSE1 <- sum((t_max_registered_1 - t_max_observed_1)^2)
	MSE2 <- sum((t_max_registered_2 - t_max_observed_2)^2)
	expect_gt(MSE1, expected = MSE2)
})
