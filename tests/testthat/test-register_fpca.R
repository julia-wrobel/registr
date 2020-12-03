context("register_fpca")


# for all tests implement for both gaussian and binary data
test_that("code only accepts supported family of distributions", {
	Y       = simulate_unregistered_curves(seed = 2020)
	
	expect_warning(register_fpca(Y, family = "binomial", max_iterations = 3),
								 "Convergence not reached. Try increasing max_iterations.")
	expect_warning(register_fpca(Y, family = "gaussian", max_iterations = 3),
								 "Convergence not reached. Try increasing max_iterations.")
								 
	Y$value = Y$value + 1 # make y values strictly positive for gamma family
	expect_warning(register_fpca(Y, family = "gamma", max_iterations = 3, fpca_index_relevantDigits = 2),
								 "Convergence not reached. Try increasing max_iterations.")
	
	expect_error(register_fpca(Y, family = "poisson"),
							 "Package currently handles only families 'gaussian', 'binomial' and 'gamma'.")
	expect_error(register_fpca(Y, family = 25),
							 "Package currently handles only families 'gaussian', 'binomial' and 'gamma'.")
})

test_that("registering binary data throws no errors",{
	Y = simulate_unregistered_curves(seed = 10001)
	expect_warning(register_fpca(Y, family = "binomial", fpca_type = "variationalEM", max_iterations = 3),
								 "Convergence not reached. Try increasing max_iterations.")
	expect_warning(register_fpca(Y, family = "binomial", fpca_type = "two-step", max_iterations = 3),
								 "Convergence not reached. Try increasing max_iterations.")
	
	dat = Y$index
	expect_error(register_fpca(dat, family = "binomial"), 
							 "Input dataset must have variables 'id', 'index', and 'value'.")
})

test_that("registering gaussian data throws no errors",{
	Y = simulate_unregistered_curves(seed = 10001)
	expect_warning(register_fpca(Y, family = "gaussian", fpca_type = "variationalEM", max_iterations = 3),
								 "Convergence not reached. Try increasing max_iterations.")
	expect_warning(register_fpca(Y, family = "gaussian", fpca_type = "two-step", max_iterations = 3, fpca_index_relevantDigits = 2),
								 "Convergence not reached. Try increasing max_iterations.")
})

test_that("registering gamma data throws no errors",{
	Y       = simulate_unregistered_curves(seed = 10001)
	Y$value = Y$value + 1 # make y values strictly positive for gamma family
	expect_warning(register_fpca(Y, family = "gamma", max_iterations = 3, fpca_index_relevantDigits = 2),
								 "Convergence not reached. Try increasing max_iterations.")
})

test_that("register_fpca output is a list with non-null items and class registration",{
	Y = simulate_unregistered_curves(seed = 1001)
	expect_warning({
		registr_object1 = register_fpca(Y, family = "binomial", fpca_type = "variationalEM", max_iterations = 3)
	}, "Convergence not reached. Try increasing max_iterations.")
	expect_warning({
		registr_object2 = register_fpca(Y, family = "binomial", fpca_type = "two-step", max_iterations = 3)
	}, "Convergence not reached. Try increasing max_iterations.")
	
	expect_equal(class(registr_object1), "registration")
	expect_equal(class(registr_object2), "registration")
	
	expect_false(any(is.na(registr_object1$loss)))
	expect_false(any(is.na(registr_object2$loss)))
	expect_false(any(is.na(registr_object1$time_warps)))
	expect_false(any(is.na(registr_object2$time_warps)))
	expect_false(any(is.na(registr_object1$Y$t_hat)))
	expect_false(any(is.na(registr_object2$Y$t_hat)))
})

test_that("register_fpca function forces gradient = FALSE in some situations",{
	Y       = simulate_unregistered_curves(seed = 2020)
	data    = data_clean(Y)
	Y       = data$Y
	
	expect_warning(register_fpca(Y, family = "binomial", warping = "piecewise_linear2", gradient = TRUE, max_iterations = 3), 
								 "gradient = TRUE is only available for warping = nonparametric. Setting gradient = FALSE.")
	expect_warning(register_fpca(Y, family = "binomial", periodic = TRUE, gradient = TRUE, max_iterations = 3), 
								 "gradient = TRUE is only available for periodic = FALSE. Setting gradient = FALSE.")
	
	Y$value = Y$value + 1 # make y values strictly positive for gamma family
	expect_warning(register_fpca(Y, family = "gamma", max_iterations = 3, fpca_index_relevantDigits = 2), 
								 "gradient = TRUE is only available for families 'gaussian' and 'binomial'. Setting gradient = FALSE.")
})

test_that("register_fpca function forces fpca_type = 'two-step' for families other than gaussian or binomial",{
	Y = registr::growth_incomplete
	Y$value = Y$value + 1 # make y values strictly positive for gamma family
	
	expect_warning(register_fpca(Y, family = "gamma", fpca_type = "variationalEM", max_iterations = 3, fpca_index_relevantDigits = 2),
								 "fpca_type = 'variationalEM' is only available for families 'gaussian' and 'binomial'. Setting fpca_type = 'two-step'.")
})

test_that("register_fpca function priors must be specified only when warping = piecewise_linear2 and family = binomial",{
	Y = simulate_unregistered_curves(seed = 2020)
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
	Y = simulate_unregistered_curves(seed = 2020)
	data = data_clean(Y)
	Y = data$Y
	
	expect_error(register_fpca(Y = Y, family = "binomial", warping = "piecewise_linear2", 
														 gradient = FALSE, priors = TRUE, prior_sd = NULL),
							 "priors = TRUE but no prior_sd supplied.")

	expect_warning(register_fpca(Y = Y, family = "binomial", warping = "piecewise_linear2", 
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
	Y = simulate_unregistered_curves(seed = 2020)
	data = data_clean(Y)
	Y = data$Y
	
	expect_error(register_fpca(Y, family = "binomial", periodic = TRUE,  warping = "nonparametric",     gradient = FALSE), NA)
	expect_error(register_fpca(Y, family = "binomial", periodic = FALSE, warping = "piecewise_linear2", gradient = FALSE), NA)
	expect_error(register_fpca(Y, family = "binomial", periodic = TRUE,  warping = "piecewise_linear2", gradient = FALSE), NA)
})

test_that("register_fpca with incompleteness = NULL leads to starting points and endpoints on the diagonal",{
	Y   = registr::growth_incomplete
	expect_warning({
		reg = register_fpca(Y = Y, family = "gaussian", incompleteness = NULL, max_iterations = 2)
	}, "Convergence not reached. Try increasing max_iterations.")
	
	t_min_observed   = tapply(X = reg$Y$tstar, INDEX = reg$Y$id, FUN = min)
	t_max_observed   = tapply(X = reg$Y$tstar, INDEX = reg$Y$id, FUN = max)
	t_min_registered = tapply(X = reg$Y$index, INDEX = reg$Y$id, FUN = min)
	t_max_registered = tapply(X = reg$Y$index, INDEX = reg$Y$id, FUN = max)
	expect_equal(t_min_registered, expected = t_min_observed)
	expect_equal(t_max_registered, expected = t_max_observed)
})

test_that("register_fpca with incompleteness != NULL and lambda_inc = 0: warping functions do not exceed the overall time domain",{
	Y   = registr::growth_incomplete
	expect_warning({
		reg1 = register_fpca(Y = Y, family = "gaussian", incompleteness = "leading", lambda_inc = 0, max_iterations = 2)
	}, "Convergence not reached. Try increasing max_iterations.")
	expect_warning({
		reg2 = register_fpca(Y = Y, family = "gaussian", incompleteness = "trailing", lambda_inc = 0, max_iterations = 2)
	}, "Convergence not reached. Try increasing max_iterations.")
	expect_warning({
		reg3 = register_fpca(Y = Y, family = "gaussian", incompleteness = "full", lambda_inc = 0, max_iterations = 2)
	}, "Convergence not reached. Try increasing max_iterations.")
	
	t_min = min(Y$index)
	t_max = max(Y$index)
	# round the registered times to prevent numerical / machine number errors
	expect_gte(round(min(reg1$Y$t_hat), 5), expected = t_min)
	expect_gte(round(min(reg2$Y$t_hat), 5), expected = t_min)
	expect_gte(round(min(reg3$Y$t_hat), 5), expected = t_min)
	expect_lte(round(max(reg1$Y$t_hat), 5), expected = t_max)
	expect_lte(round(max(reg2$Y$t_hat), 5), expected = t_max)
	expect_lte(round(max(reg3$Y$t_hat), 5), expected = t_max)
})

test_that("register_fpca with incompleteness = 'full': higher lambda_inc values lead to starting points and endpoints closer to the diagonal",{
	Y = registr::growth_incomplete
	expect_warning({
		reg1 = register_fpca(Y = Y, family = "gaussian",
												 incompleteness = "full", lambda_inc = 1, max_iterations = 2)
	}, "Convergence not reached. Try increasing max_iterations.")
	expect_warning({
		reg2 = register_fpca(Y = Y, family = "gaussian",
												 incompleteness = "full", lambda_inc = 10, max_iterations = 2)
	}, "Convergence not reached. Try increasing max_iterations.")
	
	# compare the MSE values of the warping function starting points to the diagonal
	t_min_observed_1   = tapply(X = reg1$Y$tstar, INDEX = reg1$Y$id, FUN = min)
	t_min_registered_1 = tapply(X = reg1$Y$t_hat, INDEX = reg1$Y$id, FUN = min)
	t_min_observed_2   = tapply(X = reg2$Y$tstar, INDEX = reg2$Y$id, FUN = min)
	t_min_registered_2 = tapply(X = reg2$Y$t_hat, INDEX = reg2$Y$id, FUN = min)
	MSE_min1 = sum((t_min_registered_1 - t_min_observed_1)^2)
	MSE_min2 = sum((t_min_registered_2 - t_min_observed_2)^2)
	expect_gt(MSE_min1, expected = MSE_min2)

	# compare the MSE values of the warping function endpoints to the diagonal
	t_max_observed_1   = tapply(X = reg1$Y$tstar, INDEX = reg1$Y$id, FUN = max)
	t_max_registered_1 = tapply(X = reg1$Y$t_hat, INDEX = reg1$Y$id, FUN = max)
	t_max_observed_2   = tapply(X = reg2$Y$tstar, INDEX = reg2$Y$id, FUN = max)
	t_max_registered_2 = tapply(X = reg2$Y$t_hat, INDEX = reg2$Y$id, FUN = max)
	MSE_max1 = sum((t_max_registered_1 - t_max_observed_1)^2)
	MSE_max2 = sum((t_max_registered_2 - t_max_observed_2)^2)
	expect_gt(MSE_max1, expected = MSE_max2)
})

test_that("register_fpca for gamma data throws no error",{
	Y       = registr::growth_incomplete
	Y$value = Y$value + 1 # make y values strictly positive for gamma family
	
	expect_warning({
		reg = register_fpca(Y, family = "gamma",fpca_type = "two-step", max_iterations = 3, fpca_index_relevantDigits = 2, gradient = FALSE)
	}, "Convergence not reached. Try increasing max_iterations.")
	expect_identical(class(reg), "registration")
})

test_that("register_fpca only accepts Y_template if it has the correct format.",{
	Y = registr::growth_incomplete
	
	template_ids1 = "girl18"
	Y_template1   = Y[Y$id %in% template_ids1,]
	template_ids2 = "girl36"
	Y_template2   = Y[Y$id %in% template_ids2,]
	template_ids3 = c("boy01","boy04","boy29","boy30","boy31","boy34","boy36")
	Y_template3   = Y[Y$id %in% template_ids3,]
	
	expect_error(register_fpca(Y = Y, family = "gaussian", Y_template = Y_template1$value),
							 "Y_template must have variables 'id', 'index', and 'value'.")
	expect_error(register_fpca(Y = Y, family = "gaussian", Y_template = Y_template1),
							 "The range of 'index' must be equal for Y_template and Y.")
	expect_error(register_fpca(Y = Y, family = "gaussian", Y_template = Y_template2, max_iterations = 2),
							 NA)
	expect_error(register_fpca(Y = Y, family = "gaussian", Y_template = Y_template3, max_iterations = 2),
							 NA)
})
