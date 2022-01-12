context("register_fpca")


test_that("code only accepts supported family of distributions", {
	Y = simulate_unregistered_curves(seed = 2020, I = 10, D = 50)
	
	register_fpca(Y, npc = 1, family = "binomial", max_iterations = 5)
	register_fpca(Y, npc = 1, family = "gaussian", max_iterations = 5)
								 
	Y$value = Y$value + 1 # make y values strictly positive for gamma and poisson family
	expect_warning(register_fpca(Y, npc = 1, family = "gamma", fpca_type = "two-step", gradient = FALSE,
	                             max_iterations = 1, fpca_index_significantDigits = 2),
								 "Convergence not reached. Try increasing max_iterations.")
	expect_warning(register_fpca(Y, npc = 1, family = "poisson", fpca_type = "two-step", gradient = FALSE,
	                             max_iterations = 1, fpca_index_significantDigits = 2),
								 "Convergence not reached. Try increasing max_iterations.")
	
	expect_error(register_fpca(Y, npc = 1, family = "exponential"),
							 "Package currently handles only families 'gaussian', 'binomial', 'gamma' and 'poisson'.")
	expect_error(register_fpca(Y, npc = 1, family = 25),
							 "Package currently handles only families 'gaussian', 'binomial', 'gamma' and 'poisson'.")
})

test_that("registering binary data throws no errors",{
	Y = simulate_unregistered_curves(seed = 10001)
	register_fpca(Y, npc = 1, family = "binomial", fpca_type = "variationalEM", max_iterations = 2)
	expect_warning(register_fpca(Y, npc = 1, family = "binomial", fpca_type = "two-step",
	                             max_iterations = 1, fpca_index_significantDigits = 2),
	               "Convergence not reached. Try increasing max_iterations.")
	
	dat = Y$index
	expect_error(register_fpca(dat, npc = 1, family = "binomial"), 
							 "Input dataset must have variables 'id', 'index', and 'value'.")
})

test_that("registering gaussian data throws no errors",{
	Y = simulate_unregistered_curves(seed = 10001)
	register_fpca(Y, npc = 1, family = "gaussian", fpca_type = "variationalEM", max_iterations = 2)
	expect_warning(register_fpca(Y, npc = 1, family = "gaussian", fpca_type = "two-step",
	                             max_iterations = 2, fpca_index_significantDigits = 2),
								 "Convergence not reached. Try increasing max_iterations.")
})

test_that("registering gamma data throws no errors",{
	Y       = simulate_unregistered_curves(seed = 10001)
	Y$value = Y$value + 1 # make y values strictly positive for gamma family
	expect_warning(register_fpca(Y, npc = 1, family = "gamma", fpca_type = "two-step", gradient = FALSE,
	                             max_iterations = 1, fpca_index_significantDigits = 2),
								 "Convergence not reached. Try increasing max_iterations.")
})

test_that("registering poisson data throws no errors",{
	Y       = simulate_unregistered_curves(seed = 10001)
	Y$value = Y$value + 1 # make y values nonnegative for poisson family
	expect_warning(register_fpca(Y, npc = 1, family = "poisson", fpca_type = "two-step", gradient = FALSE,
	                             max_iterations = 1, fpca_index_significantDigits = 2),
								 "Convergence not reached. Try increasing max_iterations.")
})

test_that("register_fpca output is a list with non-null items and class registration",{
  skip_on_cran()
	Y = simulate_unregistered_curves(seed = 1001)
	registr_object1  = register_fpca(Y, npc = 1, family = "binomial", fpca_type = "variationalEM", max_iterations = 2)
	registr_object1a = register_fpca(Y, npc_criterion = 0.5, family = "binomial", fpca_type = "variationalEM", max_iterations = 2)
	expect_warning({
	  registr_object2 = register_fpca(Y, npc = 1, family = "binomial", fpca_type = "two-step", max_iterations = 1, fpca_index_significantDigits = 2)
	}, "Convergence not reached. Try increasing max_iterations.")
	expect_warning({
	  registr_object2a = register_fpca(Y, npc_criterion = c(0.5,0.02), family = "binomial", fpca_type = "two-step", max_iterations = 1, fpca_index_significantDigits = 2)
	}, "Convergence not reached. Try increasing max_iterations.")
	
	expect_equal(class(registr_object1), "registration")
	expect_equal(class(registr_object1a), "registration")
	expect_equal(class(registr_object2), "registration")
	expect_equal(class(registr_object2a), "registration")
	
	expect_false(any(is.na(registr_object1$loss)))
	expect_false(any(is.na(registr_object1a$loss)))
	expect_false(any(is.na(registr_object2$loss)))
	expect_false(any(is.na(registr_object2a$loss)))
	expect_false(any(is.na(registr_object1$time_warps)))
	expect_false(any(is.na(registr_object1a$time_warps)))
	expect_false(any(is.na(registr_object2$time_warps)))
	expect_false(any(is.na(registr_object2a$time_warps)))
	expect_false(any(is.na(registr_object1$Y$t_hat)))
	expect_false(any(is.na(registr_object1a$Y$t_hat)))
	expect_false(any(is.na(registr_object2$Y$t_hat)))
	expect_false(any(is.na(registr_object2a$Y$t_hat)))
})

test_that("register_fpca function forces gradient = FALSE in some situations",{
  skip_on_cran()
	Y       = simulate_unregistered_curves(seed = 2020)
	data    = data_clean(Y)
	Y       = data$Y
	
	expect_warning(register_fpca(Y, npc = 1, family = "binomial", warping = "piecewise_linear2", gradient = TRUE, max_iterations = 1), 
								 "gradient = TRUE is only available for warping = nonparametric. Setting gradient = FALSE.")
	expect_warning(register_fpca(Y, npc = 1, family = "binomial", periodic = TRUE, gradient = TRUE, max_iterations = 1), 
								 "gradient = TRUE is only available for periodic = FALSE. Setting gradient = FALSE.")
	
	Y$value = Y$value + 1 # make y values strictly positive for gamma and poisson family
	expect_warning(register_fpca(Y, npc = 1, family = "gamma", fpca_type = "two-step", gradient = TRUE, max_iterations = 1, fpca_index_significantDigits = 2), 
								 "gradient = TRUE is only available for families 'gaussian' and 'binomial'. Setting gradient = FALSE.")
	expect_warning(register_fpca(Y, npc = 1, family = "poisson", fpca_type = "two-step", gradient = TRUE, max_iterations = 1, fpca_index_significantDigits = 2), 
								 "gradient = TRUE is only available for families 'gaussian' and 'binomial'. Setting gradient = FALSE.")
})

test_that("register_fpca function throws a warning when fpca_type = 'variationalEM' is called for families other than gaussian or binomial",{
  skip_on_cran()
	Y = registr::growth_incomplete
	Y$value = Y$value + 1 # make y values strictly positive for gamma and poisson family
	
	expect_warning(register_fpca(Y, npc = 1, family = "gamma", fpca_type = "variationalEM", gradient = FALSE,
	                             max_iterations = 1, fpca_index_significantDigits = 2),
								 "fpca_type = 'variationalEM' is only available for families 'gaussian' and 'binomial'. Calling variationalEM for 'gaussian' family.")
	expect_warning(register_fpca(Y, npc = 1, family = "poisson", fpca_type = "variationalEM", gradient = FALSE,
	                             max_iterations = 1, fpca_index_significantDigits = 2),
								 "fpca_type = 'variationalEM' is only available for families 'gaussian' and 'binomial'. Calling variationalEM for 'gaussian' family.")
})

test_that("register_fpca function priors must be specified only when warping = piecewise_linear2 and family = binomial",{
  skip_on_cran()
  Y = simulate_unregistered_curves(seed = 2020)
	data = data_clean(Y)
	Y = data$Y
	
	expect_error(register_fpca(Y = Y, npc = 1, family = "binomial", warping = "nonparametric", 
														 priors = TRUE, prior_sd = 1),
							 "priors are only available for warping = piecewise_linear2")
	
	expect_error(register_fpca(Y = Y, npc = 1, family = "gaussian", warping = "piecewise_linear2", gradient = FALSE, 
														 priors = TRUE, prior_sd = 1),
							 "priors are only available for family = binomial")
})

test_that("register_fpca function with priors on the piecewise_linear2 warping functions works as expected",{
  skip_on_cran()
  Y = simulate_unregistered_curves(seed = 2020)
	data = data_clean(Y)
	Y = data$Y
	
	expect_error(register_fpca(Y = Y, npc = 1, family = "binomial", warping = "piecewise_linear2", 
														 gradient = FALSE, priors = TRUE, prior_sd = NULL,
														 max_iterations = 1),
							 "priors = TRUE but no prior_sd supplied.")

	expect_warning(register_fpca(Y = Y, npc = 1, family = "binomial", warping = "piecewise_linear2", 
															 gradient = FALSE, priors = FALSE, prior_sd = 1,
															 max_iterations = 1),
								 "prior_sd supplied but priors = FALSE. No priors included.")
	
	expect_warning({test1 = register_fpca(Y = Y, npc = 1, family = "binomial", warping = "piecewise_linear2",
																				gradient = FALSE, priors = TRUE, prior_sd = 1,
																				max_iterations = 1)},
								 "Convergence not reached. Try increasing max_iterations.")
	expect_warning({test2 = register_fpca(Y = Y, npc = 1, family = "binomial", warping = "piecewise_linear2",
																				gradient = FALSE, priors = TRUE, prior_sd = 0.1,
																				max_iterations = 1)},
								 "Convergence not reached. Try increasing max_iterations.")

	expect_true(all(rowMeans(abs(test2$hinv_beta[c(1:2),] - 0.25)) < rowMeans(abs(test1$hinv_beta[c(1:2),] - 0.25))))
	expect_true(all(rowMeans(abs(test2$hinv_beta[c(3:4),] - 0.75)) < rowMeans(abs(test1$hinv_beta[c(3:4),] - 0.75))))
})


test_that("register_fpca with incompleteness = NULL leads to starting points and endpoints on the diagonal",{
	Y   = registr::growth_incomplete
	expect_warning({
		reg = register_fpca(Y = Y, npc = 1, family = "gaussian", incompleteness = NULL, max_iterations = 1)
	}, "Convergence not reached. Try increasing max_iterations.")
	
	t_min_observed   = tapply(X = reg$Y$tstar, INDEX = reg$Y$id, FUN = min)
	t_max_observed   = tapply(X = reg$Y$tstar, INDEX = reg$Y$id, FUN = max)
	t_min_registered = tapply(X = reg$Y$index, INDEX = reg$Y$id, FUN = min)
	t_max_registered = tapply(X = reg$Y$index, INDEX = reg$Y$id, FUN = max)
	expect_equal(t_min_registered, expected = t_min_observed)
	expect_equal(t_max_registered, expected = t_max_observed)
})

test_that("register_fpca with incompleteness != NULL and lambda_inc = 0: warping functions do not exceed the overall time domain",{
  skip_on_cran()
  Y   = registr::growth_incomplete
	expect_warning({
		reg1 = register_fpca(Y = Y, npc = 1, family = "gaussian", incompleteness = "leading", lambda_inc = 0, max_iterations = 1)
	}, "Convergence not reached. Try increasing max_iterations.")
	expect_warning({
		reg2 = register_fpca(Y = Y, npc = 1, family = "gaussian", incompleteness = "trailing", lambda_inc = 0, max_iterations = 1)
	}, "Convergence not reached. Try increasing max_iterations.")
	expect_warning({
		reg3 = register_fpca(Y = Y, npc = 1, family = "gaussian", incompleteness = "full", lambda_inc = 0, max_iterations = 1)
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

test_that("register_fpca with incompleteness = 'full': higher lambda_inc values lead to registered domain lengths closer to the observed domain lengths",{
	Y = registr::growth_incomplete
	expect_warning({
		reg1 = register_fpca(Y = Y, npc = 1, family = "gaussian",
												 incompleteness = "full", lambda_inc = 1, max_iterations = 1)
	}, "Convergence not reached. Try increasing max_iterations.")
	expect_warning({
		reg2 = register_fpca(Y = Y, npc = 1, family = "gaussian",
												 incompleteness = "full", lambda_inc = 10, max_iterations = 1)
	}, "Convergence not reached. Try increasing max_iterations.")
	
	# compare the average squared difference between registered and observed domain lengths
	lengths_obs  = Y      %>% group_by(id) %>% summarize(length = diff(range(index))) %>% pull(length)
	lengths_reg1 = reg1$Y %>% group_by(id) %>% summarize(length = diff(range(t_hat))) %>% pull(length)
	lengths_reg2 = reg2$Y %>% group_by(id) %>% summarize(length = diff(range(t_hat))) %>% pull(length)
	MSE_reg1     = mean((lengths_reg1 - lengths_obs)^2) 
	MSE_reg2     = mean((lengths_reg2 - lengths_obs)^2) 
	expect_gt(MSE_reg1, expected = MSE_reg2)
})


test_that("register_fpca only accepts Y_template if it has the correct format.",{
	Y = registr::growth_incomplete
	
	template_ids1 = "girl18"
	Y_template1   = Y[Y$id %in% template_ids1,]
	template_ids2 = "girl12"
	Y_template2   = Y[Y$id %in% template_ids2,]
	template_ids3 = c("girl12","boy01","boy04")
	Y_template3   = Y[Y$id %in% template_ids3,]
	
	expect_error(register_fpca(Y = Y, npc = 1, family = "gaussian", Y_template = Y_template1$value),
							 "Y_template must have variables 'id', 'index', and 'value'.")
	expect_error(register_fpca(Y = Y, npc = 1, family = "gaussian", Y_template = Y_template1),
							 "The range of 'index' must be equal for Y_template and Y.")
	expect_warning(register_fpca(Y = Y, npc = 1, family = "gaussian", Y_template = Y_template2, max_iterations = 1),
								 "Convergence not reached. Try increasing max_iterations.")
	expect_warning(register_fpca(Y = Y, npc = 1, family = "gaussian", Y_template = Y_template3, max_iterations = 1),
								 "Convergence not reached. Try increasing max_iterations.")
})
