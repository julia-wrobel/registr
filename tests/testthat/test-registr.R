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

test_that("registr function with preserve_domain = TRUE leads to endpoints on the diagonal",{
	Y   = registr::growth_incomplete
	reg = registr(Y = Y, family = "gaussian", preserve_domain = TRUE)
	
	t_max_observed   = tapply(X = reg$Y$tstar, INDEX = reg$Y$id, FUN = max)
	t_max_registered = tapply(X = reg$Y$index, INDEX = reg$Y$id, FUN = max)
	expect_equal(t_max_registered, expected = t_max_observed)
})

test_that("registr function with preserve_domain = FALSE and lambda_endpoint = 0: warping functions do not exceed the overall time domain",{
	Y   = registr::growth_incomplete
	reg = registr(Y = Y, family = "gaussian", preserve_domain = TRUE)
	
	t_min = min(Y$index)
	t_max = max(Y$index)
	# round the registered times to prevent numerical / machine number errors
	expect_gte(round(min(reg$Y$index), 10), expected = t_min)
	expect_lte(round(max(reg$Y$index), 10), expected = t_max)
})

test_that("registr function with preserve_domain = FALSE: higher lambda_endpoint values lead to endpoints closer to the diagonal",{
	Y = registr::growth_incomplete
	reg1 = registr(Y = Y, family = "gaussian", preserve_domain = FALSE,
								 lambda_endpoint = 1)
	reg2 = registr(Y = Y, family = "gaussian", preserve_domain = FALSE,
								 lambda_endpoint = 10)
	
	# compare the MSE values of the warping function endpoints to the diagonal
	t_max_observed_1   = tapply(X = reg1$Y$tstar, INDEX = reg1$Y$id, FUN = max)
	t_max_registered_1 = tapply(X = reg1$Y$index, INDEX = reg1$Y$id, FUN = max)
	t_max_observed_2   = tapply(X = reg2$Y$tstar, INDEX = reg2$Y$id, FUN = max)
	t_max_registered_2 = tapply(X = reg2$Y$index, INDEX = reg2$Y$id, FUN = max)
	MSE1 = sum((t_max_registered_1 - t_max_observed_1)^2)
	MSE2 = sum((t_max_registered_2 - t_max_observed_2)^2)
	expect_gt(MSE1, expected = MSE2)
})

test_that("registr function throws an error when family = 'gamma' is applied to non-strictly positive data",{
	Y = registr::growth_incomplete
	
	expect_error(registr(Y = Y, family = "gamma"),
							 "family = 'gamma' can only be applied to strictly positive data.")
})

test_that("registr function only accepts Y_template if it has the correct format.",{
	Y = registr::growth_incomplete
	
	template_ids1 = "girl18"
	Y_template1   = Y[Y$id %in% template_ids1,]
	template_ids2 = "boy30"
	Y_template2   = Y[Y$id %in% template_ids2,]
	template_ids3 = c("boy01","boy29","boy30","boy31","boy34","boy36")
	Y_template3   = Y[Y$id %in% template_ids3,]
	
	expect_error(registr(Y = Y, family = "gaussian", Y_template = Y_template1$value),
							 "Y_template must have variables 'id', 'index', and 'value'.")
	expect_error(registr(Y = Y, family = "gaussian", Y_template = Y_template1),
							 "The range of 'index' must be equal for Y_template and Y.")
	expect_error(registr(Y = Y, family = "gaussian", Y_template = Y_template2),
							 NA)
	expect_error(registr(Y = Y, family = "gaussian", Y_template = Y_template3),
							 NA)
})
