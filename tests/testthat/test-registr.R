context("registr")

## do both binary and gaussian checks for the following
test_that("registr function returns items of expected length", {
	I = 50
	Kt = 4
	Kh = 4
	Y = simulate_functional_data(I = I)$Y
	registr_object = registr(Y = Y, family = "binomial", Kt = Kt, Kh = Kh)
	
	expect_equal(length(registr_object$hinv_innerKnots), I)
	expect_equal(dim(registr_object$hinv_beta), c(Kh, I))
	
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

test_that("registr function only takes valid options for argument 'incompleteness'",{
	Y = registr::growth_incomplete
	expect_error(registr(Y = Y, family = "gaussian", incompleteness = "Obacht!"),
							 "'incompleteness' must be either 'leading', 'trailing' or 'full'.")
})

test_that("registr function with incompleteness = NULL leads to starting points and endpoints on the diagonal",{
	Y   = registr::growth_incomplete
	reg = registr(Y = Y, family = "gaussian", incompleteness = NULL)
	
	t_min_observed   = tapply(X = reg$Y$tstar, INDEX = reg$Y$id, FUN = min)
	t_max_observed   = tapply(X = reg$Y$tstar, INDEX = reg$Y$id, FUN = max)
	t_min_registered = tapply(X = reg$Y$index, INDEX = reg$Y$id, FUN = min)
	t_max_registered = tapply(X = reg$Y$index, INDEX = reg$Y$id, FUN = max)
	expect_equal(t_min_registered, expected = t_min_observed)
	expect_equal(t_max_registered, expected = t_max_observed)
})

test_that("registr function with incompleteness = 'leading' leads to endpoints on the diagonal",{
	Y   = registr::growth_incomplete
	reg = registr(Y = Y, family = "gaussian", incompleteness = "leading", lambda_inc = 0)
	
	t_max_observed   = tapply(X = reg$Y$tstar, INDEX = reg$Y$id, FUN = max)
	t_max_registered = tapply(X = reg$Y$index, INDEX = reg$Y$id, FUN = max)
	expect_equal(t_max_registered, expected = t_max_observed)
})

test_that("registr function with incompleteness = 'trailing' leads to starting points on the diagonal",{
	Y   = registr::growth_incomplete
	reg = registr(Y = Y, family = "gaussian", incompleteness = "trailing", lambda_inc = 0)
	
	t_min_observed   = tapply(X = reg$Y$tstar, INDEX = reg$Y$id, FUN = min)
	t_min_registered = tapply(X = reg$Y$index, INDEX = reg$Y$id, FUN = min)
	expect_equal(t_min_registered, expected = t_min_observed)
})

test_that("registr function with incompleteness != NULL and lambda_inc = 0: warping functions do not exceed the overall time domain",{
	Y    = registr::growth_incomplete
	reg1 = registr(Y = Y, family = "gaussian", incompleteness = "leading",  lambda_inc = 0)
	reg2 = registr(Y = Y, family = "gaussian", incompleteness = "trailing", lambda_inc = 0)
	reg3 = registr(Y = Y, family = "gaussian", incompleteness = "full",     lambda_inc = 0)
	
	t_min = min(Y$index)
	t_max = max(Y$index)
	# round the registered times to prevent numerical / machine number errors
	expect_gte(round(min(reg1$Y$index), 5), expected = t_min)
	expect_gte(round(min(reg2$Y$index), 5), expected = t_min)
	expect_gte(round(min(reg3$Y$index), 5), expected = t_min)
	expect_lte(round(max(reg1$Y$index), 5), expected = t_max)
	expect_lte(round(max(reg2$Y$index), 5), expected = t_max)
	expect_lte(round(max(reg3$Y$index), 5), expected = t_max)
})

test_that("registr function with incompleteness = 'leading': higher lambda_inc values lead to starting points closer to the diagonal",{
	Y = registr::growth_incomplete
	reg1 = registr(Y = Y, family = "gaussian",
								 incompleteness = "leading", lambda_inc = 1)
	reg2 = registr(Y = Y, family = "gaussian",
								 incompleteness = "leading", lambda_inc = 10)
	
	# compare the MSE values of the warping function starting points to the diagonal
	t_min_observed_1   = tapply(X = reg1$Y$tstar, INDEX = reg1$Y$id, FUN = min)
	t_min_registered_1 = tapply(X = reg1$Y$index, INDEX = reg1$Y$id, FUN = min)
	t_min_observed_2   = tapply(X = reg2$Y$tstar, INDEX = reg2$Y$id, FUN = min)
	t_min_registered_2 = tapply(X = reg2$Y$index, INDEX = reg2$Y$id, FUN = min)
	MSE_min1 = sum((t_min_registered_1 - t_min_observed_1)^2)
	MSE_min2 = sum((t_min_registered_2 - t_min_observed_2)^2)
	expect_gt(MSE_min1, expected = MSE_min2)
})

test_that("registr function with incompleteness = 'trailing': higher lambda_inc values lead to endpoints closer to the diagonal",{
	Y = registr::growth_incomplete
	reg1 = registr(Y = Y, family = "gaussian",
								 incompleteness = "trailing", lambda_inc = 1)
	reg2 = registr(Y = Y, family = "gaussian",
								 incompleteness = "trailing", lambda_inc = 10)
	
	# compare the MSE values of the warping function endpoints to the diagonal
	t_max_observed_1   = tapply(X = reg1$Y$tstar, INDEX = reg1$Y$id, FUN = max)
	t_max_registered_1 = tapply(X = reg1$Y$index, INDEX = reg1$Y$id, FUN = max)
	t_max_observed_2   = tapply(X = reg2$Y$tstar, INDEX = reg2$Y$id, FUN = max)
	t_max_registered_2 = tapply(X = reg2$Y$index, INDEX = reg2$Y$id, FUN = max)
	MSE_max1 = sum((t_max_registered_1 - t_max_observed_1)^2)
	MSE_max2 = sum((t_max_registered_2 - t_max_observed_2)^2)
	expect_gt(MSE_max1, expected = MSE_max2)
})

test_that("registr function throws an error when family = 'gamma' is applied to non-strictly positive data",{
	Y = registr::growth_incomplete
	
	expect_error(registr(Y = Y, family = "gamma"),
							 "family = 'gamma' can only be applied to strictly positive data.")
})

test_that("registr function throws an error when family = 'poisson' is applied to negative data",{
	Y = registr::growth_incomplete
	
	expect_error(registr(Y = Y, family = "poisson"),
							 "family = 'poisson' can only be applied to nonnegative data.")
})

test_that("registr function only accepts Y_template if it has the correct format.",{
	Y = registr::growth_incomplete
	
	template_ids1 = "girl18"
	Y_template1   = Y[Y$id %in% template_ids1,]
	template_ids2 = "girl12"
	Y_template2   = Y[Y$id %in% template_ids2,]
	template_ids3 = c("girl12","boy04","boy29")
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
