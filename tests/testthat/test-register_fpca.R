context("register_fpca")


# for all tests implement for both gaussian and binary data
test_that("code only accepts supported family of distributions", {
	Y = simulate_unregistered_curves()
	registr_object = register_fpca(Y, family = "binomial", iterations = 3)
	
	#expect_error(register_fpca(Y, family = "gaussian", iterations = 3))
	expect_error(register_fpca(Y, family = "poisson"),
							 "Package currently handles only 'binomial' or 'gaussian' families.")
	

})

test_that("registering binary data works",{
	Y = simulate_unregistered_curves()
	registr_object = register_fpca(Y, family = "binomial", iterations = 20)
	# these are just going to be overall tests that stuff works
	# 
	#registration_object = register_fpca(Y, Kt = 8, Kh = 4, family = "binomial", npc = 2)
	
	
})


test_that("loss function is generally decreasing", {
	## use diff function and losses that are returned
})