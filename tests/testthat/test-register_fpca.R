context("register_fpca")

# for all tests implement for both gaussian and binary data
test_that("registering binary data works",{
	Y = simulate_unregistered_curves()
	# these are just going to be overall tests that stuff works
	# 
	#registration_object = register_fpca(Y, Kt = 8, Kh = 4, family = "binomial", npc = 2)
	
	
})


test_that("loss function is generally decreasing", {
	## use diff function and losses that are returned
})