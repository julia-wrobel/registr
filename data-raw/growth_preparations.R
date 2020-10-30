
# code for preparing the Berkeley growth study dataset with simulated incompleteness

library(dplyr)
library(mgcv)  # for smoothing the raw curves



# basic data preparation --------------------------------------------------
dat_raw <- fda::growth

# transform to long dataset (and ensure appropriate sorting)
dat <- data.frame(id    = factor(rep(c(colnames(dat_raw$hgtm), colnames(dat_raw$hgtf)),
																		 each = length(dat_raw$age))),
									index = rep(dat_raw$age, times = ncol(dat_raw$hgtm) + ncol(dat_raw$hgtf)),
									value = c(as.vector(dat_raw$hgtm), as.vector(dat_raw$hgtf))) %>%
	arrange(id, index)

# slightly smooth the raw curves
smooth_list <- lapply(unique(dat$id), function(curve_id) {
	d <- dat %>% filter(id == curve_id)
	m <- mgcv::gam(value ~ s(index, bs = "cr", k = 15), data = d)
	d$value <- unname(mgcv::predict.gam(m))
	return(d)
})
dat_smooth <- dplyr::bind_rows(smooth_list)



# take the first derivative -----------------------------------------------
deriv_list <- lapply(levels(dat$id), function(curve_id) {
	d <- dat_smooth %>% filter(id == curve_id)
	data.frame(id               = curve_id,
						 index            = d$index[2:nrow(d)],
						 value            = diff(d$value) / diff(d$index),
						 stringsAsFactors = FALSE)
})
dat_deriv <- dplyr::bind_rows(deriv_list) %>%
	mutate(id = factor(id))



# simulate incompleteness -------------------------------------------------
index_raw <- unique(dat_deriv$index)
ids       <- unique(dat_deriv$id)

# for each curve, draw a random end time point in the second half of the domain
index_end <- index_raw[index_raw > (max(index_raw) / 2)]
set.seed(2020)
endpoints <- sample(index_end,
										size    = length(ids),
										replace = TRUE)
growth_incomplete <- dat_deriv
for (i in 1:length(ids)) {
	cut_rows <- which(growth_incomplete$id == ids[i] & growth_incomplete$index > endpoints[i])
	if (length(cut_rows) > 0)
		growth_incomplete <- growth_incomplete %>% slice(-cut_rows)
}



# save the prepared dataset -----------------------------------------------
usethis::use_data(growth_incomplete, overwrite = TRUE)


