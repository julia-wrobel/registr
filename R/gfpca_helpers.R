#' Coarsen an index vector to a given resolution
#' 
#' Reduce the resolution of a numeric vector by specifying the number of
#' \code{significant_digits} to which the numbers should be rounded. \cr \cr
#' Internal function used to coarsen the index vector before estimating the
#' two-step GFPCA with \code{\link{gfpca_twoStep}}.
#' 
#' @param index Numeric vector of index values.
#' @param significant_digits Positive integer value.
#' 
#' @return Numeric vector of rounded index values.
#' 
#' @author Alexander Bauer \email{alexander.bauer@@stat.uni-muenchen.de}
#' 
#' @examples
#' index_vector = c(0.7892, 0.2984, 0.328)
#' registr:::coarsen_index(index_vector, 1)
#' registr:::coarsen_index(index_vector, 3)
#' 
#' index_vector2 = c(2803, -7639, 13)
#' registr:::coarsen_index(index_vector2, 1)
#' registr:::coarsen_index(index_vector2, 3)
#' 
coarsen_index = function(index, significant_digits) {
	
	if (significant_digits < 1) {
		stop("'significant_digits' must be a positive integer.")
	}
	
	# remove potential signs
	signs     = ifelse(index > 0, +1, -1)
	index_abs = abs(index)
	
	digits_preDecimal = max(nchar(round(index_abs, 0)))
	# if all values are between 0 and 1, set digits_preDecimal to zero
	if ((digits_preDecimal == 1) & all(substr(index_abs, 1, 1) == 0))
		digits_preDecimal = 0
	
	index_abs_scaled         = index_abs / 10^(digits_preDecimal - 1)
	index_abs_scaled_rounded = round(index_abs_scaled, significant_digits - 1)
	index_abs_rounded        = index_abs_scaled_rounded * 10^(digits_preDecimal - 1)
	
	# add the signs again
	index_rounded = signs * index_abs_rounded
	
	return(index_rounded)
}


#' Estimate the derivative of the logit function
#' 
#' Compute the derivative of the logit function for a given point \code{x}.
#' 
#' @param x Value at which the derivative is computed
#' 
#' @author Alexander Bauer \email{alexander.bauer@@stat.uni-muenchen.de}
#' 
deriv.inv.logit = function(x) {
  exp(x) / ( (1 + exp(x))^2 )
}
