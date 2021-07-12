#' Convert data to a \code{refund} object
#' 
#' Function used for data cleaning.
#'
#' @param data Dataframe. Should have values id, value, index.
#' 
#' @importFrom dplyr mutate group_by filter select ungroup arrange
#' @importFrom tidyr spread
#' @importFrom magrittr %>%
#' 
#' @return An list containing:
#' \item{Y}{The original data sorted by id and index.}
#' \item{Y_rows}{A dataframe containing the first and last row for each subject.}
#' 
#' @export
#' 
data_clean = function(data){
  
	## NULLify global values called by tidyverse functions
	value = index = NULL
	
	# create a stop here if names index, id, value do not exist
	if(!all(c("id", "index", "value") %in% names(data)) ){
		stop("Input dataset must have variables 'id', 'index', and 'value'.")
	}
	
	# check if data are sorted by id and index; if not do so and return message
	if (!identical(data, data %>% arrange(id, index))) {
		message("Data have been sorted by id and index; all output will be in this format")
		data = data %>% arrange(id, index)
	}
	
  # get row numbers for each subject
  data_rows = data %>% 
  	select(id, index, value) %>%
  	mutate(row = row_number()) %>% 
  	group_by(id) %>% 
    filter(row_number() == 1 | row_number() == n()) %>% 
  	mutate(index = c("first", "last")) %>%
    select(-value) %>% 
  	tidyr::spread(index, row) %>% 
  	ungroup() %>% 
  	mutate(subject = row_number())
  
  data = data %>%
  	group_by(id) %>%
  	mutate(index_scaled = (index - min(index))/ (max(index) - min(index)) ) %>%
  	ungroup()
  
  colnames(data_rows) = c("id", "first_row", "last_row", "subject")
  I = dim(data_rows)[1]
  
  return(list(Y = data, I = I, Y_rows = data_rows))
  
}
