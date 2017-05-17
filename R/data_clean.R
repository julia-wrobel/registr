#' convert data to a refund object
#' 
#' Function used for data cleaning
#'
#' @param data input data. Can be long format dataframe or wide format matrix. If you provide a matrix
#' @param family \code{gaussian} or \code{binomial}.
#' 
#' @importFrom dplyr mutate group_by filter select ungroup
#' @importFrom tidyr spread
#' @importFrom magrittr %>%
#' @importFrom refund.shiny as_refundObj
#' 
#' @export

data_clean <- function(data, family = "gaussian"){
  # eventually want to create this function with no extra packages
  # that way we do not have to include them as recquired
  
  # if matrix in wide format convert to refund object
  if (is.matrix(data)) {
    index = seq(0, 1, length.out = dim(data)[2])
    data = as_refundObj(data, index = index)
  }
  
  # insert check for matrix names
  # create a stop here if names index, id, value do not exist
  # create a stop here if value is not binary when family is binomial
  
  # get row numbers for each subject
  data_rows = data %>% mutate(row = row_number()) %>% group_by(id) %>% 
    filter(row_number()==1 | row_number()==n()) %>% mutate(index = ifelse(round(index) == 0, "first", "last")) %>%
    select(-value )  %>% spread(index, row) %>% ungroup() %>% mutate(subject = row_number())
  
  colnames(data_rows) <- c("id", "first_row", "last_row", "subject")
  I = dim(data_rows)[1]
  
  return(list(Y = data, I = I, Y_rows = data_rows))
  
}