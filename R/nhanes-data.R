#' NHANES activity data
#'
#' Subset of 24 hours of activity data for 50 subjects from 2003-2004
#' National Health and Nutrition Examination Survey (NHANES). 
#' Each subject is observed over 24 hours on a Sunday and wore the 
#' activity collection device for a minimum of 10 hours. Activity is measured each minute over 24 hours.
#'
#' @docType data
#'
#' @usage data(nhanes)
#'
#' @format A dataframe made up of \describe{
#'  \item{id}{A unique subject identifier;}
#'  \item{age}{Age of survey participant;}
#'  \item{gender}{Gender of survey participant;}
#'  \item{index}{Observed time of activity measurement. Integers from 1 to 1440, indicating minutes
#'  from midnight to midnight;}
#'  \item{value}{Binary value of zero or one indicating inactivity or activity;}
#'  \item{raw_activity}{Raw activity count.}
#' }
#'
#' @keywords datasets
"nhanes"