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



#' Berkeley Growth Study data with simulated incompleteness
#' 
#' This dataset from the Berkeley Growth Study comprises the height
#' development of 39 boys and 54 girls between ages 1 and 18.
#' It is based on the dataset \code{fda::growth} and focuses not on the observed
#' heights, but on the first derivatives of the curves. Before taking the
#' first derivative, the curves were slightly smoothed. \cr \cr
#' To showcase the functionality of the \code{registr} package regarding the
#' analysis of incomplete curves, the growth curves were artificially made
#' incomplete. For each child, leading incompleteness was simulated by drawing
#' a random initial age in the first quarter of the domain.
#' Also, trailing incompleteness was simulated by drawing
#' a random cut-off age in the second half of the domain.
#' 
#' @docType data
#' 
#' @usage data(growth_incomplete)
#' 
#' @format A dataframe made up of \describe{
#'  \item{id}{A unique subject identifier;}
#'  \item{index}{Observed age of the child's height;}
#'  \item{value}{First derivative of the height development in the given age.}
#' }
#' 
#' @references Ramsay, J. O., and Silverman, B. W. (2006),
#' \emph{Functional Data Analysis, 2nd ed.}, Springer, New York.
#' 
#' Tuddenham, R. D., and Snyder, M. M. (1954).
#' Physical growth of California boys and girls from birth to age 18.
#' \emph{University of California Publications in Child Development}, 1, 183-364.
#' 
#' 
#' @keywords datasets
"growth_incomplete"