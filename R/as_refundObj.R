#' Convert data to refund objects for use in functional data analyses
#'
#' Very experimental function, primarily used to convert matrices storing functional data
#' to data.frames with specific variable names.
#'
#' @param obj Object to be converted. Currently supports class \code{matrix}, formatted
#' so that rows contain functional observations on subjects.
#' @param ... additional arguments to be passed to methods.
#'
#' @author Jeff Goldsmith \email{jeff.goldsmith@@columbia.edu}
#'
#' @return An object of classes \code{data.frame} and \code{refund.object}, the latter of
#' which is so far not used. Columns are \code{id} (taken from the rownames of \code{obj},
#' if they exist), \code{index} (with behavior described above), and \code{value} (taken
#' from entries in \code{obj}).
#'
#' @export as_refundObj
#'
#' @examples
#'
#' \dontrun{
#' library(ggplot2)
#' library(refund)
#'
#' cca_df = as_refundObj(DTI$cca)
#'
#' ggplot(cca_df, aes(x = index, y = value, group = id)) + geom_line()
#' }
#'
as_refundObj <- function(obj, ...){
	UseMethod("as_refundObj")
}
