#' Convert matrices to dataframes for use in functional data analyses
#'
#' @param obj Matrix object to be converted; rows contain functional observations on subjects.
#' @param index Time grid on which functional data are observed; defaults to \code{NULL},
#' which assumes an equally-spaced grid on [0,1].
#' @param ... additional arguments to be passed to methods (not used).
#'
#' @author Jeff Goldsmith \email{jeff.goldsmith@@columbia.edu}
#'
#' @return An object of classes \code{data.frame} and \code{refund.object}, the latter of
#' which is so far not used. Columns are \code{id} (taken from the rownames of \code{obj},
#' if they exist), \code{index} (with behavior described above), and \code{value} (taken
#' from entries in \code{obj}).
#'
#' @export
#' @importFrom stats complete.cases setNames
#' @examples
#'
#' \dontrun{
#' library(ggplot2)
#' library(refund)
#'
#' cca_df = as_refundObj(DTI$cca)
#' ggplot(cca_df, aes(x = index, y = value, group = id)) + geom_line()
#' }
#'
as_refundObj.matrix = function(obj, index = NULL, ...) {

	if (is.null(index)) {
		index = seq(0, 1, length = dim(obj)[2])
	}

	if (!is.null(rownames(obj))) {
		id = rownames(obj)
	} else {
		id = 1:dim(obj)[1]
	}

	df = data.frame(
		id = rep(id, each = dim(obj)[2]),
		index = rep(index, dim(obj)[1]),
		value = as.vector(t(obj)),
		stringsAsFactors = FALSE
	)

	df = df[complete.cases(df),]

	class(df) = c("data.frame", "refund_object")

	df
}
