% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gfpca_covHall.R
\name{crossprods_irregular}
\alias{crossprods_irregular}
\title{Crossproduct computation for highly irregular grids}
\usage{
crossprods_irregular(Y)
}
\arguments{
\item{Y}{Dataframe with the centered observations.
Should have values id, centered, index.}
}
\description{
Compute the crossproduct in a fast way for highly irregular grids
(index values are mostly unique).
Only used internally in \code{cov_hall()}.
}
