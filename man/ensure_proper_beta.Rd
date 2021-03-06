% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/constrOptim_helpers.R
\name{ensure_proper_beta}
\alias{ensure_proper_beta}
\title{Correct slightly improper parameter vectors}
\usage{
ensure_proper_beta(beta, t_min, t_max)
}
\arguments{
\item{beta}{Parameter vector.}

\item{t_min, t_max}{Minimum and maximum of the underlying time domain in the
registration step.}
}
\value{
A slightly changed parameter vector that ensures a proper solution
in the optimization of the registration step.
}
\description{
Internal function. In the joint iterations between registration and GFPCA,
the optimization with \code{constrOptim()} in the registration step sometimes
leads to slightly improper solutions, which cause the optimization to
throw an error in the following optimization step. This function corrects
the parameter vector if one of the following slight inconsistencies occurs
that can mess with the optimization of \code{constrOptim()}: \cr
- two neighboring values of the parameter vector are too similar \cr
- the initial values of the parameter vector are smaller than \code{t_min},
the minimum of the underlying time domain \cr
- the last values of the parameter vector are greater than \code{t_max},
the maximum of the underlying time domain \cr
- one parameter value is slightly greater than its following value, i.e.
the parameter vector is not monotone.
}
\examples{
beta_improper = c(0.24, 1.000047, 1.000002)
registr:::ensure_proper_beta(beta_improper, t_min = 0, t_max = 1)
}
\author{
Alexander Bauer \email{alexander.bauer@stat.uni-muenchen.de}
}
