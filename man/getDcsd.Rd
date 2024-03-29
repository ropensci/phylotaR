% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/a.R
\name{getDcsd}
\alias{getDcsd}
\title{Get extinct tips from a tree}
\usage{
getDcsd(tree, tol = 1e-08)
}
\arguments{
\item{tree}{\code{TreeMan} object}

\item{tol}{zero tolerance}
}
\description{
Return all extinct tip \code{ID}s.
}
\details{
Returns a vector.
}
\examples{

tree <- randTree(10)
(getDcsd(tree))
}
\seealso{
\code{\link{getLvng}}, \code{\link{isUltrmtrc}},
\url{https://github.com/DomBennett/treeman/wiki/get-methods}
}
