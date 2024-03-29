% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/a.R
\name{getLvng}
\alias{getLvng}
\title{Get extant tips from a tree}
\usage{
getLvng(tree, tol = 1e-08)
}
\arguments{
\item{tree}{\code{TreeMan} object}

\item{tol}{zero tolerance}
}
\description{
Return all extant tip \code{ID}s.
}
\details{
Returns a vector.
}
\examples{

tree <- randTree(10)
(getLvng(tree))
}
\seealso{
\code{\link{getDcsd}}, \code{\link{isUltrmtrc}},
\url{https://github.com/DomBennett/treeman/wiki/get-methods}
}
