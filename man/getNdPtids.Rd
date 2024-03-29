% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/a.R
\name{getNdPtids}
\alias{getNdPtids}
\title{Get post-nodes to tips}
\usage{
getNdPtids(tree, id)
}
\arguments{
\item{tree}{\code{TreeMan} object}

\item{id}{node id}
}
\description{
Return node ids for connecting \code{id} to kids.
}
\details{
Returns a vector.
}
\examples{

tree <- randTree(10)
# get all nodes from root to tip
getNdPtids(tree, id = "n1")
}
\seealso{
\code{\link{getNdsPtids}},
\code{\link{getNdPrids}},
\code{\link{getNdsPrids}},
\url{https://github.com/DomBennett/treeman/wiki/get-methods}
}
