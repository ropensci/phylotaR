% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/a.R
\name{getNdPD}
\alias{getNdPD}
\title{Get phylogenetic diversity of node}
\usage{
getNdPD(tree, id)
}
\arguments{
\item{tree}{\code{TreeMan} object}

\item{id}{node id}
}
\description{
Return summed value of all descending spns
}
\details{
Sums the lengths of all descending branches from a node.
}
\examples{

tree <- randTree(10)
getNdPD(tree, id = "n1") # return PD of n1 which in this case is for the whole tree
}
\seealso{
\code{\link{getNdsPD}},
\url{https://github.com/DomBennett/treeman/wiki/get-methods}
}
