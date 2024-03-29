% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/a.R
\name{setNdSpn}
\alias{setNdSpn}
\title{Set the branch length of a specific node}
\usage{
setNdSpn(tree, id, val)
}
\arguments{
\item{tree}{\code{TreeMan} object}

\item{id}{id of node whose preceding edge is to be changed}

\item{val}{new span}
}
\description{
Return a tree with the span of a node altered.
}
\details{
Takes a tree, a node ID and a new value for the node's preceding branch length (span).
}
\examples{

tree <- randTree(10)
tree <- setNdSpn(tree, id = "t1", val = 100)
tree <- updateSlts(tree)
summary(tree)
}
\seealso{
\code{\link{setNdsSpn}}
\url{https://github.com/DomBennett/treeman/wiki/set-methods}
}
