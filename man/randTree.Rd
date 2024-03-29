% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/a.R
\name{randTree}
\alias{randTree}
\title{Generate a random tree}
\usage{
randTree(n, wndmtrx = FALSE, parallel = FALSE)
}
\arguments{
\item{n}{number of tips, integer, must be 3 or greater}

\item{wndmtrx}{T/F add node matrix? Default FALSE.}

\item{parallel}{T/F run in parallel? Default FALSE.}
}
\description{
Returns a random \code{TreeMan} tree with \code{n}
tips.
}
\details{
Equivalent to \code{ape}'s \code{rtree()} but returns a
\code{TreeMan} tree. Tree is always rooted and bifurcating.
}
\examples{

tree <- randTree(5)
}
\seealso{
\code{\link{TreeMan-class}}, \code{\link{blncdTree}},
\code{\link{unblncdTree}}
}
