% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/a.R
\name{blncdTree}
\alias{blncdTree}
\title{Generate a balanced tree}
\usage{
blncdTree(n, wndmtrx = FALSE, parallel = FALSE)
}
\arguments{
\item{n}{number of tips, integer, must be 3 or greater}

\item{wndmtrx}{T/F add node matrix? Default FALSE.}

\item{parallel}{T/F run in parallel? Default FALSE.}
}
\description{
Returns a balanced \code{TreeMan} tree with \code{n}
tips.
}
\details{
Equivalent to \code{ape}'s \code{stree(type='balanced')} but returns a
\code{TreeMan} tree. Tree is always rooted and bifurcating.
}
\examples{

tree <- blncdTree(5)
}
\seealso{
\code{\link{TreeMan-class}}, \code{\link{randTree}},
\code{\link{unblncdTree}}
}
