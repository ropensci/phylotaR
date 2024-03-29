% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/a.R
\name{getPrnt}
\alias{getPrnt}
\title{Get parent}
\usage{
getPrnt(tree, ids)
}
\arguments{
\item{tree}{\code{TreeMan} object}

\item{ids}{vector of node ids}
}
\description{
Return parental (most recent common ancestor) node id for \code{ids}.
}
\details{
Returns a character.
}
\examples{

data(mammals)
# choosing ids from the two main branches of apes allows to find the parent for all apes
ape_id <- getPrnt(mammals, ids = c("Homo_sapiens", "Hylobates_concolor"))
}
\seealso{
\code{\link{getSubtree}},
\url{https://github.com/DomBennett/treeman/wiki/get-methods}
}
