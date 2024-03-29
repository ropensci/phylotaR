% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/a.R
\name{getAge}
\alias{getAge}
\title{Get age of tree}
\usage{
getAge(tree, parallel = FALSE)
}
\arguments{
\item{tree}{\code{TreeMan} object}

\item{parallel}{logical, make parallel?}
}
\description{
Returns age, numeric, of tree
}
\details{
Calculates the age of a tree, determined as the maximum tip to root
distance.
}
\examples{

tree <- randTree(10)
(getAge(tree))
}
\seealso{
\code{\link{updateSlts}},
\url{https://github.com/DomBennett/treeman/wiki/get-methods}
}
