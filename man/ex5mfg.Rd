% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DTMC.R
\name{ex5mfg}
\alias{ex5mfg}
\title{Manufacturing Systems}
\usage{
ex5mfg(A, B, a1, a2)
}
\arguments{
\item{A}{size of the bin for machine 1;}

\item{B}{size of the bin for machine 2;}

\item{a1}{prob(non-defective) for machine 1;}

\item{a2}{prob(non-defective) for machine 2.}
}
\value{
P = transitrion probability matrix for the manufacturing system. (See Example 5.7.)
}
\description{
Manufacturing Systems
}
\examples{
ex5mfg(5,3,0.5,0.3)
}
