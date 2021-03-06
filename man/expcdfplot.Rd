% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distr.R
\name{expcdfplot}
\alias{expcdfplot}
\title{Plot of Exponential Cumulative Distribution Function}
\usage{
expcdfplot(l = 1, x = c(0.4, 0.3, 0.2, 0.1))
}
\arguments{
\item{l}{>= 0}

\item{x}{a non negative row vector}
}
\value{
Plot of the cdf of exp(l) (y vs. x).
}
\description{
Plot of Exponential Cumulative Distribution Function
}
\examples{
expcdfplot(1,c(0.4,0.3,0.2,0.1))
}
