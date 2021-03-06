% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/angle-functions.R
\name{cellPairs}
\alias{cellPairs}
\title{Find Pairs of Tracks}
\usage{
cellPairs(X)
}
\arguments{
\item{X}{a tracks object}
}
\value{
A dataframe with two columns: one for each of the track ids in the pair.
Each row represents a pair.
}
\description{
Get all unique combinations of two track ids.
}
\examples{
## Find all pairs of cells in the T cell data
pairs <- cellPairs( TCells )
}
