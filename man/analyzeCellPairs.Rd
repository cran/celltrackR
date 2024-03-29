% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/angle-functions.R
\name{analyzeCellPairs}
\alias{analyzeCellPairs}
\title{Find Distances and Angles for all Pairs of Tracks}
\usage{
analyzeCellPairs(X, searchRadius = Inf, quietly = FALSE, ...)
}
\arguments{
\item{X}{a tracks object}

\item{searchRadius}{if specified, only return analysis for pairs of cells that
are within distance searchRadius from each other at least at one point in time.}

\item{quietly}{(default FALSE) if TRUE, suppress warnings}

\item{...}{further arguments passed on to \code{angleCells}}
}
\value{
A dataframe with four columns: two for the indices of cellpairs,
one for the distance between them, and one for their angle. Note that the 
distance will be NA for pairs of tracks that do not share time points, but
their angle will still be computed.
}
\description{
Find all pairs of cells and return the shortest distance between them at any
point in time (if they share any time points), as well as the angle between 
their overall displacement vectors.
}
\details{
Analyzing track angles at different distances can be useful to detect
directional bias or local crowding effects; see (Beltman et al, 2009).

Internally, the function uses \code{\link{cellPairs}}, \code{\link{angleCells}},
and \code{\link{distanceCells}}.
}
\examples{
## Plot distance versus angle for all cell pairs. Sample T-cell data here for speed.
pairs <- analyzeCellPairs( sample( TCells, 100 ) )
scatter.smooth( pairs$dist, pairs$angle )
}
\references{
Joost B. Beltman, Athanasius F.M. Maree and Rob. J. de Boer (2009),
Analysing immune cell migration. \emph{Nature Reviews Immunology} \bold{9},
789--798. doi:10.1038/nri2638
}
\seealso{
\code{\link{analyzeStepPairs}} to do something similar for single steps
rather than entire tracks.
}
