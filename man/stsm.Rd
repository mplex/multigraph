\name{stsm}
\alias{stsm}
\title{
Stress majorization layout
}
\description{
A function to compute the graph coordinated system with a stress majorization layout algorithm
}
\usage{
stsm(net, seed = seed, maxiter = 40, drp, jitter, method, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{net}{
	an array representing the network relations
}
  \item{seed}{
        (mandatory) the seed of the initial layout (see details)
}
  \item{maxiter}{
        (optional) the maximum number of iterations 
}
  \item{drp}{
	(optional) for \code{valued} networks, drop values less than specified
}
  \item{jitter}{
	(optional) jitter in the layout
}
  \item{method}{
	(optional) initial distance method (default \code{binary})
}
  \item{\dots}{
	Additional argument items 
}
}
\details{
This function is meant as an internal routine for graph visualization. However, it can be used with the \code{coord} option both in \code{multigraph} and in \code{bmgraph} where \code{NULL} in \code{seed} implies a random seed based on the clock watch of the computer. 
}
\value{
A data frame with a coordinated system with two columns representing the abscissa and the ordinate in a two-dimensional rectangular Cartesian coordinate system.
}
\references{
Gansner, E.R., Koren, Y., & North, S. \emph{Graph drawing by stress majorization}. In Graph Drawing: 12th International Symposium, gd 2004, New York, NY, USA, September 29 - October 2, 2004, revised selected papers. Berlin Heidelberg: Springer. pp. 239-250. 2005.
}
\author{
Antonio Rivero Ostoic
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{multigraph}}, \code{\link{bmgraph}}, \code{\link{frcd}}, \code{\link{conc}}
}
\examples{
## Create the data: two binary relations among three elements
arr <- round( replace( array(runif(18), c(3,3,2)), array(runif(18),
       c(3,3,2))>.5, 3 ) )

## Coordinates for the stress majorization layout with random start
coord <- stsm(arr, seed = NULL)

## Plot multigraph with costumized coordenates
multigraph(arr, coord = coord)

}
\keyword{ graphics }
\keyword{ manip }
