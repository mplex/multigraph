\name{conc}
\alias{conc}
\title{
Concentric layout
}
\description{
A function to compute the graph coordinated system with a concentric layout
}
\usage{
conc(net, nr, irot, inv, flip, mirror=c("N","X","Y","D","L"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{net}{
	an array representing the network relations
}
  \item{nr}{
	a scalar with the number of radii, or a vector with the clustering of the actors.
}
  \item{irot}{
	a scalar or vector with the ``internal rotation'' for each circle from closer to the center point to further away
}
  \item{inv}{
	(optional and logical) should the circles be with an inverted ordering?
}
  \item{flip}{
	(optional and logical) should the alternating circles be flipped ?
}
  \item{mirror}{
	mirror transformation
       \code{N} identity (default)

       \code{X} reflection through the vertical center line

       \code{Y} reflection through the horizontal center line

       \code{D} reflection across diagonal Y=X

       \code{L} reflection across diagonal Y=-X
}
  \item{\dots}{
	Additional argument items 
}
}
\details{
In a Euclidean plane computes the coordinated system with a concentric layout with at least two radii (unless \eqn{n = 1}). In case \code{nr} is not specified, approx. half of the vertices are located at one radius and half in another one.

The clustering of the actors may be used to establish the location of the vertices in different radii as a numerical, character, or a factor vector.
}
\value{
A data frame with the coordinated system with two columns representing the abscissa and the ordinate in a two-dimensional rectangular Cartesian coordinate system.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Antonio Rivero Ostoic
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{multigraph}}, \code{\link{bmgraph}}, \code{\link{frcd}}, \code{\link{stsm}}
}
\examples{
## Create the data: two binary relations among three elements
arr <- round( replace( array(runif(18), c(3,3,2)), array(runif(18),
       c(3,3,2))>.5, 3 ) )

## Coordinates for the concentric layout with two radii
coord <- conc(arr, nr = 2)

## Plot multigraph with costumized coordenates
multigraph(arr, coord = coord)

}
\keyword{ graphics }
\keyword{ manip }
