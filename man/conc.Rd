\name{conc}
\alias{conc}
\title{
Concentric Layout
}
\description{
A function to compute the graph coordinated system with a concentric layout.
}
\usage{
conc(net, nr, irot, inv, iflip, mirror = c("N", "X", "Y", "D", "L"), ...)
}
\arguments{
  \item{net}{Array representing the network relations.}
  \item{nr}{Scalar with number of radii, or a vector with clustering of actors.}
  \item{irot}{Scalar or vector with \dQuote{internal rotation} for each circle from closer to center point to further away.}
  \item{inv}{(optional and logical) Invert circles ordering?}
  \item{iflip}{(optional and logical) Alternate circles flipping?}
  \item{mirror}{
     Mirror transformation:
      \itemize{ 
       \item \code{N} Identity (default)
       \item \code{X} Reflection through the vertical center line
       \item \code{Y} Reflection through the horizontal center line
       \item \code{D} Reflection across diagonal \eqn{Y=X}
       \item \code{L} Reflection across diagonal \eqn{Y=-X}
     }}
  \item{\dots}{Additional argument items.}
}
\details{
An \dQuote{Euclidean plane} computes the coordinated system with a concentric layout with at least two radii (unless \eqn{n = 1}). 
In case that the number of radii is not specified in \code{nr}, approx. half of the vertices are located at one radius and half in another one.

The clustering of the actors may be used to establish the location of the vertices in different radii as a numerical, character, or a factor vector.
}
\value{
A data frame with a coordinated system with two columns representing the abscissa and the ordinate in a two-dimensional rectangular Cartesian coordinate system.
}
%\references{} ~
%\author{} ~
%\note{} ~
%\section{Warning }{....} ~
\seealso{
\code{\link{multigraph}}, \code{\link{bmgraph}}, \code{\link{frcd}}, \code{\link{stsm}}
}
\examples{
# Two binary relations among three elements
arr <- round( replace( array(runif(18), c(3,3,2)), array(runif(18),
       c(3,3,2))>.5, 3 ) )

# Coordinates for the concentric layout with two radii
coord <- conc(arr, nr = 2)

# Plot multigraph with customized coordinates
multigraph(arr, coord = coord)

}
\keyword{graphics}
\keyword{manip}
