\name{frcd}
\alias{frcd}
\title{
Force directed layout
}
\description{
A function to compute the graph coordinated system with a force directed layout algorithm
}
\usage{
frcd(net, seed = seed, maxiter, drp, scl, mov, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{net}{
	an array representing the network relations
}
  \item{seed}{
        (mandatory) the seed of the initial layout (see \emph{details})
}
  \item{maxiter}{
        (optional) the maximum number of iterations 
}
  \item{\dots}{
	Additional argument items 
}
  \item{scl}{
	(optional and experimental) numerical scalar (\eqn{x} and \eqn{y}) or vector (\eqn{x}, \eqn{y}) of the graph's scale 
}
  \item{mov}{
	(optional and experimental)  numerical scalar (\eqn{x} and \eqn{y}) or vector (\eqn{x}, \eqn{y}) to move the graph
}
  \item{drp}{
	(optional) for \code{valued} networks, drop values less than specified
}
}
\details{
This function is meant as an internal routine for graph visualization with a force-directed layout procedure. 
However, it can be used to set the coordinate system with the \code{coord} option in functions \code{\link{multigraph}} and in \code{\link{bmgraph}}. 
In such case, the coordinate system of the graph starts with a random displacement of nodes where \code{NULL} in the \code{seed} argument implies 
an initial seed based on the computer clock watch, and the number of iterations in \code{maxiter} is \eqn{60+n}.
}
\value{
A data frame with a coordinated system with two columns representing the abscissa and the ordinate in a two-dimensional rectangular Cartesian coordinate system.
}
\references{
Fruchterman, T.M.J., & Reingold, E.M. Graph drawing by force-directed placement. \emph{Software-Practice & Experience}, 21(11), 1129-1164. 1991.
}
\author{
Antonio Rivero Ostoic
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{multigraph}}, \code{\link{bmgraph}}, \code{\link{stsm}}, \code{\link{conc}}
}
\examples{
## Create the data: two binary relations among three elements
arr <- round( replace( array(runif(18), c(3,3,2)), array(runif(18),
       c(3,3,2))>.5, 3 ) )

## Coordinates for the force directed layout with random start
coord <- frcd(arr, seed = NULL)

## Plot multigraph with customized coordinates
multigraph(arr, coord = coord)

}
\keyword{ graphics }
\keyword{ manip }
