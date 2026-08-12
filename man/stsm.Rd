\name{stsm}
\alias{stsm}
\title{
Stress Majorization Layout
}
\description{
A function to compute the graph coordinated system with a stress majorization layout algorithm.
}
\usage{
stsm(net, seed = seed, maxiter = 40, drp, jitter, method, ...)
}
\arguments{
  \item{net}{Array for network relations.}
  \item{seed}{(mandatory) the seed of the initial layout (see \emph{details}).}
  \item{maxiter}{(optional) Maximum number of iterations.}
  \item{drp}{(optional, for \code{valued} networks) Drop values less than specified.}
  \item{jitter}{(optional) The amount of jitter to add on coordinates.}
  \item{method}{(optional, default \code{binary}) Initial distance method for graph vertices.}
  \item{\dots}{Additional argument items.}
}
\details{
Like the function \code{\link{frcd}}, this routine serves as an internal tool for graph visualization. 
It is also designed to establish the coordinate system using the \code{coord} option within the \code{\link{multigraph}} and \code{\link{bmgraph}} functions. 

In this scenario, the graph's coordinate system commences with nodes randomly positioned, and if \code{NULL} is entered in the \code{seed} argument.
An initial seed will be generated based on the computer clock watch where the number of iterations in \code{maxiter} is \eqn{40}.
}
\value{
A data frame with a coordinated system with two columns representing the abscissa and the ordinate in a two-dimensional rectangular Cartesian coordinate system.
}
\references{
Gansner, E.R., Koren, Y., & North, S. \emph{Graph drawing by stress majorization}. In Graph Drawing 2004: 12th International Symposium, New York, NY, USA, revised selected papers. Springer. pp. 239-250. 2005.
}
%\author{} ~
%\note{} ~
%\section{Warning }{....} ~
\seealso{
\code{\link{multigraph}}, \code{\link{bmgraph}}, \code{\link{frcd}}, \code{\link{conc}}
}
\examples{
# Two binary relations among three elements
arr <- round( replace( array(runif(18), c(3,3,2)), array(runif(18),
       c(3,3,2))>.5, 3 ) )

# Coordinates for stress majorization layout with random start
coord <- stsm(arr, seed = NULL)

# Plot multigraph with stress majorization coordinates
multigraph(arr, coord = coord)
}
\keyword{graphics}
\keyword{manip}
