\name{ccgraph}
\alias{ccgraph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cayley colour graph
}
\description{
A function to create and manipulate bipartite Cayley colour graphs
}
\usage{
ccgraph(x, main = NULL, seed = 0, maxiter = 100, alpha = c(1, 1, 1), scope, loops, 
        collRecip, undRecip, showLbs, cex.main, conc, coord, clu, cex, lwd, pch, 
        lty, bwd, att, bg, mar, pos, asp, ecol, vcol, vcol0, lbs, col, lbat, swp, 
        swp2, scl, mirrorX, mirrorY, mirrorD, mirrorL, mirrorV, mirrorH, rot, hds, 
        vedist, ffamily, fstyle, fsize, fcol, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
an algebraic structure, typically a \code{"Semigroup"} object class
}
  \item{main}{
	(optional) title of the plot
}
  \item{seed}{
	(optional) the random seed number for the vertices' initial coordinates. Ignored except for \code{force}, \code{stress} and \code{rand}
}
  \item{maxiter}{
	(optional) the maximum number of iterations in layout algorithms. Ignored except for \code{force}, \code{stress} and \code{rand}
}
  \item{alpha}{
	vector (vertex, edge, \code{bg}) with the alpha color transparecy
}
  \item{scope}{
   (optional) the scope of the graph (see details)
}
  \item{loops}{
	(optional, logical, and experimental) plot graph loops?
}
  \item{collRecip}{
	(optional and logical) whether or not collapse reciprocated edges (for undirected graphs)
}
  \item{undRecip}{
	(optional and logical) whether or not plot reciprocated edges as undirected graph
}
  \item{showLbs}{
	(optional and logical) whether or not show the vertex labels when dimnames available
}
  \item{cex.main}{
	(optional) the size of the plot's title
}
  \item{conc}{
	(optional and logical) whether the layout is concentric or not
}
  \item{coord}{
	(optional) data frame with the coordinates of the vertices. If coordinates are given then the \code{layout} option is ignored
}
  \item{clu}{
	(optional) the clustering of the vertices (see \emph{details})
}
  \item{cex}{
	(optional) the size of the vertices
}
  \item{lwd}{
	(optional) the width of the edges. Ignored if \code{weighted} is set to \code{TRUE}
}
  \item{pch}{
	(optional) the symbol representing the vertices
}
  \item{lty}{
	(optional) the shape of the edges
}
  \item{bwd}{
	(optional) the width of the bundle edges. Ranges from \code{0} (edges collapsed) to the default \code{1} (depending on the vertices' size). For \code{weighted} a value greater than one is possible
}
  \item{att}{
	(optional) a vector or an array representing the vertex attributes
}
  \item{bg}{
	(optional) the background color of the plot
}
  \item{mar}{
	(optional) the margins of the plot
}
  \item{pos}{
	(optional) the position of the vertices' labels (\code{0} means ``at the center of the vertex'')
}
  \item{asp}{
	(optional) the aspect ratio of the plot
}
  \item{ecol}{
	(optional) the color of the edges
}
  \item{vcol}{
	(optional) the color of the vertices
}
  \item{vcol0}{
	(optional) the color of the vertices' contour (only works for \code{pch 21} through \code{25}
}
  \item{lbs}{
	(optional) the vertex labels
}
  \item{col}{
	(optional) alias for \code{vcol}
}
  \item{lbat}{
	(optional) the labels for the vertex attributes
}
  \item{swp}{
	(optional and logical) whether or not swap the bundle patterns
}
  \item{swp2}{
	(optional and logical) whether or not swap reciprocals
}
  \item{scl}{
	(optional and experimental) numerical scalar (x and y) or vector (x, y) of the graph's scale 
}
  \item{mirrorX}{
	(optional) mirror of the X axis
}
  \item{mirrorY}{
	(optional) mirror of the Y axis
}
  \item{mirrorD}{
	(optional) mirror reflection across diagonal Y=X
}
  \item{mirrorL}{
	(optional) mirror reflection across diagonal Y=-X
}
  \item{mirrorV}{
	same as \code{mirrorX}
}
  \item{mirrorH}{
	same as \code{mirrorY}
}
  \item{rot}{
	(optional) clockwise rotation of the graph in degrees 
}
  \item{hds}{
	(optional and experimental) arcs' head scale
}
  \item{vedist}{
	(optional and experimental) a real number with vertex - edge distance
}
  \item{ffamily}{
	the font family
}
  \item{fstyle}{
	the font style
}
  \item{fsize}{
	the font size
}
  \item{fcol}{
	the font color
}
  \item{\dots}{
	Additional argument items (see e.g. \code{\link[graphics:par]{par}})
}
}
%
\details{
The Cayley colour graph is a graphical representation of the relationships among relations in the relational structure of a given multiplex network. Both nodes and directed edges represent string relations and each shape (and color) corresponds to a specific generator relation of the semigroup structure.

}
\value{
A plot of the semigroup or group structure.
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

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[multiplex:semigroup]{semigroup}}, \code{\link{multigraph}}, \code{\link{frcd}}, \code{\link{conc}}
}
\examples{
## Create an abstract semigroup
arr <- round( replace( array(runif(18), c(3,3,2)), array(runif(18),
       c(3,3,2))>.5, 1 ) )
S <- semigroup(arr)

## plot its Cayley graph
ccgraph(S)
}
\keyword{ graphics }
\keyword{ math }
