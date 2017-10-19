\name{multigraph}
\alias{multigraph}
\title{
Multigraphs and weighted multigraphs
}
\description{
A function to create and manipulate multigraphs and weighted multigraphs with different layout options
}
\usage{
multigraph(net, layout = c("circ", "force", "stress", "conc", "rand"), outline, 
    directed = TRUE, main = NULL, seed = NULL, maxiter = 100, alpha = c(1, 1, 1), 
    collRecip, showLbs, showAtts, weighted, cex.main, coord, clu, cex, lwd, pch, 
    lty, bwd, tcol, tcex, att, bg, mar, pos, asp, ecol, vcol, vcol0, hds, vedist, 
    rot, mirrorX, mirrorY, col, lbat, drp, swp, loops, swp2, signed, scl, add, 
    mirrorD, mirrorL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{net}{
	an array; usually with three dimensions of stacked matrices where the multiple relations are placed.
}
  \item{layout}{
    the visualization layout: 

       \code{circ} circular 

       \code{force} force-directed   

       \code{stress} stress-majorization  

       \code{conc} concentric 

       \code{rand} random 
}
  \item{outline}{
   (optional) the outline of the graph (see details)
}
  \item{directed}{
	(logical) whether or not the graph is directed or unidrected
}
  \item{main}{
	(optional) title of the plot
}
  \item{seed}{
	(optional) the random seed number for the vertices' initial coordinates. Ignored for \code{circ} and \code{conc}
}
  \item{maxiter}{
	(optional) the maximum number of iterations in layout algorithms. Ignored for \code{circ} and \code{conc}
}
  \item{alpha}{
	vector (vertex, edge, \code{bg}) with the alpha color transparecy
}
  \item{collRecip}{
	(optional and logical) whether or not collapse reciprocated edges in the unidrected graph
}
  \item{showLbs}{
	(optional and logical) whether or not show the vertex labels
}
  \item{showAtts}{
	(optional and logical) whether or not show the vertex attribute labels
}
  \item{weighted}{
	(optional and logical) whether or not the graph is weighted or dichotomous
}
  \item{cex.main}{
	(optional) the size of the plot's title
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
  \item{tcol}{
	(optional) the color of the vertices' labels
}
  \item{tcex}{
	size of the vertex labels. If \code{NULL}, it depends on the value of \code{cex}
}
  \item{att}{
	(optional) a vector or an array representing the vertices' attributes
}
  \item{bg}{
	(optional) the background color of the plot
}
  \item{mar}{
	(optional) the margins of the plot
}
  \item{pos}{
	(optional) the position of the vertices' labels (\code{0} means ``in the middle of the vertex'')
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
  \item{hds}{
	(optional and experimental) arcs' head scale
}
  \item{vedist}{
	(optional and experimental) a real number with vertex - edge distance
}
  \item{rot}{
	(optional) clockwise rotation of the graph in degrees 
}
  \item{mirrorX}{
	(optional) mirror of the X axis
}
  \item{mirrorY}{
	(optional) mirror of the Y axis
}
  \item{col}{
	(optional) alias for \code{vcol}
}
  \item{lbat}{
	(optional) the labels for the vertices' attributes
}
  \item{drp}{
	(optional) for \code{weighted} networks, drop values less than specified
}
  \item{\dots}{
	Additional argument items (see e.g. \code{\link[graphics:par]{par}})
}
  \item{swp}{
	(optional and logical) whether or not swap the bundle patterns
}
  \item{loops}{
	(optional, logical, and experimental) plot graph loops?
}
  \item{swp2}{
	(optional and logical) whether or not swap reciprocals
}
  \item{signed}{
	(optional and logical) whether or not the graph is a signed structure
}
  \item{scl}{
	(optional and experimental) numerical scalar (x and y) or vector (x, y) of the graph's scale 
}
  \item{add}{
	(optional) nodes to add to the graph
}
  \item{mirrorD}{
	(optional) mirror reflection across diagonal Y=X
}
  \item{mirrorL}{
	(optional) mirror reflection across diagonal Y=-X
}
}
\details{
Multigraph are graph having parallel edges depicting different types of relations in a network. By default a circular layout is applied where each type of tie has a distinctive shape and gray color scale. For a better visualization, undirected multigraphs automatically collapse the reciprocal ties, and there as an argument to prevent this to happen. It is possible to combine the symbols and color of vertices by assigning a class to each network member in the clustering option. Vertices can also have different sizes by specifying the argument with a vector with a length size similar to the network order.

Since this function can handle a large number of arguments, these can be stored as a list object that is passed through the \code{outline} option. In this case a vector made of lists and scalars or combinations of these are accepted.
}
\value{
A plot of the network as a multigraph or a weighted multigraph.
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
\code{\link{bmgraph}}, \code{\link{frcd}}, \code{\link{stsm}}, \code{\link{conc}}
}
\examples{
## Create the data: two binary relations among three elements
arr <- round( replace( array(runif(18), c(3,3,2)), array(runif(18),
       c(3,3,2))>.5, 3 ) )

## Plot the multigraph of this network
multigraph(arr)

## Now with a force directed algorithm
multigraph(arr, layout = "force")

## As weighted graph
multigraph(arr, weighted = TRUE)

## As signed graph
multigraph(arr, signed = TRUE)

## With loops and a costumized vertex size
multigraph(arr, cex = 3, loops = TRUE)

}

\keyword{ graphics }
