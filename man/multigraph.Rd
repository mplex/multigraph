\name{multigraph}
\alias{multigraph}
\title{
Multigraph
}
\description{
A function to create and manipulate multigraphs with different layouts
}
\usage{
multigraph(net, layout = c("circ", "stress", "rand"), directed = TRUE, 
           collRecip = TRUE, coord = NULL, tcol = 1, bwd = 1, clu = NULL, 
           cex = NULL, tcex = NULL, showLbs = TRUE, showAtts = TRUE, att = NULL, 
           lbat = "1", seed = NULL, maxiter = 100, alpha = c(1, 1, 1), 
           main = NULL, cex.main, bg, mar, pos, asp, lwd, pch, lty, ecol, vcol, 
           vcol0, hds, vedist, rot, mirrorX, mirrorY, col, ...)
}
\arguments{
  \item{net}{
	an array; usually with three dimensions of stacked matrices where the multiple relations are placed.
}
  \item{layout}{
	the layout for the visualization:

       \code{circ} circular

       \code{stress} ``stress majorization'' force-directed algorithm 

       \code{rand} random
}
  \item{directed}{
	(logical) whether or not graph is directed or unidrected
}
  \item{collRecip}{
	(logical) whether or not collapse reciprocated edges in the unidrected graph
}
  \item{coord}{
	(optional) data frame with the coordinates of the vertices. If given \code{layout} is ignored
}
  \item{tcol}{
	(optional) the color of the vertices' labels
}
  \item{bwd}{
	(optional) the width of the bundle edges. Ranges from \code{0} (edges collapsed) to the default \code{1} (depending on the vertices' size)
}
  \item{clu}{
	(optional) the clustering of the vertices (see \emph{details})
}
  \item{cex}{
	(optional) the size of the vertices
}
  \item{tcex}{
	size of the vertex labels. If \code{NULL}, it depends on the value of \code{cex}
}
  \item{showLbs}{
	(logical) whether or not show the vertex labels
}
  \item{showAtts}{
	(logical) whether or not show the vertex attribute labels
}
  \item{att}{
	(optional) a vector or an array representing the vertices' attributes
}
  \item{lbat}{
	(optional) the labels for the vertices' attributes
}
  \item{seed}{
	(optional) the random seed number for the vertices' initial coordinates. Only for \code{stress} and \code{rand}
}
  \item{maxiter}{
	(optional) the maximum number of iterations in the stress majorization (default \code{100})
}
  \item{alpha}{
	vector (vertex, edge, \code{bg}) with the alpha color transparecy
}
  \item{main}{
	(optional) title of the plot
}
  \item{cex.main}{
	(optional) the size of the plot's title
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
  \item{lwd}{
	(optional) the width of the edges
}
  \item{pch}{
	(optional) the symbol representing the vertices
}
  \item{lty}{
	(optional) the shape of the edges
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
	(optional, experimental) arcs' head scale
}
  \item{vedist}{
	(optional, experimental) a real with vertex - edge distance
}
  \item{rot}{
	(optional) rotation of the graph in degrees
}
  \item{mirrorX}{
	(optional) mirror of the X axis
}
  \item{mirrorY}{
	(optional) mirror of the Y axis
}
  \item{col}{
	(optional) same as \code{vcol}
}
  \item{\dots}{
	Additional items
}
}
\details{
Multigraph are graph having parallel edges depicting different types of relations in a network. By default a circular layout is applied where each type of tie has a distinctive shape and gray color scale. For a better visualization, undirected multigraphs automatically collapse the reciprocal ties, and there as an argument to prevent this to happen. It is possible to combine the symbols and color of vertices by assigning a class to each network member in the clustering option. Vertices can also have different sizes by specifying the argument with a vector with a length size similar to the network order.
}
\value{
A plot of the network as a multigraph.
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
\code{\link{bmgraph}}
}
\examples{
## Create the data: two binary relations among three elements
arr <- round( replace( array(runif(18), c(3,3,2)), array(runif(18),
       c(3,3,2))>.5, 3 ) )

## Plot the multigraph of this network
multigraph(arr)

## Now with a force directed algorithm
multigraph(arr, layout = "stress")

}
\keyword{ graphics }

