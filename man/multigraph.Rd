\name{multigraph}
\alias{multigraph}
\title{
Multigraphs and valued multigraphs
}
\description{
A function to create and manipulate multigraphs and valued multigraphs with different layout options
}
\usage{
multigraph(net, layout=c("circ", "force", "stress", "conc", "rand"), scope, directed=TRUE, 
           loops, signed, valued, values, lbs, showLbs, att, lbat, showAtts, main=NULL, 
           cex.main, col.main, font.main, coord, collRecip, undRecip, seed=NULL, 
           maxiter=100, clu, cex, cex2, pch, lwd, lty, vcol, vcol0, col, ecol, bwd, bwd2, 
           pos, bg, bg2, asp, drp, add, swp, swp2, alpha=c(1, 1, 1, 1), rot, mirrorX, 
           mirrorY, mirrorD, mirrorL, mirrorV, mirrorH, scl, hds, vedist, mar, ffamily, 
           fstyle, fsize, fsize2, fcol, fcol2, lclu, sel, new, mai, lscl, rm.isol, ...)
}
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
  \item{scope}{
   (optional) the scope of the graph (see details)
}
  \item{directed}{
	(logical) whether or not the graph is directed or unidrected
}
  \item{loops}{
	(optional, logical, and experimental) plot graph loops?
}
  \item{signed}{
	(optional and logical) whether or not the graph is a signed structure
}
  \item{valued}{
	(optional and logical) whether the graph is depicyed as valued or not
}
  \item{values}{
	(optional and logical) print the values of the bonds in edges?
}
  \item{lbs}{
	(optional) the vertices labels
}
  \item{showLbs}{
	(optional and logical) whether or not show the vertex labels
}
  \item{att}{
	(optional) a vector or an array representing the vertex attributes
}
  \item{lbat}{
	(optional) the labels for the vertices' attributes
}
  \item{showAtts}{
	(optional and logical) whether or not show the vertex attribute labels
}
  \item{main}{
	(optional) title of the plot
}
  \item{cex.main}{
	(optional) the size of the plot's title
}
  \item{col.main}{
%%     ~~Describe \code{col.main} here~~
}
  \item{font.main}{
%%     ~~Describe \code{font.main} here~~
}
  \item{coord}{
	(optional) data frame with the coordinates of the vertices. If coordinates are given then the \code{layout} option is ignored
}
  \item{collRecip}{
	(optional and logical) whether or not collapse reciprocated edges in the unidrected graph
}
  \item{undRecip}{
	(optional and logical) whether or not plot reciprocated edges as undirected
}
  \item{seed}{
	(optional) the random seed number for the vertices' initial coordinates. Ignored for \code{circ} and \code{conc}
}
  \item{maxiter}{
	(optional) the maximum number of iterations in layout algorithms. Only for \code{force}, \code{stress}, and \code{rand}
}
  \item{clu}{
	(optional) the clustering of the vertices (see \emph{details})
}
  \item{cex}{
	(optional) the size of the vertices
}
  \item{cex2}{
	the size of the background for the values with the \code{valued} option
}
  \item{pch}{
	(optional) the symbol representing the vertices
}
  \item{lwd}{
	(optional) the width of the edges. Ignored if \code{weighted} is set to \code{TRUE}
}
  \item{lty}{
	(optional) the shape of the edges
}
  \item{vcol}{
	(optional) the color of the vertices
}
  \item{vcol0}{
	(optional) the color of the vertices' contour (only works for \code{pch 21} through \code{25}
}
  \item{col}{
	(optional) alias for \code{vcol}
}
  \item{ecol}{
	(optional) the color of the edges
}
  \item{bwd}{
	(optional) the width of the bundle edges. 
}
  \item{bwd2}{
	(optional) the width of the bundle loop edges. 
}
  \item{pos}{
	(optional) the position of the vertices' labels (\code{0} means ``in middle of vertex'')
}
  \item{bg}{
	(optional) the background color of the plot
}
  \item{bg2}{
	(optional) the background color for \code{values}
}
  \item{asp}{
	(optional) the aspect ratio of the plot
}
  \item{drp}{
	(optional) for \code{valued} networks, drop values less than the specified
}
  \item{add}{
	(optional) nodes to add to the graph
}
  \item{swp}{
	(optional and logical) whether or not swap the bundle patterns
}
  \item{swp2}{
	(optional and logical) whether or not swap reciprocals
}
  \item{alpha}{
	vector (vertex, edge, \code{bg}) with the alpha color transparecy
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
  \item{scl}{
	(optional and experimental) numerical scalar (x and y) or vector (x, y) of the graph's scale 
}
  \item{hds}{
	(optional and experimental) arcs' head scale
}
  \item{vedist}{
	(optional and experimental) a real number with vertex - edge distance
}
  \item{mar}{
	(optional) the margins of the plot
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
  \item{fsize2}{
	the font size for \code{values}
}
  \item{fcol}{
	the font color
}
  \item{fcol2}{
	the font color for \code{values}
}
  \item{lclu}{
	(optional, vector) ``levels'' in \code{clu} (see Details)
}
  \item{sel}{
	(optional, vector) selection of node's labels to plot
}
  \item{new}{
	(optional, logical) new graph on an existing plot?
}
  \item{mai}{
	(optional, vector) plot inner margins
}
  \item{lscl}{
	(optional for valued graphs) loop scale
}
  \item{rm.isol}{
	(optional) remove isolated vertices?
}
  \item{\dots}{
	Additional argument items (see e.g. \code{\link[graphics:par]{par}})
}
}
\details{
Multigraphs are graphs having parallel edges depicting different types of relations in a network. By default, a circular layout is applied where each type of tie has a distinctive shape and gray color scale. 
For better visualization, undirected multigraphs automatically collapse the reciprocal relations, and there is an argument to prevent this from happening. It is possible to combine the symbols and colors of vertices by assigning a class to each network member in the clustering option. Vertices can also have different sizes by specifying the argument with a vector with a length size similar to the network order.

Since this function can handle a large number of arguments, these can be stored as a list object that is passed through the \code{scope} option. In this case, a vector made of lists and scalars or combinations of these is accepted for describing characteristics.

The bundle width specified by \code{bwd} (and \code{bwd2} for loops) ranges from \code{0} (edges collapsed) to the default \code{1} (depending on the vertices' size). For the \code{valued} option, numbers higher than one are possible. 
Use \code{vedist} to adjust vertex--edge distance for large and dense networks.

In some cases, such as when working with dynamic networks, it is needed to specify the ordering of the ``levels'' of the clustering information given in \code{clu}, and this is done in argument \code{lclu}. 

When using \code{new} for plotting the graph with a background image, the previous plot(s), however, can require having an equivalent command to \code{graphics::plot.new()} (cf. e.g. \code{sdam::plot.map()} function).
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
\code{\link{bmgraph}}, \code{\link{ccgraph}}, \code{\link{frcd}}, \code{\link{stsm}}, \code{\link{conc}}
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
