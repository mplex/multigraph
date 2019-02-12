\name{bmgraph}
\alias{bmgraph}
\title{
Bipartite multigraph
}
\description{
A function to create and manipulate bipartite multigraphs
}
\usage{
bmgraph(net, layout = c("bip", "bip3", "bip3e", "bipc", "force", "rand", "circ", 
       "stress", "CA", "circ2"), scope, coord, alpha = c(1, 1, 1), showLbs, showAtts, 
       att = NULL, lbat = "1", main = NULL, cex.main, bg, mar, directed, weighted, 
       collRecip, cex, pos, lwd, lty, col, ecol, vcol, vcol0, asp, seed = NULL, 
       maxiter = 100, bwd, clu, pch, fcol, rot, mirrorX, mirrorY, mirrorV, mirrorH, 
       hds, vedist, jitter, sort, add, adc, perm, ffamily, fstyle, fsize, ...)
}
\arguments{
  \item{net}{
	data frame or array representing the two-mode network (see \emph{details})
}
  \item{layout}{
    the visualization layout:   

       \code{bip} (default) bipartite graph   

       \code{bip3} bipartite graph with three columns   

       \code{bip3e} bipartite graph with three columns for events   

       \code{bipc} ``clustered'' bipartite graph   

       \code{force} force-directed algorithm    

       \code{rand} random   

       \code{circ} circular   

       \code{stress} stress-majorization algorithm    

       \code{CA} correspondence analysis   

       \code{circ2} two semi-circles  
}
  \item{scope}{
   (optional) the scope of the graph (see details)
}
  \item{coord}{
	(optional) data frame with the coordinates of the vertices. If coordinates are given then the \code{layout} option is ignored
}
  \item{alpha}{
	vector (vertex, edge, \code{bg}) with the alpha color transparecy
}
  \item{showLbs}{
	(optional and logical) whether or not show the vertex labels when dimnames available
}
  \item{showAtts}{
	(optional and logical) whether or not show the vertex attribute labels
}
  \item{att}{
	(optional) a vector or an array representing the vertex attributes
}
  \item{lbat}{
	(optional) the labels for the vertex attributes
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
  \item{directed}{
	(optional and logical) whether or not the graph is directed or unidrected
}
  \item{weighted}{
	(optional and logical) whether or not the graph is weighted or dichotomous
}
  \item{collRecip}{
	(optional and logical) whether or not collapse reciprocated edges in the unidrected graph
}
  \item{cex}{
	(optional) the size of the vertices
}
  \item{pos}{
	(optional) the position of the vertices' labels (\code{0} means ``at the center of the vertex'')
}
  \item{lwd}{
	(optional) the width of the edges. Ignored if \code{weighted} is set to \code{TRUE}
}
  \item{lty}{
	(optional) the shape of the edges
}
  \item{col}{
	(optional) alias for \code{vcol}
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
  \item{asp}{
	(optional) the aspect ratio of the plot
}
  \item{seed}{
	(optional) the random seed number for the vertices' initial coordinates. Ignored except for \code{force}, \code{stress} and \code{rand}
}
  \item{maxiter}{
	(optional) the maximum number of iterations in layout algorithms. Ignored except for \code{force}, \code{stress} and \code{rand}
}
  \item{bwd}{
	(optional) the width of the bundle edges. Ranges from \code{0} (edges collapsed) to the default \code{1} (depending on the vertices' size). For \code{weighted} a value greater than one is possible
}
  \item{clu}{
	(optional) the clustering of the vertices (see \emph{details})
}
  \item{pch}{
	(optional) the symbol representing the vertices
}
  \item{fcol}{
	the font color
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
  \item{mirrorV}{
	same as \code{mirrorX}
}
  \item{mirrorH}{
	same as \code{mirrorY}
}
  \item{hds}{
	(optional and experimental) arcs' head scale
}
  \item{vedist}{
	(optional and experimental) a real number with vertex - edge distance
}
  \item{jitter}{
	(optional) the jitter in \code{stress} or \code{CA}
}
  \item{sort}{
	(optional and logical) sort the vertex labels
}
  \item{add}{
	add nodes to the graph's domain
}
  \item{adc}{
	add nodes to the graph's codomain
}
  \item{perm}{
	(optional) a list of vectors for the permutation of network members in both the domain and codomain
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
  \item{\dots}{
	Additional argument items (see e.g. \code{\link[graphics:par]{par}})
}
}
\details{
Bipartite graphs are visualization devices for two-mode networks. Although this type of data would typically record as a data frame, it is possible to use even three dimensional arrays where each level corresponds to a particular type of tie, and thus the bipartite graphs in this case will be depicted with parallel edges. Besides, it is possible to obtain a graph of the bipartite network using the binomial approach to two-mode data, and plot it with a forced directed algorithm.

The \code{clu} argument serves to identify classes of actors and events that is used in the column disposition with the \code{"bipc"} layout option.
}
\value{
A plot of the two-mode networks as bipartite graph or multigraph
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
\code{\link{multigraph}}, \code{\link{frcd}}, \code{\link{stsm}}, \code{\link{conc}}
}
\examples{
## Create the data: two binary relations among three elements
arr <- round( replace( array(runif(18), c(3,3,2)), array(runif(18),
       c(3,3,2))>.5, 3 ) )

## Plot this network as Bipartite graph
bmgraph(arr)

## Now with a Force Directed algorithm
bmgraph(arr, layout = "force")

## And with a Correspondence Analysis method
bmgraph(arr, layout = "CA", asp = NA)

}
\keyword{ graphics }
