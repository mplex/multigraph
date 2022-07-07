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
          att = NULL, lbat = "1", main = NULL, cex.main, bg, mar, directed, valued, 
          collRecip, cex, pos, lwd, lty, col, ecol, vcol, vcol0, asp, seed = NULL, 
          maxiter = 100, bwd, clu, pch, rot, mirrorX, mirrorY, mirrorV, mirrorH, hds, 
          vedist, jitter, sort, add, adc, perm, ffamily, fstyle, fsize, fcol, vclu, ...)
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
   (optional) scope of the graph (see details)
}
  \item{coord}{
	(optional) data frame with the coordinates of the vertices. If coordinates are given then the \code{layout} option is ignored
}
  \item{alpha}{
	vector (vertex, edge, \code{bg}) with the alpha color transparency
}
  \item{showLbs}{
	(optional and logical) whether or not to show the vertex labels when dimnames available
}
  \item{showAtts}{
	(optional and logical) whether or not to show the vertex attribute labels
}
  \item{att}{
	(optional) a vector or an array representing the vertex attributes
}
  \item{lbat}{
	(optional) labels for the vertex attributes
}
  \item{main}{
	(optional) title of the plot
}
  \item{cex.main}{
	(optional) size of the plot's title
}
  \item{bg}{
	(optional) background color of the plot
}
  \item{mar}{
	(optional) margins of the plot
}
  \item{directed}{
	(optional and logical) whether or not the graph is directed or undirected
}
  \item{valued}{
	(optional and logical) whether or not the graph is valued or with dichotomous data
}
  \item{collRecip}{
	(optional and logical) whether or not collapse reciprocated edges in the undirected graph
}
  \item{cex}{
	(optional) size of the vertices
}
  \item{pos}{
	(optional) position of the vertices' labels (\code{0} means ``at the center of the vertex'')
}
  \item{lwd}{
	(optional) width of the edges. Ignored if \code{valued} is set to \code{TRUE}
}
  \item{lty}{
	(optional) shape of the edges
}
  \item{col}{
	(optional) alias for \code{vcol}
}
  \item{ecol}{
	(optional) color of the edges
}
  \item{vcol}{
	(optional) color of the vertices
}
  \item{vcol0}{
	(optional) color of the vertices' contour (only works for \code{pch 21} through \code{25}
}
  \item{asp}{
	(optional) aspect ratio of the plot
}
  \item{seed}{
	(optional) random seed number for the vertices' initial coordinates. Ignored except for \code{force}, \code{stress} and \code{rand}
}
  \item{maxiter}{
	(optional) maximum number of iterations in layout algorithms. Ignored except for \code{force}, \code{stress} and \code{rand}
}
  \item{bwd}{
	(optional) width of the bundle edges. Ranges from \code{0} (edges collapsed) to the default \code{1} (depending on the vertices' size). For \code{valued} a value greater than one is possible
}
  \item{clu}{
	(optional) clustering of the vertices (see \emph{details})
}
  \item{pch}{
	(optional) symbol representing the vertices
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
	\emph{same as \code{mirrorX}}
}
  \item{mirrorH}{
	\emph{same as \code{mirrorY}}
}
  \item{hds}{
	(optional and experimental) arcs' head scale
}
  \item{vedist}{
	(optional and experimental) a real number with vertex - edge distance
}
  \item{jitter}{
	(optional) jitter in \code{stress} or \code{CA}
}
  \item{sort}{
	(optional and logical) sort the vertex labels
}
  \item{add}{
	(optional) add nodes to the graph's domain
}
  \item{adc}{
	(optional) add nodes to the graph's codomain
}
  \item{perm}{
	(optional) a list of vectors for the permutation of network members in both the domain and codomain
}
  \item{ffamily}{
	(optional) font family
}
  \item{fstyle}{
	(optional) font style
}
  \item{fsize}{
	(optional) font size
}
  \item{fcol}{
	(optional) font color
}
  \item{vclu}{
	(optional) a list of vectors the clustering information in both the domain and the codomain
%%     ~~Describe \code{vclu} here~~ \code{cluc} was 
}
  \item{\dots}{
	Additional argument items (see e.g. \code{\link[graphics:par]{par}})
}
}
\details{
Bipartite graphs are visualization devices for two-mode networks. Although this type of data would typically record as a data frame, it is possible to use even three-dimensional arrays where each level corresponds to a particular type of tie. Thus the bipartite graphs, in this case, will be depicted with parallel edges. 
Besides, it is possible to obtain a figure of the bipartite network using the binomial approach to two-mode data and plot it with a force-directed algorithm. 
}
\value{
A plot of the two-mode network as a bipartite graph or multigraph with a projection
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
