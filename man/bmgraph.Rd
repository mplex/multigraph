\name{bmgraph}
\alias{bmgraph}
\title{
Bipartite multigraph
}
\description{
A function to create and manipulate bipartite multigraphs with different layouts
}
\usage{
bmgraph(net, layout = c("bip", "bip3", "bip3e", "bip4", "stress", "rand", "circ"), 
        coord = NULL, tcex = NULL, alpha = c(1, 1, 1), showLbs = TRUE, 
        showAtts = TRUE, att = NULL, lbat = "1", main = NULL, cex.main, bg, mar, 
        cex, pos, lwd, lty, ecol, vcol, vcol0, asp, directed, collRecip, 
        seed = NULL, maxiter, bwd, clu, pch, tcol, hds, vedist, rot, 
        mirrorX, mirrorY, col, ...)
}
\arguments{
  \item{net}{
	data frame or array representing the two-mode network (see \emph{details})
}
  \item{layout}{
	the layout for the visualization:

       \code{bip} (default) bipartite graph

       \code{bip3} bipartite graph with three columns

       \code{bip3e} bipartite graph with three columns for events

       \code{bip4} ``square'' bipartite graph

       \code{stress} ``stress majorization'' force-directed algorithm 

       \code{rand} random

       \code{circ} circular

}
  \item{coord}{
	(optional) data frame with the coordinates of the vertices.  If given \code{layout} is ignored
}
  \item{tcex}{
	size of the vertex labels. If \code{NULL}, it depends on the value of \code{cex}
}
  \item{alpha}{
	vector (vertex, edge, \code{bg}) with the alpha color transparecy
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
  \item{cex}{
	(optional) the size of the vertices
}
  \item{pos}{
	(optional) the position of the vertices' labels (\code{0} means ``in the middle of the vertex'')
}
  \item{lwd}{
	(optional) the width of the edges
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
  \item{asp}{
	(optional) the aspect ratio of the plot
}
  \item{directed}{
	(optional, logical) whether or not graph is unidrected (default)
}
  \item{collRecip}{
	(optional, logical) whether or not collapse reciprocated edges in the unidrected graph (default \code{TRUE})
}
  \item{seed}{
	(optional) the random seed number for the vertices' initial coordinates. Ignored for other than \code{stress} or \code{rand}
}
  \item{maxiter}{
	(optional) the maximum number of iterations in the stress majorization (default \code{100})
}
  \item{bwd}{
	(optional) the width of the bundle edges. Ranges from \code{0} (edges collapsed) to the default \code{1} (depending on the vertices' size)
}
  \item{clu}{
	(optional) the clustering of the vertices (see \emph{details})
}
  \item{pch}{
	(optional) the symbol representing the vertices
}
  \item{tcol}{
	(optional) the color of the vertices' labels
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
Bipartite graphs are visualization devices for two-mode networks. Although this type of data would typically record as a data frame, it is possible to use even three dimensional arrays where each level corresponds to a particular type of tie, and thus the bipartite graphs in this case will be depicted with parallel edges. Besides, it is possible to obtain a graph of the bipartite network using the binomial approach to two-mode data, and plot it with a forced directed algorithm.
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
\code{\link{multigraph}}
}
\examples{
## Create the data: two binary relations among three elements
arr <- round( replace( array(runif(18), c(3,3,2)), array(runif(18),
       c(3,3,2))>.5, 3 ) )

## Plot this network as bipartite graph
bmgraph(arr)

}
\keyword{ graphics }

