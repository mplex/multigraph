\name{bmgraph}
\alias{bmgraph}
\title{
Bipartite Graphs for Affiliation Networks
}
\description{
A function to create and manipulate bipartite graphs for affiliation or two-mode network networks.
}
\usage{
bmgraph(net, layout = c("bip", "bip3", "bip3e", "bipc", "CA", "force", "rand", "stress", 
        "circ", "conc", "circ2"), scope, coord, alpha = c(1, 1, 1), showLbs, showAtts, 
        att = NULL, lbat = "1", main = NULL, cex.main, col.main, font.main, directed, 
        valued, collRecip, cex, pos, lwd, lty, col, ecol, vcol, vcol0, bg, asp, 
        seed = NULL, maxiter = 100, bwd, clu, pch, rot, mirrorX, mirrorY, mirrorV, 
        mirrorH, hds, vedist, jitter, sort, add, adc, perm, ffamily, fstyle, fsize, fcol, 
        mar, scl, new, ...)
}
\arguments{
  \item{net}{A data frame or array representing an affiliation network incidence matrix.}
  \item{layout}{
    Visualization layout:
     \itemize{ 
       \item \code{bip} Bipartite graph (default)
       \item \code{bip3} Bipartite graph with three columns
       \item \code{bip3e} Bipartite graph with three columns for events
       \item \code{bipc} \dQuote{Clustered} bipartite graph
       \item \code{force} Apply force-directed algorithm
       \item \code{rand} Random
       \item \code{circ} Circular
       \item \code{stress} Apply stress-majorization algorithm
       \item \code{CA} Apply correspondence analysis
       \item \code{circ2} Two semi-circles.
      }}
  \item{scope}{(optional) List with scope of graph, edge and vertex characteristics.}
  \item{coord}{(optional, ignores \code{layout}) Two-column data frame of vertex coordinates in graph.}
  \item{alpha}{Vector with the alpha color transparencies for (vertex, edge, \code{bg}).}
  \item{showLbs}{(optional and logical) Show vertex labels in graph when available?}
  \item{showAtts}{(optional and logical) Show vertex attribute labels in graph?}
  \item{att}{(optional) Vector or array representing vertex attributes.}
  \item{lbat}{(optional) Labels for vertex attributes.}
  \item{main}{(optional) Main title of the plot graph.}
  \item{cex.main}{Font size of title in plot graph.}
  \item{col.main}{Font colour of title in plot graph.}
  \item{font.main}{Font family of title in plot graph.}
  \item{directed}{(optional and logical) Is the bipartite graph directed?}
  \item{valued}{ Is the bipartite graph valued?}
  \item{collRecip}{(optional and logical) Collapse reciprocated edges in (undirected) graph?}
  \item{cex}{(optional) Vertex size.}
  \item{pos}{(optional) Position of labels from vertex:
   \itemize{ 
       \item {\code{1}}: Below
       \item {\code{2}}: Left
       \item {\code{3}}: Above
       \item {\code{4}}: Right
       \item {\code{0}}: Center
      }}
  \item{lwd}{(optional, ignored for \code{valued}) Width of edges.}
  \item{lty}{(optional) Shape of edges.}
  \item{col}{(optional) Alias for \code{vcol}.}
  \item{ecol}{(optional) Colour of edges.}
  \item{vcol}{(optional) Colour of vertices.}
  \item{vcol0}{(optional, only for \code{pch=21} through \code{pch=25}) Colour of vertices contour.}
  \item{bg}{(optional) Background colour of plot graph.}
  \item{asp}{(optional, \code{asp=NA} cover whole canvas) Aspect ratio of plot.}
  \item{seed}{(optional, for \code{force}, \code{stress} and \code{rand}) Random seed number for initial coordinates of vertices.}
  \item{maxiter}{(optional, for \code{force}, \code{stress} and \code{rand}) Maximum number of iterations in layout algorithm.}
  \item{bwd}{(optional, default \code{1}) Real number for width between bundle edges.}
  \item{clu}{(optional) Clustering of vertices (see \emph{details}).}
  \item{pch}{(optional) Symbol(s) for vertices.}
  \item{rot}{(optional) Clockwise rotation in degrees on graph.}
  \item{mirrorX}{(optional) Mirror by \eqn{X} axis.}
  \item{mirrorY}{(optional) Mirror by \eqn{Y} axis.}
  \item{mirrorV}{\dQuote{Same as \code{mirrorX}}.}
  \item{mirrorH}{\dQuote{Same as \code{mirrorY}}.}
  \item{hds}{(optional, experimental) Head scale for arcs.}
  \item{vedist}{(optional, experimental) Vertex--edge distance.}
  \item{jitter}{(optional, in \code{stress} or \code{CA}) Amount of jitter to add on coordinates.}
  \item{sort}{(optional and logical) Sort vertex labels?}
  \item{add}{(optional) Nodes to add nodes to network domain.}
  \item{adc}{(optional) Nodes to add nodes to network co-domain.}
  \item{perm}{(optional) List of vectors for permutation of network nodes in domain and codomain.}
  \item{ffamily}{(optional) Font family of vertex labels in graph.}
  \item{fstyle}{(optional) Font style of vertex labels in graph with options: \code{bold}, \code{italic}, \code{bolditalic}.}
  \item{fsize}{(optional) Font size of vertex labels in graph.}
  \item{fcol}{(optional) Font colour of text in graph vertices.}
  \item{mar}{(optional) Plot margins.}
  \item{scl}{(optional, experimental) Numerical scalar (\eqn{x} and \eqn{y}) or vector (\eqn{x}, \eqn{y}) of graph scale.}
  \item{new}{(optional, logical) Plot new graph on canvas afterwards?}
  \item{\dots}{Additional argument items.}
}
\details{
Bipartite graphs serve as visual aids for affiliation or two-mode networks. 
While affiliation networks are typically represented as data frames, they can also be visualized using three-dimensional arrays. 
In such case, each level corresponds to a specific type of connection, resulting in parallel edges within the bipartite graph. 
Additionally, a bipartite network can be generated using a force-directed algorithm to create a visual representation. 

With bipartite graphs consisting of two sets of vertices, clustering information, such as vertex colours, can be stored in a list \code{vclu} with two vectors, one for each vertex set. 
It is possible to group all members of a vertex set into a single class by setting the corresponding vector to \code{NULL}.
}
\value{
A plot of the two-mode network as a bipartite graph or multigraph with a projection.
}
%\references{} ~
%\author{} ~
%\note{} ~
%\section{Warning }{....} ~
\seealso{
\code{\link{multigraph}}, \code{\link{frcd}}, \code{\link{stsm}}, \code{\link{conc}}
}
\examples{
# Two binary relations among three elements
arr <- round( replace( array(runif(18), c(3,3,2)), array(runif(18),
       c(3,3,2))>.5, 3 ) )

# Plot network as bipartite graph
bmgraph(arr)

# Apply force directed algorithm for plot
bmgraph(arr, layout = "force")

# Apply Correspondence Analysis
bmgraph(arr, layout = "CA", asp = NA)
}
\keyword{graphics}
