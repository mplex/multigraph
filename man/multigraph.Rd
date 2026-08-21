\name{multigraph}
\alias{multigraph}
\title{
Multigraphs and Valued Multigraphs
}
\description{
A function to create and manipulate multigraphs and valued multigraphs with different layout options.
}
\usage{
multigraph(net, layout = c("circ", "force", "stress", "conc", "rand"), scope, 
           directed = TRUE, loops, signed, valued, values, lbs, showLbs, att, 
           lbat, showAtts, main = NULL, cex.main, col.main, font.main, coord, 
           collRecip, undRecip, seed = NULL, maxiter = 100, clu, cex, cex2, 
           pch, lwd, lty, vcol, vcol0, col, ecol, bwd, bwd2, pos, bg, bg2, asp, 
           drp, add, swp, swp2, flip, alpha = c(1, 1, 1, 1), rot, mirrorX, 
           mirrorY, mirrorD, mirrorL, mirrorV, mirrorH, scl, hds, vedist, 
           vldist, mar, ffamily, fstyle, fsize, fsize2, fcol, fcol2, lclu, sel, 
           new, mai, lscl, rm.isol, nr, jitter, ...)
}
\arguments{
  \item{net}{Array with usually three dimensions of stacked matrices where the multiple relations are placed.}
  \item{layout}{
    Visualization layout: 
   \itemize{ 
       \item \code{circ} Circular 
       \item \code{force} Force-directed 
       \item \code{stress} Stress-majorization 
       \item \code{conc} Concentric 
       \item \code{rand} Random 
   }}
  \item{scope}{(optional) Scope of graph, edges, and vertices (see \emph{details}).}
  \item{directed}{(logical) Are graph edges directed?}
  \item{loops}{(optional and logical) Plot loops in graph?}
  \item{signed}{(optional and logical) Has the graph a signed structure?}
  \item{valued}{(optional and logical) Depict a valued graph?}
  \item{values}{(optional and logical) Print values of bonds in edges?}
  \item{lbs}{(optional) Vector of customized vertex labels.}
  \item{showLbs}{(optional and logical) Display vertex labels in graph?}
  \item{att}{(optional) Vector or array representing vertex attributes.}
  \item{lbat}{(optional) Vector of labels for vertex attributes.}
  \item{showAtts}{(optional and logical) Display vertex attribute labels?}
  \item{main}{(optional) Main title of the plot graph.}
  \item{cex.main}{Font size of title in plot graph.}
  \item{col.main}{Font colour of title in plot graph.}
  \item{font.main}{Font family of title in plot graph.}
  \item{coord}{(optional) Data frame with graph coordinates for vertices.}
  \item{collRecip}{(optional and logical) Collapse reciprocated edges in graph?}
  \item{undRecip}{(optional and logical) Plot reciprocated edges as undirected?}
  \item{seed}{(optional, for \code{force}, \code{stress} and \code{rand}) Random seed number for initial coordinates of vertices.}
  \item{maxiter}{(optional when \code{seed}) Maximum number of iterations.}
  \item{clu}{(optional) Clustering for vertices (see \emph{details}).}
  \item{cex}{(optional) Vertex size.}
  \item{cex2}{(optional, for \code{values}) Size of background in values.}
  \item{pch}{(optional) Symbol(s) for vertices.}
  \item{lwd}{(optional, ignored for \code{valued}) Width of edges.}
  \item{lty}{(optional) Shape of edges.}
  \item{vcol}{(optional) Colour of vertices.}
  \item{vcol0}{(optional, only for \code{pch=21} through \code{pch=25}) Colour of vertices contour.}
  \item{col}{(optional) Alias for \code{vcol}.}
  \item{ecol}{(optional) Colour of edges.}
  \item{bwd}{(optional, default \code{1}) Real number for width between bundle edges.}
  \item{bwd2}{(optional) Width of bundle loop edges.}
  \item{pos}{(optional) Position of labels from vertex:
   \itemize{ 
       \item {\code{1}}: Below
       \item {\code{2}}: Left
       \item {\code{3}}: Above
       \item {\code{4}}: Right
       \item {\code{0}}: Centre
      }}
  \item{bg}{(optional) Background colour of plot graph.}
  \item{bg2}{(optional, for \code{values}) Background colour.}
  \item{asp}{(optional, \code{asp=NA} cover whole canvas) Aspect ratio of plot.}
  \item{drp}{(optional) for \code{valued}) Drop values in  network less than specified.}
  \item{add}{(optional) Vector of vertices to add to graph.}
  \item{swp}{(optional and logical) Swap bundle patterns?}
  \item{swp2}{(optional and logical) Swap reciprocal arcs?}
  \item{flip}{(optional and logical) Flip loops?}
  \item{alpha}{Vector with the alpha colour transparencies for (vertex, edge, \code{bg}).}
  \item{rot}{(optional) Clockwise rotation in degrees on graph.}
  \item{mirrorX}{(optional) Mirror by \eqn{X} axis.}
  \item{mirrorY}{(optional) Mirror by \eqn{Y} axis.}
  \item{mirrorD}{(optional) Mirror reflection across diagonal \eqn{Y=X}.}
  \item{mirrorL}{(optional) Mirror reflection across diagonal \eqn{Y=-X}.}
  \item{mirrorV}{\dQuote{Same as \code{mirrorX}}.}
  \item{mirrorH}{\dQuote{Same as \code{mirrorY}}.}
  \item{scl}{(optional) Numerical scalar (\eqn{x} and \eqn{y}) or vector (\eqn{x}, \eqn{y}) of graph scale.}
  \item{hds}{(optional, for \code{directed}) Head scale of arcs.}
  \item{vedist}{(optional) Vertex--edge distance.}
  \item{vldist}{(optional) Vertex--loop distance.}
  \item{mar}{(optional) Vector with margins for plot.}
  \item{ffamily}{(optional) Font family of vertex labels in graph.}
  \item{fstyle}{(optional) Font style of vertex labels in graph with options: \code{bold}, \code{italic}, \code{bolditalic}.}
  \item{fsize}{(optional) Font size of vertex labels in graph.}
  \item{fsize2}{(optional) Font size for \code{values}.}
  \item{fcol}{(optional) Font colour of text in graph vertices.}
  \item{fcol2}{(optional) Font colour of text in \code{values}.}
  \item{lclu}{(optional) Vector with \dQuote{levels} in \code{clu} (see \emph{details}).}
  \item{sel}{(optional) Vector with selection of vertex labels to plot.}
  \item{new}{(optional, logical) Plot a new graph on existing canvas?}
  \item{mai}{(optional) Vector with inner margins for plot.}
  \item{lscl}{(optional, for \code{valued}) Loop scale in graph.}
  \item{rm.isol}{(optional and logical) Remove isolated vertices?}
  \item{nr}{(optional, for \code{conc}) Number of radii.}
  \item{jitter}{(optional, for \code{stress}) Amount of jitter to add on coordinates.}
  \item{\dots}{Additional argument items.}
}
\details{
Multigraphs contain parallel edges that represent different relation types within a network. By default, the function uses a circular layout in which tie types are distinguished by shape and a grey colour scale. 
For undirected multigraphs, reciprocal relations are collapsed automatically, although this behaviour can be disabled. 
Vertex symbols and colours can be combined by assigning each network member to a class through the clustering option. 
Vertex sizes can also be specified with a vector whose length matches the network order.

Because the function accepts many arguments, these may be collected in a list and supplied through the \code{scope} option. 
The list may include vectors, lists, scalars, or combinations of these to describe graph characteristics.

The bundle width \code{bwd}, and \code{bwd2} for loops, ranges from \code{0}, where edges are collapsed, to the default value of \code{1}, depending on vertex size. 
Values above one are allowed when the \code{valued} option is used. 
For large or dense networks, argument \code{vedist} adjusts the distance between vertices and edges.

For some applications, such as dynamic networks, the ordering of clustering \dQuote{levels} in \code{clu} must be specified; this is done with argument \code{lclu}.

If coordinates are supplied, the \code{layout} option is ignored.

When argument \code{new} is used to plot a graph with a background image, previous plots may need to be initialized with package functions equivalent to \code{\link[graphics:plot.new]{graphics::plot.new()}}.
}
\value{
A plot of the network as a multigraph or a valued multigraph.
}
%\references{} ~
%\author{} ~
%\note{} ~
%\section{Warning }{....} ~
\seealso{
\code{\link{bmgraph}}, \code{\link{ccgraph}}, \code{\link{frcd}}, \code{\link{stsm}}, \code{\link{conc}}
}
\examples{
# Two binary relations among three elements
arr <- round( replace( array(runif(18), c(3,3,2)), array(runif(18),
       c(3,3,2))>.5, 3 ) )

# Plot multigraph of this network
multigraph(arr)

# Plot with a force directed algorithm
multigraph(arr, layout = "force")

# As valued graph
multigraph(arr, valued = TRUE)

# As signed graph
multigraph(arr, signed = TRUE)

# With loops and custumized vertex size
multigraph(arr, cex = 3, loops = TRUE)
}
\keyword{graphics}
