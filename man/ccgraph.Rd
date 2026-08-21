\name{ccgraph}
\alias{ccgraph}
\title{
Cayley Colour Graph
}
\description{
A function to create and manipulate Cayley colour graphs of semigroup structures.
}
\usage{
ccgraph(x, main = NULL, seed = 0, maxiter = 100, alpha = c(1, 1, 1), 
        scope, loops, collRecip, undRecip, showLbs, cex.main, col.main, font.main, 
        conc, coord, clu, cex, lwd, pch, lty, bwd, bwd2, att, bg, mar, pos, asp, 
        ecol, vcol, vcol0, lbs, col, lbat, swp, swp2, rot, vedist, vldist, 
        mirrorX, mirrorY, mirrorD, mirrorL, mirrorV, mirrorH, ffamily, fstyle, 
        fsize, fcol, scl, hds, nr, gens, flip, ...)
}
\arguments{
  \item{x}{An algebraic structure like a \sQuote{\code{Semigroup}} class object (cf. \code{\link[multiplex:semigroup]{semigroup}}).}
  \item{main}{(optional) Main title of the plot graph.}
  \item{seed}{(optional) Random seed number for initial coordinates of vertices.}
  \item{maxiter}{(optional) Maximum number of iterations in layout algorithm.}
  \item{alpha}{Vector with the alpha colour transparencies for (vertex, edge, \code{bg}).}
  \item{scope}{(optional) List with scope of graph, edge and vertex characteristics.}
  \item{loops}{(optional and logical) Plot loops in graph?}
  \item{collRecip}{(optional and logical) Collapse reciprocated edges in graph?}
  \item{undRecip}{(optional and logical) Plot reciprocated edges as undirected?}
  \item{showLbs}{(optional and logical) Show vertex labels in graph when available?}
  \item{cex.main}{Font size of title in plot graph.}
  \item{col.main}{Font colour of title in plot graph.}
  \item{font.main}{Font family of title in plot graph.}
  \item{conc}{(optional and logical) Is layout concentric?}
  \item{coord}{(optional, ignored if \code{"conc"}) Two-column data frame of vertex coordinates in graph.}
  \item{clu}{(optional) Clustering of vertices.}
  \item{cex}{(optional) Vertex size.}
  \item{lwd}{(optional) Width of edges.}
  \item{pch}{(optional) Symbol for vertices.}
  \item{lty}{(optional) Shape of edges.}
  \item{bwd}{(optional, default \code{1}) Real number for width between bundle edges.}
  \item{bwd2}{(optional, experimental) width of the bundle loop edges.}
  \item{att}{(optional) Vector or array representing vertex attributes.}
  \item{bg}{(optional) Background colour of plot graph.}
  \item{mar}{(optional) Vector with margins for plot.}
  \item{pos}{(optional) Position of labels from vertex:
   \itemize{ 
       \item {\code{1}}: Below
       \item {\code{2}}: Left
       \item {\code{3}}: Above
       \item {\code{4}}: Right
       \item {\code{0}}: Centre
      }}
  \item{asp}{(optional, \code{asp=NA} cover whole canvas) Aspect ratio of plot.}
  \item{ecol}{(optional) Colour of edges.}
  \item{vcol}{(optional) Colour of vertices.}
  \item{vcol0}{(optional) Colour of vertices contour.}
  \item{lbs}{(optional) Vertex labels.}
  \item{col}{(optional) Alias for \code{vcol}.}
  \item{lbat}{(optional) Labels for vertex attributes.}
  \item{swp}{(optional and logical) Swap bundle patterns?}
  \item{swp2}{(optional and logical) Swap reciprocals?}
  \item{rot}{(optional) Clockwise rotation in degrees on graph.}
  \item{vedist}{(optional) Vertex--edge distance.}
  \item{vldist}{(optional) Vertex--loop distance.}
  \item{mirrorX}{(optional) Mirror by \eqn{X} axis.}
  \item{mirrorY}{(optional) Mirror by \eqn{Y} axis.}
  \item{mirrorD}{(optional) Mirror reflection across diagonal \eqn{Y=X}.}
  \item{mirrorL}{(optional) Mirror reflection across diagonal \eqn{Y=-X}.}
  \item{mirrorV}{\dQuote{Same as \code{mirrorX}}.}
  \item{mirrorH}{\dQuote{Same as \code{mirrorY}}.}
  \item{ffamily}{(optional) Font family of vertex labels in graph.}
  \item{fstyle}{(optional) Font style of vertex labels in graph with options: \code{bold}, \code{italic}, \code{bolditalic}.}
  \item{fsize}{(optional) Font size of vertex labels in graph.}
  \item{fcol}{(optional) Font colour of text in graph vertices.}
  \item{scl}{(optional) Numerical scalar (\eqn{x} and \eqn{y}) or vector (\eqn{x}, \eqn{y}) of graph scale.}
  \item{hds}{(optional) Head scale for arcs.}
  \item{nr}{(\code{conc} layout) Number of radii.}
  \item{gens}{(optional when absent) Semigroup generators in \code{x}.}
  \item{flip}{(optional and logical) Flip loops per \eqn{D}-class?}
  \item{\dots}{Additional argument items.}
}
\details{
The Cayley colour graph is a graphical representation of the relationships among relations in the relational structure of a given multiplex network. 
Both nodes and directed edges represent string relations, and each shape and colour corresponds to a specific generator relation of the semigroup structure.
}
\value{
A plot of the semigroup or group structure.
}
%\references{} ~
%\author{} ~
%\note{} ~
%\section{Warning }{....} ~
\seealso{
\code{\link[multiplex:semigroup]{semigroup}}, \code{\link{multigraph}}, \code{\link{frcd}}, \code{\link{conc}}
}
\examples{
# Create an abstract semigroup from random data
arr <- round( replace( array(runif(18), c(3,3,2)), array(runif(18),
       c(3,3,2))>.5, 1 ) )

S <- semigroup(arr)

# plot semigroup Cayley graph
ccgraph(S)
}
\keyword{graphics}
\keyword{math}
