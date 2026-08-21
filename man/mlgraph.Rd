\name{mlgraph}
\alias{mlgraph}
\title{
Multilevel Graph
}
\description{
A function to create and manipulate multilevel graphs.
}
\usage{
mlgraph(net, layout = c("circ", "force", "stress", "rand", "conc", "bip"), main = NULL, 
        seed = NULL, maxiter = 100, directed = TRUE, alpha = c(1, 1, 1), scope, collRecip, 
        lbs, showLbs, showAtts, cex.main, col.main, font.main, coord, clu, cex, lwd, pch, 
        lty, bwd, bwd2, att, bg, mar, pos, asp, ecol, vcol, vcol0, col, lbat, loops, rot, 
        mirrorX, mirrorY, mirrorD, mirrorL, mirrorV, mirrorH, hds, scl, vedist, vldist, 
        flip, ffamily, fstyle, fsize, fcol, valued, modes, elv, lng, nr, new, jitter, ...)
}
\arguments{
  \item{net}{A \dQuote{\code{Multilevel}} class object or a 3D array with clustering information.}
  \item{layout}{
    Visualization layout: 
      \itemize{ 
       \item \code{circ} Circular
       \item \code{force} Force-directed
       \item \code{stress} Stress-majorization
       \item \code{rand} Random
       \item \code{conc} Concentric
       \item \code{bip} Bipartite graph
      }
}
  \item{main}{(optional) Main title of the plot graph.}
  \item{seed}{(optional, for \code{force}, \code{stress} and \code{rand}) Random seed number for initial coordinates of vertices.}
  \item{maxiter}{(optional when \code{seed}) Maximum number of iterations.}
  \item{directed}{(logical) Are graph edges in domain directed?}
  \item{alpha}{Vector with the alpha colour transparencies for (vertex, edge, \code{bg}).}
  \item{scope}{(optional) Scope of graph, edges, and vertices.}
  \item{collRecip}{(optional, logical and experimental) Collapse reciprocated edges in graph domain?}
  \item{lbs}{(optional) List of customized vertex labels for domain and codomain.}
  \item{showLbs}{(optional and logical) Display vertex labels in graph?}
  \item{showAtts}{(optional and logical) Display vertex attribute labels?}
  \item{cex.main}{Font size of title in plot graph.}
  \item{col.main}{Font colour of title in plot graph.}
  \item{font.main}{Font family of title in plot graph.}
  \item{coord}{(optional) Data frame with graph coordinates for vertices after binomial projection.}
  \item{clu}{(optional) Clustering of vertices as a list of vectors with integers or \code{NULL} (see \emph{details}).}
  \item{cex}{(optional) Vertex size.}
  \item{lwd}{(optional, ignored for \code{valued}) Width of edges.}
  \item{pch}{(optional) Symbols for vertices in domain and in codomain.}
  \item{lty}{(optional) Shapes of edges in domain and in codomain.}
  \item{bwd}{(optional, default \code{1}) Real number for width between bundle edges.}
  \item{bwd2}{(optional) Width of bundle loop edges.}
  \item{att}{(optional) Vector or array representing vertex attributes.}
  \item{bg}{(optional) Background colour of plot graph.}
  \item{mar}{(optional) Vector with margins for plot.}
  \item{pos}{(optional) Positions of labels from vertex:
   \itemize{ 
       \item {\code{1}}: Below
       \item {\code{2}}: Left
       \item {\code{3}}: Above
       \item {\code{4}}: Right
       \item {\code{0}}: Centre
      }}
  \item{asp}{(optional, \code{asp=NA} cover whole canvas) Aspect ratio of plot.}
  \item{ecol}{(optional) Colour of edges by type in domain and codomain.}
  \item{vcol}{(optional) Colour of vertices.}
  \item{vcol0}{(optional, only for \code{pch=21} through \code{pch=25}) Colour of vertices contour.}
  \item{col}{(optional) Alias for \code{vcol}.}
  \item{lbat}{(optional) Labels for vertex attributes.}
  \item{loops}{(optional and logical) Plot loops in graph?}
  \item{rot}{(optional) Clockwise rotation in degrees on graph.}
  \item{mirrorX}{(optional) Mirror by \eqn{X} axis.}
  \item{mirrorY}{(optional) Mirror by \eqn{Y} axis.}
  \item{mirrorD}{(optional) Mirror reflection across diagonal \eqn{Y=X}.}
  \item{mirrorL}{(optional) Mirror reflection across diagonal \eqn{Y=-X}.}
  \item{mirrorV}{\dQuote{Same as \code{mirrorX}}.}
  \item{mirrorH}{\dQuote{Same as \code{mirrorY}}.}
  \item{hds}{(optional, for \code{directed}) Head scale of arcs.}
  \item{scl}{(optional) Numerical scalar (\eqn{x} and \eqn{y}) or vector (\eqn{x}, \eqn{y}) of graph scale.}
  \item{vedist}{(optional) Vertex--edge distance.}
  \item{vldist}{(optional) Vertex--loop distance.}
  \item{flip}{(optional and logical) Flip loops?}
  \item{ffamily}{(optional) Font family of vertex labels in graph.}
  \item{fstyle}{(optional) Font style of vertex labels in graph with options: \code{bold}, \code{italic}, \code{bolditalic}.}
  \item{fsize}{(optional) Font size of vertex labels in graph.}
  \item{fcol}{(optional) Font colour of text in graph vertices.}
  \item{valued}{(optional and logical) Depict a valued graph?}
  \item{modes}{(optional, for \dQuote{\code{Multilevel}} class) Vector of domains and codomains in \code{net} matrices.}
  \item{elv}{(experimental for \code{"bip"}) Elevation of intra-domain edges.}
  \item{lng}{(experimental for \code{"bip"}) Elevation control intra-domain edges.}
  \item{nr}{(optional, for \code{conc}) Integer or \code{NULL} with number of radii layout.}
  \item{new}{(optional, logical) Plot a new graph on existing canvas?}
  \item{jitter}{(optional, for \code{stress}) Amount of jitter to add on coordinates.}
  \item{\dots}{Additional argument items.}
}
\details{
Multilevel graphs represent networks with multiple \dQuote{levels}, such as distinct domains within a network structure. 
These networks may include ties both within and between domains.

Because this function accepts many arguments, they may be collected in a list and supplied through the \code{scope} option. 
The list can contain vectors, scalars, or combinations of both for rendering.

The bundle widths \code{bwd} and \code{bwd2} range from \code{0}, where edges are collapsed, to the default value of \code{1}, depending on vertex size. 
Values greater than one are allowed when \code{valued} option is used.

In multilevel structures, argument \code{clu} serves to assign network members to classes. 
Setting the corresponding vector to \code{NULL} places all members of a domain or codomain in a single class. 
Likewise, setting argument \code{nr} to \code{NULL} in the \code{conc} layout uses two radii, one for each domain.
}
\value{
A plot of the multilevel graph structure for the network.
}
%\references{} ~
%\author{} ~
%\note{} ~
%\section{Warning }{....} ~
\seealso{
\code{\link[multiplex:mlvl]{mlvl}}, \code{\link{multigraph}}, \code{\link{bmgraph}}, \code{\link{frcd}}, \code{\link{stsm}}, \code{\link{conc}}
}
%\examples{
%}
\keyword{graphics}
