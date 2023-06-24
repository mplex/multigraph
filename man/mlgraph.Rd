\name{mlgraph}
\alias{mlgraph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multilevel graph
}
\description{
A function to create and manipulate multilevel graphs
}
\usage{
mlgraph(net, layout = c("circ", "force", "stress", "rand", "conc", "bip"), main = NULL, 
        seed = NULL, maxiter = 100, directed = TRUE, alpha = c(1, 1, 1), scope, collRecip, 
        undRecip, showLbs, showAtts, cex.main, coord, clu, cex, lwd, pch, lty, bwd, bwd2, 
        att, bg, mar, pos, asp, ecol, vcol, vcol0, col, lbat, swp, loops, swp2, mirrorX, 
        mirrorY, mirrorD, mirrorL, lbs, mirrorV, mirrorH, rot, hds, scl, vedist, ffamily, 
        fstyle, fsize, fcol, valued, modes, elv, lng, nr, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{net}{
  a \code{"Multilevel"} class object or a three dimensional array with clustering information
}
  \item{layout}{
    the visualization layout: 

       \code{circ} circular

       \code{force} force-directed

       \code{stress} stress-majorization

       \code{rand} random 

       \code{conc} concentric

       \code{bip} as bipartite graph
}
  \item{main}{
	(optional) title of the plot
}
  \item{seed}{
	(optional) random seed number for the vertices' initial coordinates. Ignored except for \code{force}, \code{stress} and \code{rand}
}
  \item{maxiter}{
	(optional) maximum number of iterations in layout algorithms. Ignored except for \code{force}, \code{stress} and \code{rand}
}
  \item{directed}{
	(logical) whether or not the graph is directed or undirected
}
  \item{alpha}{
	vector (vertex, edge, \code{bg}) with the alpha color transparency
}
  \item{scope}{
   (optional) scope of the graph (see details)
}
  \item{collRecip}{
	(optional and logical) whether or not collapse reciprocated edges in the undirected graph
}
  \item{undRecip}{
	(optional and logical) whether or not plot reciprocated edges as undirected
}
  \item{showLbs}{
	(optional and logical) whether or not to show the vertex labels
}
  \item{showAtts}{
	(optional and logical) whether or not to show the vertex attribute labels
}
  \item{cex.main}{
	(optional) size of the plot's title
}
  \item{coord}{
	(optional) data frame with the coordinates of the vertices. If coordinates are given then the \code{layout} option is ignored
}
  \item{clu}{
	(optional) clustering of the vertices as a list of vectors with integers or NULL (see \emph{details})
}
  \item{cex}{
	(optional) size of the vertices
}
  \item{lwd}{
	(optional) width of the edges; ignored if \code{valued} is set to \code{TRUE}
}
  \item{pch}{
	(optional) symbol representing the vertices
}
  \item{lty}{
	(optional) shape of the edges
}
  \item{bwd}{
	(optional) width of the bundle edges. Ranges from \code{0} (edges collapsed) to the default \code{1} (depending on the vertices' size), and 
	for \code{valued} a value greater than one is possible
}
  \item{bwd2}{
	(optional) width of the bundle loop edges. 
}
  \item{att}{
	(optional) a vector or an array representing the vertex attributes
}
  \item{bg}{
	(optional) background color of the plot
}
  \item{mar}{
	(optional) margins of the plot
}
  \item{pos}{
	(optional) position of the vertices' labels (\code{0} means ``at the center of the vertex'')
}
  \item{asp}{
	(optional) aspect ratio of the plot
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
  \item{col}{
	(optional) alias for \code{vcol}
}
  \item{lbat}{
	(optional) labels for the vertex attributes
}
  \item{swp}{
	(optional and logical) whether or not to swap the bundle patterns
}
  \item{loops}{
	(optional, logical, and experimental) plot graph loops?
}
  \item{swp2}{
	(optional and logical) whether or not to swap reciprocals
}
  \item{mirrorX}{
	(optional) mirror of the \eqn{X} axis
}
  \item{mirrorY}{
	(optional) mirror of the \eqn{Y} axis
}
  \item{mirrorD}{
	(optional) mirror reflection across diagonal \eqn{Y=X}
}
  \item{mirrorL}{
	(optional) mirror reflection across diagonal \eqn{Y=-X}
}
  \item{lbs}{
	(optional) vertex labels
}
  \item{mirrorV}{
	same as \code{mirrorX}
}
  \item{mirrorH}{
	same as \code{mirrorY}
}
  \item{rot}{
	(optional) clockwise rotation of the graph in degrees 
}
  \item{hds}{
	(optional and experimental) arcs' head scale
}
  \item{scl}{
	(optional and experimental) numerical scalar (\eqn{x} and \eqn{y}) or vector (\eqn{x}, \eqn{y}) of the graph's scale 
}
  \item{vedist}{
	(optional and experimental) a real number with vertex - edge distance
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
  \item{fcol}{
	the font color
}
  \item{valued}{
	(optional and logical) whether the graph is depicyed as valued or not
}
  \item{modes}{
	(optional) a vector indicating which matrices are domains and which codomains (works only with a \code{"Multilevel"} class object)
}
  \item{elv}{
	(experimental) control loops 1
}
  \item{lng}{
	(experimental) control loops 2
}
  \item{nr}{
	integer or \code{NULL} with the number of radii for \code{conc} layout (see \emph{details})
}
  \item{\dots}{
	Additional argument items (see e.g. \code{\link[graphics:par]{par}})
}
}
\details{
Multilevel graphs serve to represent networks with different ``levels'' such as different domains in the network structure. 
A characteristic of multilevel networks is the existence of ties within and across domains.

Since this function can handle a large number of arguments, these can be stored as a list object that is passed through the \code{scope} option. 
In this case, a vector made of lists and scalars or combinations of these is accepted.

The bundle width specified by \code{bwd} and \code{bwd2} ranges from \code{0} (edges collapsed) to the default \code{1} (depending on the vertices' size). 
For the \code{valued} option, a number greater than one is possible.

In a multilevel structure, argument \code{clu} is to class network members 
and it is possible to class all members of the domain or co-domain into a single class by setting the vector to \code{NULL}. 
Similarly, \code{NULL} in argument \code{nr} for the \code{conc} layout implies the use of two radii, one for each domain.
}
\value{
A plot of the multilevel graph structure for the network
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Antonio Rivero Ostoic
}
\note{
Multilevel graphs depend on multilevel class objects
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[multiplex:mlvl]{mlvl}}, \code{\link{multigraph}}, \code{\link{bmgraph}}, \code{\link{frcd}}, \code{\link{stsm}}, \code{\link{conc}}
}
\examples{
## Not run:
# create network data as arrays
arr <- round( replace( array(runif(18), c(3,3,2)), array(runif(18),
       c(3,3,2))>.5, 3 ) )
arr2 <- round( replace( array(runif(18), c(3,3,2)), array(runif(18),
       c(3,3,2))>.5, 3 ) )

# create multilevel class object and plot multilevel graph
require(multiplex)
mlvl(arr, arr2) |> 
  mlgraph()
## End(**Not run**)
}
\keyword{ graphics }
