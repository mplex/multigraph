\name{multigraph-package}
\alias{multigraph-package}
%\alias{multigraph}
\docType{package}
\title{
Plot and Manipulate Multigraphs
}
\description{
Functions to create and network manipulate multigraphs, bipartite graphs, multilevel graphs, and semigroup Cayley colour graphs.
}
\details{
  \tabular{ll}{
     Package: \tab multigraph \cr
     Type: \tab Package \cr
     Version: \tab 1.0 \cr
     Depends: \tab \CRANpkg{multiplex} (>= 4.0) \cr
     Date: \tab 18 August 2026 \cr
     License: \tab GPL-3 \cr
  }
  The \pkg{multigraph} package provides functions for plotting complex network structures, including dichotomous and valued graphs. 
  For one-mode data, \code{multigraph} displays signed and valued multigraphs with loops, whereas \code{bmgraph} represents bipartite graphs for affiliation networks. 
  The \code{mlgraph} function visualizes multilevel graphs that \dQuote{combine} one- and two-mode network data. 
  Cayley colour graphs are drawn with \code{ccgraph} to depict relations among ties in multiplex networks encoded in the associated algebraic semigroup structure.
}
\author{
J. Antonio Rivero Ostoic

Maintainer: Antonio Rivero Ostoic <multiplex@post.com>
}
\references{
Ostoic, J.A.R. \emph{Algebraic Analysis of Social Networks: Models, Methods and Applications Using R}, Wiley, 2021 \cr
Ostoic, J.A.R. \dQuote{Algebraic Analysis of Multiple Social Networks with \code{multiplex}.} \emph{Journal of Statistical Software}, 91(11), 1-41. <doi:10.18637/jss.v092.i11>
}
\seealso{
\code{\link[multiplex:multiplex-package]{multiplex-package}}, \code{\link[multiplex:incubs]{incubs}}, \code{\link[multiplex:read.dl]{read.dl}}, 
\code{\link[multiplex:read.gml]{read.gml}}, \code{\link[multiplex:transf]{transf}}, \code{\link[multiplex:edgel]{edgel}}, \code{\link[multiplex:zbind]{zbind}}, 
\code{\link[multiplex:semigroup]{semigroup}}
}
\keyword{graphs}
\keyword{manip}
\keyword{IO}
\keyword{data}
\keyword{file}
