\name{multigraph-package}
\alias{multigraph-package}
%\alias{multigraph}
\docType{package}
\title{
Plot and Manipulate Multigraphs
}
\description{
Functions to create and manipulate multigraphs, bipartite graphs, and valued multigraphs
}
\details{
\tabular{ll}{
Package: \tab multigraph \cr
Type: \tab Package \cr
Version: \tab 0.93 (devel) \cr
Depends: \tab \CRANpkg{multiplex} (>= 2.9) \cr
Date: \tab 03 January 2021 \cr
License: \tab GPL-3 \cr
}
This package contains functions to plot diverse types of graphs representing complex network structures. 
For one-mode data, it is possible to depict signed and valued multigraphs, and bipartite graphs for two-mode data as well. 
Moreover, multilevel graphs that \emph{combine} one- and two-mode network data are represented with the latest function. 
Finally, Cayley graphs serve to depict relations among the ties with multiple edges. 

Please note that this package still under current development.
}
\author{
J. Antonio Rivero Ostoic

Maintainer: Antonio Rivero Ostoic <multiplex@post.com>
}
\references{
Ostoic, J.A.R. `Algebraic Analysis of Multiple Social Networks with \code{multiplex}.' \emph{Journal of Statistical Software}, 91(11), 1-41. <doi:10.18637/jss.v092.i11>
}

\seealso{
\code{\link[multiplex:multiplex-package]{multiplex-package}}, \code{\link[multiplex:incubs]{incubs}}, \code{\link[multiplex:zbind]{zbind}}, \code{\link[multiplex:transf]{transf}}
}

\keyword{ graphs}
\keyword{ manip }
\keyword{ IO }
\keyword{ data }
\keyword{ file }
