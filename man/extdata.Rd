\name{ucinet}
\docType{data}
\alias{ucinet}
%\alias{system.file}
\alias{davis.dl}
\alias{sampson.dl}
\alias{padgett.dl}
\alias{padgw.dl}
\title{
External Datasets in \code{dl} Format.
}
\description{
External datasets in \pkg{multigraph} package are for illustration purposes, and they have a \pkg{Ucinet} \code{dl} format. 
The four datasets need to be transformed with function \code{\link[multiplex:read.dl]{read.dl}} from \CRANpkg{multiplex} package 
for visualization according the nature of the network.
}
%\usage{
%system.file("extdata", "davis.dl", package="multigraph")
%system.file("extdata", "sampson.dl", package="multigraph")
%system.file("extdata", "padgett.dl", package="multigraph")
%system.file("extdata", "padgw.dl", package="multigraph")
%}
\format{
 \itemize{ 
     \item {\code{davis.dl} Southern woman}: (Davis, Gardner, Gardner, 1941) affiliation network \eqn{18 \times 14}.
     \item {\code{sampson.dl} Monastery novices}: (Sampson, 1969) directed, multiplex, valued, longitudinal network \eqn{18 \times 18 \times 10}.
     \item {\code{padgett.dl} Florentine families}: (Padgett and Ansell, 1993) undirected, multiplex network \eqn{16 \times 16 \times 2}.
     \item {\code{padgw.dl} Florentine families attributes}: table \eqn{16 \times 3}.
    }
}
\source{
\emph{\pkg{UCINET} 6 for windows: Software for social network analysis.}
}
%
\keyword{datasets}
