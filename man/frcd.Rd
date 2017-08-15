\name{frcd}
\alias{frcd}
\title{
Force directed layout
}
\description{
A function to compute the graph coordinated system with a force directed layout algortihm
}
\usage{
frcd(net, seed = seed, maxiter, drp, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{net}{
	an array representing the network relations
}
  \item{seed}{
        (mandatory) the seed of the initial layout (see details)
}
  \item{maxiter}{
        (optional) the maximum number of iterations 
}
  \item{drp}{
	(optional) for \code{weighted} networks, drop values less than specified
}
  \item{\dots}{
	Additional argument items 
}
}
\details{
This function was meant as an internal routine for graph visaulization; however, it can be used with the \code{coord} option both in \code{multigraph} and \code{bmgraph} where a random \code{seed} is stated by \code{NULL}
}
\value{
A data frame with the coordinated system with two columns representing the abscissa and the ordinate in a two-dimensional rectangular Cartesian coordinate system.
}
\references{
Fruchterman, T.M.J., & Reingold, E.M. Graph drawing by force-directed placement. \emph{Software-Practice & Experience}, 21(11), 1129-1164. 1991.
}
\author{
Antonio Rivero Ostoic
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{multigraph}}, \code{\link{bmgraph}}, \code{\link{stsm}}, \code{\link{conc}}
}
%\examples{

%}
\keyword{ graphics }
\keyword{ manip }
