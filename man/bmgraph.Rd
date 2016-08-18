\name{bmgraph}
\alias{bmgraph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create bipartite multigraph
}
\description{
A function to create and manipulate bipartite multigraphs with different layouts
}
\usage{
bmgraph(net, layout = c("bip", "bip3", "bip3e", "bip4", "stress", "rand", "circ"), coord = NULL, tcex = NULL, alpha = c(1, 1, 1), showLbs = TRUE, showAtts = TRUE, att = NULL, lbat = "1", main = NULL, cex.main, bg, mar, cex, pos, lwd, lty, ecol, vcol, vcol0, asp, directed, collRecip, seed = NULL, maxiter, bwd, clu, pch, tcol, hds, vedist, rot, mirrorX, mirrorY, col, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{net}{
%%     ~~Describe \code{net} here~~
}
  \item{layout}{
%%     ~~Describe \code{layout} here~~
}
  \item{coord}{
%%     ~~Describe \code{coord} here~~
}
  \item{tcex}{
%%     ~~Describe \code{tcex} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
  \item{showLbs}{
%%     ~~Describe \code{showLbs} here~~
}
  \item{showAtts}{
%%     ~~Describe \code{showAtts} here~~
}
  \item{att}{
%%     ~~Describe \code{att} here~~
}
  \item{lbat}{
%%     ~~Describe \code{lbat} here~~
}
  \item{main}{
%%     ~~Describe \code{main} here~~
}
  \item{cex.main}{
%%     ~~Describe \code{cex.main} here~~
}
  \item{bg}{
%%     ~~Describe \code{bg} here~~
}
  \item{mar}{
%%     ~~Describe \code{mar} here~~
}
  \item{cex}{
%%     ~~Describe \code{cex} here~~
}
  \item{pos}{
%%     ~~Describe \code{pos} here~~
}
  \item{lwd}{
%%     ~~Describe \code{lwd} here~~
}
  \item{lty}{
%%     ~~Describe \code{lty} here~~
}
  \item{ecol}{
%%     ~~Describe \code{ecol} here~~
}
  \item{vcol}{
%%     ~~Describe \code{vcol} here~~
}
  \item{vcol0}{
%%     ~~Describe \code{vcol0} here~~
}
  \item{asp}{
%%     ~~Describe \code{asp} here~~
}
  \item{directed}{
%%     ~~Describe \code{directed} here~~
}
  \item{collRecip}{
%%     ~~Describe \code{collRecip} here~~
}
  \item{seed}{
%%     ~~Describe \code{seed} here~~
}
  \item{maxiter}{
%%     ~~Describe \code{maxiter} here~~
}
  \item{bwd}{
%%     ~~Describe \code{bwd} here~~
}
  \item{clu}{
%%     ~~Describe \code{clu} here~~
}
  \item{pch}{
%%     ~~Describe \code{pch} here~~
}
  \item{tcol}{
%%     ~~Describe \code{tcol} here~~
}
  \item{hds}{
%%     ~~Describe \code{hds} here~~
}
  \item{vedist}{
%%     ~~Describe \code{vedist} here~~
}
  \item{rot}{
%%     ~~Describe \code{rot} here~~
}
  \item{mirrorX}{
%%     ~~Describe \code{mirrorX} here~~
}
  \item{mirrorY}{
%%     ~~Describe \code{mirrorY} here~~
}
  \item{\dots}{
Additional items
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A bipartite multigraph
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
%\examples{
%
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ graphics }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
