\name{modify.igraph}
\alias{modify.igraph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Modify \code{igraph}
}
\description{
Internal function for \code{plotiGraph} that modifies functions from the package \code{igraph}, so that plotting of different vertex circles is possible.
}
\usage{
modify.igraph()
}
\value{
no values}

\author{
G. Csardi, T. Nepusz, J. Fruth
}

\references{
Csardi, G. and Nepusz, T. (2006) The igraph software package for complex network research. \emph{InterJournal Complex Systems}, \bold{Complex Systems}, 1695.
}

\examples{
fanovaGraph:::modify.igraph()
g <- graph.ring(4)
#plot(g, vertex.shape="circle2", vertex.size=40,vertex.frame.width=c(1,3,6,9))
}
