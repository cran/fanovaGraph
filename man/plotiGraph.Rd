\name{plotiGraph}
\alias{plotiGraph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot Graph via Package \code{\link[igraph0:igraph-package]{igraph0}}
}
\description{
Plot FANOVA graphs using functions from package \code{\link[igraph0:igraph-package]{igraph0}}.}
\usage{
plotiGraph(graphlist, names = NULL, i2 = NULL, layout = NULL, ploti1 = TRUE, maxThickness=20, circleDiameter=40)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{graphlist}{a list representing the graph structure as obtained from \code{\link{estimateGraph}}
}
  \item{names}{optional character string, names of vertices, defaults to \code{1:d}
}
    \item{i2}{optional vector of second order interaction indices (thickness of inner edges)
}
    \item{ploti1}{optional boolean, if TRUE main effects are drawn in the graph by vertices thicknesses, should be set to FALSE when only total interaction indices are of interest
}
  \item{layout}{optional layout for the graph as in \code{\link[igraph0:igraph-package]{igraph0}},
  default is \code{layout.fruchterman.reingold}
}
  \item{maxThickness}{optional value for the maximal line thickness, defaults to 20
}
  \item{circleDiameter}{optional value for the circle diameter, defaults to 40
}
}
\references{
Muehlenstaedt, T.; Roustant, O.; Carraro, L.; Kuhnt, S. (2011) Data-driven Kriging models based on FANOVA-decomposition, \emph{Statistics and Computing}.

Csardi, G.; Nepusz, T. (2006) The igraph software package for complex network research, \emph{InterJournal Complex Systems}, \bold{Complex Systems}, 1695.
}
\author{J. Fruth, O. Roustant, S. Hess, S. Neumaerker
}

\seealso{
\code{\link{plotGraphChange}}
}
\examples{
op <- par(no.readonly=TRUE)
g1 <- estimateGraph(f.mat=ishigami.fun, d=3, q.arg=list(min=-pi,max=pi), N=10000)
plotiGraph(g1)
plotiGraph(g1, names=c("A","B","C"))
plotiGraph(g1, names=c("A","B","C"), ploti1 = FALSE)

# include pure second order indices
g2 <- estimateGraph(f.mat=function(x) x[,1]*x[,2]*x[,3]+x[,2]*x[,3], d=3, q.arg=list(min=-1,max=1), N=10000)

plotiGraph(g2)
plotiGraph(g2, i2 = c(0.001, 0.001, 0.05))

# equal layouts and different edge thicknesses and circle diameters
g3 <- estimateGraph(f.mat=function(x) x[,1]*x[,2]*x[,3]*x[,4]*x[,5], d=5, q.arg=list(min=-1,max=1), N=10000)

g4 <- estimateGraph(f.mat=function(x) x[,1]*x[,2]*x[,3]+x[,4]*x[,5], d=5, q.arg=list(min=-1,max=1), N=10000)

graphClassIgraph <- graph.full(n = 5, directed = FALSE)
layout <- layout.circle(graphClassIgraph)

plotiGraph(g3, maxThickness= 10, circleDiameter= 30, layout=layout)
plotiGraph(g4, maxThickness= 30, circleDiameter= 50, layout=layout)
}