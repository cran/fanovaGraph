\name{plotiGraph}
\alias{plotiGraph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot Graph via Package \code{\link[igraph:igraph-package]{igraph}}
}
\description{
More flexible function to plot FANOVA graphs than \code{\link{plotGraph}}, 
that uses functions from package \code{\link[igraph:igraph-package]{igraph}} to adjust the position of the vertices.}
\usage{
plotiGraph(E, d, names=1:d, i1=NULL, tii=NULL, i2=NULL, layout=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{E}{matrix of edges in the graph
}
  \item{d}{integer, number of input factors (vertices)
}
  \item{names}{optional character string, names of vertices, defaults to \code{1:d}
}
  \item{i1}{optional vector of length d, first order indices (thickness of vertices)
}
  \item{tii}{optional vector of length dim(E)[1], total interaction indices (thickness of edges)
}
  \item{i2}{optional vector of length dim(E)[1], second order interaction indices (thickness of inner edges)
}
  \item{layout}{optional layout for the graph as in \code{\link[igraph:igraph-package]{igraph}},
  default is \code{layout.fruchterman.reingold}
}
}
\references{
Muehlenstaedt, T.; Roustant, O.; Carraro, L.; Kuhnt, S. (2011) Data-driven Kriging models based on FANOVA-decomposition, \emph{Statistics and Computing}.

Csardi, G.; Nepusz, T. (2006) The igraph software package for complex network research, \emph{InterJournal Complex Systems}, \bold{Complex Systems}, 1695.
}
\author{J. Fruth, O. Roustant, S. Hess
}

\seealso{
\code{\link{plotGraph}}, \code{\link{plotGraphChange}}
}
\examples{
E <- matrix(c(1,1,2,1,2,3,3,4),ncol=2)
#plotiGraph(E, 4)
#plotiGraph(E, 4, names=c("A","B","C","D"))
#plotiGraph(E, 4, names=c("A","B","C","D"), i1=1:4/4, tii=1:4/8)

E <- matrix(c(1,1,2,2,2,3,4,3,5,4,5,6,5,6),7,2)
vertex.weight <- c(0.003,0.2,0.06,0.1,0.02,0.05)
edge.weight <- c(0.12,0.05,0.02,0.11, 0.04,0.08,0.03)
#plotiGraph(E, 6, i1=vertex.weight, tii= edge.weight, i2=edge.weight/2)
}