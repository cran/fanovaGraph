\name{plotGraph}
\alias{plotGraph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot Graph
}
\description{A graph is plotted, if weights for vertices and edges are provided, 
  they are represented through the thickness of lines and circles respectively, vertices can be named.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotGraph(E, d, names = 1:d, Maineffects = 1, d.bar = 1, MAIN = "Estimated Graph")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{E}{an integer valued matrix of two columns indicating the edges of the graph
}
  \item{d}{an integer indicating the number of vertices
}
  \item{names}{optional, character vector of length d indicating the names of the vertices
}
  \item{Maineffects}{optional, an vector of length d indicating weights for the vertices (here sobol indices)
}
  \item{d.bar}{optional, an vector of length dim(E)[1] indicating weights for the edges (here total interaction indices)
}
  \item{MAIN}{optional, title of the plot
}
}
\references{
Muehlenstaedt, T.; Roustant, O.; Carraro, L.; Kuhnt, S. (2011) Data-driven Kriging models based on FANOVA-decomposition, \emph{Statistics and Computing}.
}
\author{T. Muehlenstaedt, J. Fruth, O. Roustant
}
\note{
less flexible then \code{\link{plotiGraph}}, which is based on package \code{\link[igraph:igraph-package]{igraph}}
}

\seealso{
\code{\link{plotiGraph}}
}
\examples{
E <- matrix(c(1,1,2,1,2,3,3,4),ncol=2)
plotGraph(E, 4)
plotGraph(E,4, names=c("A","B","C","D"))
plotGraph(E,4, names=c("A","B","C","D"), Maineffects=1:4/4, d.bar=1:4/8)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
