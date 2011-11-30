\name{threshold}
\alias{threshold}

\title{Treshold Step of Graph Creation
}
\description{Takes total interaction indices, Sobol indices and delta and creates graph with cliques by thresholding on delta.
}
\usage{
threshold(delta, totalInt, d, dall)
}

\arguments{
  \item{delta}{numeric threshold between 0 and 1, edges with index > delta are kept in the graph
}
  \item{totalInt}{a matrix of three rows as output of \code{\link{estimateGraph}}, the third  row contains the unnormalized total interactions, the first two rows indicate the corresponding vertices
}
  \item{d}{integer, number of input factors (vertices)
}
\item{dall}{value of the overall variance
}
}
\value{
list with components
\item{CL}{list of cliques}
\item{E}{matrix with two columns containing the edges of the graph}
\item{tii.scaled}{vector containing the total indices of the edges of the graph scaled by the overall variance \code{dall}}
}
\references{
Muehlenstaedt, T.; Roustant, O.; Carraro, L.; Kuhnt, S. (2011) Data-driven Kriging models based on FANOVA-decomposition, \emph{Statistics and Computing}.
}
\author{
T. Muehlenstaedt, O. Roustant, J. Fruth
}
\examples{
f <- function(x) x[,1]*x[,3]*x[,4]
soverall <- fast99(model = f, factors = 4, 
   n = 1000, q = "qunif", q.arg = list(min = 0, max = 1))
totalInt <- estimateGraph(f=f, d=4, nMC=20, nfast99=500)
threshold(0.01, totalInt=totalInt, d=4, dall=mean(soverall$V))
}
