\name{plotGraphChange, plotTk, plotManipulate}
\alias{plotGraphChange}
\alias{plotTk}
\alias{plotManipulate}
\title{Plot Graph as It Changes with Delta
}
\description{Graphs are plottet depending on a change on delta, the threshold for edges to appear in the graph, to enable a visual decision for delta by graph behavior.
}
\usage{
plotGraphChange(totalInt, d, dall, fix.layout = TRUE, delta.layout = 0.01)
plotTk(totalInt, d, dall, delta.layout=0.01)
plotManipulate(totalInt, d, dall, delta.layout=0.01)
}
\arguments{
  \item{totalInt}{a matrix of three rows as output of \code{\link{estimateGraph}}, the third  row contains the unnormalized total interactions, the first two rows indicate the corresponding vertices
}
  \item{d}{integer, number of input factors (vertices)
}
\item{dall}{value of the overall variance
}
  \item{fix.layout}{logical, if TRUE (default) the position of the vertices is fixed for all plots such that the positions are optimal for delta = delta.layout
}
  \item{delta.layout}{optional value between 0 and 1, see \code{fix.layout}, defaults to 0.01
}
}

\author{J. Fruth, O. Roustant
}
\note{\code{plotGraphChange} shows the changing of the graph step by step by changing plots as in \code{demo}, \code{plotTk} is an interactive version using \code{\link[tcltk:tcltk-package]{tcltk}}, \code{plotManipulate} is an interactive version using \code{manipulate}
}

\seealso{
\code{\link{plotDeltaSteps}}, \code{\link{plotiGraph}}
}
\examples{
### example direct for ishigami instead of metamodel
soverall1 <- fast99(model = ishigami.fun, factors = 3, 
   n = 2000, q = "qunif", q.arg = list(min = -pi, max = pi))
totalInt1 <- estimateGraph(f.mat=ishigami.fun, d=3, nMC=20, nfast99=500, 
q = "qunif", q.arg = list(min = -pi, max = pi), method="fixed")

#plotGraphChange(totalInt=totalInt1, d=3, dall=mean(soverall1$V), fix.layout=TRUE, delta.layout=0.1)

#library(tcltk)
#plotTk(totalInt=totalInt1, d=3, dall=mean(soverall1$V), delta.layout=0.1)

#plotManipulate(totalInt=totalInt1, d=3, dall=mean(soverall1$V), delta.layout=0.1)
}