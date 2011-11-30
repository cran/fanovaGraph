\name{plotDeltaSteps}
\alias{plotDeltaSteps}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Delta Step Plot
}
\description{Threshold discision plot. \code{plotDeltaSteps} plots the threshold steps (the values of delta at which the graph changes) equidistant against the number of cliques and the values of delta on the real axis.
}
\usage{
plotDeltaSteps(totalInt, d, dall, interval = c(0, 1), meanCliqueSize = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{totalInt}{a matrix of three rows as output of \code{\link{estimateGraph}}, the third row
   contains the unnormalized total interactions, the first two rows indicate the corresponding vertices
}
  \item{d}{integer, number of input factors (vertices)
}
\item{dall}{value of the overall variance
}
  \item{interval}{an optional vector of size 2, interval for the delta steps to be shown, defaults to c(0,1)
}
  \item{meanCliqueSize}{logical, if TRUE (default) an additional line is drawn representing 
        the mean of the number of vertices in the cliques
}
}
\details{
The plots shall give help in the choice for the treshold. In the first plot a small number of cliques might be preferable in order to have less parameters to estimate. If several values result in the same number of cliques the ones with higher mean clique size are possibly preferable. 

In the second plot a sharp increase in the curve indicates a point of big distance between two successive edge indices and thus a clear change in the graph structure. Those points are probably good choices for the tresholdhold.
}

\author{
J. Fruth, O. Roustant
}

\seealso{
\code{\link{plotGraphChange}}
}
\examples{
totalInt <- rbind(rep(1:5,5:1),c(2:6,3:6,4:6,5:6,6),
c(0.061,0.037,0.0011,0.0012,0.0013,0.041,0.002,0.0014,0.0015,0.002,0.002,0.0016,0.019,0.018,0.021))
dall <- 0.8
plotDeltaSteps(totalInt=totalInt, d=6,dall=dall)
plotDeltaSteps(totalInt=totalInt, d=6,dall=dall, meanCliqueSize = TRUE)

totalInt2 <- rbind(c(1,1,2),c(2,3,3),c(0.58,3.42,0.34))
dall2 <- 13.7
plotDeltaSteps(totalInt=totalInt2, d=3, dall=dall2)
}
