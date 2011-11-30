\name{fanovaGraph-package}
\alias{fanovaGraph}
\docType{package}
\title{Building Kriging models with FANOVA graphs
}
\description{
Estimates and plots the FANOVA graph of a function to identify its interaction structure and fits a kriging model modified by the identified structure
}
\details{
Important functions:
       \tabular{ll}{
         \code{estimateGraph} \tab Estimates the total interaction indices of a function \cr
         \code{threshold} \tab Takes the total interaction indices and creates the graph structure, including a threshold cut for indices lower than a given value \cr
         \code{plotiGraph} \tab Plots a given graph structure \cr
         \code{plotDeltaSteps} \tab Provides plots for the choice of the threshold \cr
         \code{MLoptimConstrained} \tab Kriging model estimation with the kernel modified by the graph \cr
         \code{yhat} \tab Prediction function for the modified model
       }
}
\author{
J. Fruth, T. Muehlenstaedt, O. Roustant
}
\references{
Muehlenstaedt, T.; Roustant, O.; Carraro, L.; Kuhnt, S. (2011) Data-driven Kriging models based on FANOVA-decomposition, \emph{Statistics and Computing}.

Fruth, J.; Roustant, O.; Kuhnt, S. (2011+) Total interaction indices for the decomposition of functions with high complexity.

Sobol', I. M. (1993) Sensitivity estimates for nonlinear mathematical models, \emph{Mathematical Modeling and Computational Experiment}, \bold{1}, 407-414.

Mara, T.A (2009) Extension of the RBD-FAST method to the computation of global sensitivity indices, \emph{Reliability Engineering & System Safety}, \bold{94} no. 8, 1274-1281.
}
\keyword{package}
\seealso{
\code{\link[DiceKriging:DiceKriging-package]{DiceKriging}}, 
\code{\link[sensitivity:sensitivity-package]{sensitivity}},
\code{\link[igraph:igraph-package]{igraph}}
}
\examples{
#demo(ExampleIshigami)
#demo(Example6D)
#demo(Estimation)
#demo(Threshold)
}