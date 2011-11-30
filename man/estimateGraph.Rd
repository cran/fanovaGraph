\name{estimateGraph}
\alias{estimateGraph}
\alias{estimateGraphFix}
\alias{estimateGraphRBD}
\alias{estimateGraphSob}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimation of the total interaction indices
}
\description{Estimation of the total interaction indices of a function for the FANOVA graph edges.
}
\usage{
estimateGraph(f.mat, d, q, q.arg, N, method = "fixed", nMC, nfast99 = 500, L, M = 6, Nsobol, ...)
}
\arguments{
  \item{f.mat}{vectorized function of which total interaction indices (weights for the graph edges) shall be estimated
}
  \item{d}{integer, number of input factors (vertices)
}
  \item{q}{a vector of character strings of quantile functions corresponding to the factors distributions, it can be a single character string meaning same distribution for all, if not specified \code{"qunif"} is taken
}
  \item{q.arg}{a list of lists of quantile functions parameters of the distributions in \code{q}, it can be a single list meaning same parameters for all, if not specified the default values of the respective distributions are taken
}
  \item{N}{optional integer, total number of function evaluations, special parameters (\code{nMC}, \code{L} and \code{Nsobol} depending on \code{method}) can be specified instead}
  \item{method}{character string specifying the estimation method of the total interaction indices, to be chosen between \code{"fixed"}, \code{"RBD"} and \code{"Sobol"}, defaults to "fixed"}
  \item{nMC}{optional integer, only if \code{method="fixed"}, number of Monte Carlo simulations for the expectation in fixing method}
  \item{nfast99}{optional integer, only if \code{method="fixed"}, number of design points for FAST algorithm, defaults to 500}
  \item{L}{optional integer, only if \code{method="RBD"}, parameter L in RBD-FAST method}
  \item{M}{optional integer, only if \code{method="RBD"}, parameter M in RBD-FAST method}
  \item{Nsobol}{optional integer, only if \code{method="Sobol"}, number of Monte Carlo simulations in Sobol method}
  \item{...}{additional arguments to be passed to the function \code{f.mat}}
}
\value{a matrix containing the unnormalized total interaction indices}

\author{
J. Fruth, T. Muehlenstaedt
}

\references{
Fruth, J.; Roustant, O.; Kuhnt, S. (2011+) Total interaction indices for the decomposition of functions with high complexity.

Muehlenstaedt, T.; Roustant, O.; Carraro, L.; Kuhnt, S. (2011) Data-driven Kriging models based on FANOVA-decomposition, \emph{Statistics and Computing}.

Sobol', I. M. (1993) Sensitivity estimates for nonlinear mathematical models, \emph{Mathematical Modeling and Computational Experiment}, \bold{1}, 407-414.

Mara, T.A (2009) Extension of the RBD-FAST method to the computation of global sensitivity indices, \emph{Reliability Engineering & System Safety}, \bold{94} no. 8, 1274-1281.
}
\examples{
# Ishigami function, true analytical values: D12 = D23 = 0, D13 =~ 3.374
estimateGraph(f.mat=ishigami.fun, d=3, q.arg=list(min=-pi,max=pi), N=10000, method="fixed")
estimateGraph(f.mat=ishigami.fun, d=3, q.arg=list(min=-pi,max=pi), N=10000, method="RBD")
estimateGraph(f.mat=ishigami.fun, d=3, q.arg=list(min=-pi,max=pi), N=10000, method="Sobol")
}