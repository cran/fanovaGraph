\name{estimateGraph}
\alias{estimateGraph}
\alias{estimateGraphFixLO}
\alias{estimateGraphFixFast}
\alias{estimateGraphRBD}
\alias{estimateGraphSob}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{FANOVA graph estimation.
}
\description{Estimates the structure of the FANOVA graph by estimating the total interaction indices for the graph edges (a particular case of superset importance introduced by Liu and Owen, 2006), the main effect indices for the graph vertices and the overall variance for normalizing the indices (both by function \code{\link{fast99}} from package \code{\link[sensitivity:sensitivity-package]{sensitivity}}) and finding the clique structure of the estimated graph.
}
\usage{
estimateGraph(f.mat, d, q = NULL, q.arg = NULL, N = NULL, method = "FixLO", nLO = NULL, nMC = NULL, nfast99 = 500, L = NULL, M = 6, Nsobol = NULL, confInt = TRUE, ...)
}
\arguments{
  \item{f.mat}{vectorized function for which the FANOVA graph shall be estimated
}
  \item{d}{integer, number of input factors (vertices)
}
  \item{q}{a vector of character strings of quantile functions corresponding to the factors distributions, it can be a single character string meaning same distribution for all, if not specified \code{"qunif"} is taken
}
  \item{q.arg}{a list of lists of quantile functions parameters of the distributions in \code{q}, it can be a single list meaning same parameters for all, if not specified the default values of the respective distributions are taken
}
  \item{N}{optional integer, total number of function evaluations, instead of \code{N} method related parameters (\code{nLO}, \code{nMC}, \code{L} or \code{Nsobol}) can be provided}
  \item{method}{character string specifying the estimation method of the total interaction indices, to be chosen between \code{"FixLO"}, \code{"FixFast"}, \code{"RBD"} and \code{"Sobol"}, defaults to "FixLO", see references for further details}
  \item{nLO}{optional integer, only if \code{method="FixLO"}, number of Monte Carlo simulations in fixing method of Liu and Owen}
  \item{nMC}{optional integer, only if \code{method="FixFast"}, number of Monte Carlo simulations for the expectation in fixing method using FAST}
  \item{nfast99}{optional integer, only if \code{method="fixed"}, number of design points for FAST algorithm, defaults to 500}
  \item{L}{optional integer, only if \code{method="RBD"}, parameter L in RBD-FAST method}
  \item{M}{optional integer, only if \code{method="RBD"}, parameter M in RBD-FAST method}
  \item{Nsobol}{optional integer, only if \code{method="Sobol"}, number of Monte Carlo simulations in Sobol method}
  \item{confInt}{optional boolean, if \code{TRUE}, standard error and 95\% confidence intervals of the indices are computed additionally for \code{method="FixLO"}, defaults to TRUE}
  \item{...}{additional arguments to be passed to the function \code{f.mat}}
}
\value{a list containing a graph structure that includes
       
       \item{d}{number of input factors}
\item{tii}{matrix containing the unscaled total interaction indices and if \code{confInt = TRUE} their standard error and lower and upper confidence limits}
         \item{i1}{matrix containing the unscaled main effect indices}
         \item{V}{overall variance}
         \item{tii.scaled}{matrix containing the scaled total interaction indices}
         \item{cliques}{list of cliques}

 }


\author{
J. Fruth, T. Muehlenstaedt
}

\references{
Fruth, J.; Roustant, O.; Kuhnt, S. (2012+) Total interaction index: A variance-based sensitivity index for interaction screening.

Janon, A.; Klein, T.; Lagnoux-Renaudie, A.; Nodet, M.; Prieur, C. (2012+) Asymptotic normality and efficiency of two Sobol index estimators.

Liu, R.; Owen, A.B. (2006) Estimating mean dimensionality of analysis of variance decompositions, \emph{Journal of the American Statistical Association}, \bold{101} 474, 712-721.

Mara, T.A (2009) Extension of the RBD-FAST method to the computation of global sensitivity indices, \emph{Reliability Engineering & System Safety}, \bold{94} no. 8, 1274-1281.

Muehlenstaedt, T.; Roustant, O.; Carraro, L.; Kuhnt, S. (2011) Data-driven Kriging models based on FANOVA-decomposition, \emph{Statistics and Computing}.

Sobol', I. M. (1993) Sensitivity estimates for nonlinear mathematical models, \emph{Mathematical Modeling and Computational Experiment}, \bold{1}, 407-414.
}
\examples{
# Ishigami function, true analytical values: D12 = D23 = 0, D13 =~ 3.374

estimateGraph(f.mat=ishigami.fun, d=3, q.arg=list(min=-pi,max=pi), N=10000, method="FixLO")
estimateGraph(f.mat=ishigami.fun, d=3, q.arg=list(min=-pi,max=pi), N=10000, method="FixFast")
estimateGraph(f.mat=ishigami.fun, d=3, q.arg=list(min=-pi,max=pi), N=10000, method="RBD")
estimateGraph(f.mat=ishigami.fun, d=3, q.arg=list(min=-pi,max=pi), N=10000, method="Sobol")
}