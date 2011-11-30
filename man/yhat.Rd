\name{yhat}
\alias{yhat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Prediction Function with Modified Kernel
}
\description{Standard kriging prediction function for the modified correlation functions.
}
\usage{
yhat(x, y, parameter, COR = "gauss", L, eps.R = 1e-08, Cl, iso = FALSE)
}
\arguments{
  \item{x}{matrix containing the points where to perform predictions
}
  \item{y}{vector of output data
}
  \item{parameter}{(by \code{\link{MLoptimConstrained}} estimated) kriging parameters
}
  \item{COR}{an optional character string specifying the covariance structure to be used, 
    to be chosen between "gauss", "matern5_2", "matern3_2", "exp" or "powexp" (see DiceKriging), defaults to "gauss"
}
  \item{L}{matrix of input data
}
  \item{eps.R}{small positive number indicating the nugget effect added to the covariance matrix diagonalk, defaults to \code{eps.R = 1e-08}
}
  \item{Cl}{list of cliques
}
  \item{iso}{boolean vector indicating for each clique if it is isotropic (TRUE) or anisotropic (FALSE) with isotropic cliques at the end, defaults to \code{iso = FALSE} (all cliques anisotropic)
}
}
\value{matrix of predicted values (kriging mean) at \code{x} and kriging standard deviation at \code{x}
}
\references{
Muehlenstaedt, T.; Roustant, O.; Carraro, L.; Kuhnt, S. (2011) Data-driven Kriging models based on FANOVA-decomposition, \emph{Statistics and Computing}.
}
\author{
T. Muehlenstaedt, O. Roustant, J. Fruth
}

\seealso{
\code{\link{MLoptimConstrained}}
}
\examples{
# example for ishigami function with cliques {1,3} and {2}
d <- 3
L <- matrix(runif(100*d,-pi,pi),nc=d)
y <- ishigami.fun(L)

Cl <- list(c(1,3),c(2))

# constrained ML optimation with kernel defined by the cliques
parameter <- MLoptimConstrained(seed = 1, L, y, Cl = Cl)

# prediction with the new model
xpred <- matrix(runif(500 * d,-pi,pi), ncol = d)
ypred <- yhat(xpred, y, parameter, L=L, Cl=Cl)
yexact <- ishigami.fun(xpred)

# rmse
sqrt(mean((ypred[,1]- yexact)^2))

# scatterplot
par(mfrow=c(1,1))
plot(yexact, ypred[,1], asp = 1)
abline(0, 1)
}
