\name{MLoptimConstrained}
\alias{MLoptimConstrained}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Constrained MLE Optimization
}
\description{Constrained MLE optimization for kernels defined by cliques using \code{\link{constrOptim}}
}
\usage{
MLoptimConstrained(x, y, n.initialtries = 50, limits = NULL, eps.R = 1e-08, Cl, covtype = "gauss", eps.Var = 1e-06, MAXIT = 1000, iso = FALSE)
}

\arguments{
  \item{x}{a design matrix of input variables, number of columns should be number of variables
}
  \item{y}{a vector of output variables of the same length as the columns of L
}
  \item{n.initialtries}{number of random initial parameters for optimization, defaults to 50
}
  \item{limits}{a list with items lower, upper containing boundaries for the covariance parameter vector theta, if \code{NULL} suitable bounds are computed from the range of \code{x}
}
  \item{eps.R}{small positive number indicating the nugget effect added to the covariance matrix diagonalk, defaults to \code{eps.R = 1e-08}
}
  \item{Cl}{list of cliques, can be obtained by function \code{\link{threshold}}
}
  \item{covtype}{an optional character string specifying the covariance structure to be used, 
    to be chosen between "gauss", "matern5_2", "matern3_2", "exp" or "powexp" (see \code{\link[DiceKriging:DiceKriging-package]{DiceKriging}}), defaults to "gauss"
}
  \item{eps.Var}{small positive number providing the limits for the alpha parameters in order to guarantee strict inequalities (0+eps.Var <= alpha <= 1-esp.Var), defaults to \code{eps.Var = 1e-06}
}
  \item{MAXIT}{maximum number of iterations for optimization, defaults to \code{MAXIT=1000}
}
  \item{iso}{boolean vector indicating for each clique if it is isotropic (TRUE) or anisotropic (FALSE) with isotropic cliques at the end, defaults to \code{iso = FALSE} (all cliques anisotropic)
}
}
\value{
list of estimated parameter 'alpha' and 'theta' corresponding to the clique structure in 'Cl'}
\references{
Muehlenstaedt, T.; Roustant, O.; Carraro, L.; Kuhnt, S. (2011) Data-driven Kriging models based on FANOVA-decomposition, \emph{Statistics and Computing}.
}
\author{
T. Muehlenstadt, O. Roustant, J. Fruth
}
\seealso{
\code{\link{yhat}}
}
\examples{
### example for ishigami function with cliques {1,3} and {2}
d <- 3
x <- matrix(runif(100*d,-pi,pi),nc=d)
y <- ishigami.fun(x)

Cl <- list(c(2), c(1,3))

# constrained ML optimation with kernel defined by the cliques
parameter <- MLoptimConstrained(x, y, Cl = Cl)

# prediction with the new model
xpred <- matrix(runif(500 * d,-pi,pi), ncol = d)
ypred <- yhat(xpred, x, y, parameter, Cl=Cl)
yexact <- ishigami.fun(xpred)

# rmse
sqrt(mean((ypred[,1]- yexact)^2))

# scatterplot
par(mfrow=c(1,1))
plot(yexact, ypred[,1], asp = 1)
abline(0, 1)

### compare to one single clique {1,2,3}
Cl <- list(c(1,2,3))

# constrained ML optimation with kernel defined by the cliques
parameter <- MLoptimConstrained(x, y, Cl = Cl)

# prediction with the new model
ypred <- yhat(xpred, x, y, parameter, Cl=Cl)

# rmse
sqrt(mean((ypred$mean- yexact)^2))

# scatterplot
par(mfrow=c(1,1))
plot(yexact, ypred$mean, asp = 1)
abline(0, 1)

### isotropic cliques

Cl <- list(c(2),c(1,3))
parameter <- MLoptimConstrained(x, y, Cl = Cl, iso=c(FALSE,TRUE))
ypred <- yhat(xpred, x, y, parameter, Cl=Cl, iso=c(FALSE,TRUE))
sqrt(mean((ypred$mean- yexact)^2))

# the same since first clique has length 1
parameter <- MLoptimConstrained(x, y, Cl = Cl, iso=c(TRUE,TRUE))
ypred <- yhat(xpred, x, y, parameter, Cl=Cl, iso=c(TRUE,TRUE))
sqrt(mean((ypred$mean- yexact)^2))
}
