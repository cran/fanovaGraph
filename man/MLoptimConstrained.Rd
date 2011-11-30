\name{MLoptimConstrained}
\alias{MLoptimConstrained}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Constrained MLE Optimization
}
\description{Constrained MLE optimization for kernels defined by cliques using \code{\link{constrOptim}}
}
\usage{
MLoptimConstrained(seed = 1, L, y, n.initialtries = 50, limits = c(0.001, 4), eps.R = 1e-08, Cl, COR = "gauss", eps.Var = 1e-06, MAXIT = 1000, iso = FALSE)
}

\arguments{
  \item{seed}{a numerical number for set.seed command at the beginning
}
  \item{L}{a design matrix of input variables, number of columns should be number of variables
}
  \item{y}{a vector of output variables of the same length as the columns of L
}
  \item{n.initialtries}{number of random initial parameters for optimization, defaults to 50
}
  \item{limits}{vector of length 2 indicating the constraining limits for theta parameters, defaults to \code{limits=c(0.001,4)}
}
  \item{eps.R}{small positive number indicating the nugget effect added to the covariance matrix diagonalk, defaults to \code{eps.R = 1e-08}
}
  \item{Cl}{list of cliques, can be obtained by function \code{\link{threshold}}
}
  \item{COR}{an optional character string specifying the covariance structure to be used, 
    to be chosen between "gauss", "matern5_2", "matern3_2", "exp" or "powexp" (see DiceKriging), defaults to "gauss"
}
  \item{eps.Var}{small positive number providing the limits for the alpha parameters in order to guarantee strict inequalities (0+eps.Var <= alpha <= 1-esp.Var), defaults to \code{eps.Var = 1e-06}
}
  \item{MAXIT}{maximum number of iterations for optimization, defaults to \code{MAXIT=1000}

}
  \item{iso}{boolean vector indicating for each clique if it is isotropic (TRUE) or anisotropic (FALSE) with isotropic cliques at the end, defaults to \code{iso = FALSE} (all cliques anisotropic)
}
}
\value{
vector of estimated parameters for alpha and theta}
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
# example for ishigami function with cliques {1,3} and {2}
set.seed(1)
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
