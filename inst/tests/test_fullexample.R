#################################################################
# One full application example of the package 'fanovaGraph'
#################################################################

test_that("the full example works as before", {
  
  
d <- 6     
domain <- c(-1,1)

Bsp <- function(x){
  beta <- c(-0.8, -1.1, 1.1, 1)
  gamma <- c(-0.5, 0.9, 1, -1.1)
  result <- cos(cbind(1, x[,c(1,5,3)]) %*% beta) + 
    sin(cbind(1, x[,c(4,2,6)]) %*% gamma)   # function a
  return(result)
}

### maximin design with package "lhs" 

library(lhs)
set.seed(1)
L01 <- maximinLHS(100, d)
L <- L01 * (domain[2]-domain[1]) + domain[1]

### kriging model with package "DiceKriging"

y <- Bsp(L)
covtype <- "matern5_2"
set.seed(1)
KM <- km( ~ 1, design = data.frame(L), response = y, covtype = covtype,)

### prediction function (the new metamodel!)

krigingMean <- function(Xnew) 
  predict(object = KM, newdata = Xnew, type = "UK", se.compute = FALSE, checkNames=FALSE)$mean

### model validation

plot(KM)
par(mfrow = c(1, 1))
npred <- 1000
xpred <- matrix(runif(d * npred,domain[1],domain[2]), ncol = d) 
yexact <- Bsp(xpred)
yKM <- krigingMean(xpred)
plot(yexact, yKM, asp = 1)
abline(0, 1)

expect_equal(0.13046535, sqrt(mean((yKM- yexact)^2)))

### standard Sobol indices with package "sensitivity"
### and function "krigingMean"

soverall <- fast99(model = krigingMean, factors = d, n = 2000, 
                   q = "qunif", q.arg = list(min = domain[1], max = domain[2]))
expect_equal(0.7772051, mean(soverall$V))

#############################################################
### estimate total interaction indices for the graph edges
#############################################################

### estimation by function "estimateGraph" with fixing method
set.seed(1)
totalInt <- estimateGraph(f.mat=krigingMean, d=d, N=10000,
                          q.arg=list(min=domain[1],max=domain[2]), method="FixFast")
expect_equivalent(c(0.0033, 0.0095), totalInt$tii[1:2,1])  

totalInt <- estimateGraph(f.mat=krigingMean, d=d, N=10000,
                          q.arg=list(min=domain[1],max=domain[2]), method="FixLO")
expect_equivalent(c(0.0012, 0.0401), totalInt$tii[1:2,1])


##########################################
### Graph plotting
##########################################

### creating the graph with cut at indices > 0.01 by function "threshold"

graph <- threshold(totalInt, delta = 0.01)

Cl <- graph$cliques
expect_equal(list(c(1,3,5),c(2,4,6)),Cl)


##############################
### estimating the new model
##############################

## new kernel estimation by function "MLoptimConstrained"
## and prediction by function "yhat"

Cl <- list(c(1, 3, 5),c(2, 4, 6))
eps.R <-  1e-04
eps.Var <- 1e-4
nMaxit <- 100
nInitial <- 20

set.seed(1)
parameter <- MLoptimConstrained(x=L, y, n.initialtries = nInitial, eps.R = eps.R, 
      Cl=Cl, covtype = covtype, eps.Var = eps.Var, MAXIT = nMaxit, iso = FALSE)
expect_equal(0.492900404, parameter[[1]]$alpha)
ypred <- yhat(xpred, x=L, y, parameter, covtype = covtype, eps.R = eps.R, Cl)
expect_equal(0.0429707179, sqrt(mean((ypred$mean - yexact)^2)))

### isotropic clique
set.seed(1)
parameter <- MLoptimConstrained(x=L, y, n.initialtries = nInitial, eps.R = eps.R, 
                                Cl=Cl, covtype = covtype, eps.Var = eps.Var, MAXIT = nMaxit, iso = c(FALSE,TRUE))
expect_equal(1.90368446, parameter[[2]]$theta)
ypred <- yhat(xpred, x=L, y, parameter, covtype = covtype, eps.R = eps.R, Cl, iso=c(FALSE,TRUE))
expect_equal(0.0424629352, sqrt(mean((ypred$mean- yexact)^2)))

### two isotropic cliques
set.seed(1)
parameter <- MLoptimConstrained(x=L, y, n.initialtries = nInitial, 
                                eps.R = eps.R, Cl=Cl, covtype = covtype, eps.Var = eps.Var, 
                                MAXIT = nMaxit, iso = c(TRUE,TRUE))
expect_equal(1.90394986, parameter[[2]]$theta)
ypred <- yhat(xpred, x=L, y, parameter, covtype = covtype, eps.R = eps.R, 
              Cl, iso=c(TRUE,TRUE))
expect_equal(0.0424650539, sqrt(mean((ypred$mean- yexact)^2)))

### single clique (km is used for estimation and prediction)
Cl <- list(1:6)
set.seed(1)
parameter <- MLoptimConstrained(x=L, y, Cl=Cl, covtype=covtype)
expect_equal(parameter@covariance@sd2, KM@covariance@sd2)
ypred <- yhat(xpred, x=L, y, parameter, covtype = covtype, Cl=Cl)
expect_true(all(ypred$mean == yKM))

### simulation from the new kernel
Cl <- list(c(1, 3, 5),c(2, 4, 6))
parameter <- list(list(alpha=0.53,theta=c(1.7,1.8,1.8)),
                  list(alpha=0.47,theta=c(1.9,1.9,1.7)))
set.seed(1)
xsimu <- matrix(runif(d*3,domain[1],domain[2]),3,d)

expect_equal(c(-0.305388389,  0.792219843,  0.313896841),
       simAdd(xsimu, mu=0, parameter, covtype, Cl, iso=FALSE))    

parameter <- list(list(alpha=0.53,theta=c(1.7,1.8,1.8)),
                  list(alpha=0.47,theta=c(1.1)))
set.seed(1)
xsimu <- matrix(runif(d*3,domain[1],domain[2]),3,d)
expect_equal(c(-0.305388389, 0.812151639, 0.288423824),
     simAdd(xsimu, mu=0, parameter, covtype, Cl, iso=c(FALSE,TRUE)))
# [1] -0.3053884  0.8121516  0.2884238
})
