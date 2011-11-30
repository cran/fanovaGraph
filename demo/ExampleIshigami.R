#################################
# Ishigami function application #
#################################

require(fanovaGraph)

### function definition:

d <- 3
domain <- c(-pi, pi)
fun <- ishigami.fun

### maximin design via package 'lhs'

library(lhs)
L01 <- maximinLHS(100, d)
L <- L01 * (domain[2] - domain[1]) + domain[1]

### kriging model via package 'DiceKriging'

y <- fun(L)
KM <- km(~1, design = data.frame(L), response = y)

### prediction function (the new metamodel!)

krigingMean <- function(Xnew) predict.km(object = KM, newdata = Xnew, 
    type = "UK", se.compute = FALSE)$mean

### standard Sobol indices with package 'sensitivity'

i1 <- fast99(model = krigingMean, factors = d, n = 2000, q = "qunif", 
    q.arg = list(min = domain[1], max = domain[2]))
plot(i1)


### estimation of total interaction indices via fixing method

totalInt <- estimateGraph(f.mat = krigingMean, d = d, N = 40000, q.arg = 
  list(min = domain[1], max = domain[2]), method = "fixed") 

### plotting the thresholded graph without and with weights

graph <- threshold(delta = 0.1, totalInt = totalInt, d = d, dall = mean(i1$V))
plotiGraph(graph$E, d)
plotiGraph(E = graph$E, d = d, i1 = i1$D1/mean(i1$V), tii = graph$tii.scaled)

### estimate new model

Cliques <- graph[[1]]
parameter <- MLoptimConstrained(seed = 1, L, y, Cl = Cliques)

### comparison to standard kriging

xpred <- matrix(runif(d * 1000, domain[1], domain[2]), ncol = d)
y_new <- yhat(xpred, y, parameter, L = L, Cl = Cliques)
y_old <- krigingMean(xpred)
y_exact <- fun(xpred)

op <- par("mfrow")
par(mfrow = c(1, 2))
plot(y_exact, y_old, asp = 1, xlab="y, exact", ylab="y, predicted", main="Standard Kernel")
abline(0, 1)
plot(y_exact, y_new[, 1], asp = 1, xlab="y, exact", ylab="y, predicted", main="Modified Kernel")
abline(0, 1)
par(mfrow = op)

sqrt(mean((y_old - y_exact)^2))
sqrt(mean((y_new[, 1] - y_exact)^2))