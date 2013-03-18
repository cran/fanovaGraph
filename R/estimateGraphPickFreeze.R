estimateGraphPickFreeze <- function(f.mat, d, q, q.arg, n.pf, print.loop.index, ...) {
  JK <- combn(d, 2)
  p <- choose(d, 2)
  X.01 <- matrix(runif(n.pf * d), ncol = d)
  Z.01 <- matrix(runif(n.pf * d), ncol = d)
  X <- matrix(nrow = n.pf, ncol = d)
  for (j in 1:d) X[, j] <- do.call(q[j], c(list(p = X.01[, j]), 
                                           q.arg[[j]]))
  Z <- matrix(nrow = n.pf, ncol = d)
  for (j in 1:d) Z[, j] <- do.call(q[j], c(list(p = Z.01[, j]), 
                                           q.arg[[j]]))
  y.all <- f.mat(X) # ... missing!
  y.but.i <- matrix(,n.pf,d)
  for (i in 1:d){
    X2 <- X
    X2[,i] <- Z[,i]
    y.but.i[,i] <- f.mat(X2) # ... missing!
  }
  D <- var(c(y.all, y.but.i))*(n.pf-1)/n.pf
  DC.but.i <- numeric(d)
  for (i in 1:d){
    if(print.loop.index) cat("index =",i,"\n")
    mu <- 1/2*mean(y.all+y.but.i[,i])
    DC.but.i[i] <- mean(y.but.i[,i]*y.all) - mu^2
  }
  DC.but.ij <- numeric(p)
  for (r in 1:p){
    i <- JK[,r][1]
    j <- JK[,r][2]
    if(print.loop.index) cat("index =",i,j ,"\n")
    mu <- 1/2*mean(y.but.i[,i]+y.but.i[,j])
    DC.but.ij[r] <- mean(y.but.i[,i]*y.but.i[,j]) - mu^2
  }
  totalInt <- numeric(p)
  for (r in 1:p) totalInt[r] <- D + DC.but.ij[r] - sum(DC.but.i[JK[,r]])
  inter <- paste("X",JK[1,],"*","X",JK[2,], sep="")
  totalInt <- as.matrix(totalInt)
  rownames(totalInt) <- inter
  return(totalInt)
}

# test
# f.mat <- ishigami.fun
# d <- 3
# q <- rep("qunif",d)
# q.arg <- rep(list(list(min=-pi, max=pi)),d)
# n.pf = 10000
# 
# n.sim <- 500
# x <- matrix(,nrow=n.sim, ncol=3)
# for (i in 1:n.sim){
#   x[i,] <- estimateGraphPickFreezeSaltelliNdplus1(f.mat=f.mat, d=d, q=q, q.arg=q.arg, n.pf=n.pf)
# }
# apply(x, 2, var) #  0.01372603 0.03981987 0.02505435
# # 0.01279237 0.02826905 0.01406206 | 
