estimateGraphFixLO <- function(f.mat, d, q, q.arg, n.lo, confint, ...) {
  ij.set <- combn(d, 2)
  p <- choose(d, 2)
  ZX01 <- matrix(runif(n.lo * 2*d), ncol = 2*d)
  ZX <- matrix(,n.lo,2*d)
  for (j in 1:(2*d)) ZX[, j] <- do.call(rep(q,2)[j], 
                                        c(list(p = ZX01[, j]),rep(q.arg,2)[[j]]))
  Z <- ZX[,1:d]
  X <- ZX[,-(1:d)]
  
  totalInt <- numeric(p)
  if (confint) {
    varInt <- numeric(p)
  }
  for (index in 1:p)
  {
    message(paste("index = ", paste(ij.set[,index], collapse="")))
    i <- ij.set[1,index]; j <- ij.set[2,index]
    X12 <- Z
    X12[,c(i,j)] <- X[,c(i,j)]
    X1 <- Z
    X1[,i] <- X[,i]
    X2 <- Z
    X2[,j] <- X[,j]
    totalInt[index] <- 
      1/(4*n.lo) * sum((f.mat(X12,...)-f.mat(X1,...)-f.mat(X2,...)+f.mat(Z,...))^2)
    if (confint) {
      varInt[index] <- var((f.mat(X12,...)-f.mat(X1,...)-f.mat(X2,...)+f.mat(Z,...))^2) / 16
    }
  }    
  inter <- paste("X", ij.set[1, ], "*", "X", ij.set[2,], sep = "")
  res <- as.matrix(totalInt)
  if (confint){ 
    Std.Error <- sqrt(varInt)/sqrt(n.lo)
    lower <- totalInt - qnorm(0.975)*Std.Error
    upper <- totalInt + qnorm(0.975)*Std.Error
    res <- cbind(totalInt, Std.Error,lower,upper)
  }
  rownames(res) <- inter
  return(res)
}