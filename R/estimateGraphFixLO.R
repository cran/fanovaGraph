estimateGraphFixLO <- function(f.mat, d, q, q.arg, nLO, confInt, ...) {
  ij.set <- combn(d, 2)
  p <- choose(d, 2)
  ZX01 <- matrix(runif(nLO * 2*d), ncol = 2*d)
  ZX <- matrix(,nLO,2*d)
  for (j in 1:(2*d)) ZX[, j] <- do.call(rep(q,2)[j], 
                                        c(list(p = ZX01[, j]),rep(q.arg,2)[[j]]))
  Z <- ZX[,1:d]
  X <- ZX[,-(1:d)]
  
  totalInt <- numeric(p)
  if (confInt) {
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
      1/(4*nLO) * sum((f.mat(X12,...)-f.mat(X1,...)-f.mat(X2,...)+f.mat(Z,...))^2)
    if (confInt) {
      varInt[index] <- var((f.mat(X12,...)-f.mat(X1,...)-f.mat(X2,...)+f.mat(Z,...))^2) / 16
    }
  }    
  inter <- paste("X", ij.set[1, ], "*", "X", ij.set[2,], sep = "")
  res <- as.matrix(round(totalInt,4))
  if (confInt){ 
    Std.Error <- sqrt(varInt)/sqrt(nLO)
    lower <- totalInt - qnorm(0.975)*Std.Error
    upper <- totalInt + qnorm(0.975)*Std.Error
    res <- round(cbind(totalInt, Std.Error,lower,upper),4)
  }
  rownames(res) <- inter
  return(res)
}