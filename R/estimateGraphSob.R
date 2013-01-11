estimateGraphSob <- function(f.mat, d, q, q.arg, n.sobol, 
    ...) {
    JK <- combn(d, 2)
    p <- choose(d, 2)
    X.01 <- matrix(runif(n.sobol * d), ncol = d)
    X1 <- matrix(nrow = n.sobol, ncol = d)
    for (j in 1:d) X1[, j] <- do.call(q[j], c(list(p = X.01[, j]), 
        q.arg[[j]]))
    y1 <- f.mat(X1, ...)
    DTij <- numeric(p)
    DTi <- numeric(d)
    v.01 <- matrix(runif(n.sobol * 2), ncol = 2)
    for (i in 1:p) {
        message(paste("index = ", paste(JK[,i], collapse="")))
      
        y <- (1:d)[-JK[, i]]
        X2 <- matrix(, n.sobol, d)
        X2[, y] <- X.01[, y]
        X2[, -y] <- v.01
        for (j in 1:d) X2[, j] <- do.call(q[j], c(list(p = X2[, j]), 
            q.arg[[j]]))
        y2 <- f.mat(X2,...)
        D1 <- mean(y1 * y2) - mean( (y1+y2)/2 )^2 #closed index {-ij}
        DTij[i] <- mean( (y1^2+y2^2)/2) - mean( (y1+y2)/2 )^2 - D1
    }
    v.01 <- matrix(runif(n.sobol), ncol = 1)
    for (i in 1:d) {
        message(paste("index = ", i))
        (y <- (1:d)[-i])
        X2 <- matrix(, n.sobol, d)
        X2[, y] <- X.01[, y]
        X2[, -y] <- v.01
        for (j in 1:d) X2[, j] <- do.call(q[j], c(list(p = X2[, j]), 
            q.arg[[j]]))
        y2 <- f.mat(X2,...)
        D1 <- mean(y1 * y2) - mean( (y1+y2)/2 )^2
        DTi[i] <- mean( (y1^2+y2^2)/2) - mean( (y1+y2)/2 )^2 - D1
    }
    totalInt <- numeric(p)
    for (i in 1:p) totalInt[i] <- sum(DTi[JK[, i]]) - DTij[i]
    inter <- paste("X",JK[1,],"*","X",JK[2,], sep="")
    totalInt <- as.matrix(totalInt)
    rownames(totalInt) <- inter
    return(totalInt)
}
