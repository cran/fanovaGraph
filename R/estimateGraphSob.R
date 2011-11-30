estimateGraphSob <- function(f.mat, d, q, q.arg, Nsobol, 
    ...) {
    JK <- combn(d, 2)
    p <- choose(d, 2)
    X.01 <- matrix(runif(Nsobol * d), nc = d)
    X1 <- matrix(nr = Nsobol, nc = d)
    for (j in 1:d) X1[, j] <- do.call(q[j], c(list(p = X.01[, j]), 
        q.arg[[j]]))
    y1 <- f.mat(X1, ...)
    f0 <- mean(y1)
    D <- mean(y1^2) - f0^2
    DTij <- numeric(p)
    DTi <- numeric(d)
    v.01 <- matrix(runif(Nsobol * 2), nc = 2)
    for (i in 1:p) {
        print(JK[, i])
        y <- (1:d)[-JK[, i]]
        X2 <- matrix(, Nsobol, d)
        X2[, y] <- X.01[, y]
        X2[, -y] <- v.01
        for (j in 1:d) X2[, j] <- do.call(q[j], c(list(p = X2[, j]), 
            q.arg[[j]]))
        D1 <- mean(y1 * f.mat(X2, ...)) - f0^2
        DTij[i] <- D - D1
    }
    v.01 <- matrix(runif(Nsobol), nc = 1)
    for (i in 1:d) {
        print(i)
        (y <- (1:d)[-i])
        X2 <- matrix(, Nsobol, d)
        X2[, y] <- X.01[, y]
        X2[, -y] <- v.01
        for (j in 1:d) X2[, j] <- do.call(q[j], c(list(p = X2[, j]), 
            q.arg[[j]]))
        D1 <- mean(y1 * f.mat(X2, ...)) - f0^2
        DTi[i] <- D - D1
    }
    totalInt <- numeric(p)
    for (i in 1:p) totalInt[i] <- sum(DTi[JK[, i]]) - DTij[i]
    res <- rbind(JK, totalInt)
    rownames(res) <- NULL
    res
}