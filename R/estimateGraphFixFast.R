estimateGraphFixFast <- function(f.mat, d, q, q.arg, nMC, nfast99, ...) {
    # FAST frequencies
    w <- c(11, 35)
    if (nfast99 < 2 * 6 * max(w)) 
        stop("nfast99 too small to guarantee positive values")
    fast <- function(f, ...) {
        Y <- drop(f(X, ...))
        Var <- 2 * sum(colMeans(Y * sinsop)^2 + colMeans(Y * cossop)^2)
        D1 <- 2 * sum(colMeans(Y * sinsop[, harm1])^2 + colMeans(Y * 
            cossop[, harm1])^2)
        D2 <- 2 * sum(colMeans(Y * sinsop[, harm2])^2 + colMeans(Y * 
            cossop[, harm2])^2)
        return(Var - D1 - D2)
    }
    ## the fixed values, only d-2 points needed, but programming much
    #   shorter with d
    SampleFixed <- matrix(runif(nMC * d), ncol = d)
    for (j in 1:d) SampleFixed[, j] <- do.call(q[j], c(list(p = SampleFixed[, 
        j]), q.arg[[j]]))
    # function f.mat but depending only on Xjk, other variables fixed
    #   by xfixed
    fjk.mat <- function(Xjk, jk, xfixed, ...) {
        X <- rep(1, nrow(Xjk)) %*% t(xfixed)
        X[, jk[1]] <- Xjk[, 1]
        X[, jk[2]] <- Xjk[, 2]
        return(f.mat(X, ...))
    }
    JK <- t(combn(1:d, 2))
    DintJK <- numeric(dim(JK)[1])
    # Design matrix for fast:
    s <- seq(-pi, pi, len = nfast99)
    p <- 1:210
    sop <- s %o% p
    sinsop <- sin(sop)
    cossop <- cos(sop)
    harm1 <- 1:6 * w[1]
    harm2 <- 1:6 * w[2]
    X01 <- matrix(, nrow = nfast99, ncol = 2)
    for (i in 1:2)
      X01[, i] <- 1/2 + 1/pi * asin(sin(w[i] * s))
    X <- matrix(, nrow = nfast99, ncol = 2)
    for (j in 1:nrow(JK)) {
        for (i in 1:2) {
           X[, i] <- do.call(q[JK[j,i]], c(list(p = X01[, i]), q.arg[[JK[j,i]]]))
        }
      Dint <- numeric(nMC)
        print(paste("index = ", paste(JK[j, ], collapse="")))
        for (m in (1:nMC)) {
            Dint[m] <- fast(fjk.mat, jk = JK[j, ], xfixed = SampleFixed[m, ], ...)
        }
        DintJK[j] <- mean(Dint)
    }
    return(rbind(j = JK[, 1], k = JK[, 2], DTjk = DintJK))
} 
