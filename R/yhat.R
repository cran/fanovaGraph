####################################################################
#### correlation functions
####################################################################
#### Correlation functions including the additive structure,
#### employ DICEKriging correlation functions.
####################################################################

Rfunc <- function(x, theta, alpha, COR, n, n.Cl, Cl, iso) {
    R <- matrix(0, ncol = n, nrow = n)
    theta.list <- Cl
    ntemp <- 0
    n.Cl.ani <- n.Cl - sum(iso)
    for (j in (1:n.Cl.ani)) {
        theta.list[[j]] <- theta[(ntemp + 1):(ntemp + length(Cl[[j]]))]
        ntemp <- ntemp + length(Cl[[j]])
    }
    if (n.Cl.ani < n.Cl) {
        for (j in ((n.Cl.ani + 1):n.Cl)) {
            theta.list[[j]] <- theta[ntemp]
            ntemp <- ntemp + 1
        }
    }
    for (j in 1:n.Cl) {
        COV <- covStruct.create(covtype = COR, d = length(Cl[[j]]), 
            var.names = NULL, known.covparam = "None", coef.cov = theta.list[[j]], 
            coef.var = alpha[j], iso = iso[j])
        R <- R + covMatrix(object = COV, X = as.matrix(x[, Cl[[j]]]))[[1]]
    }
    return(R)
}

kleinr <- function(x, L, theta, alpha, COR, n.Cl, Cl, iso) {
    theta.list <- Cl
    ntemp <- 0
    n.Cl.ani <- n.Cl - sum(iso)
    for (j in (1:n.Cl.ani)) {
        theta.list[[j]] <- theta[(ntemp + 1):(ntemp + length(Cl[[j]]))]
        ntemp <- ntemp + length(Cl[[j]])
    }
    if (n.Cl.ani < n.Cl) {
        for (j in ((n.Cl.ani + 1):n.Cl)) {
            theta.list[[j]] <- theta[ntemp]
            ntemp <- ntemp + 1
        }
    }
    kleinr <- numeric(nrow(L))
    for (j in 1:n.Cl) {
        COV <- covStruct.create(covtype = COR, d = length(Cl[[j]]), 
            var.names = NULL, known.covparam = "None", coef.cov = theta.list[[j]], 
            coef.var = alpha[j], iso = iso[j])
        kleinr <- kleinr + covMat1Mat2(object = COV, X1 = as.matrix(L[, 
            Cl[[j]]]), X2 = t(x[Cl[[j]]]))
    }
    return(kleinr)
}

####################################################################
#### prediction function
####################################################################
#### standard kriging prediction function for the modified
#   correlation functions
####################################################################

yhat <- function(x, y, parameter, COR = "gauss", L, eps.R = 1e-08, 
    Cl, iso = FALSE) {
    n.Cl <- length(Cl)
    if (identical(iso, FALSE)) {
        iso <- rep(FALSE, n.Cl)
    }
    if (n.Cl == 1 & class(parameter) == "km") {
        warning("since only one clique is given, DiceKriging:::predict.km is used")
        pred <- predict.km(parameter, newdata = x, type = "UK")
        result <- data.frame(pred$mean, pred$sd)
    } else {
        if (n.Cl == 1 & class(parameter) != "km") 
            stop("in case of only one clique a kriging model of class km is required for 'parameter'")
        L <- as.matrix(L)
        alpha <- parameter[(1:n.Cl)]
        thetastern <- parameter[-(1:n.Cl)]
        n <- length(y)
        DM <- matrix(1, ncol = 1, nrow = n)
        p <- 1  #  number of trend parameters
        R <- Rfunc(as.matrix(L), thetastern, alpha, COR, n, n.Cl, Cl, 
            iso) + diag(eps.R, ncol = n, nrow = n)
        Rinvs <- solve(R)
        betadach <- solve(t(DM) %*% Rinvs %*% DM) %*% t(DM) %*% Rinvs %*% 
            y
        Faktor2 <- Rinvs %*% (y - DM %*% betadach)
        sigmadachhoch2 <- (1/(n - p)) * t((y - DM %*% betadach)) %*% 
            Rinvs %*% (y - DM %*% betadach)
        ydach <- numeric(length(x[, 1]))
        sigmadachx <- numeric(length(x[, 1]))
        for (i in (1:length(x[, 1]))) {
            kleinrx <- as.matrix(kleinr(x[i, ], L, thetastern, alpha, 
                COR, n.Cl, Cl, iso))
            ydach[i] <- betadach + t(kleinrx) %*% Faktor2
            sigmadachx[i] <- sigmadachhoch2 * (1 - t(kleinrx) %*% Rinvs %*% 
                kleinrx + (1 - t(kleinrx) %*% Rinvs %*% DM) * (1/(t(DM) %*% 
                Rinvs %*% DM)) * (1 - t(kleinrx) %*% Rinvs %*% DM))
        }
        result <- cbind(ydach, sigmadachx)
        colnames(result) <- c("mean", "sd")
    }
    return(result)
} 