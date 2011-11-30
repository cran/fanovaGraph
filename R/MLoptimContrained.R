LogLConstrained <- function(parameter, y, x, DM, n, COR, 
    eps.R = 1e-12, n.Cl, Cl, LogLik = "ML", iso = FALSE) {
    ### Likelihood for constrained optimization
    parameter <- c(parameter[1:(n.Cl - 1)], 1 - sum(parameter[1:(n.Cl - 
        1)]), parameter[(n.Cl):length(parameter)])
    alpha <- parameter[1:n.Cl]
    theta <- parameter[-(1:n.Cl)]
    R <- Rfunc(as.matrix(x), theta, alpha, COR, n, n.Cl, Cl, iso) + 
        diag(eps.R, ncol = n, nrow = n)
    F <- matrix(rep(1, n), ncol = 1)
    T <- chol(R)
    M <- backsolve(t(T), F, upper.tri = FALSE)
    Tinv_y <- backsolve(t(T), matrix(y, ncol = 1), upper.tri = FALSE)
    Q <- qr.Q(qr(M))
    H <- Q %*% t(Q)
    z <- Tinv_y - H %*% Tinv_y  # linear regression properties
    v <- t(z) %*% z/n
    loglik <- -0.5 * (n * log(2 * pi * v) + 2 * sum(log(diag(T))) + 
        n)
    return(loglik)
}

MLoptimConstrained <- function(seed = 1, L, y, n.initialtries = 50, 
    limits = c(0.001, 4), eps.R = 1e-08, Cl, COR = "gauss", eps.Var = 1e-06, 
    MAXIT = 1000, iso = FALSE) {
    ### based on  LogLCl4
    set.seed(seed)
    n.Cl <- length(Cl)
    if (identical(iso, FALSE)) {
        iso <- rep(FALSE, n.Cl)
    }
    if (n.Cl == 1) {
        warning("since only one clique is given, DiceKriging:::km is used and\nan object of class km returned")
        parameter <- km(~1, design = data.frame(L), response = y, covtype = COR, 
            iso = iso)
    } else {
        n.Cl.ani <- n.Cl - sum(iso)
        n.Cl.iso <- sum(iso)
        
        d <- ncol(L)
        n <- length(y)
        theta.n <- sum(as.numeric(lapply(Cl[which(iso == FALSE)], length)))
        theta.n <- theta.n + sum(iso)
        DM <- matrix(1, ncol = 1, nrow = n)
        
        Ui <- c(rep(-1, n.Cl - 1), rep(0, theta.n))
        Ui <- rbind(Ui, cbind(diag(rep(1, n.Cl - 1)), matrix(0, ncol = theta.n, 
            nrow = n.Cl - 1)))
        Ui <- rbind(Ui, cbind(matrix(0, ncol = n.Cl - 1, nrow = theta.n), 
            diag(rep(1, theta.n))))
        Ui <- rbind(Ui, cbind(matrix(0, ncol = n.Cl - 1, nrow = theta.n), 
            diag(rep(-1, theta.n))))
        Ci <- c(-(1 - eps.Var), rep(eps.Var, n.Cl - 1), rep(limits[1], 
            theta.n), rep(-limits[2], theta.n))
        ## choosing initial points
        alphastart <- matrix(runif(n.initialtries * n.Cl), ncol = n.Cl)
        alphastart <- (alphastart/apply(alphastart, 1, sum))[, -n.Cl]
        thetastart <- matrix(runif(n.initialtries * theta.n, limits[1], 
            limits[2]), nrow = n.initialtries)  #  choosing a random start vector
        parameterstart <- cbind(alphastart, thetastart)
        
        LLinitial <- apply(parameterstart, 1, LogLConstrained, y = y, 
            x = L, DM = DM, n = n, COR = COR, eps.R = eps.R, n.Cl = n.Cl, 
            Cl = Cl, iso = iso)
        parameter <- parameterstart[which.max(LLinitial), ]
        
        ### optimization the constrained LL
        test <- try(constrOptim(theta = parameter, f = LogLConstrained, 
            grad = NULL, ui = Ui, ci = Ci, y = y, x = L, DM = DM, n = n, 
            COR = COR, eps.R = eps.R, n.Cl = n.Cl, Cl = Cl, iso = iso, 
            outer.iterations = 10, outer.eps = 1e-05, control = list(fnscale = -1, 
                maxit = MAXIT, trace = 0)))
        parameter <- test$par
        parameter <- c(parameter[1:(n.Cl - 1)], 1 - sum(parameter[1:(n.Cl - 
            1)]), parameter[n.Cl:length(parameter)])
        names(parameter) <- c(paste("alpha", 1:n.Cl, sep = ""), paste("theta", 
            1:theta.n, sep = ""))
    }
    return(parameter)
} 