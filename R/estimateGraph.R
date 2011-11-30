estimateGraph <- function(f.mat, d, q = NULL, q.arg = NULL, 
    N = NULL, method = "fixed", nMC = NULL, nfast99 = 500, L = NULL, 
    M = 6, Nsobol = NULL, ...) {
    p <- choose(d, 2)
    
    if (is.null(q)) {
        q <- rep("qunif", d)
    } else if (length(q) == 1) {
        q <- rep(q, d)
    }
    if (is.null(q.arg)) {
        q.arg <- rep(list(list()), d)
    } else if (FALSE %in% sapply(q.arg, is.list)) {
        q.arg <- rep(list(q.arg), d)
    }
    
    if (!(method %in% c("fixed", "RBD", "Sobol"))) 
        stop("method must be set to 'fixed', 'RBD' or 'Sobol'")
    
    if (method == "RBD") {
        if (!is.null(L) & !is.null(N)) {
            warning("L will be omitted since N is specified")
            L <- round(N/(2 * (p + d)) - 6 * d, 0)
        }
        if (is.null(L) & !is.null(N)) 
            L <- round(N/(2 * (p + d)) - 6 * d, 0)
        if (is.null(L) & is.null(N)) 
            stop("either N or L must be specified")
        return <- estimateGraphRBD(f.mat, d, q, q.arg, L, M, ...)
    }
    
    if (method == "Sobol") {
        if (!is.null(Nsobol) & !is.null(N)) {
            warning("Nsobol will be omitted since N is specified")
            Nsobol <- round(N/(p + d + 1), 0)
        }
        if (is.null(Nsobol) & !is.null(N)) 
            Nsobol <- round(N/(p + d + 1), 0)
        if (is.null(Nsobol) & is.null(N)) 
            stop("either N or Nsobol must be specified")
        return <- estimateGraphSob(f.mat, d, q, q.arg, Nsobol, ...)
    }
    
    if (method == "fixed") {
        if (!is.null(nMC) & !is.null(N)) {
            warning("nMC will be omitted since N is specified")
            nMC <- round(N/(p * nfast99), 0)
        }
        if (is.null(nMC) & !is.null(N)) 
            nMC <- round(N/(p * nfast99), 0)
        if (is.null(nMC) & is.null(N)) 
            stop("either N or nMC must be specified")
        return <- estimateGraphFix(f.mat, d, q, q.arg, nMC, nfast99, 
            ...)
    }
    
    return(return)
}