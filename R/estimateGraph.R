estimateGraph <- function(f.mat, d, q = NULL, q.arg = NULL, 
    N = NULL, method = "FixLO", nLO = NULL, nMC = NULL, 
    nfast99 = 500, L = NULL, M = 6, Nsobol = NULL, confInt = TRUE, ...) {
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
    
    if (!(method %in% c("FixLO", "FixFast", "RBD", "Sobol"))) 
        stop("method must be set to 'FixLO', 'FixFast' ,'RBD' or 'Sobol'")

    if (method == "RBD") {
        if (!is.null(L) & !is.null(N)) {
            warning("L will be omitted since N is specified")
            L <- round(N/(2 * (p + d)) - 6 * d, 0)
        }
        if (is.null(L) & !is.null(N)) 
            L <- round(N/(2 * (p + d)) - 6 * d, 0)
        if (is.null(L) & is.null(N)) 
            stop("either N or L must be specified")
        tii <- estimateGraphRBD(f.mat, d, q, q.arg, L, M, ...)
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
        tii <- estimateGraphSob(f.mat, d, q, q.arg, Nsobol, ...)
    }
    
    if (method == "FixLO") {
      if (!is.null(nLO) & !is.null(N)) {
        warning("nLO will be omitted since N is specified")
        nLO <- round(N/(4*p), 0)
      }
      if (is.null(nLO) & !is.null(N)) 
        nLO <- round(N/(4*p), 0)
      if (is.null(nLO) & is.null(N)) 
        stop("either N or nLO must be specified")
      tii <- estimateGraphFixLO(f.mat, d, q, q.arg, nLO, confInt, ...)
    }
        
    if (method == "FixFast") {
        if (!is.null(nMC) & !is.null(N)) {
            warning("nMC will be omitted since N is specified")
            nMC <- round(N/(p * nfast99), 0)
        }
        if (is.null(nMC) & !is.null(N)) 
            nMC <- round(N/(p * nfast99), 0)
        if (is.null(nMC) & is.null(N)) 
            stop("either N or nMC must be specified")
        tii <- estimateGraphFixFast(f.mat, d, q, q.arg, nMC, nfast99, 
            ...)
    }
    soverall <- fast99(model=f.mat, factors=d, n=1000, q = q, q.arg = q.arg, ...)
    V <- round(mean(soverall$V),4)
    i1 <- as.matrix(round(soverall$D1,4))
    rownames(i1) <- paste("X",1:d,sep="")
    # estimate cliques
      E <- t(combn(d,2)[,tii[,1] > 0])
      cliques <- maximal.cliques(graph(as.vector(t(E)), d + 1, FALSE))[-1]
      tii.scaled <- tii[,1,drop=FALSE] / V
    return(list(d=d, tii=tii, i1=i1, V = V, tii.scaled = tii.scaled, cliques = cliques))
}