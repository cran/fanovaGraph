estimateGraph <- function(f.mat, d, q = NULL, q.arg = NULL, 
    n.tot = NULL, method = "FixLO", n.lo = NULL, n.mc = NULL, 
    n.fast = 500, L = NULL, M = 6, n.sobol = NULL, confint = TRUE, ...) {
  p <- choose(d, 2)
  if (is.null(q)) {
      q <- rep("qunif", d)
  } else if (length(q) == 1) {
      q <- rep(q, d)
  }
  if (is.null(q.arg)) {
      q.arg <- rep(list(list()), d)
  } else if (all(sapply(q.arg, function(x) is.numeric(x) & length(x)==1))) {
      q.arg <- rep(list(q.arg), d)
  }
#   input test
  if (!(all(sapply(q.arg, is.list)) & length(q.arg)==d))
    stop("q.arg must be a list of lists of quantile functions parameters")
  
  Xtest <- matrix(0.5,2,d)
   for (j in (1:d))
     Xtest[, j] <- do.call(rep(q,2)[j],c(list(p = Xtest[, j]),rep(q.arg,2)[[j]]))
   if (length(f.mat(Xtest,...))!=2)
     stop("f.mat must be a vectorized function (matrix input must be possible)")
  
  if (!(method %in% c("FixLO", "FixFast", "RBD", "Sobol"))) 
      stop("method must be set to 'FixLO', 'FixFast' ,'RBD' or 'Sobol'")
  
  if (method == "RBD") {
      if (!is.null(L) & !is.null(n.tot)) {
          warning("L will be omitted since n.tot is specified")
          L <- round(n.tot/(2 * (p + d)) - 6 * d, 0)
      }
      if (is.null(L) & !is.null(n.tot)) 
          L <- round(n.tot/(2 * (p + d)) - 6 * d, 0)
      if (is.null(L) & is.null(n.tot)) 
          stop("either n.tot or L must be specified")
      tii <- estimateGraphRBD(f.mat, d, q, q.arg, L, M, ...)
  }
  
  if (method == "Sobol") {
      if (!is.null(n.sobol) & !is.null(n.tot)) {
          warning("n.sobol will be omitted since n.tot is specified")
          n.sobol <- round(n.tot/(p + d + 1), 0)
      }
      if (is.null(n.sobol) & !is.null(n.tot)) 
          n.sobol <- round(n.tot/(p + d + 1), 0)
      if (is.null(n.sobol) & is.null(n.tot)) 
          stop("either n.tot or n.sobol must be specified")
      tii <- estimateGraphSob(f.mat, d, q, q.arg, n.sobol, ...)
  }
  
  if (method == "FixLO") {
    if (!is.null(n.lo) & !is.null(n.tot)) {
      warning("n.lo will be omitted since n.tot is specified")
      n.lo <- round(n.tot/(4*p), 0)
    }
    if (is.null(n.lo) & !is.null(n.tot)) 
      n.lo <- round(n.tot/(4*p), 0)
    if (is.null(n.lo) & is.null(n.tot)) 
      stop("either n.tot or n.lo must be specified")
    tii <- estimateGraphFixLO(f.mat, d, q, q.arg, n.lo, confint, ...)
  }
      
  if (method == "FixFast") {
      if (!is.null(n.mc) & !is.null(n.tot)) {
          warning("n.mc will be omitted since n.tot is specified")
          n.mc <- round(n.tot/(p * n.fast), 0)
      }
      if (is.null(n.mc) & !is.null(n.tot)) 
          n.mc <- round(n.tot/(p * n.fast), 0)
      if (is.null(n.mc) & is.null(n.tot)) 
          stop("either n.tot or n.mc must be specified")
      tii <- estimateGraphFixFast(f.mat, d, q, q.arg, n.mc, n.fast, 
          ...)
  }
  soverall <- fast99(model=f.mat, factors=d, n=1000, q = q, q.arg = q.arg,...)
  V <- mean(soverall$V)
  i1 <- round(as.matrix(soverall$D1), 27)
  rownames(i1) <- paste("X",1:d,sep="")
  tii <- round(tii, 27)
  # estimate cliques
    E <- t(combn(d,2)[,tii[,1] > 0])
    cliques <- maximal.cliques(graph(as.vector(t(E)), d + 1, FALSE))[-1]
    tii.scaled <- round(tii[,1,drop=FALSE] / V,27)
    res <- list(d=d, tii=tii, i1=i1, V = V, 
       tii.scaled = tii.scaled, cliques = cliques)
    class(res) <- "graphlist" # define S3 object
  return(res)
}