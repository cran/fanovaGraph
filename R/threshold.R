threshold <- function(delta, totalInt, d, dall) {
    ## second part of the graph estimation procedure
    ## former EstG2c
    col.index <- which(totalInt[3, ]/dall > delta)
    E <- t(totalInt[-3, col.index])  ##  cut in delta
    tii.scaled <- totalInt[3, col.index]/dall
    E.graph <- graph(as.vector(t(E)), n = d + 1, directed = FALSE)
    CL <- maximal.cliques(E.graph)[-1]
    return(list(CL = CL, E = E, tii.scaled = tii.scaled))
} 