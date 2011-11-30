plotGraphChange <- function(totalInt, d, dall, fix.layout = TRUE, 
    delta.layout = 0.01) {
    oldPar <- par(c("usr", "xpd"))
    par(usr = c(-1.161, 1.161, -1.17826, 1.17826), xpd = TRUE)
    graph0 <- threshold(delta = 0, totalInt = totalInt, d = d, dall = dall)
    delta <- c(0, sort(graph0$tii.scaled), 1)
    delta.max <- delta[length(delta) - 1]
    n.CL <- c()
    
    # position of the vertices (defined by igraphs at the graph on
    #   delta.layout)
    E.layout <- threshold(delta = delta.layout, totalInt = totalInt, 
        d = d, dall = dall)$E
    g.layout <- graph(as.vector(t(E.layout)) - 1, n = d, directed = FALSE)
    layout <- layout.fruchterman.reingold(g.layout)
    
    devAskNewPage(ask = TRUE)
    for (i in 1:(length(delta) - 1)) {
        graph <- threshold(delta[i], totalInt, d, dall)
        if (fix.layout == TRUE) {
            plotiGraph(graph$E, d, tii = graph$tii.scaled, layout = layout)
        } else {
            plotiGraph(graph$E, d = d, tii = graph$tii.scaled)
        }
        E.graph <- graph(as.vector(t(graph$E)), n = d + 1, directed = FALSE)
        CL <- maximal.cliques(E.graph)[-1]
        n.CL[i] <- length(CL)
        text(0, -2, paste("delta = [", round(delta[i], 4), ",", round(delta[i + 
            1], 4), ")"))
        title(main = paste("number of cliques =", length(CL)))
        rect(delta[i] * 2/delta.max - 1, -1.55, delta[i + 1] * 2/delta.max - 
            1, -1.45, density = NA, col = "blue")
        lines(c(-1, 1.1), c(-1.5, -1.5))
        points(c(-1, 0.01 * 2/delta.max - 1, 1), rep(-1.5, 3), pch = "I")
        text(c(-1, 0.01 * 2/delta.max - 1, 1), c(-1.6, -1.6), c(delta[1], 
            0.01, round(delta.max, 4)), pos = 1, cex = 0.6)
        if (delta.max > 0.1) {
            points(0.1 * 2/delta.max - 1, -1.5, pch = "I")
            text(0.1 * 2/delta.max - 1, -1.6, 0.1, pos = 1, cex = 0.6)
        }
    }
    devAskNewPage(ask = FALSE)
    par(oldPar)
} 