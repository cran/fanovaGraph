plotiGraph <- function(E, d, names = 1:d, i1 = NULL, tii = NULL, 
    i2 = NULL, layout = NULL) {
    oldPar <- par(no.readonly = TRUE)  
      ## igraph changes settings at least for abline
    fanovaGraph:::modify.igraph()
    g <- graph(as.vector(t(E)) - 1, n = d, directed = FALSE)
    if (is.null(layout)) {
        layout <- layout.fruchterman.reingold(g)
    }
    max.frame.width <- 20  
      ## thickness of the greates value
    V(g)$size <- 40  
      ## vertex size (circle diameter)
    v.col <- "black"
    e.col <- "black"
      ## scaling
    max <- max(c(i1, tii, 0))  # 0 for the case i1 and tii are not given
    if (is.null(i1)) 
        vertex.weight.scale <- rep(4, d) else {
        vertex.weight.scale <- i1 * max.frame.width/max
        v.col <- "darkgreen"
    }
    if (is.null(tii)) 
        edge.weight.scale <- rep(4, dim(E)[1]) else {
        edge.weight.scale <- tii * max.frame.width/max
        e.col <- "darkgreen"
    }
    #### plotting
    plot(g, layout = layout, vertex.shape = "circle2", vertex.frame.width = vertex.weight.scale, 
        edge.width = edge.weight.scale, vertex.frame.color = v.col, 
        vertex.color = "white", edge.color = e.col, vertex.label = names, 
        vertex.label.color = "black")
    
    if (!is.null(i2)) {
        # if i2 given additionally
        if (any(i2 > tii)) {
            warning("some interaction indices larger then corresponding total interaction indices, omitted")
            i2 <- pmin(i2, tii)
        }
        edge.weight.scale2 <- i2 * max.frame.width/max
        par(new = TRUE)
        plot(g, layout = layout, vertex.shape = "circle2", vertex.frame.width = vertex.weight.scale, 
            edge.width = edge.weight.scale2, vertex.frame.color = v.col, 
            vertex.color = "white", edge.color = "lightgreen", vertex.label = names, 
            vertex.label.color = "black")
    }
    par(oldPar)
} 