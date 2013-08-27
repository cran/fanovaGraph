plot.graphlist <- function(x, names = NULL, i2 = NULL, layout = NULL, 
                        max.thickness=20, circle.diameter=40, ...) {
    oldPar <- par(no.readonly = TRUE)
      ## igraph changes settings at least for abline
    graphlist <- x
   
    d <- graphlist$d
    V <- graphlist$V
    i1 <- graphlist$i1[,1]
    tii <- graphlist$tii[,1]
    active <- which(tii > 0)
    tii <- tii[active]
    E <- t(combn(d,2))
    E <- E[active,]
    
    g <- graph(as.vector(t(E)) , n = d, directed = FALSE)
    if (is.null(layout)) {
        layout <- layout.fruchterman.reingold(g)
    }
    
    if (is.null(names)) {
      names <- 1:d
    }
    max.frame.width <- max.thickness
      ## thickness of the greates value
    V(g)$size <- circle.diameter
      ## vertex size (circle diameter)
  ###################################
   
  ###################################
    ## scaling
    max <- max(c(tii,0.000001))  # 0.000001 for the case everything is zero
   
    
    
    
    edge.weight.scale <- tii * max.frame.width/max
    e.col <- "darkgreen"
    
  
    v.col <- "black"
    
    
    #### plotting
    plot(g, layout = layout, vertex.shape = "circle", 
        edge.width = edge.weight.scale, vertex.frame.color = v.col,
        vertex.color = "white", edge.color = e.col, vertex.label = names,
        vertex.label.color = "black")
    
    if (!is.null(i2)) {
        if (any(i2 > tii)) {
            warning("some interaction indices larger then corresponding total interaction indices, omitted")
            i2 <- pmin(i2, tii)
        }
        edge.weight.scale2 <- i2 * max.frame.width/max
        par(new = TRUE)
        plot(g, layout = layout, vertex.shape = "circle", 
            edge.width = edge.weight.scale2, vertex.frame.color = v.col,
            vertex.color = "white", edge.color = "lightgreen", vertex.label = names,
            vertex.label.color = "black")
    }
    par(oldPar)
}