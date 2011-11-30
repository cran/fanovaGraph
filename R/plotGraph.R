plotGraph <- function(E, d, names = 1:d, Maineffects = 1, 
    d.bar = 1, MAIN = "Estimated Graph") {
    ### Maineffects are the first order Sobol indices, the line width
    #   of the circles is adjusted according
    ### to them if available
    if (length(Maineffects) > 1) {
        if (length(d.bar) != dim(E)[1]) 
            d.bar <- rep(1, dim(E)[1])  # if d.bar is not provided
        MIN <- min(c(d.bar, Maineffects))
        MAX <- max(c(d.bar, Maineffects))
        ## for plotting, Maineffects and d.bar will be scaled to the
        #   interval (1, 21)
        MaineffectsS <- (Maineffects - MIN)/(MAX - MIN)  # scaling on (0,1)
        d.barS <- (d.bar - MIN)/(MAX - MIN)
        MaineffectsS <- MaineffectsS * 20 + 1
        d.barS <- d.barS * 20 + 1
    } else {
        MaineffectsS <- rep(5, d)
        d.barS <- rep(5, nrow(E))
    }
    
    positionen <- seq(0, 1, length.out = d + 1) * 2 * pi
    positionen <- positionen[-(d + 1)]
    plot(1, type = "n", xlab = "", ylab = "", xlim = c(-1.3, 1.3), 
        ylim = c(-1.3, 1.3), axes = FALSE, main = MAIN)
    
    for (i in (1:nrow(E))) {
        lines(x = c(sin(positionen[E[i, 1]]), sin(positionen[E[i, 2]])), 
            y = c(cos(positionen[E[i, 1]]), cos(positionen[E[i, 2]])), 
            lwd = d.barS[i])
    }
    points(sin(positionen), cos(positionen), type = "p", pch = 21, 
        lwd = MaineffectsS, cex = 10, bg = "white")
    text(sin(positionen), cos(positionen), labels = names, col = "black", 
        cex = 1.5)
} 