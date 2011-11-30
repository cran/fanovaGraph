plotTk <- function(totalInt, d, dall, delta.layout = 0.01) {
    require(tcltk)
    E.layout <- threshold(delta.layout, totalInt, d, dall)$E
    g.layout <- graph(as.vector(t(E.layout)) - 1, n = d, directed = FALSE)
    layout <- layout.fruchterman.reingold(g.layout)
    max.delta <- max(totalInt[3, ]/dall + 0.001)
    
    plotting <- function(delta) {
        graph <- threshold(delta, totalInt, d, dall)
        plotiGraph(graph$E, d, tii = graph$tii.scaled, layout = layout)
        E.graph <- graph(as.vector(t(graph$E)), n = d + 1, directed = FALSE)
        CL <- maximal.cliques(E.graph)[-1]
        n.CL <- length(CL)
        title(main = paste("delta =", round(delta, 5)))
        title(sub = paste("number of cliques =", length(CL)))
    }
    
    variable <- tclVar(delta.layout)  # start value
    refresh <- function(...) {
        # function for every change of input
        delta <- as.numeric(tclvalue(variable))
        plotting(delta)
    }
    m <- tktoplevel()  # Erstellen des Eingabefensters
    tkwm.title(m, "input window")
    cutFrame <- tkframe(m)
    cutSlider <- tkscale(cutFrame, command = refresh, from = 0, to = max.delta, 
        resolution = 5e-04, orient = "horiz", variable = variable)
    tkpack(tklabel(cutFrame, text = "Delta:"), side = "left")
    tkpack(cutFrame, cutSlider, side = "bottom")
}


plotManipulate <- function(totalInt, d, dall, delta.layout = 0.01) {
    require(manipulate)
    E.layout <- threshold(delta.layout, totalInt, d = d, dall = dall)$E
    g.layout <- graph(as.vector(t(E.layout)) - 1, n = d, directed = FALSE)
    layout <- layout.fruchterman.reingold(g.layout)
    max.delta <- max(totalInt[3, ]/dall + 0.001)
    
    plotting <- function(delta) {
        graph <- threshold(delta, totalInt, d, dall)
        plotiGraph(graph$E, d, tii = graph$tii.scaled, layout = layout)
        E.graph <- graph(as.vector(t(graph$E)), n = d + 1, directed = FALSE)
        CL <- maximal.cliques(E.graph)[-1]
        n.CL <- length(CL)
        title(main = paste("delta =", round(delta, 5)))
        title(sub = paste("number of cliques =", length(CL)))
    }
    manipulate(plotting(delta), delta = slider(0, max.delta, initial = delta.layout, 
        step = 5e-04))
} 