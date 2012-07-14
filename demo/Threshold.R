#################################################################
# Several graphs to support threshold decision
#################################################################

require(fanovaGraph)

### possible graph structure:
g <- list(d = 6, 
          tii = matrix(c(0.0018, 0.0265, 0.0017, 0.0277, 0.0018, 0.001, 0.028, 0.0013, 0.0212, 0.002, 0.0372, 0.0024, 0.0022, 0.0157, 0.003)),
          i1 = matrix(c(0.0901, 0.1288, 0.0683, 0.0979, 0.0882, 0.1572)),
          V = 0.8,
          cliques = list(1:6))

### plots for threshold decision:

### examine full graph

plotiGraph(g)

### Delta Jump Plot
plotDeltaJumps(g)
plotDeltaJumps(g, meanCliqueSize=TRUE)

### see graph changing as delta chages

plotGraphChange(g, fix.layout = TRUE)

### see effect of delta on graph interactively with library tcltk
plotTk(g)

### the same with library manipulate
plotManipulate(g)

