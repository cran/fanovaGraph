#################################################################
# Several graphs to support threshold decision
#################################################################

require(fanovaGraph)

### example fanovaGraph

fun <- function(x) {x[,1]*x[,2]*x[,3] + x[,4]*x[,5]*x[,6] + 
  0.7*x[,1]*x[,2]*x[,3]*x[,4]*x[,5]*x[,6]}

g <- estimateGraph(fun, d=6, q.arg=list(min=-1,max=1), n.tot=20000, method="LiuOwen")

### plots for threshold decision:

### examine full graph

plot(g)

### Delta Jump Plot
plotDeltaJumps(g)
plotDeltaJumps(g, mean.clique.size=TRUE)

### see graph changing as delta chages

plotGraphChange(g, fix.layout = TRUE)

### see effect of delta on graph interactively with library tcltk
plotTk(g)

### the same with library manipulate
plotManipulate(g)

