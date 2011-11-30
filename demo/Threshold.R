#################################################################
# Several graphs to support threshold decision
#################################################################

require(fanovaGraph)

### possible estimates for function A:
d <- 6
totalInt <- rbind(rep(1:5, 5:1), c(2:6, 3:6, 4:6, 5:6, 6), 
    c(0.0018, 0.0265, 0.0017, 0.0277, 0.0018, 0.001, 0.028, 0.0013, 0.0212, 
        0.002, 0.0372, 0.0024, 0.0022, 0.0157, 0.003))
dall <- 0.75

### plots for threshold decision:

### graph without thresholding

graph0 <- threshold(delta = 0, totalInt = totalInt, d = d, 
    dall = dall)
plotiGraph(E = graph0$E, d = d, tii = graph0$tii.scaled)

### see graph changing as delta chages

plotGraphChange(totalInt = totalInt, d = d, dall = dall, fix.layout = TRUE)

### see effect of delta on graph interactively with library tcltk
plotTk(totalInt, d, dall = dall)

### the same with library manipulate
plotManipulate(totalInt, d, dall = dall)

### Delta Step Plot
plotDeltaSteps(totalInt = totalInt, d = d, dall = dall, interval = c(0, 
    1), meanCliqueSize = FALSE)

plotDeltaSteps(totalInt = totalInt, d = d, dall = dall, interval = c(0, 
    1), meanCliqueSize = TRUE) 