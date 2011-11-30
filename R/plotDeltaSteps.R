plotDeltaSteps <- function(totalInt, d, dall, interval = c(0, 
    1), meanCliqueSize = FALSE) {
    op <- par(no.readonly = TRUE)
    par(mfrow = c(2, 1), mar = c(0, 4, 1, 0), oma = c(6, 0.5, 2, 1), 
        mgp = c(2.5, 1, 0))
    o <- order(totalInt[3, ], decreasing = FALSE)
    totalIntSort <- totalInt[, o]
    totalIntSortNorm <- rbind(totalIntSort[1:2, ], totalIntSort[3, 
        ]/dall)
    delta <- c(0, totalIntSortNorm[3, ], 1)
    delta <- delta[delta >= interval[1] & delta <= interval[2]]
    n.CL <- c()
    s.CL <- c()
    for (i in 1:length(delta)) {
        E <- t(totalIntSortNorm[-3, which(totalIntSortNorm[3, ] > delta[i])])
        E.graph <- graph(as.vector(t(E)), n = d + 1, directed = FALSE)
        CL <- maximal.cliques(E.graph)[-1]
        n.CL[i] <- length(CL)
        s.CL[i] <- mean(sapply(CL, length))
    }
    plot(1:length(delta), n.CL, type = "s", ylab = "number of cliques", 
        xaxt = "n", ylim = c(0, max(n.CL)))
    if (meanCliqueSize == TRUE) {
        lines(1:length(delta), s.CL, type = "s", lty = 3, lwd = 2)
        legend("topright", c("number of cliques", "mean clique size"), 
            lty = c(1, 3))
    }
    deltaCut <- c(delta[-length(delta)], 4/3 * delta[length(delta) - 
        1])
    plot(1:length(delta), deltaCut, type = "s", ylab = "delta", xaxt = "n")
    axis(1, at = 1:length(delta), label = round(delta, 4), las = 2, 
        outer = TRUE)
    title(xlab = "delta steps", outer = TRUE, mgp = c(4, 1, 0))
    title("Delta Step Plot", outer = TRUE)
    par(op)
} 