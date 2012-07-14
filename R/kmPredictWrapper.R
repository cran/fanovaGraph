kmPredictWrapper <- function(Xnew, km.object) 
  predict(object = km.object, newdata = Xnew, type = "UK", 
          se.compute = FALSE, checkNames = FALSE)$mean
