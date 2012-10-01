test_that("estimateGraph works for different N and d for all methods",{
  set.seed(1)
  expect_equivalent(c(0, 3.2314, 0), 
                    estimateGraph(f.mat=ishigami.fun, d=3, q.arg=list(min=-pi,max=pi), 
                             N=100, method="FixLO")$tii[,1])
  fun <- function(x) x[,1]*x[,2]
  set.seed(1)
  expect_equivalent(0.0069, estimateGraph(f.mat=fun, d=2, N=10000, method="FixLO")$tii[,1])
  
  set.seed(1)
  expect_equivalent(c(0.0027, 3.2912, 0.1045), 
                    estimateGraph(f.mat=ishigami.fun, d=3, q.arg=list(min=-pi,max=pi), 
                                  N=1000, method="FixFast")$tii[,1])
  set.seed(1)
  expect_equivalent(0.007, estimateGraph(f.mat=fun, d=2, N=10000, method="FixFast")$tii[,1])
  
  set.seed(1)
  expect_equivalent(c(0.7786, -0.4707, -4.6113), 
                    estimateGraph(f.mat=ishigami.fun, d=3, q.arg=list(min=-pi,max=pi), 
                                  N=100, method="Sobol")$tii[,1])
  
  set.seed(1)
  expect_equivalent(c(-37.7240, -60.9438, -43.1977), 
                    estimateGraph(f.mat=ishigami.fun, d=3, q.arg=list(min=-pi,max=pi), 
                                  N=100, method="RBD")$tii[,1])
  set.seed(1)
  expect_equivalent(0.0088, estimateGraph(f.mat=fun, d=2, N=10000, method="RBD")$tii[,1])
})



