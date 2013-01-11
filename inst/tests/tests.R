test_that("estimateGraph works for different n.tot and d for all methods",{
  set.seed(1)
  expect_equivalent(c(0, 3.231419477, 0), 
      estimateGraph(f.mat=ishigami.fun, d=3, q.arg=list(min=-pi,max=pi), 
                             n.tot=100, method="FixLO")$tii[,1])
  fun <- function(x) x[,1]*x[,2]

  set.seed(1)
  expect_equivalent(0.006885859902, 
                    estimateGraph(f.mat=fun, d=2, n.tot=10000, method="FixLO")$tii[,1])
  
  set.seed(1)
  expect_equivalent(3.291240347 , 
                    estimateGraph(f.mat=ishigami.fun, d=3, q.arg=list(min=-pi,max=pi), 
                                  n.tot=1000, method="FixFast")$tii[2,1])
  set.seed(1)
  expect_equivalent(0.007030367989 , 
                    estimateGraph(f.mat=fun, d=2, n.tot=10000, method="FixFast")$tii[,1])
  
  set.seed(1)
  expect_equivalent(-0.4706963569, 
                    estimateGraph(f.mat=ishigami.fun, d=3, q.arg=list(min=-pi,max=pi), 
                                  n.tot=100, method="Sobol")$tii[2,1])
  
  set.seed(1)
  expect_equivalent(-60.94382090, 
                    estimateGraph(f.mat=ishigami.fun, d=3, q.arg=list(min=-pi,max=pi), 
                                  n.tot=100, method="RBD")$tii[2,1])
  set.seed(1)
  expect_equivalent(0.008837557797, 
                    estimateGraph(f.mat=fun, d=2, n.tot=10000, method="RBD")$tii[,1])
})

test_that("estimateGraph works for very small values",{
  fun <- function(x) ishigami.fun(x)/1000
  set.seed(1)
  expect_equivalent(0.000003231419477, 
                    estimateGraph(f.mat=fun, d=3, q.arg=list(min=-pi,max=pi), 
                                  n.tot=100, method="FixLO")$tii[2,1])
})

test_that("totalIndex works",{
  set.seed(1)
  t1 <- totalIndex(f.mat=ishigami.fun, d=3, q.arg=list(min=-pi,max=pi), n.mc=1000)
  expect_equivalent(7.129535241, t1[1])
  set.seed(1)
  t2 <- totalIndex(f.mat=sobol.fun, d=8, q.arg=list(min=0,max=1), n.mc=1000)
  expect_equivalent(0.3707765176, t2[1])
})

test_that("confint works for estimateGraph",{
  set.seed(1)
  g <- estimateGraph(f.mat=ishigami.fun, d=3, q.arg=list(min=-pi,max=pi), 
                n.tot=100, method="FixLO")
  expect_equal(c(3,4), dim(g$tii))
  expect_equivalent(2.705853788, g$tii[2,2])
})
