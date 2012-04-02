simAdd <- function(newdata, mu, parameter, covtype, Cl, iso=FALSE, eps.R = 1e-08){
  newdata <- as.matrix(newdata)
  n <- nrow(newdata)
  R <- matrix(0, ncol = n, nrow = n)
  nCl <- length(Cl)
  if (length(iso)==1) 
    iso <- rep(iso,nCl)
  parameter <- paramList2Vect(parameter, Cl, iso)
  alpha <- parameter[1:nCl]
  theta <- parameter[-(1:nCl)]
  thetalist <- Cl
  # building thetalist:
  ntemp <- 0
  nClani <- nCl - sum(iso)
  for (j in (1:nClani)) {
    thetalist[[j]] <- theta[(ntemp + 1):(ntemp + length(Cl[[j]]))]
    ntemp <- ntemp + length(Cl[[j]])
  }
  if (nClani < nCl) {
    for (j in ((nClani + 1):nCl)) {
      thetalist[[j]] <- theta[ntemp]
      ntemp <- ntemp + 1
    }
  }
  #building R
  for (j in 1:nCl) { 
    cor.str <- DiceKriging:::covStruct.create(covtype = covtype, d = length(Cl[[j]]), 
            var.names = NULL, known.covparam = "All", coef.cov = thetalist[[j]], 
            coef.var = alpha[j], iso = iso[j])
    R <- R + DiceKriging:::covMatrix(object = cor.str, X = as.matrix(newdata[, Cl[[j]]]))[[1]]
  }
  R <- R + diag(eps.R, n, n)
  whitenoise <- rnorm(n)
  cR <- chol(R)
  y <- mu + t(cR) %*% whitenoise   # t() because unlike Rasmussen cR'cR=R
  return(drop(y))
}
