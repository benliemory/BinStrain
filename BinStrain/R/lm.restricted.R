lm.restricted <-
function(y,X) {
  ## X must have full column rank. deal with duplicated columns in X                                                                  
  ix=!duplicated(t(X))
  X2=X[,ix,drop=F]
  
  
  n=nrow(X2); p=ncol(X2)
  Dmat=t(X2) %*% X2
  dvec=t(X2) %*% y
  bvec=rep(0, p)
  Amat=diag(p)
  
  
  a=solve.QP(Dmat, dvec, Amat, bvec, meq=0)
  res=rep(NA, ncol(X))
  res[ix]= a[[1]]
  res
}
