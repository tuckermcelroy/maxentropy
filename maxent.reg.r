maxent.reg <- function(n,ao,ls)
{
  # maxent.reg.r by Tucker McElroy
  #   computes X regression matrix corresponding to outliers
  #
  # Inputs:
  #   n: length of time series
  #   ao: vector of elements in {1,...,n} of AO times
  #   ls: vector of elements in {1,...,n} of LS times
  
  ao.mat <- diag(n)
  ls.mat <- diag(n)
  ls.mat[lower.tri(ls.mat)] <- 1
  X.mat <- NULL
  if(length(ao) > 0)
  {
    ao <- sort(ao)
    X.mat <- ao.mat[,ao]
  }
  if(length(ls) > 0)
  {
    ls <- sort(ls)
    X.mat <- cbind(X.mat,ls.mat[,ls])
  }
  if(length(union(ao,ls))>0) 
  { 
    X.mat <- as.matrix(X.mat)[,order(union(ao,ls)),drop=FALSE]
  }
  return(X.mat)
}
