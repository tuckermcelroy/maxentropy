maxent.prep <- function(x,ao,ls,d,ds)
{
  # maxent.prep by Tucker McElroy
  #   computes matrices needed for maximum entropy calculations
  #  This version presumes D > 0
  #
  # Inputs:
  #   x: time series of length n
  #   ao: vector of elements in {1,...,n} of AO times
  #   ls: vector of elements in {1,...,n} of LS times
  #   d: order of regular differencing (altogether)
  #   ds: order of seasonal aggregation
  #
  # Requires: polymult.r, maxent.reg.r
  
  n <- length(x)
  s <- frequency(x)
  deltaS <- 1
  if(ds==1) deltaS <- rep(1,s)
  deltaT <- 1
  if(d==1) deltaT <- c(1,-1)
  if(d==2) deltaT <- c(1,-2,1)
  delta <- polymult(deltaS,deltaT)
  D <- length(delta) - 1
  
  X.mat <- maxent.reg(n,ao,ls)
  if(D > 0)
  {
    t <- 1
    t.star <- n+1
    while(t <= (n-D+1))
    {
      if(sum(X.mat[t:(t+D-1),]^2) == 0) 
      { 
        t.star <- t-1 
        t <- n
      }
      t <- t+1
    }
  }
  if(t.star <= n)
  {
    p.perm <- seq(t.star+1,t.star+D)
  }
  P.mat <- array(0,c(D,n))
  P.mat[,p.perm] <- diag(D)
  Delta.mat <- toeplitz(c(delta,rep(0,n-D-1)))
  Delta.mat[lower.tri(Delta.mat)] <- 0
  Delta.mat <- Delta.mat[1:(n-D),]
  Delta.mat <- rbind(P.mat,Delta.mat)
  Delta.inv <- solve(Delta.mat)
  
  exts <- sort(union(ao,ls))
  regs <- setdiff(seq(1,n),sort(union(ao,ls)))
  J.mat <- diag(n)[,exts]
  K.mat <- diag(n)[,regs]
  pi.mat <- t(cbind(K.mat,J.mat))
#  Xtilde.mat <- cbind(K.mat,X.mat) %*% solve(t(pi.mat))
  Xtilde.inv <- t(pi.mat) %*% solve(cbind(K.mat,X.mat))
  
  t.hash <- t.star + sum(exts <= t.star + 1)
  q.perm <- seq(1,D)
  if(t.hash > 0) { q.perm <-  c(seq(D+1,D+t.hash),seq(1,D)) }
  Q.mat <- diag(n-length(exts))
  Q.mat[1:(D+t.hash),1:(D+t.hash)] <- diag(D+t.hash)[,q.perm]
    
  Lambda.mat <- t(K.mat) %*% Xtilde.inv 
  Omega.mat <- t(J.mat) %*% Xtilde.inv 
  temp <- Q.mat %*% Lambda.mat %*% Delta.inv
  A.mat <- temp[-seq(1,D),1:D,drop=FALSE]
  B.mat <- temp[-seq(1,D),(D+1):n]
  x.diff <- cbind(-A.mat,diag(n-D-length(exts))) %*% Q.mat %*% t(K.mat) %*% Xtilde.inv %*% x
  RQ.mat <- rbind(diag(n-length(exts))[1:D,,drop=FALSE],
                  cbind(-A.mat,diag(n-D-length(exts)))) %*% Q.mat
  
  return(list(x.diff,B.mat,Delta.inv,RQ.mat,Omega.mat,Lambda.mat))
    
}