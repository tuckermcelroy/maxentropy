maxent.prep <- function(datareg,ao,ls,d,ds)
{
  # maxent.prep by Tucker McElroy
  #   computes matrices needed for maximum entropy calculations
  #  This version presumes D > 0
  #
  # Inputs:
  #   datareg: time series matrix object with n rows, 
  #     first column is data vector (with NAs for missing values),
  #     subsequent columns (if any) are regressors 
  #   ao: vector of elements in {1,...,n} of AO times
  #   ls: vector of elements in {1,...,n} of LS times
  #   d: order of regular differencing (altogether)
  #   ds: order of seasonal aggregation
  # Outputs:
  #   x.diff: differenced input x, in paper it's Delta %*% x
  #   B.mat: a matrix used in the differencing operation, see paper
  #   Delta.inv: inverse of tilde{Delta}_P from paper
  #   RQ.mat: product R %*% Q from paper
  #   Omega.mat: Omega matrix from paper
  #   Lambda.mat: Lambda matrix from paper
  # Requires: polymult.r, maxent.reg.r
  
  n <- dim(datareg)[1]
  s <- frequency(datareg)
  deltaS <- 1
  if(ds==1) deltaS <- rep(1,s)
  deltaT <- 1
  if(d==1) deltaT <- c(1,-1)
  if(d==2) deltaT <- c(1,-2,1)
  delta <- polymult(deltaS,deltaT)
  D <- length(delta) - 1
  
  nas <- NULL
  for(k in 1:n) { nas <- union(nas,seq(1,n)[is.na(datareg[,1])]) }
  nas <- sort(nas)
  exists <- setdiff(seq(1,n),nas)
  
  X.mat <- maxent.reg(n,ao,ls)
  if(D > 0)
  {
    t <- 1
    t.star <- n+1
    while(t <= (n-D+1))
    {
      if( (sum(X.mat[t:(t+D-1),]^2) == 0) && (sum(is.na(datareg[t:(t+D-1),1]))==0))
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
  regs <- setdiff(seq(1,n),exts)
  J.mat <- diag(n)[,exts]
  K.mat <- diag(n)[,regs]
  pi.mat <- t(cbind(K.mat,J.mat))
  Xtilde.mat <- cbind(K.mat,X.mat) %*% solve(t(pi.mat))
#  Xtilde.inv <- t(pi.mat) %*% solve(cbind(K.mat,X.mat))
  Xtilde.inv <- solve(Xtilde.mat[exists,exists,drop=FALSE])
  
#  t.hash <- t.star + sum(exts <= t.star + 1)
  t.flat <- t.star + sum(union(exts,nas) <= t.star + 1)
  q.perm <- seq(1,D)
#  if(t.hash > 0) { q.perm <-  c(seq(D+1,D+t.hash),seq(1,D)) }
  if(t.flat > 0) { q.perm <-  c(seq(D+1,D+t.flat),seq(1,D)) }
  Q.mat <- diag(n-length(union(exts,nas)))
#  Q.mat[1:(D+t.hash),1:(D+t.hash)] <- diag(D+t.hash)[,q.perm]
  Q.mat[1:(D+t.flat),1:(D+t.flat)] <- diag(D+t.flat)[,q.perm]
  
  Jtilde.mat <- diag(n)[exists,exts,drop=FALSE]
  Ktilde.mat <- diag(n)[exists,setdiff(seq(1,n),union(exts,nas))]
#  Lambda.mat <- t(K.mat) %*% Xtilde.inv 
#  Omega.mat <- t(J.mat) %*% Xtilde.inv 
  Lambda.mat <- t(Ktilde.mat) %*% Xtilde.inv 
  Omega.mat <- t(Jtilde.mat) %*% Xtilde.inv 
  temp <- Q.mat %*% Lambda.mat %*% Delta.inv[exists,,drop=FALSE]
  A.mat <- temp[-seq(1,D),1:D,drop=FALSE]
  B.mat <- temp[-seq(1,D),(D+1):n]
  x.diff <- cbind(-A.mat,diag(n-D-length(union(exts,nas)))) %*% 
    Q.mat %*% Lambda.mat %*% datareg[exists,]
  RQ.mat <- rbind(diag(n-length(union(exts,nas)))[1:D,,drop=FALSE],
                  cbind(-A.mat,diag(n-D-length(union(exts,nas))))) %*% Q.mat
  
  return(list(x.diff,B.mat,Delta.inv,RQ.mat,Omega.mat,Lambda.mat))
    
}