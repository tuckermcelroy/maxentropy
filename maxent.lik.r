maxent.lik <- function(psi,x,s,p,q,ps,qs,B.mat,outFlag)
{
  # maxent.lik by Tucker McElroy
  #	Computes Gaussian likelihood for a SARMA model psi fitted to stationary data x
  #	The model is SARMA(p,q,ps,qs) in Box-Jenkins notation
  #   psi contains the parameters:
  #	    first p components of psi are the AR parameters 
  #	    second q components of psi are the MA parameters
  #	    third ps components of psi are the seasonal AR parameters
  #	    fourth qs components of psi are the seasonal MA parameters
  #	    the next component is the logged innovation variance
  #     the last components are for regression effects.
  #   s is the seasonal period
  #   x is a matrix object,
  #     first column is "differenced data" defined via 
  #     B.mat %*% W = Delta %*% X in paper,
  #     subsequent columns (if any) are differenced regressors defined via
  #     Delta %*% V in paper, where V is initial given regressor matrix.
  #     For a mean effect, include a column of ones in x.
  #     Note: presumes no collinearity or zero vectors!
  #   B.mat is a matrix used in the differencing operation, see paper
  #	Outputs:
  #		value of the log likelihood excluding irrelevant constants (outFlag = 1)
  #		also yields time series residuals as alternative output (outFlag = 2)
  #  Utilizes: polymult.r, ARMAauto.r
  
  psi2phi <- function(psi)
  {
    p <- length(psi)
    pacfs <- (exp(psi)-1)/(exp(psi)+1)
    if(p==1)
    {
      phi <- pacfs
    } else
    {
      phi <- as.vector(pacfs[1])
      for(j in 2:p)
      {
        A.mat <- diag(j-1) - pacfs[j]*diag(j-1)[,(j-1):1,drop=FALSE]
        phi <- A.mat %*% phi
        phi <- c(phi,pacfs[j])
      }
    }
    return(phi)
  }
  
  # error control code
  for(k in 1:length(psi))
  {
    if (psi[k] == "NaN") psi[k] <- 0
  }
  
  v <- x[,-1,drop=FALSE]
  num.reg <- dim(v)[2]
  
  # input parameters defined, preliminary calculations	
  r <- length(psi) - (1+num.reg)
  eta <- psi[-seq(1,r+1)]
  #mu <- psi[r+2]
  if (p==0) { ar <- NULL } else { ar <- psi2phi(psi[1:p]) }
  if (q==0) { ma <- NULL } else { ma <- psi2phi(psi[(p+1):(p+q)]) }
  if (ps==0) { ars <- NULL } else { ars <- psi2phi(psi[(p+q+1):(p+q+ps)]) }	
  if (qs==0) { mas <- NULL } else { mas <- psi2phi(psi[(p+q+ps+1):(p+q+ps+qs)]) }
  
  arpoly <- NULL
  if (ps > 0) arpoly <- as.vector(matrix(c(rep(0,s-1),1),nrow=s,ncol=1) %*% 
                                    matrix(ars,nrow=1,ncol=ps))
  arpoly <- c(1,-1*arpoly)
  arpoly <- polymult(c(1,-1*ar),arpoly)
  mapoly <- NULL
  if (qs > 0) mapoly <- as.vector(matrix(c(rep(0,s-1),1),nrow=s,ncol=1) %*% 
                                    matrix(mas,nrow=1,ncol=qs))
  mapoly <- c(1,-1*mapoly)
  mapoly <- polymult(c(1,-1*ma),mapoly)
  
  m <- dim(B.mat)[2]
  x.acf <- ARMAauto(ar = -1*Re(arpoly[-1]), ma = Re(mapoly[-1]),lag.max=m)[1:m]
  Gamma.mat <- toeplitz(x.acf)
  lik.mat <- B.mat %*% Gamma.mat %*% t(B.mat)
  C.mat <- t(chol(lik.mat))
  z <- solve(C.mat,x[,1,drop=FALSE]-v %*% eta)
  Q.form <- sum(z^2)
  logdet <- 2*sum(log(diag(C.mat)))
  lik <- Q.form/exp(psi[r+1]) + logdet + (dim(B.mat)[1])*psi[r+1]
  
#  print(lik)
   
  # compute time series residuals
  sqrlik <- svd(lik.mat)
  sqrlik <- sqrlik$u %*% solve(diag(sqrt(sqrlik$d))) %*% t(sqrlik$u)
  ts.resid <- sqrlik %*% (x[,1,drop=FALSE]-v %*% eta)
  
  if (outFlag == 1) out <- lik
  if (outFlag == 2) out <- ts.resid
  if (outFlag == 3) out <- Gamma.mat
  if (outFlag == 4) out <- Q.form/exp(psi[r+1])
  
  return(out)
}	