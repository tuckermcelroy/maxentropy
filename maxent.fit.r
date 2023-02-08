maxent.fit <- function(datareg,ao,ls,p,q,ps,qs,d,ds)
{
  # maxent.fit by Tucker McElroy
  #   fits SARIMA model using maximum entropy framework
  #  This version presumes D > 0
  #
  # Inputs:
  #   datareg: time series matrix object with n rows, 
  #     first column is data vector,
  #     subsequent columns (if any) are regressors 
  #   ao: vector of elements in {1,...,n} of AO times
  #   ls: vector of elements in {1,...,n} of LS times
  #   p: order of regular AR
  #   q: order of regular MA
  #   ps: order of seasonal AR
  #   qs: order of seasonal MA
  #   d: order of regular differencing (altogether)
  #   ds: order of seasonal aggregation
  # Outputs:
  #   param: contains the fitted parameters:
  #	    first p components are the AR parameters 
  #	    second q components are the MA parameters
  #	    third ps components are the seasonal AR parameters
  #	    fourth qs components are the seasonal MA parameters
  #	    the next component is the logged innovation variance
  #     the last components are for regression effects.
  #   ent.mle: value of the divergence at mle
  #   x.resid: residuals
  # Requires: maxent.lik.r, maxent.prep.r

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
  
  n <- dim(datareg)[1]
  v <- datareg[,-1,drop=FALSE]
  num.reg <- dim(v)[2]
  s <- frequency(datareg)
  prep <- maxent.prep(datareg,ao,ls,d,ds)
  x.diff <- prep[[1]]
  B.mat <- prep[[2]]
#  x <- datareg[,1]
  
  # fit the SARIMA model
  r <- p+q+ps+qs
  psi.init <- rep(0,r+1+num.reg)
  ent.mle <- nlminb(start=psi.init,objective=maxent.lik,x=x.diff,
                    s=s,p=p,q=q,ps=ps,qs=qs,B.mat=B.mat,outFlag=1)
  psi <- ent.mle$par
  
  # get SARMA parameters
  # input parameters defined, preliminary calculations	
  if (p==0) { ar <- NULL } else { ar <- psi2phi(psi[1:p]) }
  if (q==0) { ma <- NULL } else { ma <- psi2phi(psi[(p+1):(p+q)]) }
  if (ps==0) { ars <- NULL } else { ars <- psi2phi(psi[(p+q+1):(p+q+ps)]) }	
  if (qs==0) { mas <- NULL } else { mas <- psi2phi(psi[(p+q+ps+1):(p+q+ps+qs)]) }
  param <- c(ar,ma,ars,mas,exp(psi[r+1]),psi[-seq(1,r+1)])
  
  x.resid <- ts(maxent.lik(psi,x.diff,s,p,q,ps,qs,B.mat,2),names="Residual",frequency=s)
  return(list(param,ent.mle,x.resid))
  
}