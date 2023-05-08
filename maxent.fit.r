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
  #	    the initial component is the logged innovation variance; after this --
  #	    first p components are the AR parameters 
  #	    second q components are the MA parameters
  #	    third ps components are the seasonal AR parameters
  #	    fourth qs components are the seasonal MA parameters
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
  rr <- p+q+ps+qs
  psi.init <- rep(0,rr+1+num.reg)
  ent.mle <- nlminb(start=psi.init,objective=maxent.lik,x.diff=x.diff,
                    s=s,p=p,q=q,ps=ps,qs=qs,B.mat=B.mat,outFlag=1)
  psi <- ent.mle$par
  
  # get SARMA parameters
  # input parameters defined, preliminary calculations	
  if (p==0) { ar <- NULL } else { ar <- psi2phi(psi[2:(p+1)]) }
  if (q==0) { ma <- NULL } else { ma <- psi2phi(psi[(p+2):(p+q+1)]) }
  if (ps==0) { ars <- NULL } else { ars <- psi2phi(psi[(p+q+2):(p+q+ps+1)]) }	
  if (qs==0) { mas <- NULL } else { mas <- psi2phi(psi[(p+q+ps+2):(p+q+ps+qs+1)]) }
  param <- c(exp(psi[1]),ar,ma,ars,mas,psi[-seq(1,rr+1)])
  
  x.resid <- ts(maxent.lik(psi,x.diff,s,p,q,ps,qs,B.mat,2),names="Residual",frequency=s)
  return(list(param,ent.mle,x.resid))
  
}