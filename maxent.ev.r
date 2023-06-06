maxent.ev <- function(datareg,ao,ls,psi,p,q,ps,qs,d,ds,alpha)
{
  # maxent.ev by Tucker McElroy
  #   Computes maximum entropy extreme-value adjustment.
  #   Determines projection of extremes onto regular values,
  #     and gives test statistic for its magnitude.
  #   Provides shrinkage of extremes, and returns entropified data.
  # Inputs:
  #   datareg: time series matrix object with n rows, 
  #     first column is data vector,
  #     subsequent columns (if any) are regressors 
  #   ao: vector of elements in {1,...,n} of AO times
  #   ls: vector of elements in {1,...,n} of LS times
  #   psi: vector of pre-parameters for SARIMA model
  #   p: order of regular AR
  #   q: order of regular MA
  #   ps: order of seasonal AR
  #   qs: order of seasonal MA
  #   d: order of regular differencing (altogether)
  #   ds: order of seasonal aggregation
  #   alpha: controls shrinkage; 
  #     set to 1 for full shrinkage (use conditional expectation),
  #     set to 0 for no shrinkage (keep raw data)
  # Outputs:
  #   x.extreme: extreme-values from observed sample
  #   x.adjust: adjusted extreme-values
  #   x.regular: regular values from observed sample
  #   x.nas: casted missing values
  #   x.casted: sample corrected for missing values
  #   x.entropy: sample corrected for extremes and missing values
  #   mse.casted: MSE of sample corrected for missing values
  #   mse.entropy: MSE of sample corrected for extremes and missing values
  #   wald: wald statistic for extreme values
  # Requires: maxent.lik.r, maxent.prep.r
  
  n <- dim(datareg)[1]
  v <- datareg[,-1,drop=FALSE]
  num.reg <- dim(v)[2]
#  r <- length(psi) - (1+num.reg)
  rr <- p+q+ps+qs
  eta <- psi[-seq(1,rr+1)]
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
  
  exts <- sort(union(ao,ls))
  r <- length(exts)
  prep <- maxent.prep(datareg,ao,ls,d,ds)
  x.diff <- prep[[1]]
  B.mat <- prep[[2]]
  Gamma.mat <- maxent.lik(psi,x.diff,s,p,q,ps,qs,B.mat,3)
  Gamma.mat <- exp(psi[1])*Gamma.mat
  Gamma.alt <- diag(n)[,1:(n-length(union(exts,nas))),drop=FALSE]
  Gamma.alt[(D+1):n,(D+1):(n-length(union(exts,nas)))] <- Gamma.mat %*% 
    t(B.mat) %*% solve(B.mat %*% Gamma.mat %*% t(B.mat))
  x.extreme <- prep[[5]] %*% datareg[exists,1,drop=FALSE]
  x.regular <- prep[[6]] %*% datareg[exists,1,drop=FALSE]
  x.proj <- (v %*% eta + prep[[3]] %*% Gamma.alt %*% prep[[4]] %*% 
      prep[[6]] %*% (datareg[exists,1,drop=FALSE]-v[exists,,drop=FALSE] %*% eta))
  x.adjust <- prep[[5]] %*% x.proj[exists,,drop=FALSE]
  if(length(nas) > 0) { x.nas <- x.proj[nas,,drop=FALSE] } else { x.nas <- NULL }
  Gamma.mse <- matrix(0,n,n)
  Gamma.mse[(D+1):n,(D+1):n] <- Gamma.mat - Gamma.mat %*% t(B.mat) %*% 
    solve(B.mat %*% Gamma.mat %*% t(B.mat)) %*% B.mat %*% Gamma.mat
  mse.entropy <- prep[[3]] %*% Gamma.mse %*% t(prep[[3]])
  x.mse <- prep[[5]] %*% mse.entropy[exists,exists,drop=FALSE] %*% t(prep[[5]])
  
  wald <- t(x.extreme - x.adjust) %*% solve(x.mse) %*% (x.extreme - x.adjust)
  alpha <- max(alpha,1-pchisq(wald,df=r))
  kappa <- 1 - sqrt((qchisq(1-alpha,df=r))/wald[1,1])
  x.shrink <- kappa*x.adjust + (1-kappa)*x.extreme
  
  if(length(nas) > 0) { I.mat <- diag(n)[,nas,drop=FALSE] } else { I.mat <- NULL }
  H.mat <- diag(n)[,exists,drop=FALSE]
  Trans.mat <- t(cbind(I.mat,H.mat %*% t(rbind(prep[[6]],prep[[5]]))))
  Trans.inv <- solve(Trans.mat)
  x.casted <- Trans.inv %*% c(x.nas,x.regular,x.extreme)
  x.entropy <- Trans.inv %*% c(x.nas,x.regular,x.shrink)
  mse.casted <- 0*diag(n)
  if(length(nas) > 0) 
  {
    mse.casted[1:length(nas),1:length(nas)] <- mse.entropy[nas,nas,drop=FALSE]
    mse.casted <- Trans.inv %*% mse.casted %*% t(Trans.inv) 
  }
#  mse.evadjust <- 0*diag(n)
#  if(length(exts) > 0)
#  {
#    mse.evadjust[(n-length(exts)+1):n,(n-length(exts)+1):n] <- x.mse
#    mse.evadjust <- Trans.inv %*% mse.evadjust %*% t(Trans.inv) 
#  }
  
  return(list(x.extreme,x.adjust,x.regular,x.nas,x.casted,x.entropy,
              mse.casted,mse.entropy,wald))
}