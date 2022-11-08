maxent.ev <- function(x,ao,ls,psi,p,q,ps,qs,d,ds,alpha)
{
  # maxent.ev by Tucker McElroy
  #   Computes maximum entropy extreme-value adjustment.
  #   Determines projection of extremes onto regular values,
  #     and gives test statistic for its magnitude.
  #   Provides shrinkage of extremes, and returns entropified data.
  # Inputs:
  #   x: time series matrix object with n rows, 
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
  
  n <- dim(x)[1]
  v <- x[,-1,drop=FALSE]
  num.reg <- dim(v)[2]
  r <- length(psi) - (1+num.reg)
  eta <- psi[-seq(1,r+1)]
  s <- frequency(x)
  deltaS <- 1
  if(ds==1) deltaS <- rep(1,s)
  deltaT <- 1
  if(d==1) deltaT <- c(1,-1)
  if(d==2) deltaT <- c(1,-2,1)
  delta <- polymult(deltaS,deltaT)
  D <- length(delta) - 1
  
  exts <- sort(union(ao,ls))
  prep <- maxent.prep(x,ao,ls,d,ds)
  x.diff <- prep[[1]]
  B.mat <- prep[[2]]
  Gamma.mat <- maxent.lik(psi,x.diff,s,p,q,ps,qs,B.mat,3)
  Gamma.alt <- diag(n)[,1:(n-length(exts)),drop=FALSE]
  Gamma.alt[(D+1):n,(D+1):(n-length(exts))] <- Gamma.mat %*% t(B.mat) %*% 
    solve(B.mat %*% Gamma.mat %*% t(B.mat))
  x.extreme <- prep[[5]] %*% x[,1,drop=FALSE]
  x.adjust <- prep[[5]] %*% (v %*% eta + prep[[3]] %*% Gamma.alt %*% 
    prep[[4]] %*% prep[[6]] %*% (x[,1,drop=FALSE]-v %*% eta))
  x.regular <- prep[[6]] %*% x[,1,drop=FALSE]
  Gamma.mse <- matrix(0,n,n)
  Gamma.mse[(D+1):n,(D+1):n] <- Gamma.mat - Gamma.mat %*% t(B.mat) %*% 
    solve(B.mat %*% Gamma.mat %*% t(B.mat)) %*% B.mat %*% Gamma.mat
  x.mse <- prep[[5]] %*% prep[[3]] %*% Gamma.mse %*% 
    t(prep[[3]]) %*% t(prep[[5]])
  
  wald <- t(x.extreme - x.adjust) %*% solve(x.mse) %*% (x.extreme - x.adjust)
  alpha <- max(alpha,1-pchisq(wald,df=r))
  kappa <- 1 - sqrt((qchisq(1-alpha,df=r))/wald[1,1])
  x.shrink <- kappa*x.adjust + (1-kappa)*x.extreme
  x.entropy <- solve(rbind(prep[[6]],prep[[5]]),c(x.regular,x.shrink)) 
  
  return(list(x.extreme,x.adjust,x.regular,x.entropy,wald))
}