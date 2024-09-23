## Script for entropy outlier project

library(seasonal)
library(devtools)
library(Rcpp)

setwd("C:\\Users\\neide\\OneDrive\\Documents\\GitHub\\sigex")
load_all(".")

setwd("C:\\Users\\neide\\OneDrive\\Documents\\GitHub\\maxentropy")

source("hend.r")
source("x11Filter.r")
source("maxent.reg.r")
source("maxent.prep.r")
source("maxent.lik.r")
source("maxent.fit.r")
source("maxent.ev.r")
source("maxent.plot.r")
source("hpsa.r")

get_start <- function(x,lag)
{
  period <- frequency(x)
  new <- time(x)[1] + lag/period
  return(c(floor(new),1+period*(new-floor(new))))
}

ubgenerator <- function(period,trunc.len,m)
{

  ceps2wold <- function(ceps,q)
  {
    m <- length(ceps)
    if(q > m) {	ceps <- c(ceps,rep(0,q-m)) }
    wold <- 1
    wolds <- wold
    for(j in 1:q)
    {
      wold <- sum(seq(1,j)*ceps[1:j]*wolds[j:1])/j
      wolds <- c(wolds,wold)
    }
    return(wolds)
  }

  half.len <- floor(period/2)
  if(length(trunc.len)==0) { trunc.len <- half.len }
  ceps <- rep(0,m)

  for(ell in 1:m)
  {
    ceps[ell] <- -2*sum(cos(2*pi*ell*seq(1,trunc.len)/period))/ell
  }
  wolds <- ceps2wold(ceps,2*trunc.len)

  return(wolds)
}

#######################
## Load test data
# Download: 5/18/2023

claims <- read.csv(file="C:\\Users\\neide\\OneDrive\\Documents\\Research\\MaxEntOutlier\\r539cy_new.csv",
                 header=TRUE,skip=1)
weekspan <- claims[,1]
claims <- claims[,2]
#claims <- read.table("InitialClaims.dat")
begin <- c(1967,1)
#end <- c(2021,14)
period <- 52
claims <- ts(claims*10^(-6),frequency=period,start=begin)
plot(claims)

setwd("C:\\Users\\neide\\OneDrive\\Documents\\Research\\MaxEntOutlier\\Figures")

## get a figure
#pdf(file="ClaimsStream.pdf",width=5, height=4)
plot(log(claims),ylab="Log Claims",xlab="Year")
dev.off()

## sub-window for 2019 to present
dataNew.ts <- ts(claims[2714:2938],frequency=period,start=c(2019,1))
plot(dataNew.ts,ylab="Claims",xlab="Year")
points(dataNew.ts)

x <- dataNew.ts
n <- length(x)

#####################
##  Identify outliers

# partial specification
#ao <- c(65,66,67,68,69,71,72,81,82)
#ls <- c(64,73,83)

# full specification
ao <- seq(65,82)
ls <- c(64,83)

r <- length(union(ao,ls))
data.ao <- rep(NA,n)
data.ao[ao] <- x[ao]
data.ls <- rep(NA,n)
data.ls[ls] <- x[ls]

## get a figure
#pdf(file="Claims2019.pdf",width=5, height=4)
plot(dataNew.ts,ylab="Claims",xlab="Year")
points(ts(data.ao,frequency=period,start=c(2019,1)),col=2)
points(ts(data.ls,frequency=period,start=c(2019,1)),col=3)
dev.off()

########################
## Specify and fit model

p <- 4
q <- 3
ps <- 1
qs <- 0
d <- 1
ds <- 0

datareg <- ts(cbind(x,seq(1,n)),start=start(x),frequency=period)
fit.mle <- maxent.fit(datareg,ao,ls,p,q,ps,qs,d,ds)
par.mle <- fit.mle[[1]]
psi.mle <- fit.mle[[2]]$par
ts.resid <- ts(c(rep(NA,r+d+ds*period),fit.mle[[3]]),
               frequency=period,start=c(2019,1))
plot(ts.resid)

## get a figure
#pdf(file="AcfResid-Claims.pdf",width=5, height=4)
acf(ts.resid[-seq(1,r+d+ds*period)],lag.max = period,main="Residual")
dev.off()

######################
## Obtain Shrinkage EV

# Modify to get forecast and aftcast
H <- 20
x.ext <- ts(c(rep(NA,H),x,rep(NA,H)),start=get_start(x,-H),frequency=period)
datareg <- ts(cbind(x.ext,seq(1,n+2*H)-H),start=start(x.ext),frequency=period)
ao_mod <- ao+H
ls_mod <- ls+H

# First do full shrinkage
alpha <- 1
out <- maxent.ev(datareg,ao_mod,ls_mod,psi.mle,p,q,ps,qs,d,ds,alpha)
x.casted <- ts(out[[5]],start=start(x.ext),frequency=period)
mse.casted <- out[[7]]
x.entropy <- ts(out[[6]],start=start(x.ext),frequency=period)
mse.entropy <- out[[8]]
1-pchisq(out[[9]],df=r)
kappa <- 1 - sqrt((qchisq(1-alpha,df=r))/out[[9]])

## get a figure
#pdf(file="MaxentFull.pdf",width=5,height=4)
#maxent.plot(x.ext,x.casted,mse.casted,option=2,top.bound=7,bot.bound=-3)
maxent.plot(x.ext,x.entropy,mse.entropy,option=2,top.bound=6,bot.bound=-2)
#plot(x.entropy,ylim=c(-2,6),ylab="Claims",xlab="Year",col=4)
#lines(x.ext,col=1)
#points(ts(x - x.entropy,frequency=period,start=c(2019,1)),col=6)
dev.off()

# Second do partial shrinkage
alpha <- .0001
out <- maxent.ev(datareg,ao_mod,ls_mod,psi.mle,p,q,ps,qs,d,ds,alpha)
x.casted <- ts(out[[5]],start=start(x.ext),frequency=period)
mse.casted <- out[[7]]
x.entropy <- ts(out[[6]],start=start(x.ext),frequency=period)
mse.entropy <- out[[8]]
1-pchisq(out[[9]],df=r)
kappa <- 1 - sqrt((qchisq(1-alpha,df=r))/out[[9]])

## get a figure
#pdf(file="MaxentHalf.pdf",width=5,height=4)
#maxent.plot(x.ext,x.casted,mse.casted,option=2,top.bound=7,bot.bound=-3)
maxent.plot(x.ext,x.entropy,mse.entropy,option=2,top.bound=6,bot.bound=-2)
#plot(x.entropy,ylim=c(-2,6),ylab="Claims",xlab="Year",col=4)
#lines(x.ext,col=1)
#points(ts(x - x.entropy,frequency=period,start=c(2019,1)),col=6)
dev.off()

##########################
## Get Seasonal Adjustment

## define trend and SA weekly filters
week.period <- 365.25/7
half.len <- floor(week.period/2)
p.seas <- 1
trend.filter <- ubgenerator(week.period,NULL,1000)
trend.filter <- trend.filter/sum(trend.filter)
#plot.ts(trend.filter)
detrend.filter <- c(rep(0,half.len),1,rep(0,half.len)) - trend.filter
seas.filter <- 0
for(j in 1:p.seas)
{
  week.periodj <- j*week.period
  half.lenj <- floor(week.periodj/2)
  seas.filterj <- ubgenerator(week.periodj,half.lenj-1,1000)
  seas.filterj <- polymult(seas.filterj,c(1,0,-1))
  seas.filter <- c(seas.filter,rep(0,length(seas.filterj)-length(seas.filter)))
  seas.filter <- seas.filter + seas.filterj
}
seas.filter <- c(rep(0,length(seas.filter)-1),seas.filter)
seas.filter <- seas.filter + rev(seas.filter)
seas.filter <- c(rep(0,(length(seas.filter)-1)/2),1,rep(0,(length(seas.filter)-1)/2)) - seas.filter/(2*p.seas+1)
#plot.ts(seas.filter)
ss.filter <- polymult(detrend.filter,seas.filter)
shift <- (length(ss.filter)-1)/2
sa.filter <- c(1,rep(0,shift)) - rev(ss.filter[1:(shift+1)])
sa.filter <- c(rev(sa.filter),sa.filter[-1])
#plot.ts(sa.filter)
n_seas <- n + 2*shift
psiMat_seas <- matrix(c(rep(0,n-1),ss.filter),nrow=1)
psiMat_sa <- matrix(c(rep(0,n-1),sa.filter),nrow=1)
for(i in 2:n)
{
  psiMat_seas <- rbind(c(psiMat_seas[1,-1],0),psiMat_seas)
  psiMat_sa <- rbind(c(psiMat_sa[1,-1],0),psiMat_sa)
}

# Modify to get forecast and aftcast
H <- shift
x.ext <- ts(c(rep(NA,H),x,rep(NA,H)),start=get_start(x,-H),frequency=period)
datareg <- ts(cbind(x.ext,seq(1,n+2*H)-H),start=start(x.ext),frequency=period)
ao_mod <- ao+H
ls_mod <- ls+H
alpha <- 1
out <- maxent.ev(datareg,ao_mod,ls_mod,psi.mle,p,q,ps,qs,d,ds,alpha)
x.casted <- ts(out[[5]],start=start(x.ext),frequency=period)
mse.casted <- out[[7]]
x.entropy <- ts(out[[6]],start=start(x.ext),frequency=period)
mse.entropy <- out[[8]]

# Seasonally adjust
x.evadjust <- x.casted - x.entropy
x.sa_evfree <- psiMat_sa %*% x.entropy
x.sa <- x.sa_evfree + x.evadjust[(H+1):(n_seas-H)]
x.seas <- psiMat_seas %*% x.entropy
x.sa <- ts(x.sa,start=start(x),frequency=period)
x.seas <- ts(x.seas,start=start(x),frequency=period)
mse.sa <- (psiMat_sa - diag(n_seas)[(H+1):(n_seas-H),]) %*% mse.entropy %*%
  t(psiMat_sa - diag(n_seas)[(H+1):(n_seas-H),])
mse.seas <- (psiMat_seas - diag(n_seas)[(H+1):(n_seas-H),]) %*% mse.entropy %*%
  t(psiMat_seas - diag(n_seas)[(H+1):(n_seas-H),])

## get a figure
#pdf(file="ClaimsSA.pdf",width=5,height=4)
maxent.plot(x.ext,x.sa,mse.sa,option=2,top.bound=6,bot.bound=0)
#maxent.plot(x.ext,x.seas,mse.seas,option=3,top.bound=5,bot.bound=-5)
dev.off()

#spec.ar(x.sa_evfree)
#acf(diff(x.sa_evfree),lag=100)



#############################################
##  Monthly Retail and Manufacturing Analyses

####################################
# I: Retail Trade and Food Services

#Source: Monthly Retail Trade and Food Services
#722: Food Services and Drinking Places: U.S. Total
#Not Seasonally Adjusted Sales - Monthly [Millions of Dollars]
#Period: 1992 to 2023
#Data Extracted on: June 5, 2023 (4:57 pm)

setwd("C:\\Users\\neide\\OneDrive\\Documents\\GitHub\\maxentropy")

food <- read.table(file="C:\\Users\\neide\\OneDrive\\Documents\\Research\\MaxEntOutlier\\FoodDrinknew.dat")
period <- 12
food <- ts(food[-seq(1,192),1],start=c(2008,1),frequency=period)
plot(food)

## get a figure
#pdf(file="FoodStream.pdf",width=5, height=4)
plot(log(food),ylab="Log Food and Drink",xlab="Year")
dev.off()

x <- ts(log(food),start=start(food),frequency=period)
n <- length(x)

#####################
##  Identify outliers

# partial specification
#ao <- c(148, 149, 150, 151, 160, 161)
#ls <- c(147, 159)

# full specification
ao <- seq(148,158)
ls <- c(147, 159)

r <- length(union(ao,ls))
data.ao <- rep(NA,n)
data.ao[ao] <- x[ao]
data.ls <- rep(NA,n)
data.ls[ls] <- x[ls]

## get a figure
#pdf(file="FoodDrink_Plot.pdf",width=5, height=4)
plot(x,ylab="Log Food and Drink",xlab="Year")
points(ts(data.ao,frequency=period,start=start(x)),col=2)
points(ts(data.ls,frequency=period,start=start(x)),col=3)
dev.off()

########################
## Specify and fit model

p <- 4
q <- 2
ps <- 0
qs <- 1
d <- 1
ds <- 1

datareg <- ts(cbind(x,seq(1,n)^d),start=start(x),frequency=period)
fit.mle <- maxent.fit(datareg,ao,ls,p,q,ps,qs,d,ds)
par.mle <- fit.mle[[1]]
psi.mle <- fit.mle[[2]]$par
ts.resid <- ts(c(rep(NA,r+d+ds*period),fit.mle[[3]]),
               frequency=period,start=start(x))
plot(ts.resid)

## get a figure
#pdf(file="AcfResid-Food.pdf",width=5, height=4)
acf(ts.resid[-seq(1,r+d+ds*period)],lag.max = 4*period,main="Residual")
#spec.ar(ts.resid[-seq(1,r+d+ds*period)])
dev.off()

############################################
## Obtain Shrinkage EV and Seasonally Adjust

# Construct SA and seasonal filters
p1 <- 3
p2 <- 5
Hendq <- 9
x11f_seas <- x11filter(p1,p2,Hendq,period,1)
x11f_sa <- x11filter(p1,p2,Hendq,period,2)
H <- (length(x11f_seas)-1)/2
n_seas <- n + 2*H
psiMat_seas <- matrix(c(rep(0,n-1),x11f_seas),nrow=1)
psiMat_sa <- matrix(c(rep(0,n-1),x11f_sa),nrow=1)
for(i in 2:n)
{
  psiMat_seas <- rbind(c(psiMat_seas[1,-1],0),psiMat_seas)
  psiMat_sa <- rbind(c(psiMat_sa[1,-1],0),psiMat_sa)
}

# Modify to get forecast and aftcast
x.ext <- ts(c(rep(NA,H),x,rep(NA,H)),start=get_start(x,-H),frequency=period)
datareg <- ts(cbind(x.ext,(seq(1,n+2*H)-H)^d),start=start(x.ext),frequency=period)
ao_mod <- ao+H
ls_mod <- ls+H

# Do full shrinkage
alpha <- 1
out <- maxent.ev(datareg,ao_mod,ls_mod,psi.mle,p,q,ps,qs,d,ds,alpha)
x.casted <- ts(out[[5]],start=start(x.ext),frequency=period)
mse.casted <- out[[7]]
x.entropy <- ts(out[[6]],start=start(x.ext),frequency=period)
mse.entropy <- out[[8]]
1-pchisq(out[[9]],df=r)
kappa <- 1 - sqrt((qchisq(1-alpha,df=r))/out[[9]])

## get a figure
#pdf(file="Maxent-Food.pdf",width=5,height=4)
#maxent.plot(x.ext,x.casted,mse.casted,prop=.05,2)
maxent.plot(x.ext,x.entropy,mse.entropy,prop=.05,2)
dev.off()

# Seasonally adjust
x.evadjust <- x.casted - x.entropy
x.sa_evfree <- psiMat_sa %*% x.entropy
x.sa <- x.sa_evfree + x.evadjust[(H+1):(n_seas-H)]
x.seas <- psiMat_seas %*% x.entropy
x.sa <- ts(x.sa,start=start(x),frequency=period)
x.seas <- ts(x.seas,start=start(x),frequency=period)
mse.sa <- (psiMat_sa - diag(n_seas)[(H+1):(n_seas-H),]) %*% mse.entropy %*%
  t(psiMat_sa - diag(n_seas)[(H+1):(n_seas-H),])
mse.seas <- (psiMat_seas - diag(n_seas)[(H+1):(n_seas-H),]) %*% mse.entropy %*%
  t(psiMat_seas - diag(n_seas)[(H+1):(n_seas-H),])

## get a figure
#pdf(file="FoodDrinkSA.pdf",width=5,height=4)
maxent.plot(x.ext,x.sa,mse.sa,prop=.05,2)
#maxent.plot(x.ext,x.seas,mse.seas,prop=1,3)
dev.off()

#spec.ar(x.sa_evfree)
#acf(diff(x.sa_evfree),lag=100)



#######################################################
# II: Manufacturers' Orders, Shipments, and Inventories

#Source: Manufacturers' Shipments, Inventories, and Orders
#Total Manufacturing: U.S. Total
#Not Seasonally Adjusted Value of Shipments [Millions of Dollars]
#Period: 2001 to 2021
#Data Extracted on: November 16, 2021 (10:27 am)

setwd("C:\\Users\\neide\\OneDrive\\Documents\\GitHub\\maxentropy")

manu <- read.table("manusio.dat")
period <- 12
x <- ts(manu,start=c(2001,1),frequency=period)
n <- length(x)



#####################
# III: Grocery Stores

#Source: Monthly Retail Trade and Food Services
#4451: Grocery Stores: U.S. Total
#Not Seasonally Adjusted Sales - Monthly [Millions of Dollars]
#Period: 1992 to 2021
#Data Extracted on: November 22, 2021 (10:35 pm)

setwd("C:\\Users\\neide\\OneDrive\\Documents\\GitHub\\maxentropy")

grocery <- read.table("grocery.dat")
period <- 12
x <- ts(grocery,c(2001,1),frequency=period)
n <- length(x)

###################
# IV: Durable Goods

# Source: Manufacturers' Shipments, Inventories, and Orders
# Durable Goods: U.S. Total
##Not Seasonally Adjusted Value of Shipments [Millions of Dollars]
#Period: 2001 to 2021
#Data Extracted on: November 22, 2021 (10:27 pm)

setwd("C:\\Users\\neide\\OneDrive\\Documents\\GitHub\\maxentropy")

durable <- read.table("durable.dat")
period <- 12
x <- ts(durable,c(2001,1),frequency=period)
n <- length(x)

############################################
# V: Building Materials and Garden Equipment

#Source: Monthly Retail Trade and Food Services
#444: Building Mat. and Garden Equip. and Supplies Dealers: U.S. Total
#Not Seasonally Adjusted Sales - Monthly [Millions of Dollars]
#Period: 2001 to 2021
#Data Extracted on: November 22, 2021 (11:13 pm)

setwd("C:\\Users\\neide\\OneDrive\\Documents\\GitHub\\maxentropy")

bldg <- read.table("building.dat")
period <- 12
x <- ts(bldg,c(2001,1),frequency=period)
n <- length(x)




#############################################
## Simulation Studies

setwd("C:\\Users\\neide\\OneDrive\\Documents\\Research\\MaxEntOutlier\\Numerical")

### Sim settings
monte <- 1000
n <- 120
transform <- "none"
aggregate <- FALSE
subseries <- 1
range <- NULL
data.ts <- sigex.prep(matrix(1,nrow=n,ncol=1),
                      transform,aggregate,subseries,range,TRUE)

### Outlier Timing: Set t0 to 60 or 108
#t0 <- 60
t0 <- 108
t1 <- 72

### Gaussian Airline
period <- 12
delta_air <- c(1,-1,rep(0,period-2),-1,1)
p <- 0
q <- 1
ps <- 0
qs <- 1
d <- 2
ds <- 1
mdl <- NULL
mdl <- sigex.add(mdl,seq(1,1),"sarma",c(p,q,ps,qs,period),NULL,"process",delta_air)
mdl <- sigex.meaninit(mdl,data.ts,0)
param <- matrix(c(.6,.6),ncol=2)
zeta <- sigex.par2zeta(param,mdl[[2]][[1]])
psi <- c(0,zeta,0)
ma_poly <- polymult(c(1,-1*zeta[1]),c(1,rep(0,11),-1*zeta[2]))
gamma <- ARMAauto(NULL,-1*ma_poly[-1],14)

### get differencing polynomial
deltaS <- 1
if(ds==1) deltaS <- rep(1,period)
deltaT <- 1
if(d==1) deltaT <- c(1,-1)
if(d==2) deltaT <- c(1,-2,1)
delta <- polymult(deltaS,deltaT)
D <- length(delta) - 1

### Outlier Scaling: Set v = 0, 1, 2, or 3 for stochastic outlier,
###   v = 0, .25, .5, or 1 for deterministic outlier.
###   Set u = 3 for stochastic outlier, u = 1 for deterministic outlier
v <- 3
u <- 3
beta_scale <- v*sqrt(gamma[1])
nu_scale <- u*sqrt(gamma[1])

### further settings
burnin <- 60
dof <- Inf
init <- rep(0,period+1)
ext.mat <- diag(n)
beta <- beta_scale
nu <- nu_scale
alpha <- 1

### Set extreme = TRUE for stochastic outlier, extreme = FALSE for deterministic outlier
extreme <- TRUE
#extreme <- FALSE

### Set levelshift = TRUE for LS, levelshift = FALSE for AO
levelshift <- TRUE
#levelshift <- FALSE

### Set misspec = TRUE to add an AO at time 72, else misspec = FALSE
misspec <- TRUE
#misspec <- FALSE

if(levelshift) { ao <- NULL; ls <- t0 } else { ao <- t0; ls <- NULL }
r <- length(union(ao,ls))
stats <- NULL
for(i in 1:monte)
{
  y <- ts(sigex.sim(psi,mdl,n,burnin,dof,init),start=1,frequency=period)
  if(extreme) nu <- rt(1,df=6)*nu_scale/sqrt(3/2)
  if(misspec) y <- y + nu*ext.mat[,t1]
  if(extreme) beta <- rt(1,df=6)*beta_scale/sqrt(3/2)
  if(levelshift) ext.mat[lower.tri(ext.mat)] <- 1
  x <- ts(y + beta*ext.mat[,t0],start=1,frequency=period)
  plot(x)
  lines(y,col=2)

  # maxent estimation
  datareg <- ts(cbind(x,seq(1,n)^d),start=start(x),frequency=period)
  fit.mle <- maxent.fit(datareg,ao,ls,p,q,ps,qs,d,ds)
  par.mle <- fit.mle[[1]]
  psi.mle <- fit.mle[[2]]$par
  out.mle <- maxent.ev(datareg,ao,ls,psi.mle,p,q,ps,qs,d,ds,alpha)
  out.true <- maxent.ev(datareg,ao,ls,psi,p,q,ps,qs,d,ds,alpha)
  stat <- c(out.mle[[9]],out.true[[9]])

  # regarima estimation
  X.mat <- maxent.reg(n,ao,ls)
  datareg <- ts(cbind(cbind(x,seq(1,n)^d),X.mat),start=start(x),frequency=period)
  fit.mle <- maxent.fit(datareg,NULL,NULL,p,q,ps,qs,d,ds)
  par.mle <- fit.mle[[1]]
  psi.mle <- fit.mle[[2]]$par
  beta.reg <- psi.mle[(length(psi.mle)-dim(X.mat)[2]+1):length(psi.mle)]
  Xdiff.mat <- filter(X.mat,delta,method="convolution",sides=1)[-seq(1,D),]
  prep <- maxent.prep(datareg,NULL,NULL,d,ds)
  x.diff <- prep[[1]]
  B.mat <- prep[[2]]
  Gamma.mle <- maxent.lik(psi.mle,x.diff,period,p,q,ps,qs,B.mat,3)
  # Gamma.true does not depend on the random beta, but we need to put
  #  some vector of length beta in the call...
  Gamma.true <- maxent.lik(c(psi,beta),x.diff,period,p,q,ps,qs,B.mat,3)

  fstat.mle <- exp(-1*psi.mle[1])*t(beta.reg) %*% t(Xdiff.mat) %*% solve(Gamma.mle) %*% Xdiff.mat %*% beta.reg
  fstat.true <- exp(-1*psi[1])*t(beta.reg) %*% t(Xdiff.mat) %*% solve(Gamma.true) %*% Xdiff.mat %*% beta.reg
  stat <- c(stat,fstat.mle,fstat.true)

  stats <- rbind(stats,stat)
  if(i %% 100 == 0) print(i)
}

sum(stats[,1] >= qchisq(.95,df=1))/monte
sum(stats[,2] >= qchisq(.95,df=1))/monte
sum(stats[,3] >= qchisq(.95,df=1))/monte
sum(stats[,4] >= qchisq(.95,df=1))/monte

# write output
#write(t(stats),file="statsM_ao_fixed_middle_scale1.txt",ncol=4)
#write(t(stats),file="statsM_ls_fixed_middle_scale1.txt",ncol=4)
#write(t(stats),file="statsM_ao_random_middle_scale3.txt",ncol=4)
#write(t(stats),file="statsM_ls_random_middle_scale3.txt",ncol=4)
#write(t(stats),file="statsM_ao_fixed_end_scale1.txt",ncol=4)
#write(t(stats),file="statsM_ls_fixed_end_scale0.txt",ncol=4)
#write(t(stats),file="statsM_ao_random_end_scale3.txt",ncol=4)
#write(t(stats),file="statsM_ls_random_end_scale3.txt",ncol=4)

