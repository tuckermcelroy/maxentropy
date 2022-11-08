## Script for entropy outlier project

library(seasonal)

setwd("C:\\Users\\neide\\OneDrive\\Documents\\GitHub\\maxentropy")

source("ARMAauto.r")
source("polymult.r")
source("maxent.reg.r")
source("maxent.prep.r")
source("maxent.lik.r")
source("maxent.fit.r")
source("maxent.ev.r")
source("hpsa.r")

#######################
## Load test data
# Download: 6/11/2021 at 3:44 ET

setwd("C:\\Users\\neide\\OneDrive\\Documents\\GitHub\\maxentropy")

claims <- read.table("InitialClaims.dat")
begin <- c(1967,1) 
end <- c(2021,19)
period <- 52
claims <- ts(claims*10^(-6),frequency=period,start=begin)
plot(claims)

setwd("C:\\Users\\neide\\OneDrive\\Documents\\Research\\MaxEntOutlier\\Figures")

## get a figure
pdf(file="ClaimsStream.pdf",width=5, height=4)
plot(log(claims),ylab="Log Claims",xlab="Year")
dev.off()

## sub-window for 2020-2021
dataNew.ts <- ts(claims[2714:2836,1],frequency=period,start=c(2019,1))
plot(dataNew.ts,ylab="Claims",xlab="Year")  
points(dataNew.ts)

x <- dataNew.ts
n <- length(x)

#####################
##  Identify outliers

#ao <- c(65,66,67,68,69,70,71,72)
ao <- c(65,66,67,68,69,71,72)
ls <- c(64,73)
r <- length(union(ao,ls))
data.ao <- rep(NA,n)
data.ao[ao] <- x[ao]
data.ls <- rep(NA,n)
data.ls[ls] <- x[ls]

## get a figure
pdf(file="Claims2021.pdf",width=5, height=4)
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
 
fit.mle <- maxent.fit(x,ao,ls,p,q,ps,qs,d,ds)
par.mle <- fit.mle[[1]]
psi.mle <- fit.mle[[2]]$par
ts.resid <- ts(c(rep(NA,r+d+ds*period),fit.mle[[3]]),frequency=period,start=c(2019,1))
plot(ts.resid)

## get a figure
pdf(file="AcfResid.pdf",width=5, height=4)
acf(ts.resid[-seq(1,r+d+ds*period)],lag.max = period,main="Residual")
dev.off()

######################
## Obtain Shrinkage EV

# First do full shrinkage
alpha <- 1
out <- maxent.ev(x,ao,ls,psi.mle,p,q,ps,qs,d,ds,alpha)
1-pchisq(out[[5]],df=r)
kappa <- 1 - sqrt((qchisq(1-alpha,df=r))/out[[5]])
x.entropy <- out[[4]]

## get a figure
pdf(file="MaxentFull.pdf",width=5,height=4)
plot(x,ylim=c(-2,6),ylab="Claims",xlab="Year")
lines(ts(x.entropy,frequency=period,start=c(2019,1)),col=4)
#points(ts(x - x.entropy,frequency=period,start=c(2019,1)),col=6)
dev.off() 

# Second do partial shrinkage
alpha <- .99
out <- maxent.ev(x,ao,ls,psi.mle,p,q,ps,qs,d,ds,alpha)
1-pchisq(out[[5]],df=r)
kappa <- 1 - sqrt((qchisq(1-alpha,df=r))/out[[5]])
x.entropy <- out[[4]]

## get a figure
pdf(file="MaxentHalf.pdf",width=5,height=4)
plot(x,ylim=c(-2,6),ylab="Claims",xlab="Year")
lines(ts(x.entropy,frequency=period,start=c(2019,1)),col=4)
#points(ts(x - x.entropy,frequency=period,start=c(2019,1)),col=6)
dev.off() 
  
#############################################
##  Monthly Retail and Manufacturing Analyses

####################################
# I: Retail Trade and Food Services

#Source: Monthly Retail Trade and Food Services	
#722: Food Services and Drinking Places: U.S. Total	
#Not Seasonally Adjusted Sales - Monthly [Millions of Dollars]	
#Period: 2001 to 2021
#Data Extracted on: November 22, 2021 (10:32 pm)	

setwd("C:\\Users\\neide\\OneDrive\\Documents\\GitHub\\maxentropy")

food <- read.table("fooddrink.dat")
period <- 12
x <- ts(food,start=c(2001,1),frequency=period)
n <- length(x)
plot(x)
x.sub <- ts(x[1:228],start=start(x),frequency=period)

#m0 <- seas(x = x.sub, transform.function = "auto", x11 = "",outlier.types = "none")
#out(m0)


data.ts <- log(x)
ao <- c(seq(232,235),244,245)
ls <- c(231,243)
r <- length(union(ao,ls))
data.ao <- rep(NA,n)
data.ao[ao] <- data.ts[ao]
data.ls <- rep(NA,n)
data.ls[ls] <- data.ts[ls]

plot(data.ts,ylab="Food and Drink",xlab="Year")  
points(ts(data.ao,frequency=period,start=start(data.ts)),col=2)
points(ts(data.ls,frequency=period,start=start(data.ts)),col=3)

p <- 2
q <- 1
ps <- 0
qs <- 1
d <- 2
ds <- 1

datareg.ts <- ts(cbind(data.ts,rep(1,n)),start=start(x),frequency=period)
fit.mle <- maxent.fit(datareg.ts,ao,ls,p,q,ps,qs,d,ds)
par.mle <- fit.mle[[1]]
psi.mle <- fit.mle[[2]]$par
ts.resid <- ts(c(rep(NA,r+d+ds*period),fit.mle[[3]]),
               frequency=period,start=start(data.ts))
plot(ts.resid)
acf(ts.resid[-seq(1,r+d+ds*period)],lag.max = 4*period,main="Residual")
spec.ar(ts.resid[-seq(1,r+d+ds*period)])

alpha <- 1
out <- maxent.ev(data.ts,ao,ls,psi.mle,p,q,ps,qs,d,ds,alpha)
1-pchisq(out[[5]],df=r)
kappa <- 1 - sqrt((qchisq(1-alpha,df=r))/out[[5]])
x.entropy <- out[[4]]
filts <- hpsa(n,period,1/100,1)
x.sa1 <- x.entropy - filts[[2]] %*% x.entropy
x.sa2 <- data.ts - filts[[2]] %*% x.entropy

setwd("C:\\Users\\neide\\OneDrive\\Documents\\Research\\MaxEntOutlier\\Figures")

#pdf(file="FoodDrinkSA.pdf",width=5,height=4)
plot(data.ts,ylim=c(10,11.5),ylab="",xlab="Year")
lines(ts(x.entropy,frequency=period,start=start(data.ts)),col=4)
lines(ts(x.sa1,frequency=period,start=start(data.ts)),col=3)
lines(ts(x.sa2,frequency=period,start=start(data.ts)),col=2)
dev.off()


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


