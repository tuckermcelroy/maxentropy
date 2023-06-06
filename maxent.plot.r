maxent.plot <- function(x.ext,x.est,mse.est,prop=0,option=1,
                        top.bound = NULL, bot.bound = NULL)
{
  # maxent.plot by Tucker McElroy
  #   Plots maximum entropy extreme-value adjustment components.
  # Inputs:
  #   x.ext: data with original NAs, and extended with NAs fore and aft
  #   x.est: modified x.ext, such as:
  #     x.casted: sample corrected for missing values
  #     x.entropy: sample corrected for extremes and missing values
  #     x.sa: sample filtered to removed seasonality, and extremes added back
  #   mse.est: MSE of x.est, such as:
  #     mse.casted: MSE of sample corrected for missing values
  #     mse.entropy: MSE of sample corrected for extremes and missing values
  #     mse.sa: MSE of seasonal adjustment
  #   prop: scaling proportion to y axis to help with visualizing 
  #   option: various plotting options
  #     1: plot x.ext alone
  #     2: plot x.ext overlaid with x.est with uncert bands
  #     3: plot x.est alone with uncert bands
  #   top.bound, bot.bound: top and bottom bounds for the ylimit
  
  period <- frequency(x)
  eps <- 10e-10
 
  if(option==1)
  {
    plot(x.ext,ylab="",xlab="Year",lwd=2)
  }
  if(option==2)
  {
#    top.bound <- (1+sign(max(x.ext,na.rm=TRUE))*prop)*max(x.ext,na.rm=TRUE)
#    bot.bound <- (1-sign(min(x.ext,na.rm=TRUE))*prop)*min(x.ext,na.rm=TRUE)
    if(length(top.bound)==0)
    { top.bound <- (1+sign(max(x.est,na.rm=TRUE))*prop)*max(x.est,na.rm=TRUE) }
    if(length(bot.bound)==0)
    { bot.bound <- (1-sign(min(x.est,na.rm=TRUE))*prop)*min(x.est,na.rm=TRUE) }
    x.est_low <- x.est - 2*sqrt(diag(mse.est)+eps)
    x.est_hi <- x.est + 2*sqrt(diag(mse.est)+eps)
    plot(x.est,col=4,ylim=c(bot.bound,top.bound),ylab="",xlab="Year")
    lines(x.ext,lwd=2)
    polygon(c(time(x.est),rev(time(x.est))),
            c(pmax(bot.bound,x.est_low),rev(pmin(top.bound,x.est_hi))),
              col="#0000FF40",border=NA)
  }
  if(option==3)
  {
    if(length(top.bound)==0)
    { top.bound <- (1+sign(max(x.est,na.rm=TRUE))*prop)*max(x.est,na.rm=TRUE) }
    if(length(bot.bound)==0)
    { bot.bound <- (1-sign(min(x.est,na.rm=TRUE))*prop)*min(x.est,na.rm=TRUE) }
    x.est_low <- x.est - 2*sqrt(diag(mse.est)+eps)
    x.est_hi <- x.est + 2*sqrt(diag(mse.est)+eps)
    plot(x.est,col=2,ylim=c(bot.bound,top.bound),ylab="",xlab="Year")
    polygon(c(time(x.est),rev(time(x.est))),
            c(pmax(bot.bound,x.est_low),rev(pmin(top.bound,x.est_hi))),
            col="#FF000040",border=NA)
  }
   
}


