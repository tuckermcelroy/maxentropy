maxent.plot <- function(x.ext,x.casted,x.entropy,mse.casted,mse.entropy,option=1)
{
  # maxent.plot by Tucker McElroy
  #   Plots maximum entropy extreme-value adjustment components.
  # Inputs:
  #   x.ext: data with original NAs, and extended with NAs fore and aft
  #   x.casted: sample corrected for missing values
  #   x.entropy: sample corrected for extremes and missing values
  #   mse.casted: MSE of sample corrected for missing values
  #   mse.entropy: MSE of sample corrected for extremes and missing values
  #   option: various plotting options
  #     1: plot x.ext alone
  #     2: plot x.ext overlaid with x.casted with uncert bands
  #     3: plot x.ext overlaid with x.entropy with uncert bands
  
  period <- frequency(x)
  eps <- 10e-10
 
  if(option==1)
  {
    plot(x.ext,ylab="",xlab="Year",lwd=2)
  }
  if(option==2)
  {
    x.casted_low <- x.casted - 2*sqrt(diag(mse.casted)+eps)
    x.casted_hi <- x.casted + 2*sqrt(diag(mse.casted)+eps)
    plot(x.casted,col=4,ylim=c(min(min(x.casted_low),min(x.ext,na.rm=TRUE)),
                               max(max(x.casted_hi),max(x.ext,na.rm=TRUE))))
    lines(x.ext)
    polygon(c(time(x.casted),rev(time(x.casted))),
            c(x.casted_low,rev(x.casted_hi)),
              col="#0000FF40",border=NA)
  }
  if(option==3)
  {
    x.entropy_low <- x.entropy - 2*sqrt(diag(mse.entropy)+eps)
    x.entropy_hi <- x.entropy + 2*sqrt(diag(mse.entropy)+eps)
    plot(x.entropy,col=2,ylim=c(min(min(x.entropy_low),min(x.ext,na.rm=TRUE)),
                                max(max(x.entropy_hi),max(x.ext,na.rm=TRUE))))
    lines(x.ext)
    polygon(c(time(x.entropy),rev(time(x.entropy))),
            c(x.entropy_low,rev(x.entropy_hi)),
            col="#FF000040",border=NA)
    
  }
   
}


