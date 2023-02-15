polyLagMatrix <- function(poly,dim)
{
  # forms a matrix that is polynomial in the lag matrix
  l <- length(poly)
  p <- l-1  
  L <- matrix(0,nrow = dim, ncol=dim)
  if(dim==1) polymat <- matrix(poly[1],nrow=1,ncol=1)
  if(dim > 1)
  {
    for(i in 2:dim)	{ L[i,(i-1)] <- 1 }
    Lpow <- diag(rep(1,dim))
    polymat <- 0
    for(j in 1:l)
    {
      polymat <- polymat + poly[j] * Lpow
      Lpow <- Lpow %*% L
    }
  }
  return(polymat)
}	