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
