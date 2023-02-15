x11filter <- function(p1,p2,q,per,comp)
{
	#  Function to generate filter coefficients for x-11 filters
	#   output is "two-sided" since all the filters are symmetric;
	# p1 is the p in first seasonal filter
	# p2 is the p in second seasonal filter
	# q is the order of Henderson trend
	# per is the seasonal period, either 4 or 12
	# comp is a flag for desired component filter: 1 for seasonal, 2 for nonseasonal,
	#  3 for trend, and 4 for irregular
	# requires function polymult.r and hend.r

  nu3 <- c(1,rep(0,per-1),1,rep(0,per-1),1)/3
  nu5 <- c(1,rep(0,per-1),1,rep(0,per-1),1,rep(0,per-1),1,rep(0,per-1),1)/5
  nu7 <- c(1,rep(0,per-1),1,rep(0,per-1),1,rep(0,per-1),1,rep(0,per-1),1,rep(0,per-1),1,rep(0,per-1),1)/7
  nu9 <- c(1,rep(0,per-1),1,rep(0,per-1),1,rep(0,per-1),1,rep(0,per-1),1,rep(0,per-1),1,rep(0,per-1),1,rep(0,per-1),1,rep(0,per-1),1)/9
  null <- c(1,rep(0,per-1),1,rep(0,per-1),1,rep(0,per-1),1,rep(0,per-1),1,rep(0,per-1),1,rep(0,per-1),1,rep(0,per-1),1,rep(0,per-1),1,rep(0,per-1),1,rep(0,per-1),1)/11
  nu13 <- c(1,rep(0,per-1),1,rep(0,per-1),1,rep(0,per-1),1,rep(0,per-1),1,rep(0,per-1),1,rep(0,per-1),1,rep(0,per-1),1,rep(0,per-1),1,rep(0,per-1),1,rep(0,per-1),1,rep(0,per-1),1,rep(0,per-1),1)/13
  nu15 <- c(1,rep(0,per-1),1,rep(0,per-1),1,rep(0,per-1),1,rep(0,per-1),1,rep(0,per-1),1,rep(0,per-1),1,rep(0,per-1),1,rep(0,per-1),1,rep(0,per-1),1,rep(0,per-1),1,rep(0,per-1),1,rep(0,per-1),1,rep(0,per-1),1,rep(0,per-1),1)/15
  if (p1 == 3) { lamp1 <- polymult(nu3,nu3) }
	if (p1 == 5) { lamp1 <- polymult(nu3,nu5) }
  if (p1 == 7) { lamp1 <- polymult(nu3,nu7) }
	if (p1 == 9) { lamp1 <- polymult(nu3,nu9) }
	if (p1 == 11) {	lamp1 <- polymult(nu3,nu11) }
	if (p1 == 13) { lamp1 <- polymult(nu3,nu13) }
	if (p1 == 15) {	lamp1 <- polymult(nu3,nu15) }
	if (p2 == 3) { lamp2 <- polymult(nu3,nu3) }
	if (p2 == 5) { lamp2 <- polymult(nu3,nu5) }
  if (p2 == 7) { lamp2 <- polymult(nu3,nu7) }
	if (p2 == 9) { lamp2 <- polymult(nu3,nu9) }
	if (p2 == 11) {	lamp2 <- polymult(nu3,nu11) }
	if (p2 == 13) {	lamp2 <- polymult(nu3,nu13) }
	if (p2 == 15) {	lamp2 <- polymult(nu3,nu15) }
	hq <- hend((q-1)/2)
	mu <- c(1,rep(2,per-1),1)/(2*per)
	muc <- c(rep(0,per/2),1,rep(0,per/2)) - mu

	omegaS <- polymult(lamp1,polymult(muc,muc))	
	m <- length(omegaS)
	omegaS <- polymult(hq,c(rep(0,(m-1)/2),1,rep(0,(m-1)/2)) - omegaS)
	m <- length(omegaS)
	omegaS <- polymult(muc,polymult(lamp2,c(rep(0,(m-1)/2),1,rep(0,(m-1)/2)) - omegaS))
  m <- length(omegaS)
	omegaN <- c(rep(0,(m-1)/2),1,rep(0,(m-1)/2)) - omegaS
  omegaT <- polymult(hq,omegaN)
	m <- length(hq)	
	omegaI <- polymult(omegaN,c(rep(0,(m-1)/2),1,rep(0,(m-1)/2))-hq)

	if(comp == 1) { out <- omegaS } 
	if(comp == 2) { out <- omegaN } 
	if(comp == 3) { out <- omegaT } 
	if(comp == 4) { out <- omegaI } 		
	    
	return(out)
}