hend <- function(q)
{
	# Function to compute coefficients of Henderson filter
	#  just one-sided output, starting with zeroth coefficient
	# Note that the order of the Henderson is 2q+1 
	#   (the number of coefficients)
	coeff <- ((q+1)^2)*((q+2)^2)*((q+3)^2)
	for(j in 1:q)
	{
		coeff <- c(coeff,(((q+1)^2-j^2)*((q+2)^2-j^2)*((q+3)^2-j^2)))
	}
	beta <- 0
	gamma <- 0
	for(j in 1:q)
	{
		beta <- beta + j^2*coeff[j+1]
		gamma <- gamma + j^4*coeff[j+1]
	}
	alpha <- 2*sum(coeff) - coeff[1]
	beta <- 2*beta
	gamma <- 2*gamma
	a <- -gamma/(beta^2 - alpha*gamma)
	b <- beta/(beta^2 - alpha*gamma)
	for(j in 1:(q+1))
	{
		coeff[j] <- coeff[j]*(a + b*(j-1)^2)
	}
	coeff <- c(rev(coeff),coeff[-1])
	return(coeff)
}