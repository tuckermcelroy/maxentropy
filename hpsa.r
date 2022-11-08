hpsa <- function(n,period,q,r)
{
	# hpsa
	#     gives an HP filter for seasonal data
	#	presumes trend+seas+irreg structure
	#		trend is integrated rw
	#		seas is seasonal rw
	#		irreg is wn
	#	q is snr for trend to irreg
	#	r is snr for seas to irreg

# define trend differencing matrix

delta.mat <- diag(n)
temp.mat <- 0*diag(n)
temp.mat[-1,-n] <- -2*diag(n-1)
delta.mat <- delta.mat + temp.mat
temp.mat <- 0*diag(n)
temp.mat[c(-1,-2),c(-n,-n+1)] <- 1*diag(n-2)
delta.mat <- delta.mat + temp.mat
diff.mat <- delta.mat[3:n,]

# define seasonal differencing matrix

delta.mat <- diag(n)
temp.mat <- 0*diag(n)
inds <- 0
for(t in 1:(period-1))
{
	temp.mat <- 0*diag(n)
	temp.mat[-(1+inds),-(n-inds)] <- 1*diag(n-t)
	delta.mat <- delta.mat + temp.mat
	inds <- c(inds,t)
}
sum.mat <- delta.mat[period:n,]

# define two-comp sig ex matrices

#trend.mat <- solve(diag(n) + t(diff.mat) %*% diff.mat/q)
#seas.mat <- solve(diag(n) + t(sum.mat) %*% sum.mat/r)
trend.mat <- diag(n) - t(diff.mat) %*% solve(q*diag(n-2) + diff.mat %*% 
	t(diff.mat)) %*% diff.mat
seas.mat <- diag(n) - t(sum.mat) %*% solve(r*diag(n-period+1) + sum.mat %*% 
	t(sum.mat)) %*% sum.mat

# define three-comp sig ex matrices

trend.filter <- solve(diag(n) - trend.mat %*% seas.mat) %*%
	trend.mat %*% (diag(n) - seas.mat)
seas.filter <- solve(diag(n) - seas.mat %*% trend.mat) %*%
	seas.mat %*% (diag(n) - trend.mat)
irreg.filter <- diag(n) - (trend.filter + seas.filter)

filters <- list(trend.filter,seas.filter,irreg.filter)
return(filters)
}



