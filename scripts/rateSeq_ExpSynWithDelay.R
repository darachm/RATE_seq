ExpSynWithDelay <- function(t, alpha, t_d, n){
	#	alpha - decay rate constant (sum from degradation and dilution)
	#	t_d    - lag time (delay, akin to, but not exactly the same as a dead time)
	# 	n       - number of delay steps before RNA synthesis
	delta <- n/t_d	#	reciprocal of "delay" per step
	eta <- delta/(delta-alpha)
	H <- matrix(0, nrow=length(t), ncol=n)
	for(k in 0:(n-1)){
		H[,k+1] <- (delta*t)^k/gamma(k+1)
	}
	J=0:(n-1)
	w <- 1 - eta^n*exp(-alpha*t) - rowSums(     H%*%diag(  1-eta^(n-J)  )     ) * exp(-delta*t)
	w
}