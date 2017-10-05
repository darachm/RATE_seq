 LLexpKINETICS<-function(x, y, n_delay, nu, aOFmHillParams, nH, t, grad=TRUE, H="A"){
 # LLexpKINETICS<-function(x, y, n_delay, nu, aOFmHillParams, nH, t, H="A")
 # t - time points of measurments
 # x is the vector of logs of the parameters in the first-order synthesis model with delay:
 # 																											z_ss, alpha, delta
 # z_ss - maximum amol, i.e., maximum of z(t).  m is synononmous with z_ss (see below).
 # alpha - decay rate constant
 #	delta <- n/t_d	#	reciprocal of "delay" per step
 # n_delay - number of delay stage, chosen and fixed, not fitted
 # nu is  the library-specific scale factor that convrs from amol to counts 
 # y is one row of the original count/read matrix, called mRNAmat, for example
 # aOFmHillCoefs is the vector of parameters in the
 # model for a(m), where m and z_ss are synonomous:
 # 	a(z) = a0 +(aInf-a0)* ( z^nH/(z_half^nH + z^nH) ), where
 # nH is the Hill coefficient
 #
 ##########MATHEMATICAL NOTES
 # dz_l/dt = gamma*h(t) - alpha*z_l, where z_l is labeled trancript, is a delayed
 #                                                      step-like function
 # (1/alpha) dz_l/dt = (gamma/alpha)*h(t)  - z_l
 # (1/alpha) dz_l/dt = z_ss*h(t)  - z_l, where z_ss is the steady-state value of z
 #########END MATHEMATICAL NOTES 
 # nu is  as defined in notes
 	y<-as.numeric(y)
	z_ss <- exp(x[1]); alpha <- exp(x[2]); delta <- exp(x[3])
	m <- z_ss
	a0    <- aOFmHillParams["a0"]
	aInf  <- aOFmHillParams["aInf"]
	zHalf <- aOFmHillParams["zHalf"]
	a <- a0 +( aInf - a0 )* ( z_ss^nH/(zHalf^nH + z_ss^nH) )

		w <- ExpSynWithDelay(t, alpha, t_d=n_delay/delta, n_delay)
		eta <- delta/(delta-alpha)
		n<-n_delay
		B <- matrix(0, nrow=length(t), ncol=n_delay)
		C <- B
		for(k in 0:(n-1)){
			B[,k+1] <-                (delta*t)^k/gamma(k+1)
			C[,k+1] <- ( (delta-alpha)*t  )^k/gamma(k+1)
		}

	LL <- lgamma(a+y) - lgamma(a) - lgamma(y+1) +
		      y * ( log(m*nu*w) - log(a+m*nu*w) ) +
		      a * ( log(a)      - log(a+m*nu*w) )
		
	dLL.dw <- ( a/(a+m*nu*w) ) * ( y - m*nu*w )/w
	
	J<-0:(n-1)

	dw.dalpha <-  eta^n*( t*exp(-alpha*t) )  -  ( n*eta^(n-1)*delta/(delta-alpha)^2  )*exp(-alpha*t) -
							exp(-delta*t)*eta^n*(1/(delta-alpha))*rowSums(  C%*%diag( J-n )  )
	
	dLL.dalpha <- dLL.dw*dw.dalpha
	
	dw.ddelta <-                n*eta^(n-1)*(alpha/(delta-alpha)^2)*exp(-alpha*t) +
	             exp(-delta*t)*rowSums(  B%*%diag( (n-J)*eta^(n-J-1) *
	             (-alpha/(delta-alpha)^2) - (J/delta)*(1-eta^(n-J)) )  )  +
	             t*exp(-delta*t)*rowSums(  B%*%diag( 1-eta^(n-J) )  )
	dLL.ddelta <- dLL.dw*dw.ddelta                     
		
	da.dm <- (aInf-a0)*nH* m^(nH-1) * zHalf^nH / (zHalf^nH + m^nH)^2

	dLL.dm <- ( (y-m*nu*w)/(a+m*nu*w) )*( a/m - da.dm ) +
	          ( digamma(a+y) - digamma(a) + log(a) - log(a+nu*m*w) )*da.dm
	
	   
	dLL.dx1         <- dLL.dm * m
	dLL.dx1.Sum <- sum(dLL.dx1)
	dLL.dx2          <- dLL.dalpha*alpha
	dLL.dx2.Sum <- sum(dLL.dx2)
	dLL.dx3         <- dLL.ddelta*delta
	dLL.dx3.Sum  <- sum(dLL.dx3)
	gradLL.x <- c(dLL.dx1.Sum, dLL.dx2.Sum, dLL.dx3.Sum)
	
	gradLL      <- gradLL.x
	LL.sum       <- sum(LL)
	
	res<- -LL.sum
	if(grad==TRUE){
			attr(res,"gradient") <- -gradLL
	}
	res
		           
}