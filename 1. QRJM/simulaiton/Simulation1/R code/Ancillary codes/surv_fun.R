sim_Ti = function(x, alpha1, alpha2, delta, gamma){
	library(mvtnorm)
	# random effects
	Sigma=matrix(c(0.09, 0.09*0.16, 0.09*0.16, 0.09), nrow=2, byrow=T)
	u=rmvnorm(n=1,c(0,0),Sigma)
	# fixed effects
	H = rnorm(2, 1, 1)
	W = rnorm(2, 2, 4)

	# define the survival function
	surv = function(t) {
		if(alpha1==0 & alpha2==0) {res=t*exp(W%*%gamma)}
		else{
			res = ((exp(alpha2*u[1]+alpha2*t*u[2]+alpha1*delta[1]*H[1]+alpha1*delta[2]*H[2]*t+W%*%gamma)-exp(alpha1*H[1]*delta[1]+W%*%gamma+alpha2*u[1]))/(alpha2*u[2]+alpha1*H[2]*delta[2]))
		}	
		return(exp(-res))
	}
	# simulate event time
	surv(t=x)
}

sim_Ti(x=seq(0, 10), alpha1=1, alpha2=1, delta=c(2,-1), gamma=c(3, -1))