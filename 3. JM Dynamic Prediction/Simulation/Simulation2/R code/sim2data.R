## Simulation study 2: to check the accuracy of predictions
## steps:
# - simulate data and fit the joint model
# - took a random sample of patients
# - these patients have measurements up to t
# - using their available measurements calculate the random effects
# - then calculate predictions for u > t, with u \in (4, 6, 8, 10, 14, 16,
# 20, 22, 26)



########################################################################
# --------------------------- 1. simulate data ------------------------#
########################################################################

###############################################################
#----------- function to simulate survival time --------------#
###############################################################
# survival function is given by: S(t)= exp(- exp(B) * (exp(A*t) - 1) ) / A)
sim_Ti = function(n=500, alpha1, alpha2, delta, gamma){
	library(mvtnorm)
	# random effects
	Sigma=matrix(c(0.09, 0.09*0.16, 0.09*0.16, 0.09), nrow=2, byrow=T)
	u=rmvnorm(n,c(0,0),Sigma)
	# fixed effects
	H = matrix(rnorm(2*n, 1, 1),n,2)
	W = matrix(rnorm(2*n, 2, 4),n,2)

	# define the survival function
	surv = function(t) {
		if(alpha1==0 & alpha2==0) {res=t*exp(W%*%gamma)}
		else{
			res = ((exp(alpha2*u[,1]+alpha2*t*u[,2]+alpha1*delta[1]*H[,1]+alpha1*delta[2]*H[,2]*t+W%*%gamma)-exp(alpha1*H[,1]*delta[1]+W%*%gamma+alpha2*u[,1]))/(alpha2*u[,2]+alpha1*H[,2]*delta[2]))
		}
		return(exp(-res))
	}
	# simulate event time
	rnd =runif(n)
	Ti = rep(Inf,n)
	w = which(surv(1e100)-rnd < 0)
	# w = which(1-surv(101)-rnd > 0)
	Ti[w] = sapply(w,function(j) uniroot(function(x) surv(x)[j]-rnd[j],lower=0,upper=1e100)$root)
	Ti[Ti==0] = min(Ti[Ti>0])/2
	Ci = rexp(n, 0.01)
	Delta = as.numeric(Ti < Ci)
	Ti2 = pmin(Ti, Ci)

	list(Ti=matrix(Ti2, ncol=1), event=matrix(Delta, ncol=1), H=H, W=W, U=u)
}

# test_Ti = sim_Ti(alpha1=1, alpha2=1, delta=c(-2,-1), gamma=c(-3,-1))

# test_Ti[[1]]
# mean(test_Ti[[2]])
# sum(test_Ti[[1]]>=2)
# sum(test_Ti[[1]]>=4)
# sum(test_Ti[[1]]>=8)
# sum(test_Ti[[1]]>=16)
# summary(test_Ti[[1]])

###############################################################
#----- function to simulate longitudinal data --------------- #
###############################################################

sim_longitudinal_data = function(survival_data=surdata, n=500, time=c(0, 0.25, 0.5, 1, 1.5, 2, 3), tau, sigma=1, beta, delta){
	# survival_data - data simulated from survival model
	# n - # of subjects
	# time - time points of observations
	# tau - quantile
	# sigma - scale parameter

	# random generation of ALD(0,sigma,p)
	rald = function(n, location=0, scale, p){
		u = rnorm(n)
		z = rexp(n)
		v = scale*z
		theta = (1-2*p)/(p*(1-p))
		tau = sqrt(2/(p*(1-p)))
		samples = theta * z + tau * sqrt(z) * u
		samples
	}

	y = matrix(NA, nrow=n, ncol=length(time)) # wide format
	Ti = survival_data$Ti
	U = survival_data$U # random effects
	H = survival_data$H
	X = cbind(1, rnorm(n))
	count = sapply(Ti, function(x) sum(x > time)) # number of observations after drop-outs

	for (i in 1:n){
		for (j in 1:count[i]){
			location = beta %*% X[i, ] + delta %*% c(H[i,1], H[i,2]*time[j]) + U[i,] %*% c(1, time[j])
			y[i,j] = location + rald(1, scale=sigma, p=tau)
		}
	}
	list(y = y, X = X, J=matrix(count, ncol=1))
}

###############################################################
#----- function to simulate multiple joint data sets ---------#
###############################################################
sim_multiple_data = function(N, n, sur_fun=sim_Ti, longi_fun=sim_longitudinal_data, alpha, tau, beta, delta, gamma){
	# N - number of data sets to generate
	# n - sample size in each data set
	# sur_fun - function to simulate survival data
	# longi_fun - function to simulate longitudinal data
	# alpha - association mechanism for JM
	# tau - quantile

	outdata = vector(mode='list', N)
	for (i in 1:N){
		sur_data = sur_fun(n=n, alpha1=alpha[1], alpha2=alpha[2], delta=delta, gamma=gamma)
		longi_data = longi_fun(n=n, sur_data, tau=tau, beta=beta, delta=delta)
		outdata[[i]] = list(survival_data=sur_data, longitudinal_data=longi_data)
	}
	outdata
}

# sim_data = sim_multiple_data(3, alpha=c(1,1),tau=0.5, beta=c(1,-1), delta=c(-2,1), gamma=c(-3,-1))




# sim_data[[1]][[1]]$Ti[1:30]
# sim_data[[1]][[2]]$y[1:30,]
# summary(sim_data[[1]][[1]]$Ti)
# range(sim_data[[1]][[1]]$Ti)
# mean(sim_data[[1]][[1]]$event)

