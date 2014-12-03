##### R code to do simulation study of longitudianl quantile regression with informative drop-outs using JM approach ####

#### date: 2014/10/22 (happy birthday to Lan!)####
#### author: ming yang ####
rm(list=ls())
########################################################################
---------------------------- 1. simulate data --------------------------
########################################################################
library(LaplacesDemon)
library(MASS)

###############################################################
############ function to simulate survival time ###############
###############################################################
# survival function is given by: S(t)= exp(- exp(B) * (exp(A*t) - 1) ) / A)
sim_Ti = function(n=500, alpha, delta=c(1,1), gamma=c(1,1)){
	Time = numeric(n)
	S = runif(n) # survival probability
	H = matrix(rnorm(2*n), ncol=2)
	U = mvrnorm(n, mu=c(0,0), Sigma=matrix(c(0.09, 0.09*0.16, 0.09*0.16, 0.09), nrow=2, byrow=T))
	attributes(U)[[2]]=NULL # remove 'dimnames' attribute
	W = matrix(rnorm(n*2), ncol=2)
	
	if (alpha[1]==0 & alpha[2]==0){
		for (i in 1:n){
			Time[i] = - log(S[i]) / exp(gamma %*% W[i,])
		}
	}
	
	else{
		for (i in 1:n){
			B = alpha[1] * delta[1] * H[i,1] + alpha[2] * U[i,1] + gamma %*% W[i,]
			A = alpha[2] * U[i,2] + alpha[1] * delta[2] * H[i,2]
			Time[i] = log(1-log(S[i])*A/exp(B))/A
		}
	}	
	Ti_id = which(!is.na(Time))
	Time=Time[Ti_id][1:250]
	U=U[Ti_id, ][1:250, ]
	H=H[Ti_id, ][1:250, ]
	W=W[Ti_id, ][1:250, ]
		
	list(Ti=Time, H=H, U=U, W=W)	
}

####### try it #######
surdata = sim_Ti(alpha=c(1,1))
# when A is negative there may be no solution for t, so remove those observations


###############################################################
###### function to simulate longitudinal data #################
###############################################################
sim_longitudinal_data = function(survival_data=surdata, n=250, time=c(0, 0.25, 0.5, 0.75, 1, 3), tau, sigma=1, beta=c(1,1), delta=c(1,1)){
	# survival_data - data simulated from survival model
	# n - # of subjects
	# time - time points of observations
	# tau - quantile
	# sigma - scale parameter
	time = time # at most # = length(time) observations per patient
	y = matrix(NA, nrow=n, ncol=length(time)) # wide format
	Ti = survival_data$Ti
	U = survival_data$U # random effects
	H = survival_data$H
	X = cbind(1, rnorm(n))
	count = sapply(Ti, function(x) sum(x > time)) # number of observations after drop-outs
	event = as.numeric(count != 6) # 1 for event, 0 for censoring

	for (i in 1:n){
		for (j in 1:count[i]){
			location = beta %*% X[i, ] + delta %*% c(H[i,1], H[i,2]*time[j]) + U[i,] %*% c(1, time[j])
			y[i,j] = ralaplace(1, location, scale=sigma, kappa=tau)		
		}	
	}
	
	list(y = y, X = X, J=count, event=event)		
}

longidata = sim_longitudinal_data(surdata, time=c(0, 0.25, 0.5, 0.75, 1, 3), tau=0.25)


###############################################################
###### function to simulate multiple joint data sets ##########
###############################################################
sim_multiple_data = function(N, sur_fun=sim_Ti, longi_fun=sim_longitudinal_data, alpha, tau){
	# N - number of data sets to generate
	# sur_fun - function to simulate survival data
	# longi_fun - function to sumualte longitudinal data
	# association mechanism for JM
	outdata = vector(mode='list', N)
	for (i in 1:N){
		sur_data = sur_fun(alpha=alpha)
		longi_data = longi_fun(sur_data, tau=tau)
		outdata[[i]] = list(survival_data=sur_data, longitudinal_data=longi_data)		
	}
	outdata	
}

## try it out
t1=Sys.time()
testdata = sim_multiple_data(N=500, alpha=c(1,1), tau=0.25)
t2=Sys.time()
t2-t1

check = function(data, ind1, ind2, var, fun, ...){
	# function used to check the simulated data
	var = data[[ind1]][[ind2]][[var]]
	fun(var, ...)
}

# check(testdata, 8, 2, 'event', fun=table)


########################################################################
---------------------------- 2. run the model --------------------------
########################################################################
library(R2jags)

QRJM_jags<-function(data,tau,I=250){
	setwd("/Users/askming/Documents/github/Projects/1.\ QRJM/simulaiton/model\ files")
	model = "/Users/askming/Documents/github/Projects/1.\ QRJM/simulaiton/model\ files/QRJM.txt"
	# file.show(model)
	y = data[[2]][['y']]
	X = data[[2]][['X']]
	J = data[[2]][['J']]
	event =  data[[2]][['event']]

	Ti = data[[1]][['Ti']]
	H = data[[1]][['H']]
	W = data[[1]][['W']]

	jags.data = list(y=y, X=X, H=H, W=W, Ti=Ti, event=event, qt=tau, t=c(0, 1/4, 1/2, 3/4, 1, 3), I=I, J=J)
 	jags.params = c("beta", "delta","gamma","alpha1","alpha2","sigma")
 	jags.inits = function(){	list(beta=c(0.1,0.1), delta=c(0.1, 0.1), gamma=c(0.1, 0.1), sigma=0.1, alpha1=0.1, alpha2=0.1)	
  	}	

  jags(data=jags.data, inits=jags.inits, jags.params, n.iter=20000, n.burnin=10000, model.file=model)	
}



testfit = QRJM_jags(testdata[[1]], tau=0.25)

print(testfit)
traceplot(testfit)
xyplot(as.mcmc(testift))




