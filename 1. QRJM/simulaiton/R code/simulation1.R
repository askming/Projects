##### R code to do simulation study of longitudianl quantile regression with informative drop-outs using JM approach ####

#### date: 2014/10/22 (happy birthday to Lan!)####
#### author: ming yang ####

########################################################################
---------------------------- 1. simulate data --------------------------
########################################################################
library(LaplacesDemon)
library(MASS)

# function to simulate longitudinal data
sim_longitudinal_data = function(n=250, t=6, tau, sigma=1, beta=c(1,1), delta=c(1,1)){
	# n - # of subjects
	# t - # of measures per subjects without drop-outs
	# tau - quantile
	# sigma - scale parameter
	time = c(0, 0.25, 0.5, 0.75, 1, 3)
	y = matrix(NA, nrow=n, ncol=t) # wide format
	X = H = Z = array(NA, dim=c(n, t, 2)) # array, wide format
	 # leng = n
	X[ , , 1] = Z [, , 1]= 1  
	X[ , , 2] = matrix(rep(rnorm(n), each=t), ncol=t, byrow=T)
	Z[ , , 2] = matrix(rep(time, n), nrow=n, byrow=T)
	H[ , , 1] = matrix(rep(rnorm(n), each=t), ncol=t, byrow=T)
	H[ , , 2] = matrix(rep(rnorm(n), each=t), nrow=n, byrow=T)
	U = NULL # random effects
	for (i in 1:n){
		u = mvrnorm(1, mu = c(0,0), Sigma = matrix(c(0.09, 0.09*0.16, 0.09*0.16, 0.09), nrow=2, byrow=T))
		for (j in 1:t){
			U = rbind(U, u)
			location = beta %*% X[i,j, ] + delta %*% c(H[i,j,1], H[i,j,2]*time[j]) + u %*% Z[i,j, ]
			y[i,j] = ralaplace(1, location, scale = sigma, kappa = tau)
		
		}	
	}
	
	list(y = y, X = X, Z = Z, H = H, U = U)		
}
# try it
# testdata = sim_longitudinal_data(tau=0.25, beta=c(1,1), delta=c(1,1))

# function to simulate survival time
sim_Ti = function(longidata, n=250, alpha, delta=c(1,1), gamma=c(1,1)){
	Time = numeric(250)
	S = runif(n) # survival probability
	Z = longidata$Z
	H = longidata$H
	U = longidata$U[seq(1,n*6,6),]
	W = matrix(rnorm(n*2), ncol=2)
	
	if (alpha[1]==0 & alpha[2]==0){
		for (i in 1:n){
			Time[i] = - log(S[i]) / exp(gamma %*% W[i,])
		}
	}
	
	else{
		for (i in 1:n){
			B = exp(alpha[1] * delta[1] * H[i,1,1] + alpha[2] * U[i,1] + gamma %*% W[i,])
			A = alpha[2] * U[i,2] + alpha[1] * delta[2] * H[i,1,2]
			Time[i] = log(1-log(S[i])*A/exp(B))/A
		}
	}		
	return(Time)	
}

# try it
# Ti = sim_Ti(testdata, alpha=c(0,0), delta=c(1,1), gamma=c(1,1))
# Ti = sim_Ti(testdata, alpha=c(1,1), delta=c(1,1), gamma=c(1,1))
# when A is negative there may be no solution for t


C = rbeta(250, 4, 1) *5

sum(Ti > 3)/250



### 2. run the model ###
