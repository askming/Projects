##### R code to do simulation study of longitudianl quantile regression with informative drop-outs using JM approach ####

#### date: 2014/10/22 (happy birthday to Lan!)####
#### author: ming yang ####


### 1. simulate data ###
library(LaplacesDemon)
# ralaplace(n, location, scale)
library(MASS)

sim_longitudinal_data = function(n=250, t=6, tau, sigma=1, beta, delta, re){
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


calculate_Ti = function(longidata, n=250, alpha, delta, gamma){
	Time[i] = numeric(250)
	S = runif(n) # survival probability
	Z = longidata$Z
	H = longidata$H
	U = longidata$U
	W = matrix(rnorm(n*2), ncol=2)
	
	if (alpha == c(0, 0)){
		for (i in 1:n){
			fun = function(gamma, t, W){
				exp(-t*exp(gamma%*% W[i,])) - S[i]
			}
		}
	}
	
	else{
		for (i in 1:n){
			fun = function(n, alpha, delta, i, t, W){	
				# survival function formula	
				exp(- (exp(alpha[1]*(delta %*% c(H[i,1,1], H[i,1,2]*t)) + alpha[2]*(U[i,]%*%c(Z[i,1,1], t)) + gamma%*% W[i,]) -exp(alpha[1]*delta[1]*H[i,1,1] + alpha[2]*U[i,1])+gamma%*% W[i,])/(alpha[2]*U[i,2]+alpha[1]*delta[2]*H[i,1,2])) - S[i]
		
			}
			Time[i] = optimize(f=fun, interval=c(-1e10, 1e10))		
		}
	}		
	return(Time)	
}





### 2. run the model ###
