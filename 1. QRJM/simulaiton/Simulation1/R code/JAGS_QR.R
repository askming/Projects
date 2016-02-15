# library(LaplacesDemon)
library(MASS)

###############################################################
################# functions related with ALD  #################
###############################################################
# desity function of ALD(mu, sigma,p)
# dald = function(x, p, mu, sigma){
# 	# density functiol of ALD(mu, sigma, tau) based on the definition
# 	rho = (abs(x-mu) + (2*p -1)*(x-mu)) / (2*sigma)
# 	den = p*(1-p)/sigma*exp(-rho)
# }

# r_ald = function(n, location=0, sigma, p){
# 	# random variable simulation based on normal+exponetial mixture
# 		u = rnorm(n)
# 		z = rexp(n)
# 		v = sigma*z
# 		theta = (1-2*p)/(p*(1-p))
# 		tau = sqrt(2/(p*(1-p)))
# 		sample = theta * z + tau * sqrt(z) * u
# 	}

# # check the validity of above functions
# x = seq(-20,20,0.01)
# plot(x, dald(x, sigma=1, mu=0, p=0.25), ylim=c(0,0.45), main='', ylab="density",type="l", lwd=2)
# lines(x, dald(x, p=0.75, mu=0, sigma=1), col="blue", lty=2, lwd=2)
# lines(x, dald(x, p=0.5, mu=0, sigma=1), col="red",lty=3, lwd=2)
# lines(x, dnorm(x, 0, 1),col="green",lty=4, lwd=2)
# legend("topright", legend=c("N(0,1)", "LD(0, 1)", "ALD(0, 1, 0.25)", "ALD(0, 1, 0.75)"), lty=c(4, 3, 2, 1), col=c("green", "red", "blue", "black"), lwd=c(2,2,2,2), cex=0.75)


###############################################################
###### function to simulate longitudinal data #################
###############################################################
sim_longitudinal_data = function(n=250, time=c(0, 0.25, 0.5, 0.75, 1, 3), tau, sigma=1, beta=c(1,1), delta=c(1,1)){
	# random generation of ALD(0,sigma,p)
	rald = function(n, location=0, scale, p){
		u = rnorm(n)
		z = rexp(n)
		v = scale*z
		theta = (1-2*p)/(p*(1-p))
		tau = sqrt(2/(p*(1-p)))
		sample = theta * z + tau * sqrt(z) * u
	}
	# survival_data - data simulated from survival model
	# n - # of subjects
	# time - time points of observations
	# tau - quantile
	# sigma - scale parameter
	time = time # at most # = length(time) observations per patient
	y = matrix(NA, nrow=n, ncol=length(time)) # wide format
	# Ti = survival_data$Ti
	U = mvrnorm(n, mu=c(0,0), Sigma=matrix(c(0.09, 0.09*0.16, 0.09*0.16, 0.09), nrow=2, byrow=T))
	attributes(U)[[2]]=NULL # remove 'dimnames' attribute
	H = matrix(rnorm(2*n), ncol=2)
	X = cbind(1, rnorm(n))
	count = sample.int(6, n, replace=TRUE) # number of observations after drop-outs
	# J = length(time)

	for (i in 1:n){
		for (j in 1:count[i]){
			location = beta %*% X[i, ] + delta %*% c(H[i,1], H[i,2]*time[j]) + U[i, ] %*% c(1, time[j])
			y[i,j] = location + rald(1, scale=sigma, p=tau)		
		}	
	}
	
	list(y = y, X = X, H = H, J=count)		
}

ald_data = sim_longitudinal_data(tau=0.25)


########################################################################
---------------------------- 2. run the model --------------------------
########################################################################
library(R2jags)

QRJM_jags = function(data,tau,I=250){
	setwd("/Users/askming/Documents/github/Projects/1.\ QRJM/simulaiton/model\ files")
	model = "/Users/askming/Documents/github/Projects/1.\ QRJM/simulaiton/model\ files/QRJM_QR.txt"
	# file.show(model)
	# load longitudinal data
	y = data[['y']]
	X = data[['X']]
	H = data[['H']]
	J = data[['J']]
	
	jags.data = list(y=y, X=X, H=H, qt=tau, t=c(0, 1/4, 1/2, 3/4, 1, 3), I=I, J=J)
 	jags.params = c("beta", "delta","sigma")
 	jags.inits = function(){	list(beta=c(0.1,.1), delta=c(.1, .1), sigma=.1)	
  	}	
  jags(data=jags.data, inits=jags.inits, jags.params, n.iter=15000, n.burnin=5000, model.file=model)	
}

qrfit = QRJM_jags(ald_data, tau=0.25)
print(qrfit)

### seems to work




