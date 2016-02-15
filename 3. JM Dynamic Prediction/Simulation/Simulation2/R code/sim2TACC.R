rm(list=ls())

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
	H = matrix(rnorm(2*n, 0, 1),n,2)
	W = matrix(rnorm(2*n, 0, 1),n,2)

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
	w = which(surv(100)-rnd < 0)
	# print(w)
	# w = which(1-surv(101)-rnd > 0)
	Ti[w] = sapply(w,function(j) uniroot(function(x) surv(x)[j]-rnd[j],lower=0,upper=100)$root)
	Ti[Ti==0] = min(Ti[Ti>0])/2
	Ci = rexp(n, 0.01)
	Delta = as.numeric(Ti < Ci)
	Ti2 = pmin(Ti, Ci)

	list(Ti=Ti2, event=Delta, H=H, W=W, U=u)
}

# test_Ti = sim_Ti(alpha1=1, alpha2=1, delta=c(-3,-3), gamma=c(-2,-3))

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
	list(y = y, X = X, J=count)
}

###############################################################
#----- function to simulate multiple joint data sets ---------#
###############################################################
sim_multiple_data = function(N, sur_fun=sim_Ti, longi_fun=sim_longitudinal_data, alpha, tau,beta, delta, gamma){
	# N - number of data sets to generate
	# sur_fun - function to simulate survival data
	# longi_fun - function to sumualte longitudinal data
	# alpha - association mechanism for JM
	# tau - quantile

	outdata = vector(mode='list', N)
	for (i in 1:N){
		sur_data = sur_fun(alpha1=alpha[1], alpha2=alpha[2], delta=delta, gamma=gamma)
		longi_data = longi_fun(sur_data, tau=tau, beta=beta, delta=delta)
		outdata[[i]] = list(survival_data=sur_data, longitudinal_data=longi_data)
	}
	outdata
}

## simulate multiple data sets
sim_data = sim_multiple_data(30, alpha=c(1,1),tau=0.25, beta=c(1, 1), delta=c(-1,1), gamma=c(-1,1))

# sim_data = sim_multiple_data(30, alpha=c(2,2),tau=0.25, beta=c(2,-2), delta=c(-4,2), gamma=c(-3,-1))
# define function to run the model
QRJM_jags = function(data,tau,I=500){
	library(R2jags)
	setwd("/work/02784/myang3/QRJM/model")
	model = "/work/02784/myang3/QRJM/model/QRJM_Cholesky.txt"
	# file.show(model)
	# load longitudinal data
	y = data[[2]][['y']]
	X = data[[2]][['X']]
	J = data[[2]][['J']]

	# load survival data
	Ti = data[[1]][['Ti']]
	H = data[[1]][['H']]
	W = data[[1]][['W']]
	event =  data[[1]][['event']]

	# alpha2 = 1
	zeros = rep(0, I)
	zero = c(0,0)

	jags.data = list(y=y, X=X, H=H, W=W, Ti=Ti, event=event, qt=tau, t=c(0, 0.25, 0.5, 1, 1.5, 2, 3), I=I, J=J, zeros=zeros, zero=zero)
 	jags.params = c("beta", "delta","gamma","alpha1","alpha2","sigma", "c", "w11","w21","w22")
 	jags.inits = function(){	list(beta=c(0.1,0.1), delta=c(0.1, 0.1), gamma=c(0.1, 0.1), sigma=0.1, alpha1=0.1, alpha2=0.1, c=0.1)
  	}
  jags(data=jags.data, inits=jags.inits, jags.params, n.iter=15000, n.burnin=5000, model.file=model)
}



### run parallely on TACC
# setting up parallel computing in cluster
library(foreach)
library(doParallel)
jEOM = F #job exec on multiplenodes
raw_PE <- switch(jEOM+1, F = "localhost.ls4.tacc.utexas.edu", T = system("cat $PE_HOSTFILE",intern=T))  # or use $HOSTNAME ,`hostname -s`
catt('raw_PE is',raw_PE)
nodes_name <- sub("^([\\w-]+)\\..*$","\\1",raw_PE,perl=T)
nodes_slots <- switch(jEOM+1, F = as.numeric(system("cat /proc/cpuinfo | grep processor|wc -l",intern=T)),T = as.numeric( sub("^([\\w-]+)\\..* (\\d+?) .*$","\\2",raw_PE,perl=T) ) )
#---------------------------------------#
######################################
## code confirmed correct.
socketHosts_name <- rep(nodes_name,times=nodes_slots)
cl <- makePSOCKcluster(socketHosts_name,homogeneous=TRUE)
clusterCall(cl, function() Sys.info()[c("nodename","machine")])
registerDoParallel(cl)
clusterExport(cl=cl,ls())



### 3 run the model and save the output
t1 = Sys.time()
sim2_25_0520_originafit2 <- foreach(i=1:30,.packages=c("rjags","R2jags","foreach"),.verbose=T )%dopar%{
	set.seed(123)
  	temp <- QRJM_jags(sim_data[[i]], tau=0.25)
  	temp$BUGSoutput$summary
}
Sys.time()-t1

setwd("/work/02784/myang3/QRJM_sim2/results")
save(sim2_25_0520_originafit2,file="sim2_25_0520_originafit2.Rdata")



### try to use LMJM to fit the data
t1 = Sys.time()
sim2_25_0520_LMJM <- foreach(i=1:30,.packages=c("rjags","R2jags","foreach"),.verbose=T )%dopar%{
	set.seed(123)
  	temp <- NormJM_jags(sim_data[[i]])
  	temp$BUGSoutput$summary
}
Sys.time()-t1

setwd("/work/02784/myang3/QRJM_sim2/results")
save(sim2_25_0520_LMJM,file="sim2_25_500_0520_LMJM.Rdata")


