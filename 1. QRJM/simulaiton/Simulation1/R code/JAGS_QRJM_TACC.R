rm(list=ls())

###############################################################
#----------- function to simulate survival time --------------#
###############################################################
# survival function is given by: S(t)= exp(- exp(B) * (exp(A*t) - 1) ) / A)
sim_Ti = function(n=800, alpha1, alpha2, delta=c(-1,1), gamma=c(-1,1)){
	library(mvtnorm)
	# random effects
	Sigma=matrix(c(0.09, 0.09*0.16, 0.09*0.16, 0.09), nrow=2, byrow=T)
	u = rmvnorm(n,c(0,0),Sigma)
	H = matrix(rnorm(2*n),n,2)
	W = matrix(rnorm(2*n),n,2)

	# define the survival function
	surv = function(t) {
		if(alpha1==0 & alpha2==0) {res=t*exp(W%*%gamma)}
		else{
			res = ((exp(alpha2*u[,1]+alpha2*t*u[,2]+alpha1*delta[1]*H[,1]+alpha1*delta[2]*H[,2]*t+W%*%gamma)-exp(alpha1*H[,1]*delta[1]+W%*%gamma+alpha2*u[,1]))/(alpha2*u[,2]+alpha1*H[,2]*delta[2]))
		}

		return(exp(-res))
	}

	rnd =runif(n)
	Ti = rep(Inf,n)
	# w = which(surv(1e8)-rnd < 0)
	w = which(1-surv(101)-rnd > 0)
	Ti[w] = sapply(w,function(j) uniroot(function(x) 1-surv(x)[j]-rnd[j],lower=0,upper=101)$root)
	Ti[Ti==0] = min(Ti[Ti>0])/2
	Ci = 5*rbeta(n,4,1)
	Delta = as.numeric(Ti < Ci)
	Ti2 = pmin(Ti, Ci)

	list(Ti=Ti2, event=Delta, H=H, W=W, U=u)
}

######################################################################
# sim_multiple_surdata = function(N, sur_fun=sim_Ti, alpha){
# 	# N - number of data sets to generate
# 	# sur_fun - function to simulate survival data
# 	# longi_fun - function to sumualte longitudinal data
# 	# alpha - association mechanism for JM
# 	# tau - quantile

# 	outdata = vector(mode='list', N)
# 	for (i in 1:N){
# 		sur_data = sur_fun(alpha1=alpha[1], alpha2=alpha[2])
# 		#longi_data = longi_fun(sur_data, tau=tau)
# 		outdata[[i]] = sur_data
# 	}
# 	outdata
# }

# sur_data = sim_multiple_surdata(30, alpha=c(0,0))




###############################################################
#----- function to simulate longitudinal data independetly ----------------#
###############################################################
sim_longitudinal_data = function(survival_data=surdata, n=800, time=c(0, 0.25, 0.5, 0.75, 1, 3), tau, sigma=1, beta=c(1,1), delta=c(-1,1)){
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

	# random variable generation from ALD(0,sigma,p)
	rald = function(n, location=0, scale, p){
		u = rnorm(n)
		z = rexp(n)
		v = scale*z
		theta = (1-2*p)/(p*(1-p))
		tau = sqrt(2/(p*(1-p)))
		sample = theta * z + tau * sqrt(z) * u
	}

	for (i in 1:n){
		for (j in 1:count[i]){
			location = beta %*% X[i, ] + delta %*% c(H[i,1], H[i,2]*time[j]) + U[i,] %*% c(1, time[j])
			y[i,j] = location + rald(1, scale=sigma, p=tau)
		}
	}

	list(y = y, X = X, J=count)
}


######################################################################
# sim_multiple_QRdata = function(N, QR_fun=sim_longitudinal_data, tau){
# 	# N - number of data sets to generate
# 	# sur_fun - function to simulate survival data
# 	# longi_fun - function to sumualte longitudinal data
# 	# alpha - association mechanism for JM
# 	# tau - quantile

# 	outdata = vector(mode='list', N)
# 	for (i in 1:N){
# 		QR_data = QR_fun(tau=tau)
# 		#longi_data = longi_fun(sur_data, tau=tau)
# 		outdata[[i]] = QR_data
# 	}
# 	outdata
# }

# QR_data = sim_multiple_QRdata(30, tau=0.25)




###############################################################
#----- function to simulate multiple joint data sets ---------#
###############################################################
sim_multiple_data = function(N, sur_fun=sim_Ti, longi_fun=sim_longitudinal_data, alpha, tau){
	# N - number of data sets to generate
	# sur_fun - function to simulate survival data
	# longi_fun - function to sumualte longitudinal data
	# alpha - association mechanism for JM
	# tau - quantile

	outdata = vector(mode='list', N)
	for (i in 1:N){
		print(i)
		sur_data = sur_fun(alpha1=alpha[1], alpha2=alpha[2])
		longi_data = longi_fun(sur_data, tau=tau)
		outdata[[i]] = list(survival_data=sur_data, longitudinal_data=longi_data)
	}
	outdata
}

# single survival model
# QRJM_jags = function(data, I=250){
# 	setwd("/work/02784/myang3/QRJM/model")
# 	model = "/work/02784/myang3/QRJM/model/QRJM_survival00.txt"

# 	# load survival data
# 	Ti = data[['Ti']]
# 	# H = data[['H']]
# 	W = data[['W']]
# 	event =  data[['event']]
# 	zeros = rep(0,I)

# 	jags.data = list(W=W, Ti=Ti, event=event, I=I, zeros=zeros)
#  	jags.params = c("gamma","c")
#  	jags.inits = function(){ list(gamma=c(1,1),c=1) }
#   	jags(data=jags.data, inits=jags.inits, jags.params, n.iter=20000, n.burnin=10000, model.file=model, DIC=FALSE)
# }

# single quantile model
# QRJM_jags = function(data,tau,I=250){
# 	setwd("/work/02784/myang3/QRJM/model")
# 	model = "/work/02784/myang3/QRJM/model/QRJM_QR.txt"
# 	# file.show(model)
# 	# load longitudinal data
# 	y = data[['y']]
# 	X = data[['X']]
# 	H = data[['H']]
# 	J = data[['J']]

# 	jags.data = list(y=y, X=X, H=H, qt=tau, t=c(0, 1/4, 1/2, 3/4, 1, 3), I=I, J=J)
#  	jags.params = c("beta", "delta","sigma")
#  	jags.inits = function(){	list(beta=c(0.1,.1), delta=c(.1, .1), sigma=.1)
#   	}
#   jags(data=jags.data, inits=jags.inits, jags.params, n.iter=15000, n.burnin=5000, model.file=model)
# }

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

	jags.data = list(y=y, X=X, H=H, W=W, Ti=Ti, event=event, qt=tau, t=c(0, 1/4, 1/2, 3/4, 1, 3), I=I, J=J, zeros=zeros, zero=zero)
 	jags.params = c("beta", "delta","gamma","alpha1","alpha2","sigma", "c", "w11","w21","w22")
 	jags.inits = function(){	list(beta=c(0.1,0.1), delta=c(0.1, 0.1), gamma=c(0.1, 0.1), sigma=0.1, alpha1=0.1, alpha2=0.1, c=0.1)
  	}
  jags(data=jags.data, inits=jags.inits, jags.params, n.iter=25000, n.burnin=5000, model.file=model)
}



### run parallely on TACC
# setting up parallel computing in cluster
library(foreach)
library(doParallel)
jEOM = F
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

### 1
# sur_fit <- foreach(i=1:30,.packages=c("rjags","R2jags","foreach"),.verbose=T )%dopar%{
# 	set.seed(123)
#   	temp <- QRJM_jags(sur_data[[i]])
#   	temp$BUGSoutput$summary
# }

# setwd("/work/02784/myang3/QRJM/dataout/20141217")
# save( sur_fit,file="sur_fit00.Rdata")


### 2
# QR_fit <- foreach(i=1:30,.packages=c("rjags","R2jags","foreach"),.verbose=T )%dopar%{
# 	set.seed(123)
#   	temp <- QRJM_jags(QR_data[[i]], tau=0.25)
#   	temp$BUGSoutput$summary
# }

# setwd("/work/02784/myang3/QRJM/dataout/20141217")
# save( QR_fit,file="QR_fit.Rdata")

# simulate data
multi_data = sim_multiple_data(N=50, alpha=c(1,1), tau=0.5)

### 3
t1 = Sys.time()
qt50_50data_medianfit_thetahat_newpars <- foreach(i=1:50,.packages=c("rjags","R2jags","foreach"),.verbose=T )%dopar%{
	set.seed(123)
  	temp <- QRJM_jags(multi_data[[i]], tau=0.5)
  	temp$BUGSoutput$sims.matrix
}
Sys.time()-t1

# setwd("/work/02784/myang3/QRJM/dataout/JAGS_QRJM_0118")
setwd("/work/02784/myang3/QRJM_sim2/results/20150721/thetahat")
save(qt50_50data_medianfit_thetahat_newpars,file="qt50_50data_medianfit_thetahat_newpars.Rdata")


