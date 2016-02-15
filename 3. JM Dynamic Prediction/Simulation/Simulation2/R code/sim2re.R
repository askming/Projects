# to simulate random effect given the historical values of y

# input: posterior sample of parameters, histological longitudianl outcomes
# for each $\hat_{\theta}$, simulate a series of R.E. and take the mean  




QRJM_jags_simre = function(data,par,tau){
	library(R2jags)
	setwd("/Users/askming/Documents/github/Projects/1.\ QRJM/simulaiton/model\ files")
	model = "/Users/askming/Documents/github/Projects/1.\ QRJM/simulaiton/model\ files/QRJM_Cholesky.txt"
	# file.show(model)
	# data variable should include both yi(t) (up to t) and t (survival time)
	# par is the posterior sample of parameters, a list
	# load longitudinal data
	y = data[[2]][['y']]
	X = data[[2]][['X']]
	J = length(y)
	# time = c(0,2,4,6)
	
	# load survival data
	t = data[[1]][['Ti']]
	H = data[[1]][['H']]
	W = data[[1]][['W']]

	zeros = rep(0, I)
	zero = c(0,0)

	jags.data = list(y=y, X=X, H=H, W=W, t=t, event=0, qt=tau, t=time, J=J, zeros=zeros, zero=zero)
 	# jags.params = c("beta", "delta","gamma","alpha1","alpha2","sigma", "c", "w11", "w21", "w22")
 	jags.params = c("u")
 	jags.inits = function(){ list(u=c(1,1)) }	
  jags(data=jags.data, inits=jags.inits, jags.params, n.iter=50, n.burnin=0, model.file=model)	
}



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

system.time(
  sur_fit <- foreach(i=1:1000,.packages=c("rjags","R2jags","foreach"),.verbose=T )%dopar%{ 
  	set.seed(123)
  	temp <- QRJM_jags_simre(data, par[[i]], tau=0.5)
  	temp$BUGSoutput
  	}  
)





