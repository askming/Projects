

NormJM_jags = function(data){
	library(R2jags)
	setwd("/work/02784/myang3/QRJM/model")
	model = "/work/02784/myang3/QRJM/model/NormalJM_Cholesky.txt" #JM using LMM
	# file.show(model)
	# load longitudinal data
	y = data[[2]][['y']]
	X = data[[2]][['X']]
	J = data[[2]][['J']]
	I = length(y)
	# load survival data
	Ti = data[[1]][['Ti']]
	H = data[[1]][['H']]
	W = data[[1]][['W']]
	event =  data[[1]][['event']]

	# alpha2 = 1
	zeros = rep(0, I)
	zero = c(0,0)

	jags.data = list(y=y, X=X, H=H, W=W, Ti=Ti, event=event, t=c(0, 2, 4, 6, 8, 10), I=I, J=J, zeros=zeros, zero=zero)
 	jags.params = c("beta", "delta","gamma","alpha1","alpha2","sigma", "c", "w11","w21","w22")
 	jags.inits = function(){	list(beta=c(0.1,0.1), delta=c(0.1, 0.1), gamma=c(0.1, 0.1), sigma=0.1, alpha1=0.1, alpha2=0.1, c=0.1)
  	}
  jags(data=jags.data, inits=jags.inits, jags.params, n.iter=15000, n.burnin=5000, model.file=model)
}



### run parallely on TACC
# setting up parallel computing in cluster
# library(foreach)
# library(doParallel)
# jEOM = F
# raw_PE <- switch(jEOM+1, F = "localhost.ls4.tacc.utexas.edu", T = system("cat $PE_HOSTFILE",intern=T))  # or use $HOSTNAME ,`hostname -s`
# catt('raw_PE is',raw_PE)
# nodes_name <- sub("^([\\w-]+)\\..*$","\\1",raw_PE,perl=T)
# nodes_slots <- switch(jEOM+1, F = as.numeric(system("cat /proc/cpuinfo | grep processor|wc -l",intern=T)),T = as.numeric( sub("^([\\w-]+)\\..* (\\d+?) .*$","\\2",raw_PE,perl=T) ) )
# #---------------------------------------#
# ######################################
# ## code confirmed correct.
# socketHosts_name <- rep(nodes_name,times=nodes_slots)
# cl <- makePSOCKcluster(socketHosts_name,homogeneous=TRUE)
# clusterCall(cl, function() Sys.info()[c("nodename","machine")])
# registerDoParallel(cl)
# clusterExport(cl=cl,ls())



### 3 run the model and save the output
setwd("/work/02784/myang3/QRJM_sim2/in_data")
load("in_data_25_12.Rdata")
input_data = lapply(sel_sim_data_25_12 ,`[[`, 'selected')
t1 = Sys.time()
sim2_25_norm <- foreach(i=1:12,.packages=c("rjags","R2jags","foreach"),.verbose=T )%dopar%{
	set.seed(123)
  	temp <- NormJM_jags(input_data[[i]])
  	temp$BUGSoutput$sims.matrix # there will be 1000 poterior samples for each parameter
}
Sys.time()-t1

setwd("/work/02784/myang3/QRJM_sim2/results/pred")
save(sim2_25_norm,file="sim2_25_12_thetahat_normal.Rdata")


