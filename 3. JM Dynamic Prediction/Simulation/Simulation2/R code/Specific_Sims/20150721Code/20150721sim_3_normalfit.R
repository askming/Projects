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
socketHosts_name <- rep(nodes_name,times=60)
cl <- makePSOCKcluster(socketHosts_name,homogeneous=TRUE)
clusterCall(cl, function() Sys.info()[c("nodename","machine")])
registerDoParallel(cl)
clusterExport(cl=cl,ls())


setwd('/work/02784/myang3/QRJM_sim2/code')
source('sim2TACC_functions.R')


# ############################## generate data ######################################

############################## 1. load the data ##############################
setwd('/work/02784/myang3/QRJM_sim2/data/20150721')
load(dir()[1])
input_data = lapply(sim_data_tau25_100, `[[`, 'selected')


library(R2jags)
############################## 2. fit the data using median regression for inference ##############################
t1 = Sys.time()
sim_data_tau25_100_normfit_thetahat <- foreach(i=1:100,.packages=c("rjags","R2jags","foreach"),.verbose=T )%dopar%{
    set.seed(123)
    temp <- NormJM_jags(input_data[[i]])
    temp$BUGSoutput$sims.matrix # there will be 3000 poterior samples for each parameter
}
Sys.time()-t1
setwd("/work/02784/myang3/QRJM_sim2/results/20150721/thetahat")
save(sim_data_tau25_100_normfit_thetahat,file="sim_data_tau25_100_normfit_thetahat.Rdata")


############################## 3. prepare the validation data ##############################
setwd('/work/02784/myang3/QRJM_sim2/data/20150721')
load(dir()[1])
valid_data = lapply(sim_data_tau25_100 ,`[[`, 'valid')
# choose a time window and truncate the data
time = c(0,0.25)
curt_data = lapply(valid_data, select_valid_data, time=time)

# prepare posterior samples of theta from previous JM fitting
setwd("/work/02784/myang3/QRJM_sim2/results/20150721/thetahat")
load(dir()[2])
post_par = lapply(sim_data_tau25_100_normfit_thetahat, as.data.frame)
### randomly choose 1k posterior samples out of 3k
# sample1k = sample(seq(1:3e3), 1e3, replace=FALSE)
# dput(sample1k, 'sample1k_normfit.txt')
sample1k = dget('sample1k_normfit.txt')
post_par = lapply(post_par, function(x) x[sample1k, ])
N = dim(post_par[[1]])[1] # number of posterior samples


############################## 4. make prediction of the RE for validation data ##############################
### split the data into two parts and predict the RE separately
### 1-30 part
t1 = Sys.time()
RE_all100qt25_meanfit_1_30 = vector(length = 30, mode='list')
for (i in 1:30){
	print(i)
    RE_all100qt25_meanfit_1_30[[i]] = LMJM_pred_re(curt_data[i], post_par[i])
}
Sys.time()-t1
setwd("/work/02784/myang3/QRJM_sim2/results/20150721/pred_re")
save(RE_all100qt25_meanfit_1_30, file="RE_all100qt25_meanfit_1_30.Rdata")


### ------------------------ ###
### 31-65 part
t1 = Sys.time()
RE_all100qt25_meanfit_31_65 = vector(length = 35, mode='list')
for (i in 31:65){
    RE_all100qt25_meanfit_31_65[[i-30]] = LMJM_pred_re(curt_data[i], post_par[i])
}
Sys.time()-t1
setwd("/work/02784/myang3/QRJM_sim2/results/20150721/pred_re")
save(RE_all100qt25_meanfit_31_65, file="RE_all100qt25_meanfit_31_65.Rdata")


### ------------------------ ###
### 66-100 part
t1 = Sys.time()
RE_all100qt25_meanfit_66_100 = vector(length = 35, mode='list')
for (i in 66:100){
    RE_all100qt25_meanfit_66_100[[i-65]] = LMJM_pred_re(curt_data[i], post_par[i])
}
Sys.time()-t1
setwd("/work/02784/myang3/QRJM_sim2/results/20150721/pred_re")
save(RE_all100qt25_meanfit_66_100, file="RE_all100qt25_meanfit_66_100.Rdata")





############################## 5. calculate the gold standard survival probabilities for specific time window
# gold_std_res = gold_std_pred_surv(indata=curt_data, t1=0.25, dt=0.25)
# gold_std_res_025_1 = gold_std_pred_surv(indata=curt_data, t1=0.25, dt=1)
# gold_std_res_025_2 = gold_std_pred_surv(indata=curt_data, t1=0.25, dt=2)

## load the gold standard
setwd("/work/02784/myang3/QRJM_sim2/results/20150721/pred_re")
load(dir()[1])


############################## 6. make prediction of survival probabilities for different time windows
pred_surv_LMJM = function(data, par, pred_re, t, dt, gold_std){
    predicted = model_pred_surv(indata=data, t1=t, dt=dt, post_par=par, pred_re)
    res = unlist(lapply(predicted, colMeans))
    logit = function(a){log(a/(1-a))}# was incorretly defined before
    list(mc_samples = res, mean_diff=mean(res-gold_std), mean_diff_logit = mean(logit(res)-logit(gold_std)))
}

pre_surv_LMJM_025_025 = pred_surv_LMJM(data=curt_data, par=post_par, pred_re=RE_all100qt25_meanfit, t=0.25, dt=0.25, gold_std=gold_std_res)
pre_surv_LMJM_025_1 = pred_surv_LMJM(data=curt_data, par=post_par, pred_re=RE_all100qt25_meanfit, t=0.25, dt=1, gold_std=gold_std_res_025_1)
pre_surv_LMJM_025_2 = pred_surv_LMJM(data=curt_data, par=post_par, pred_re=RE_all100qt25_meanfit, t=0.25, dt=2, gold_std=gold_std_res_025_2)

QRJM_pred_qt25data_meanfit = list(pre_surv_QRJM_025_025, pre_surv_QRJM_025_1, pre_surv_QRJM_025_2)
setwd("/work/02784/myang3/QRJM_sim2/results/20150721/pred_re")
save(QRJM_pred_qt25data_meanfit, file='QRJM_pred_qt25data_meanfit.Rdata')


