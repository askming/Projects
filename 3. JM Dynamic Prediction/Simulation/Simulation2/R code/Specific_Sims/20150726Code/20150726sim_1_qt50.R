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
socketHosts_name <- rep(nodes_name,times=24)
cl <- makePSOCKcluster(socketHosts_name,homogeneous=TRUE)
clusterCall(cl, function() Sys.info()[c("nodename","machine")])
registerDoParallel(cl)
clusterExport(cl=cl,ls())
getDoParWorkers()
# stopCluster(cl)

setwd('/work/02784/myang3/QRJM_sim2/code')
source('sim2TACC_functions.R')


############################## generate data #######################################################
### sample size 600, 500 for inference and 100 for prediction
### number of data sets = 100
# sim_data_tau25_100 = prep_sim_data(N=100, sam_size=600, alpha=c(1,1),tau=0.25, beta=c(1, 1), delta=c(-1,1), gamma=c(-1,1))
# # sim_data_tau50_100 = prep_sim_data(N=100, sam_size=600, alpha=c(1,1),tau=0.5, beta=c(1, 1), delta=c(-1,1), gamma=c(-1,1))
# # sim_data_norm_100 = prep_sim_data(N=100, sam_size=600, alpha=c(1,1), beta=c(1, 1), delta=c(-1,1), gamma=c(-1,1), normal_data=TRUE)

### save the data
# setwd('/work/02784/myang3/QRJM_sim2/data/20150721')
# save(sim_data_tau25_100, file='sim_data_tau25_100.Rdata')
# save(sim_data_tau50_100, file='sim_data_tau50_100.Rdata')
# save(sim_data_norm_100, file='sim_data_norm_100.Rdata')


############################## 1. load the data ##############################
setwd('/work/02784/myang3/QRJM_sim2/data/20150721')
load(dir()[1])
input_data = lapply(sim_data_tau50_100, `[[`, 'selected')


############################## 2. fit the data using median regression for inference ##############################
t1 = Sys.time()
sim_data_tau50_100_medianfit_thetahat <- foreach(i=1:100,.packages=c("rjags","R2jags","foreach"),.verbose=T )%dopar%{
    set.seed(123)
    temp <- QRJM_jags(input_data[[i]], tau=0.5)
    temp$BUGSoutput$sims.matrix # there will be 3000 poterior samples for each parameter
}
Sys.time()-t1
setwd("/work/02784/myang3/QRJM_sim2/results/20150721/thetahat")
save(sim_data_tau50_100_medianfit_thetahat,file="sim_data_tau50_100_medianfit_thetahat_try2.Rdata")


############################## 3. prepare the validation data ##############################
setwd('/work/02784/myang3/QRJM_sim2/data/20150721')
load(dir()[1])
valid_data = lapply(sim_data_tau50_100 ,`[[`, 'valid')
# choose a time window and truncate the data
time = c(0,0.25)
curt_data = lapply(valid_data, select_valid_data, time=time)

# prepare posterior samples of theta from previous JM fitting
setwd("/work/02784/myang3/QRJM_sim2/results/20150721/thetahat")
load(dir()[3])
post_par = lapply(sim_data_tau50_100_medianfit_thetahat, as.data.frame)
### randomly choose 1k posterior samples out of 3k
sample1k = sample(seq(1:3e3), 1e3, replace=FALSE)
dput(sample1k, 'sample1k_dataqt50_qt25fit.txt')
sample1k = dget('sample1k_dataqt50_qt25fit.txt')
post_par = lapply(post_par, function(x) x[sample1k, ])
N = dim(post_par[[1]])[1] # number of posterior samples


############################## 4. make prediction of the RE for validation data ##############################
### split the data into two parts and predict the RE separately
### 1-50 part
t1 = Sys.time()
RE_all100qt50_medianfit_1_50 = vector(length = 50, mode='list')
for (i in 1:50){
    RE_all100qt50_medianfit_1_50[[i]] = QRJM_pred_re(curt_data[i], post_par[i], tau=0.50)
}
Sys.time()-t1
setwd("/work/02784/myang3/QRJM_sim2/results/20150721/pred_re")
save(RE_all100qt50_medianfit_1_50, file="RE_all100qt50_medianfit_1_50.Rdata")


### ------------------------ ###
### 51-100 part
t1 = Sys.time()
RE_all100qt50_medianfit_51_100 = vector(length = 50, mode='list')
for (i in 51:100){
    RE_all100qt50_medianfit_51_100[[i-50]] = QRJM_pred_re(curt_data[i], post_par[i], tau=0.50)
}
Sys.time()-t1
setwd("/work/02784/myang3/QRJM_sim2/results/20150721/pred_re")
save(RE_all100qt50_medianfit_51_100, file="RE_all100qt50_medianfit_51_100.Rdata")




##############################  5. calculate the gold standard survival probabilities for specific time window
gold_std_res = gold_std_pred_surv(indata=curt_data, t1=0.25, dt=0.25)
gold_std_res_025_1 = gold_std_pred_surv(indata=curt_data, t1=0.25, dt=1)
gold_std_res_025_2 = gold_std_pred_surv(indata=curt_data, t1=0.25, dt=2)

## save the gold standard
gold_std_qt50data = list(gold_std_res, gold_std_res_025_1, gold_std_res_025_2)
setwd("/work/02784/myang3/QRJM_sim2/results/20150721/pred_re")
save(gold_std_qt50data, file='gold_std_qt50data.Rdata')


############################## 6. make prediction of survival probabilities for different time windows
pred_surv_QRJM = function(data, par, pred_re, t, dt, gold_std){
    predicted = model_pred_surv(indata=data, t1=t, dt=dt, post_par=par, pred_re)
    res = unlist(lapply(predicted, colMeans))
    logit = function(a){log(a/(1-a))} # was incorrectly defined before
    list(mc_samples = res, mean_diff=mean(res-gold_std), mean_diff_logit = mean(logit(res)-logit(gold_std)))
}

### combine the prediction of RE into one data set, a list of 100
RE_all100qt50_medianfit = vector(length=100, mode='list')
for (i in 1:50) RE_all100qt50_medianfit[[i]] = RE_all100qt50_medianfit_1_50[[i]]
for (i in 51:100) RE_all100qt50_medianfit[[i]] = RE_all100qt50_medianfit_1_50[[i-50]]
setwd('/work/02784/myang3/QRJM_sim2/results/20150721/pred_re')
save(RE_all100qt50_medianfit, file="RE_all100qt50_medianfit.Rdata")

pre_surv_QRJM_025_025 = pred_surv_QRJM(data=curt_data, par=post_par, pred_re=RE_all100qt50_medianfit, t=0.25, dt=0.25, gold_std=gold_std_res)
pre_surv_QRJM_025_1 = pred_surv_QRJM(data=curt_data, par=post_par, pred_re=RE_all100qt50_medianfit, t=0.25, dt=1, gold_std=gold_std_res_025_1)
pre_surv_QRJM_025_2 = pred_surv_QRJM(data=curt_data, par=post_par, pred_re=RE_all100qt50_medianfit, t=0.25, dt=2, gold_std=gold_std_res_025_2)

QRJM_pred_qt50data_medianfit = list(pre_surv_QRJM_025_025, pre_surv_QRJM_025_1, pre_surv_QRJM_025_2)
setwd("/work/02784/myang3/QRJM_sim2/results/20150721/pred_re")
save(QRJM_pred_qt50data_medianfit, file='QRJM_pred_qt50data_medianfit.Rdata')


