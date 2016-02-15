### code to run simulation on TACC ##
#   From Rizopoulos:
#     - simulate data and fit the joint model
#     - took a random sample of patients
#     - these patients have measurements up to t
#     - using their available measurements calculate the random effects
#     - then calculate predictions for u > t, with u \in (4, 6, 8, 10, 14, 16,
# 20, 22, 26)




setwd("/work/02784/myang3/QRJM_sim2/code")
# load pre-defined functions
source("sim2TACC_functions.R")


## create multiple nodes for parallel
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


# a=sim_Ti(n=500, alpha1=1, alpha2=1, delta=c(-2,1), gamma=c(-3,-1))
### generate data and save the data
# step to sumulate data
# t1 = Sys.time()

### Following may not work
# t1 = Sys.time()
# sim2_data_qt50_30 = prep_sim_data(N=30, alpha=c(2,2),tau=0.5, beta=c(2,-2), delta=c(-4,2), gamma=c(-3,-1))
# Sys.time()-t1
# # save the simulated data
# setwd("/work/02784/myang3/QRJM_sim2/results/20150618")
# save(sim2_data_qt25_30, file="sim2_data_qt25_30.Rdata")

### simulate quantile data
t1 = Sys.time()
sim2_data_qt50_30 = prep_sim_data(N=30, alpha=c(1,1),tau=0.5, beta=c(1, 1), delta=c(-1,1), gamma=c(-1,1))
Sys.time()-t1
# save the simulated data
setwd("/work/02784/myang3/QRJM_sim2/data/20150618")
save(sim2_data_qt50_30, file="sim2_data_qt50_30.Rdata")

input_data = lapply(sim2_data_qt50_30, `[[`, 'selected')

### simulate normal random error data
t1 = Sys.time()
sim2_normdata_30 = prep_sim_data(N=30, alpha=c(1,1), beta=c(1, 1), delta=c(-1,1), gamma=c(-1,1), normal_data=TRUE)
Sys.time()-t1
# save the simulated data
setwd("/work/02784/myang3/QRJM_sim2/data/20150618")
save(sim2_normdata_30, file="sim2_normdata_30.Rdata")

input_data = lapply(sim2_normdata_30,`[[`, 'selected')



### 1. fit model using original quantile model and input_data & save the results
# input_data is the data used to fit the model, valid_data is the validation data

t1 = Sys.time()
sim2_qt25_30_thetahat <- foreach(i=1:30,.packages=c("rjags","R2jags","foreach"),.verbose=T )%dopar%{
    set.seed(123)
    temp <- QRJM_jags(input_data[[i]], tau=0.25)
    temp$BUGSoutput$sims.matrix # there will be 3000 poterior samples for each parameter
}
Sys.time()-t1
setwd("/work/02784/myang3/QRJM_sim2/results/20150520")
save(sim2_qt25_30_thetahat,file="1.sim2_qt25_30_thetahat.Rdata")


## 2.1 fit model using median regression and input_data & save the results
t1 = Sys.time()
sim2_norm_30_50fit_thetahat <- foreach(i=1:30,.packages=c("rjags","R2jags","foreach"),.verbose=T )%dopar%{
    set.seed(123)
    temp <- QRJM_jags(input_data[[i]], tau=0.5)
    temp$BUGSoutput$sims.matrix # there will be 3000 poterior samples for each parameter
}
Sys.time()-t1
setwd("/work/02784/myang3/QRJM_sim2/results/20150618/thetahat")
save(sim2_norm_30_50fit_thetahat,file="sim2_norm_30_medianfit_theta_hat.Rdata")



## 2.2 fit model using qt=0.75 regression and input_data & save the results
t1 = Sys.time()
sim2_25_30_75fit_thetahat <- foreach(i=1:30,.packages=c("rjags","R2jags","foreach"),.verbose=T )%dopar%{
    set.seed(123)
    temp <- QRJM_jags(input_data[[i]], tau=0.75)
    temp$BUGSoutput$sims.matrix # there will be 3000 posterior samples for each parameter
}
Sys.time()-t1
setwd("/work/02784/myang3/QRJM_sim2/results/pred")
save(sim2_25_30_75fit_thetahat,file="4.sim2_25_30_75fit_theta_hat.Rdata")



## 3. fit model using mean regression and input_data & save the reuslts
t1 = Sys.time()
sim2_50_30_meanfit_theta_hat <- foreach(i=1:30,.packages=c("rjags","R2jags","foreach"),.verbose=T )%dopar%{
    set.seed(123)
    temp <- NormJM_jags(input_data[[i]])
    temp$BUGSoutput$sims.matrix # there will be 3000 poterior samples for each parameter
}
Sys.time()-t1
setwd("/work/02784/myang3/QRJM_sim2/results/20150618/thetahat")
save(sim2_50_30_meanfit_theta_hat,file="sim2_50_30_meanfit_theta_hat.Rdata")




#---------------------------------------------------------------#
#################################################################
#---------------------------------------------------------------#


### make predictions based infrence and save the results

# load valid_data
setwd("/work/02784/myang3/QRJM_sim2/data/20150618")
load(dir()[1])
valid_data = lapply(sim2_data_qt50_30 ,`[[`, 'valid')

# choose a time window and truncate the data
time = c(0,0.25)
curt_data = lapply(valid_data, select_valid_data, time=time)

# load poterior samples of theta from previous JM fitting
setwd("/work/02784/myang3/QRJM_sim2/results/20150520")
load(dir()[1:3])
post_par1= lapply(sim2_qt25_30_thetahat, as.data.frame)
post_par2= lapply(sim2_25_30_50fit_thetahat, as.data.frame)
post_par3 = lapply(sim2_25_30_75fit_thetahat, as.data.frame)
post_par4 = lapply(sim2_25_30_meanfit_thetahat, as.data.frame)
N = dim(post_par2[[1]])[1] # number of poterior samples

post_par2= lapply(sim2_norm_30_50fit_thetahat, as.data.frame)
post_par4 = lapply(sim2_norm_30_meanfit_thetahat, as.data.frame)
post_par4 = lapply(sim2_50_30_meanfit_theta_hat, as.data.frame)
#  for EACH validation data set do the following:
#  for each subject in the data set, loop over the posterior samples of theta; for each posterior sample get a prediction of the ramdom effects

# t1 = Sys.time()
# for (j in 1:dim(curt_data[[1]]$longitudina_data$y)[1]){#iterate over subjects in data
#     assign(paste('sim2_pred_q25_subject', j, sep=''), foreach(i=1:N,.packages=c("rjags","R2jags","foreach"),.verbose=T )%dopar%{ # iterate over posterior sampels
#     set.seed(123)
#     temp <- QRJM_jags_pred(curt_data[[1]], time=time, tau=0.25, par=post_par[[1]][i,], id=j)
#     temp$BUGSoutput$sims.matrix # collect predictions random effects
#         }
#     )
# }



# length(gold_std_res)
# [1] 409

#  predict random effects from posterior samples of parameters and the data

## for single data set
## QRJM tau=0.25, original one ###
t2=Sys.time()
RE1 = QRJM_pred_re(curt_data[1], post_par1[1], tau=0.25)
Sys.time()-t2

t1 = Sys.time()
RE_all30_25fit_t025 = vector(length=30, mode='list')
for (i in 1:30){
    RE_all30_25fit_t025[[i]] = QRJM_pred_re(curt_data[i], post_par1[i], tau=0.25)
}
Sys.time()-t1
setwd("/work/02784/myang3/QRJM_sim2/results/20150520/pred_re")
save(RE_all30_25fit_t025, file="RE_all30_25fit_t025.Rdata")



#### QRJM, median ####
t1 = Sys.time()
RE_all30norm_medianfit = vector(length=30, mode='list')
for (i in 1:30){
    RE_all30norm_medianfit[[i]] = QRJM_pred_re(curt_data[i], post_par2[i], tau=0.5)
}
Sys.time()-t1
## 3.97 hours
setwd("/work/02784/myang3/QRJM_sim2/results/20150618/pred_re")
save(RE_all30norm_medianfit, file="RE_all30norm_medianfit.Rdata")

#### QRJM, tau=0.75 ####

# complete data sets
t3=Sys.time()
RE_all30_75fit_t4 = vector(length=30, mode='list')
for (i in 1:30){
    RE_all30_75fit_t4[[i]] = QRJM_pred_re(curt_data[i], post_par3[i], tau=0.75)
    }
Sys.time()-t3

setwd("/work/02784/myang3/QRJM_sim2/results/pred")
save(RE_all30_75fit_t4, file="RE_all30_75fit_t4.Rdata")

#### LMJM ####
t4=Sys.time()
RE_all30_meanfit = vector(length=30, mode='list')
for (i in 1:30){
    RE_all30_meanfit[[i]] = LMJM_pred_re(curt_data[i], post_par4[i])
    }
Sys.time()-t4
setwd("/work/02784/myang3/QRJM_sim2/results/20150618/pred_re")
save(RE_all30_meanfit, file="RE_all30_meanfit.Rdata")



#  The gold standard survival probabilities for specific time window
gold_std_res = gold_std_pred_surv(indata=curt_data, t1=0.25, dt=0.25)
gold_std_res_025_1 = gold_std_pred_surv(indata=curt_data, t1=0.25, dt=1)
gold_std_res_025_2 = gold_std_pred_surv(indata=curt_data, t1=0.25, dt=2)

gold_std_qt50data = list(gold_std_res, gold_std_res_025_1, gold_std_res_025_2)
setwd("/work/02784/myang3/QRJM_sim2/results/20150618/pred_re")
save(gold_std_qt50data, file='gold_std_qt50data.Rdata')

# pred_surv(indata=curt_data, dataid=1, post_par=post_par[1], RE=REs, id=11, t=2)
#  QR model, qt=0.25. The predicted survival probabilities for specific time window
pred_surv_QRJM = function(data, par, pred_re, t, dt, gold_std){
    predicted = model_pred_surv(indata=data, t1=t, dt=dt, post_par=par, pred_re)
    res = unlist(lapply(predicted, colMeans))
    logit = function(a){log(a/(1-a))} # was incorrectly defined before
    list(mc_samples = res, mean_diff=mean(res-gold_std), mean_diff_logit = mean(logit(res)-logit(gold_std)))
}

pre_surv_QRJM_025_025 = pred_surv_QRJM(data=curt_data, par=post_par2, pred_re=RE_all30norm_medianfit, t=0.25, dt=0.25, gold_std=gold_std_res)
pre_surv_QRJM_025_1 = pred_surv_QRJM(data=curt_data, par=post_par2, pred_re=RE_all30norm_medianfit, t=0.25, dt=1, gold_std=gold_std_res_025_1)
pre_surv_QRJM_025_2 = pred_surv_QRJM(data=curt_data, par=post_par2, pred_re=RE_all30norm_medianfit, t=0.25, dt=2, gold_std=gold_std_res_025_2)

QRJM_pred_normdata_medianfit = list(pre_surv_QRJM_025_025, pre_surv_QRJM_025_1, pre_surv_QRJM_025_2)
setwd("/work/02784/myang3/QRJM_sim2/results/20150618/pred_re")
save(QRJM_pred_normdata_medianfit, file='QRJM_pred_surv_prob_normdata_medianfit.Rdata')

# pre_surv_data30_25 = model_pred_surv(indata=curt_data, t1=2, dt=2, post_par=post_par, pred_re=RE_all30)
# pre_surv_data30_25_means = lapply(pre_surv_data30_25, colMeans)
# pre_surv_25 = unlist(pre_surv_data30_25_means)

# #  QR model, qt=0.5 The predicted survival probabilities for specific time window
# pre_surv_data30_median = model_pred_surv(indata=curt_data, t1=2, dt=2, post_par=post_par2, pred_re=RE_all30)
# pre_surv_data30_median_means = lapply(pre_surv_data30_median, colMeans)
# pre_surv_median = unlist(pre_surv_data30_median_means)

## mean fit, to predicvt survival probabilities
pred_surv_LMJM = function(data, par, pred_re, t, dt, gold_std){
    predicted = model_pred_surv(indata=data, t1=t, dt=dt, post_par=par, pred_re)
    res = unlist(lapply(predicted, colMeans))
    logit = function(a){log(a/(1-a))}# was incorretly defined before
    list(mc_samples = res, mean_diff=mean(res-gold_std), mean_diff_logit = mean(logit(res)-logit(gold_std)))
}

pre_surv_LMJM_025_025 = pred_surv_LMJM(data=curt_data, par=post_par4, pred_re=RE_all30_meanfit, t=0.25, dt=0.25, gold_std=gold_std_res)
pre_surv_LMJM_025_1 = pred_surv_LMJM(data=curt_data, par=post_par4, pred_re=RE_all30_meanfit, t=0.25, dt=1, gold_std=gold_std_res_025_1)
pre_surv_LMJM_025_2 = pred_surv_LMJM(data=curt_data, par=post_par4, pred_re=RE_all30_meanfit, t=0.25, dt=2, gold_std=gold_std_res_025_2)

LMJM_pred_qt50data = list(pre_surv_LMJM_025_025, pre_surv_LMJM_025_1, pre_surv_LMJM_025_2)
setwd("/work/02784/myang3/QRJM_sim2/results/20150618/pred_re")
save(LMJM_pred_qt50data, file='LMJM_pred_surv_prob_qt50_meanfit.Rdata')



