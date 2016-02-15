#### set up parallel R
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

### Load R functions for following procedure
setwd('/work/02784/myang3/QRJM_sim2/code')
source('sim2TACC_functions.R')

#### try to redefine the prediction function and see if it can be faster
#### in the output, each list stands for a data, in which each list stands for one RE for each subject in the list
# QRJM_pred_re = function(indata, post_par, tau){
#     # indata - validation data sets with qualified subjects only, a list, each element in the list is a list as well that contains two lists of longitudinal and survival data
#     # post_par - posterior samples of the parameters, a list
#     N = length(indata) # number of multiple simulated data
#     K = dim(post_par[[1]])[1] # number of posterior samples
#     predicted_RE = vector(length=N, mode="list")

#     ### functions to collect predicted RE ###
#     collect_re = function(pred_data){
#         re = unlist(lapply(pred_data, colMeans))
#         u1 = re[names(re)=='u[1]']
#         u2 = re[names(re)=='u[2]']
#         cbind(u1, u2) # output a matrix of K times 2 MC samples of REs
#     }

#     for (i in 1:N){
#         for (k in 1:K){
#             current_re = foreach(j = 1:length(indata[[i]][['longitudina_data']]$id),.packages=c("rjags","R2jags","foreach"),.verbose=T, .export=c("QRJM_jags_pred"))%dopar%{ # iterate over posterior sampels
#                     set.seed(123)
#                     temp <- QRJM_jags_pred(indata[[i]], tau=tau, par=post_par[[i]][k,], id=j)
#                     temp$BUGSoutput$sims.matrix# collect predictions random effects
#             }
#         current_re = collect_re(current_re) # a matrix of K times 2
#         predicted_RE[[i]][[k]]= current_re
#         }
#     }
#     return(predicted_RE)
# }







# ############################## generate data #######################################################
# ### sample size 600, 500 for inference and 100 for prediction
# ### number of data sets = 100
# sim_data_tau25_100 = prep_sim_data(N=100, sam_size=600, alpha=c(1,1),tau=0.25, beta=c(1, 1), delta=c(-1,1), gamma=c(-1,1))

# # sim_data_tau50_100 = prep_sim_data(N=100, sam_size=600, alpha=c(1,1),tau=0.5, beta=c(1, 1), delta=c(-1,1), gamma=c(-1,1))

# # sim_data_norm_100 = prep_sim_data(N=100, sam_size=600, alpha=c(1,1), beta=c(1, 1), delta=c(-1,1), gamma=c(-1,1), normal_data=TRUE)

# ### save the data
# setwd('/work/02784/myang3/QRJM_sim2/data/20150721')
# save(sim_data_tau25_100, file='sim_data_tau25_100.Rdata')
# # save(sim_data_tau50_100, file='sim_data_tau50_100.Rdata')
# # save(sim_data_norm_100, file='sim_data_norm_100.Rdata')

## 1. load the data
setwd('/work/02784/myang3/QRJM_sim2/data/20150721')
load(dir()[1])
input_data = lapply(sim_data_tau25_100, `[[`, 'selected')


## 2. fit the data using median regression for inference
t1 = Sys.time()
sim_data_tau25_100_qt25fit_thetahat_fewpost <- foreach(i=1:100,.packages=c("rjags","R2jags","foreach"),.verbose=T )%dopar%{
    set.seed(123)
    temp <- QRJM_jags(input_data[[i]], tau=0.25)
    temp$BUGSoutput$sims.matrix # there will be 900 poterior samples for each parameter
}
Sys.time()-t1
setwd("/work/02784/myang3/QRJM_sim2/results/20150721/thetahat")
save(sim_data_tau25_100_qt25fit_thetahat_fewpost,file="sim_data_tau25_100_qt25fit_thetahat_fewpost.Rdata")


### 3. prepare the validation data
setwd('/work/02784/myang3/QRJM_sim2/data/20150721')
load(dir()[1])
valid_data = lapply(sim_data_tau25_100 ,`[[`, 'valid')
# choose a time window and truncate the data
time = c(0,0.25)
curt_data = lapply(valid_data, select_valid_data, time=time)

# prepare posterior samples of theta from previous JM fitting
setwd("/work/02784/myang3/QRJM_sim2/results/20150721/thetahat")
load(dir()[3])
post_par = lapply(sim_data_tau25_100_qt25fit_thetahat, as.data.frame)
### randomly choose 1k posterior samples out of 3k
# sample1k = sample(seq(1:3e3), 1e3, replace=FALSE)
# dput(sample1k, 'sample1k_qt25.txt')
sample1k = dget('sample1k_qt25.txt')
post_par = lapply(post_par, function(x) x[sample1k, ])
N = dim(post_par[[1]])[1] # number of posterior samples


### 4. make prediction of the RE for validation data
### split the data into two parts and predict the RE separately
### 1-30 part
t1 = Sys.time()
RE_all100qt25_qt25fit_1_30 = vector(length = 30, mode='list')
for (i in 1:30){
    RE_all100qt25_qt25fit_1_30[[i]] = QRJM_pred_re(curt_data[i], post_par[i], tau=0.25)
}
Sys.time()-t1
setwd("/work/02784/myang3/QRJM_sim2/results/20150721/pred_re")
save(RE_all100qt25_qt25fit_1_30, file="RE_all100qt25_qt25fit_1_30.Rdata")


### ------------------------ ###
### 31-65 part
RE_all100qt25_qt25fit_31_65 = vector(length = 35, mode='list')
for (i in 31:65){
    RE_all100qt25_qt25fit_31_65[[i-30]] = QRJM_pred_re(curt_data[i], post_par[i], tau=0.25)
}
Sys.time()-t1
setwd("/work/02784/myang3/QRJM_sim2/results/20150721/pred_re")
save(RE_all100qt25_qt25fit_31_65, file="RE_all100qt25_qt25fit_31_65.Rdata")


### ------------------------ ###
### 66-100 part
RE_all100qt25_qt25fit_66_100 = vector(length = 35, mode='list')
for (i in 66:100){
    RE_all100qt25_qt25fit_66_100[[i-65]] = QRJM_pred_re(curt_data[i], post_par[i], tau=0.25)
}
Sys.time()-t1
setwd("/work/02784/myang3/QRJM_sim2/results/20150721/pred_re")
save(RE_all100qt25_qt25fit_66_100, file="RE_all100qt25_qt25fit_66_100.Rdata")




### 5. calculate the gold standard survival probabilities for specific time window
gold_std_res = gold_std_pred_surv(indata=curt_data, t1=0.25, dt=0.25)
gold_std_res_025_1 = gold_std_pred_surv(indata=curt_data, t1=0.25, dt=1)
gold_std_res_025_2 = gold_std_pred_surv(indata=curt_data, t1=0.25, dt=2)

## save the gold standard
# gold_std_qt25data = list(gold_std_res, gold_std_res_025_1, gold_std_res_025_2)
# setwd("/work/02784/myang3/QRJM_sim2/results/20150721/pred_re")
# save(gold_std_qt25data, file='gold_std_qt25data.Rdata')


### 6. make prediction of survival probabilities for different time windows
pred_surv_QRJM = function(data, par, pred_re, t, dt, gold_std){
    predicted = model_pred_surv(indata=data, t1=t, dt=dt, post_par=par, pred_re)
    res = unlist(lapply(predicted, colMeans))
    logit = function(a){log(a/(1-a))} # was incorrectly defined before
    list(mc_samples = res, mean_diff=mean(res-gold_std), mean_diff_logit = mean(logit(res)-logit(gold_std)))
}

pre_surv_QRJM_025_025 = pred_surv_QRJM(data=curt_data, par=post_par, pred_re=RE_all100qt25_qt25fit, t=0.25, dt=0.25, gold_std=gold_std_res)
pre_surv_QRJM_025_1 = pred_surv_QRJM(data=curt_data, par=post_par, pred_re=RE_all100qt25_qt25fit, t=0.25, dt=1, gold_std=gold_std_res_025_1)
pre_surv_QRJM_025_2 = pred_surv_QRJM(data=curt_data, par=post_par, pred_re=RE_all100qt25_qt25fit, t=0.25, dt=2, gold_std=gold_std_res_025_2)

QRJM_pred_qt25data_qt25fit = list(pre_surv_QRJM_025_025, pre_surv_QRJM_025_1, pre_surv_QRJM_025_2)
setwd("/work/02784/myang3/QRJM_sim2/results/20150721/pred_re")
save(QRJM_pred_qt25data_qt25fit, file='QRJM_pred_qt25data_qt25fit.Rdata')


