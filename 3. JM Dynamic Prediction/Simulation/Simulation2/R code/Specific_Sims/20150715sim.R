
## simulation: fit qt=0.25 data with median regression to make the simulation study more complete
setwd('/work/02784/myang3/QRJM_sim2/code')
source('sim2TACC_functions.R')

## 1. load the data
setwd('/work/02784/myang3/QRJM_sim2/data/20150520')
load(dir()[1])
input_data = lapply(sim2_data_qt25_30, `[[`, 'selected')


## 2. fit the data using median regression for inference
t1 = Sys.time()
sim2_qt25_30_medianfit_thetahat <- foreach(i=1:30,.packages=c("rjags","R2jags","foreach"),.verbose=T )%dopar%{
    set.seed(123)
    temp <- QRJM_jags(input_data[[i]], tau=0.5)
    temp$BUGSoutput$sims.matrix # there will be 3000 poterior samples for each parameter
}
Sys.time()-t1
setwd("/work/02784/myang3/QRJM_sim2/results/20150520/thetahat")
save(sim2_qt25_30_medianfit_thetahat,file="sim2_qt25_30_medianfit_thetahat.Rdata")


### 3. prepare the validation data
setwd('/work/02784/myang3/QRJM_sim2/data/20150520')
load(dir()[1])
valid_data = lapply(sim2_data_qt25_30 ,`[[`, 'valid')
# choose a time window and truncate the data
time = c(0,0.25)
curt_data = lapply(valid_data, select_valid_data, time=time)
# prepare posterior samples of theta from previous JM fitting
post_par = lapply(sim2_qt25_30_medianfit_thetahat, as.data.frame)
N = dim(post_par[[1]])[1] # number of poterior samples


### 4. make prediction of the RE for validation data
RE_all30qt25_medianfit = vector(length = 30, mode='list')
for (i in 1:30){
    RE_all30qt25_medianfit[[i]] = QRJM_pred_re(curt_data[i], post_par[i], tau=0.5)
}
Sys.time()-t1
setwd("/work/02784/myang3/QRJM_sim2/results/20150520/pred_re")
save(RE_all30qt25_medianfit, file="RE_all30qt25_medianfit.Rdata")


### 5. calculate the gold standard survival probabilities for specific time window
gold_std_res = gold_std_pred_surv(indata=curt_data, t1=0.25, dt=0.25)
gold_std_res_025_1 = gold_std_pred_surv(indata=curt_data, t1=0.25, dt=1)
gold_std_res_025_2 = gold_std_pred_surv(indata=curt_data, t1=0.25, dt=2)

# gold_std_qt25data = list(gold_std_res, gold_std_res_025_1, gold_std_res_025_2)
# setwd("/work/02784/myang3/QRJM_sim2/results/20150520/pred_re")
# save(gold_std_qt25data, file='gold_std_qt25data.Rdata')


### 6. make prediction of survival probabilities for different time windows
pred_surv_QRJM = function(data, par, pred_re, t, dt, gold_std){
    predicted = model_pred_surv(indata=data, t1=t, dt=dt, post_par=par, pred_re)
    res = unlist(lapply(predicted, colMeans))
    logit = function(a){log(a/(1-a))} # was incorrectly defined before
    list(mc_samples = res, mean_diff=mean(res-gold_std), mean_diff_logit = mean(logit(res)-logit(gold_std)))
}

pre_surv_QRJM_025_025 = pred_surv_QRJM(data=curt_data, par=post_par, pred_re=RE_all30qt25_medianfit, t=0.25, dt=0.25, gold_std=gold_std_res)
pre_surv_QRJM_025_1 = pred_surv_QRJM(data=curt_data, par=post_par, pred_re=RE_all30qt25_medianfit, t=0.25, dt=1, gold_std=gold_std_res_025_1)
pre_surv_QRJM_025_2 = pred_surv_QRJM(data=curt_data, par=post_par, pred_re=RE_all30qt25_medianfit, t=0.25, dt=2, gold_std=gold_std_res_025_2)

QRJM_pred_qt25data_medianfit = list(pre_surv_QRJM_025_025, pre_surv_QRJM_025_1, pre_surv_QRJM_025_2)
setwd("/work/02784/myang3/QRJM_sim2/results/20150520/pred_re")
save(QRJM_pred_qt25data_medianfit, file='QRJM_pred_qt25data_medianfit.Rdata')


