rm(list=ls())
################## function to prepare data for comparison ########################
pred_data = function(gold_std, pred){
    l = length(gold_std)
    y = c(gold_std, pred)
    data = as.data.frame(y)
    data$meth=c(rep("goldstd", l), rep("pred", l))
    data$item = c(seq(1, l), seq(1, l))
    data$rep = c(rep(1,l*2))
    return(data)
}
#################################################################################


###################### function to draw BA plot #################################
plot_BA = function(data, ...){
    library(MethComp)
    data_to_plot = Meth(data)
    par(mar=c(6,4,4,5))
    BA.plot(data_to_plot,cex.point=0.5, diflim=c(-0.9,0.9), xlab="mean", ylab="difference", ...)
}
#################################################################################

####################### function to calculate MSE and bias ########################
pred_mse_bias = function(gold_std, pred){
    mse = mean((pred - gold_std)^2)
    bias = mean(pred - gold_std)
    out = cbind(mse, bias)
    colnames(out) = c('MSE', 'bias')
    return(out)
}
#################################################################################


############### Simulation results on 07/21/2015 ##############################
#### -------------- data generated from qt=0.25 error ------------------------- ###
setwd('/Users/askming/Dropbox/Research\ works/Self/Phd\ thesis/Simulation\ reports/pred_surv_qt25data')
load(dir()[1])
load(dir()[2])
load(dir()[3])
load(dir()[4])

gold_std1 = gold_std_qt25data[[1]]
gold_std2 = gold_std_qt25data[[2]]
gold_std3 = gold_std_qt25data[[3]]

# normal model
pred1 = QRJM_pred_qt25data_meanfit[[1]][[1]]
pred2 = QRJM_pred_qt25data_meanfit[[2]][[1]]
pred3 = QRJM_pred_qt25data_meanfit[[3]][[1]]

mse_bias_LMJMpred1 = pred_mse_bias(gold_std1, pred1)
mse_bias_LMJMpred2 = pred_mse_bias(gold_std2, pred2)
mse_bias_LMJMpred3 = pred_mse_bias(gold_std3, pred3)
round(rbind(mse_bias_LMJMpred1, mse_bias_LMJMpred2, mse_bias_LMJMpred3),3)
#        MSE   bias
# [1,] 0.003 -0.005
# [2,] 0.011 -0.033
# [3,] 0.024 -0.053

# qt=25 model
predQR1 = QRJM_pred_qt25data_qt25fit[[1]][[1]]
predQR2 = QRJM_pred_qt25data_qt25fit[[2]][[1]]
predQR3 = QRJM_pred_qt25data_qt25fit[[3]][[1]]

mse_bias_QRJMpred1 = pred_mse_bias(gold_std1, predQR1)
mse_bias_QRJMpred2 = pred_mse_bias(gold_std2, predQR2)
mse_bias_QRJMpred3 = pred_mse_bias(gold_std3, predQR3)
round(rbind(mse_bias_QRJMpred1, mse_bias_QRJMpred2, mse_bias_QRJMpred3),3)
#        MSE   bias
# [1,] 0.003  0.011
# [2,] 0.005  0.002
# [3,] 0.009 -0.005

# qt=50 model
predQR1_median = QRJM_pred_qt25data_medianfit[[1]][[1]]
predQR2_median = QRJM_pred_qt25data_medianfit[[2]][[1]]
predQR3_median = QRJM_pred_qt25data_medianfit[[3]][[1]]

mse_bias_QRJMpred1_median = pred_mse_bias(gold_std1, predQR1_median)
mse_bias_QRJMpred2_median = pred_mse_bias(gold_std2, predQR2_median)
mse_bias_QRJMpred3_median = pred_mse_bias(gold_std3, predQR3_median)
rbind(mse_bias_QRJMpred1_median, mse_bias_QRJMpred2_median, mse_bias_QRJMpred3_median)
#              MSE         bias
# [1,] 0.003039118  0.011072325
# [2,] 0.004717866  0.002366571
# [3,] 0.009354163 -0.003332169


###### plotting
# qt=25 model
comp_QRJMqt25_dt025 = pred_data(gold_std1, predQR1)
comp_QRJMqt25_dt1 = pred_data(gold_std2, predQR2)
comp_QRJMqt25_dt2 = pred_data(gold_std3, predQR3)

par(mfrow=c(2,2))
plot_BA(comp_QRJMqt25_dt025, main=expression(paste(Delta, 't1')))
plot_BA(comp_QRJMqt25_dt1, main=expression(paste(Delta, 't2')))
plot_BA(comp_QRJMqt25_dt2, main=expression(paste(Delta, 't3')))


# qt=50 model
comp_QRJMqt50_dt025 = pred_data(gold_std1, predQR1_median)
comp_QRJMqt50_dt1 = pred_data(gold_std2, predQR2_median)
comp_QRJMqt50_dt2 = pred_data(gold_std3, predQR3_median)

par(mfrow=c(2,2))
plot_BA(comp_QRJMqt50_dt025, main=expression(paste(Delta, 't1')))
plot_BA(comp_QRJMqt50_dt1, main=expression(paste(Delta, 't2')))
plot_BA(comp_QRJMqt50_dt2, main=expression(paste(Delta, 't3')))

# normal model
comp_LMJM_dt025 = pred_data(gold_std1, pred1)
comp_LMJM_dt1 = pred_data(gold_std2, pred2)
comp_LMJM_dt2 = pred_data(gold_std3, pred3)

par(mfrow=c(2,2))
plot_BA(comp_LMJM_dt025, main=expression(paste(Delta, 't1')))
plot_BA(comp_LMJM_dt1, main=expression(paste(Delta, 't2')))
plot_BA(comp_LMJM_dt2, main=expression(paste(Delta, 't3')))




# #### -------------- data generated from qt=0.5 error -------------------------- ###
rm(list=ls())
setwd('/Users/askming/Dropbox/Research\ works/Self/Phd\ thesis/Simulation\ reports/pred_surv_qt50data')
load(dir()[1])
load(dir()[2])
load(dir()[3])
load(dir()[4])

gold_std1 = gold_std_qt50data[[1]]
gold_std2 = gold_std_qt50data[[2]]
gold_std3 = gold_std_qt50data[[3]]

pred1 = QRJM_pred_qt50data_meanfit[[1]][[1]]
pred2 = QRJM_pred_qt50data_meanfit[[2]][[1]]
pred3 = QRJM_pred_qt50data_meanfit[[3]][[1]]

predQR1 = QRJM_pred_qt50data_medianfit[[1]][[1]]
predQR2 = QRJM_pred_qt50data_medianfit[[2]][[1]]
predQR3 = QRJM_pred_qt50data_medianfit[[3]][[1]]


predqt25QR1 = QRJM_pred_qt50data_qt25fit[[1]][[1]]
predqt25QR2 = QRJM_pred_qt50data_qt25fit[[2]][[1]]
predqt25QR3 = QRJM_pred_qt50data_qt25fit[[3]][[1]]

# ############################## Steps to do the comparison ########################
# ## compare results from LMJM with gold standard
comp_LMJM_dt025 = pred_data(gold_std1, pred1)
comp_LMJM_dt1 = pred_data(gold_std2, pred2)
comp_LMJM_dt2 = pred_data(gold_std3, pred3)

par(mfrow=c(2,2))
plot_BA(comp_LMJM_dt025)
plot_BA(comp_LMJM_dt1)
plot_BA(comp_LMJM_dt2)

mse_bias_LMJMpred1 = pred_mse_bias(gold_std1, pred1)
mse_bias_LMJMpred2 = pred_mse_bias(gold_std2, pred2)
mse_bias_LMJMpred3 = pred_mse_bias(gold_std3, pred3)
round(rbind(mse_bias_LMJMpred1, mse_bias_LMJMpred2, mse_bias_LMJMpred3),3)
#        MSE   bias
# [1,] 0.002  0.003
# [2,] 0.008 -0.023
# [3,] 0.022 -0.046

# ## compare results from QRJM (tau=0.5) with gold standard
comp_QRJM_dt025 = pred_data(gold_std1, predQR1)
comp_QRJM_dt1 = pred_data(gold_std2, predQR2)
comp_QRJM_dt2 = pred_data(gold_std3, predQR3)

par(mfrow=c(2,2))
plot_BA(comp_QRJM_dt025)
plot_BA(comp_QRJM_dt1)
plot_BA(comp_QRJM_dt2)

mse_bias_QRJMpred1 = pred_mse_bias(gold_std1, predQR1)
mse_bias_QRJMpred2 = pred_mse_bias(gold_std2, predQR2)
mse_bias_QRJMpred3 = pred_mse_bias(gold_std3, predQR3)
round(rbind(mse_bias_QRJMpred1, mse_bias_QRJMpred2, mse_bias_QRJMpred3),3)
#        MSE   bias
# [1,] 0.003  0.010
# [2,] 0.005  0.002
# [3,] 0.009 -0.004

# ## compare results from QRJM (tau=0.25) with gold standard
comp_QRJMqt25_dt025 = pred_data(gold_std1, predqt25QR1)
comp_QRJMqt25_dt1 = pred_data(gold_std2, predqt25QR2)
comp_QRJMqt25_dt2 = pred_data(gold_std3, predqt25QR3)

par(mfrow=c(2,2))
plot_BA(comp_QRJMqt25_dt025)
plot_BA(comp_QRJMqt25_dt1)
plot_BA(comp_QRJMqt25_dt2)

#### -------------- data generated from normal error -------------------------- ###
rm(list=ls())
setwd('/Users/askming/Dropbox/Research\ works/Self/Phd\ thesis/Simulation\ reports/pred_surv_normdata')
load(dir()[1])
load(dir()[2])
load(dir()[3])

######################### results from LMJM  ###########################################
gold_std1 = gold_std_normdata[[1]]
gold_std2 = gold_std_normdata[[2]]
gold_std3 = gold_std_normdata[[3]]

pred1 = QRJM_pred_normdata_meanfit[[1]][[1]]
pred2 = QRJM_pred_normdata_meanfit[[2]][[1]]
pred3 = QRJM_pred_normdata_meanfit[[3]][[1]]

predQR1 = QRJM_pred_normdata_medianfit[[1]][[1]]
predQR2 = QRJM_pred_normdata_medianfit[[2]][[1]]
predQR3 = QRJM_pred_normdata_medianfit[[3]][[1]]

# ## compare results from LMJM with gold standard
comp_LMJM_dt025 = pred_data(gold_std1, pred1)
comp_LMJM_dt1 = pred_data(gold_std2, pred2)
comp_LMJM_dt2 = pred_data(gold_std3, pred3)

par(mfrow=c(2,2))
plot_BA(comp_LMJM_dt025)
plot_BA(comp_LMJM_dt1)
plot_BA(comp_LMJM_dt2)

mse_bias_LMJMpred1 = pred_mse_bias(gold_std1, pred1)
mse_bias_LMJMpred2 = pred_mse_bias(gold_std2, pred2)
mse_bias_LMJMpred3 = pred_mse_bias(gold_std3, pred3)
round(rbind(mse_bias_LMJMpred1, mse_bias_LMJMpred2, mse_bias_LMJMpred3),3)
#        MSE   bias
# [1,] 0.004  0.013
# [2,] 0.009  0.002
# [3,] 0.022 -0.012



######################### results from QRJM (tau=0.5)  ##################################
comp_QRJM_dt025 = pred_data(gold_std1, predQR1)
comp_QRJM_dt1 = pred_data(gold_std2, predQR2)
comp_QRJM_dt2 = pred_data(gold_std3, predQR3)

par(mfrow=c(2,2))
plot_BA(comp_QRJM_dt025)
plot_BA(comp_QRJM_dt1)
plot_BA(comp_QRJM_dt2)

mse_bias_QRJMpred1 = pred_mse_bias(gold_std1, predQR1)
mse_bias_QRJMpred2 = pred_mse_bias(gold_std2, predQR2)
mse_bias_QRJMpred3 = pred_mse_bias(gold_std3, predQR3)
round(rbind(mse_bias_QRJMpred1, mse_bias_QRJMpred2, mse_bias_QRJMpred3),3)
# ## mse and bias for QRJM
#        MSE   bias
# [1,] 0.002  0.006
# [2,] 0.004 -0.001
# [3,] 0.007 -0.005