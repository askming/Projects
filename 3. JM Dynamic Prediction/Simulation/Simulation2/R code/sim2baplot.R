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
plot_BA = function(data){
    library(MethComp)
    data_to_plot = Meth(data)
    par(mar=c(6,4,4,5))
    BA.plot(data_to_plot,cex.point=0.5, diflim=c(-0.7,0.7), xlab="mean", ylab="difference")
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


############### Simulation results on 05/21/2015 ##############################
#### -------------- data generated from qt=0.25 error ------------------------- ###
setwd('/Users/askming/Documents/github/Projects/3.\ JM\ Dynamic\ Prediction/Simulation/Simulation2/Results/0520-0521_preds')
load(dir()[3])
load(dir()[1])
load(dir()[2])

gold_std1 = gold_std[[1]]
gold_std2 = gold_std[[2]]
gold_std3 = gold_std[[3]]

pred1 = LMJM_pred[[1]][[1]]
pred2 = LMJM_pred[[2]][[1]]
pred3 = LMJM_pred[[3]][[1]]

predQR1 = QRJM_pred[[1]][[1]]
predQR2 = QRJM_pred[[2]][[1]]
predQR3 = QRJM_pred[[3]][[1]]

mse_bias_LMJMpred1 = pred_mse_bias(gold_std1, pred1)
mse_bias_LMJMpred2 = pred_mse_bias(gold_std2, pred2)
mse_bias_LMJMpred3 = pred_mse_bias(gold_std3, pred3)
rbind(mse_bias_LMJMpred1, mse_bias_LMJMpred2, mse_bias_LMJMpred3)
#              MSE         bias
# [1,] 0.003351773  0.009696817
# [2,] 0.004566668 -0.004989684
# [3,] 0.009382506 -0.012867794

mse_bias_QRJMpred1 = pred_mse_bias(gold_std1, predQR1)
mse_bias_QRJMpred2 = pred_mse_bias(gold_std2, predQR2)
mse_bias_QRJMpred3 = pred_mse_bias(gold_std3, predQR3)
rbind(mse_bias_QRJMpred1, mse_bias_QRJMpred2, mse_bias_QRJMpred3)
#              MSE         bias
# [1,] 0.002514678  0.008509040
# [2,] 0.003685394 -0.003772431
# [3,] 0.008289105 -0.008087641

############### Simulation results on 06/18/2015 #################################
setwd('/Users/askming/Documents/github/Projects/3.\ JM\ Dynamic\ Prediction/Simulation/Simulation2/Results/0618_preds')

#### -------------- data generated from normal error -------------------------- ###
load(dir()[6])
load(dir()[1])
load(dir()[3])


gold_std1 = gold_std_normdata[[1]]
gold_std2 = gold_std_normdata[[2]]
gold_std3 = gold_std_normdata[[3]]

pred1 = LMJM_pred[[1]][[1]]
pred2 = LMJM_pred[[2]][[1]]
pred3 = LMJM_pred[[3]][[1]]

predQR1 = QRJM_pred_normdata_medianfit[[1]][[1]]
predQR2 = QRJM_pred_normdata_medianfit[[2]][[1]]
predQR3 = QRJM_pred_normdata_medianfit[[3]][[1]]

## mse and bias for LMJM
#              MSE         bias
# [1,] 0.002520234  0.004967799
# [2,] 0.005541016 -0.007837558
# [3,] 0.009819399 -0.013371873

## mse and bias for QRJM
#              MSE         bias
# [1,] 0.001948289  0.005029296
# [2,] 0.005527062 -0.001633113
# [3,] 0.009534088 -0.002727554

#### -------------- data generated from qt=0.5 error -------------------------- ###
load(dir()[8])
load(dir()[2])
load(dir()[4])

gold_std1 = gold_std_qt50data[[1]]
gold_std2 = gold_std_qt50data[[2]]
gold_std3 = gold_std_qt50data[[3]]

pred1 = LMJM_pred_qt50data[[1]][[1]]
pred2 = LMJM_pred_qt50data[[2]][[1]]
pred3 = LMJM_pred_qt50data[[3]][[1]]

predQR1 = QRJM_pred_qt50_medianfit[[1]][[1]]
predQR2 = QRJM_pred_qt50_medianfit[[2]][[1]]
predQR3 = QRJM_pred_qt50_medianfit[[3]][[1]]


############################## Steps to do the comparison ########################
## compare results from LMJM with gold standard
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
rbind(mse_bias_LMJMpred1, mse_bias_LMJMpred2, mse_bias_LMJMpred3)
#              MSE          bias
# [1,] 0.001591336 -0.0003156817
# [2,] 0.005952234 -0.0183922787
# [3,] 0.023398132 -0.0441487874

## compare results from QRJM with gold standard
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
rbind(mse_bias_QRJMpred1, mse_bias_QRJMpred2, mse_bias_QRJMpred3)
#              MSE        bias
# [1,] 0.001838707 0.005545511
# [2,] 0.004500276 0.004882855
# [3,] 0.008519653 0.001912865