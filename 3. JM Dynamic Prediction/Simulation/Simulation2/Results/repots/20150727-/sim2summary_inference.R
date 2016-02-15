setwd('/Users/askming/Dropbox/Research\ works/Self/Phd\ thesis/Simulation\ reports/thetahat_0721')

setwd('/Users/askming/Dropbox/Research\ works/Self/Phd\ thesis/Simulation\ reports/normdata_thetahat')



true_par = c(1,1,1,1,1,-1,1,-1,1,1)

summarize_thetahat = function(data, true_par){
	all_post_samples = NULL
	for(i in seq_along(data)){
		all_post_samples = rbind(all_post_samples, data[[i]])
	}
	post_means = colMeans(all_post_samples)[c(1:7, 9:11)]
	bias = post_means - true_par
	se = apply(all_post_samples[, c(1:7, 9:11)], 2, sd)
	MSE = sapply(seq(1:dim(all_post_samples[,c(1:7, 9:11)])[2]), function(x) mean((all_post_samples[,c(1:7, 9:11)][,x] - true_par[x])^2))

	output = cbind(bias, se, MSE)
	colnames(output) = c('bias', 'se', 'MSE')
	return(round(output, 3))
}

load(dir()[1])
load(dir()[2])
load(dir()[3])

qt25_qt25fit = summarize_thetahat(sim_data_tau25_100_qt25fit_thetahat, true_par)
qt25_qt50fit = summarize_thetahat(sim_data_tau25_100_medianfit_thetahat, true_par)
qt25_meanfit = summarize_thetahat(sim_data_tau25_100_normfit_thetahat, true_par)
qt25data_inference = cbind(qt25_qt25fit, qt25_qt50fit, qt25_meanfit)
xtable(qt25data_inference)


load(dir()[4])
load(dir()[5])
load(dir()[6])

# qt50_qt25fit = summarize_thetahat(sim_data_tau50_100_qt25fit_thetahat, true_par)
qt50_qt50fit = summarize_thetahat(sim_data_tau50_100_medianfit_thetahat, true_par)
qt50_meanfit = summarize_thetahat(sim_data_tau50_100_normfit_thetahat, true_par)
# qt50data_inference = cbind(qt50_qt25fit, qt50_qt50fit, qt50_meanfit)
qt50data_inference = cbind(qt50_qt50fit, qt50_meanfit)
xtable(qt50data_inference)



normal_qt50fit = summarize_thetahat(sim_data_norm_100_medianfit_thetahat, true_par)
normal_meanfit = summarize_thetahat(sim_data_norm_100_normfit_thetahat, true_par)
normdata_inference = cbind(normal_qt50fit, normal_meanfit)
xtable(normdata_inference)

