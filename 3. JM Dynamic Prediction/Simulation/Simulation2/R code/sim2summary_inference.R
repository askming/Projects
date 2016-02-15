setwd('Z:/thetahat/')

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


qt25_qt25fit = summarize_thetahat(sim2_qt25_30_thetahat, true_par)
qt25_qt50fit = summarize_thetahat(sim2_qt25_30_medianfit_thetahat, true_par)
qt25_meanfit = summarize_thetahat(sim2_25_30_meanfit_thetahat, true_par)
qt25_inference = cbind(qt25_qt25fit, qt50_qt50fit, qt50_meanfit)

qt50_qt50fit = summarize_thetahat(sim2_25_30_50fit_thetahat, true_par)
qt50_meanfit = summarize_thetahat(sim2_50_30_meanfit_theta_hat, true_par)

normal_qt50fit = summarize_thetahat(sim2_norm_30_50fit_thetahat, true_par)
normal_meanfit = summarize_thetahat(sim2_norm_30_meanfit_thetahat, true_par)

