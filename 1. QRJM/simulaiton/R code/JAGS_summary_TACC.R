# calculate the posterior means and sd of the estimates
post_summary = function(fit){
	len = length(fit)
	means = sapply(fit, function(x) x[,1])
	m = round(rowMeans(means), 3)
	sds = sapply(fit, function(x) x[,2])
	sd = round(rowMeans(sds),3)

	list(mean=m, sd=sd)
}

post_summary(multi_fit11_100_500_75)
post_summary(multi_fit01)