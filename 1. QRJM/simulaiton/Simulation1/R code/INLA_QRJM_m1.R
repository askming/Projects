### Joint model using INLA ###
### 2015-01-10 ###
### MODEL 1
rm(list=ls())


sim_Ti = function(n=250, alpha1=1, alpha2=1, gamma=c(1,1), delta=c(1,1)){
	u = rnorm(n)
	W = matrix(rnorm(2*n),n,2)

	rnd =runif(n) 
	Ti = -log(rnd)/exp(W%*%gamma + alpha*u)
	Ci = 5*rbeta(n,4,1)
	Delta = as.numeric(Ti < Ci)
	Ti2 = pmin(Ti, Ci)

	data.frame(Ti=Ti2, event=Delta, W=W, U=u)
}

# simulate quantile longitudinal data in long format
sim_longitudinal_data = function(survival_data=surdata, n=250, time=c(0, 0.25, 0.5, 0.75, 1, 3), tau, sigma=1, beta=c(1,1)){
	n_obs = length(time)
	time_long = rep(time, n) # at most # = length(time) observations per patient
 
	Y = rep(NA, n*n_obs) # long format
	U = rep(survival_data$U, each=n_obs) # random effects
	X = cbind(1, rep(rnorm(n), each=n_obs))
	Ti = survival_data$Ti
	count = sapply(Ti, function(x) sum(x > time)) # number of observations after drop-outs
	end = seq(6, n*n_obs, 6)

	# random variable generator: ALD(0,sigma,p)
	rald = function(n, location=0, scale, p){
		u = rnorm(n)
		z = rexp(n)
		v = scale*z
		theta = (1-2*p)/(p*(1-p))
		tau = sqrt(2/(p*(1-p)))
		sample = theta * z + tau * sqrt(z) * u
	}
	# create outcome
	Y = X %*% beta + U + rnorm(n*n_obs)	

	# create drop out missing
	for(i in seq_along(end)){
		if(count[i]!=6){
					Y[min(count[i]+6*(i-1)+1,end[i]):end[i]] = NA
		}

	}

	obs_id = rep(seq(1:6), n)
	data.frame(obs_id, Y, X, time_long)
}

sim_multi_data = function(n){
	out = vector(n, mode="list")
	for (i in 1:n){
		surdata = sim_Ti()
		long_data = sim_longitudinal_data(surdata, tau=0.25)
		out[[i]] = list(surdata, long_data)
	}
	out
}

surdata = sim_Ti()
long_data = sim_longitudinal_data(surdata)



str(long_data)

library(INLA)

run_multi_data = function(data){
	long_data = data[[2]]
	surdata = data[[1]]


	N = 250
	ng =dim(long_data)[1]
	ns = N

	## prepare the response variable
	y_long = c(long_data$Y, rep(NA, ns))
	y_surv = inla.surv(time=c(rep(NA, ng), surdata$Ti), event=c(rep(NA, ng), surdata$event))
	Yjoint = list(y_long, y_surv)

	# prepare the fixed covariate
	linear.covariate = data.frame(
                               b1.X1 = c(long_data$X1, rep(0, ns)),
                               b2.X2 = c(long_data$X2, rep(0, ns)),
                               time = c(long_data$time_long, rep(0, ns)),
                               g1.W1 = c(rep(0, ng), surdata$W.1),
                               g2.W2 = c(rep(0, ng), surdata$W.2),
                               Ti = c(rep(0, ng), surdata$Ti)
                               )


	# prepare random effects
	random.covariate <- list(U11 = c(rep(1:N, each=6),rep(NA, ns)),
                         	 U12 = c(rep(NA,ng), 1:N)
                         )

	joint.data <- c(linear.covariate,random.covariate)
	joint.data$Y <- Yjoint

	formula = Y ~ -1 + b1.X1 + b2.X2 + g1.W1 + g2.W2  +
	f(U11, model="iid", n=N) +
	f(U12, copy="U11", fixed=FALSE)


	# enable.model.likelihood.laplace=TRUE
	mod = inla(formula, family=c("gaussian", "exponential"), data = joint.data, verbose=TRUE)

	list(summary(mod)$fixed, summary(mod)$hyperpar)

}

multidata = sim_multi_data(30)
res = lapply(multidata, run_multi_data)


sry = function(res){

		fix_effects_mean = sapply(res, function(x) x[[1]][,1])
		fix_effects_sd = sapply(res, function(x) x[[1]][,2])
		re = sapply(res, function(x) x[[2]][3,1:2])
		mean_fix_mean = round(rowMeans(fix_effects_mean), 3)
		mean_fix_sd = round(rowMeans(fix_effects_sd), 3)
		re_mean = round(rowMeans(re), 3)

		list(mean_fix_mean, mean_fix_sd, alpha = re_mean)


}

sry(res)
