### Joint model using INLA ###
### 2015-01-10 ###
### MODEL 2 & 3
## longitudianl outcome can have missing values

# simulate survival data
sim_Ti = function(n=250, alpha1=1, alpha2=1, delta=c(1,1), gamma=c(1,1)){
	library(mvtnorm)
	# random effects
	Sigma=matrix(c(0.09, 0.09*0.16, 0.09*0.16, 0.09), nrow=2, byrow=T)
	u = rmvnorm(n,c(0,0),Sigma)
	H = matrix(rnorm(2*n),n,2)
	W = matrix(rnorm(2*n),n,2)

	# define the survival function
	rnd = runif(n)
	Ti = -log(rnd)/(exp(alpha1*u[,1] +  alpha2*u[,2] + W%*%gamma))
	Ci = 5*rbeta(n,4,1)
	Delta = as.numeric(Ti < Ci)
	Ti2 = pmin(Ti, Ci)

	data.frame(Ti=Ti2, event=Delta, W=W, U=u)
}

# simulate quantile longitudinal data in long format
sim_longitudinal_data = function(survival_data=surdata, n=250, time=c(0, 0.25, 0.5, 0.75, 1, 3), tau, sigma=1, beta=c(1,1), delta=c(1,1)){
	# survival_data - data simulated from survival model
	# n - # of subjects
	# time - time points of observations
	# tau - quantile
	# sigma - scale parameter
	n_obs = length(time)
	time_long = rep(time, n) # at most # = length(time) observations per patient
	U = matrix(NA, n*n_obs, 2)
	Y = rep(NA, n*n_obs) # long format
	U[,1] = rep(survival_data$U.1, each=n_obs) # random effects
	U[,2] = rep(survival_data$U.2, each=n_obs)
	# H[,1] = rep(survival_data$H.1, each=n_obs)
	# H[,2] = rep(survival_data$H.2, each=n_obs)*time_long
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
	for (i in 1:length(Y)) {
		Y[i] = X[i,] %*% beta + U[i,] %*% c(1, time_long[i]) + rald(1,scale=sigma, p=tau)	

	}

	# create drop out missing
	for(i in seq_along(end)){
		if(count[i]!=6){
					Y[min(count[i]+6*(i-1)+1,end[i]):end[i]] = NA
		}

	}

	obs_id = rep(seq(1:6), n)
	data.frame(obs_id, Y, X, time_long)
}


###############################################################
#----- function to simulate multiple joint data sets ---------#
###############################################################
sim_multiple_data = function(N, n, sur_fun=sim_Ti, longi_fun=sim_longitudinal_data, alpha, tau){
	# N - number of data sets to generate
	# sur_fun - function to simulate survival data
	# longi_fun - function to sumualte longitudinal data
	# alpha - association mechanism for JM
	# tau - quantile

	outdata = vector(mode='list', N)
	for (i in 1:N){
		sur_data = sur_fun(n=n, alpha1=alpha[1], alpha2=alpha[2])
		longi_data = longi_fun(n=n, sur_data, tau=tau)
		outdata[[i]] = list(survival_data=sur_data, longitudinal_data=longi_data)		
	}
	outdata	
}

# surdata = sim_Ti(alpha1=1, alpha2=1)
# long_data = sim_longitudinal_data(surdata, tau=0.25)

# head(long_data, 24)

run_multi_inla = function(data, tau){
	long_data = data$longitudinal_data
	surdata = data$survival_data
	enable.model.likelihood.laplace=TRUE
	## prepare data set
	N = 250
	ng =dim(long_data)[1]
	ns = N

	## prepare the response variable
	y_long = c(long_data$Y, rep(NA, ns))
	y_surv = inla.surv(time=c(rep(NA, ng), surdata$Ti), event=c(rep(NA, ng), surdata$event))
	Yjoint = list(y_long, y_surv)

	# prepare the fixed covariate
	linear.covariate = data.frame(mu = as.factor(c(rep(1, ng), rep(2, ns))),
                               b1.X1 = c(long_data$X1, rep(0, ns)),
                               b2.X2 = c(long_data$X2, rep(0, ns)),
                               time = c(long_data$time_long, rep(0, ns)),
                               g1.W1 = c(rep(0, ng), surdata$W.1),
                               g2.W2 = c(rep(0, ng), surdata$W.2),
                               Ti = c(rep(0, ng), surdata$Ti)
                               )


	# prepare random effects
	random.covariate <- list(U11 = c(rep(1:N, each=6),rep(NA, ns)),
                         U21 = c(rep(N+(1:N), each=6),rep(NA, ns)),
                         U12 = c(rep(NA,ng), 1:N),
                         U22 = c(rep(NA,ng), N+(1:N))
                         )

	joint.data <- c(linear.covariate,random.covariate)
	joint.data$Y <- Yjoint


	
	### MODEL
	formula = Y ~ b1.X1 + b2.X2 + g1.W1 + g2.W2 - 1 +
	f(U11, model="iid2d", n=2*N) +
	f(U21, time, copy="U11") +
	f(U12, copy="U11", fixed=FALSE) +
	f(U22, copy="U11", fixed=FALSE) 


	# enable.model.likelihood.laplace=TRUE
	mod = inla(formula, family=c("laplace", "exponential"), data = joint.data, verbose=TRUE, control.family=list(list(alpha=0.25), list()))

	list(summary(mod)$fixed, summary(mod)$hyperpar)
}


multi_data = sim_multiple_data(1, n=250, alpha=c(1,1), tau=0.25)

t1=Sys.time()
INLA_m3_100_25 = foreach(i=1:100,.packages=c("INLA","foreach"),.verbose=T)%dopar%{ 
	enable.model.likelihood.laplace=TRUE
	set.seed(123)
  	temp = run_multi_inla(multi_data[[i]], tau=0.25)
  	temp
}  
t2=Sys.time()
t2-t1


summ = function(result){

		fix_effects_mean = sapply(result, function(x) x[[1]][,1])
		fix_effects_sd = sapply(result, function(x) x[[1]][,2])
		re = sapply(result, function(x) x[[2]][5:6,1:2])
		mean_fix_mean = round(rowMeans(fix_effects_mean), 3)
		mean_fix_sd = round(rowMeans(fix_effects_sd), 3)
		re_mean = round(rowMeans(re), 3)

		list(mean_fix_mean, mean_fix_sd, re_mean)


}

summ(result)























