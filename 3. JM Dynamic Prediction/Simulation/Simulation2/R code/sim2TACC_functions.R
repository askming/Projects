rm(list=ls())

#### ------------------------------------------------ ####
################ Functions to generate data ##############
#### ------------------------------------------------ ####


#----------- function to simulate survival time --------------#
# survival function is given by: S(t)= exp(- exp(B) * (exp(A*t) - 1) ) / A)
sim_Ti = function(n=500, alpha1, alpha2, delta, gamma){
	library(mvtnorm)
	# random effects
	Sigma=matrix(c(0.09, 0.09*0.16, 0.09*0.16, 0.09), nrow=2, byrow=T)
	u=rmvnorm(n,c(0,0),Sigma)
	# fixed effects
	H = matrix(rnorm(2*n, 1, 1),n,2)
	W = matrix(rnorm(2*n, 2, 4),n,2)
    # original
    # H = matrix(rnorm(2*n),n,2)
    # W = matrix(rnorm(2*n),n,2)

	# define the survival function
	surv = function(t) {
		if(alpha1==0 & alpha2==0) {res=t*exp(W%*%gamma)}
		else{
			res = ((exp(alpha2*u[,1]+alpha2*t*u[,2]+alpha1*delta[1]*H[,1]+alpha1*delta[2]*H[,2]*t+W%*%gamma)-exp(alpha1*H[,1]*delta[1]+W%*%gamma+alpha2*u[,1]))/(alpha2*u[,2]+alpha1*H[,2]*delta[2]))
		}
		return(exp(-res))
	}
	# simulate event time
	rnd =runif(n)
	Ti = rep(Inf,n)
	w = which(surv(500)-rnd < 0)
	# print(w)
	# print(cbind(w, surv(500)[w], rnd[w]))
	# w = which(1-surv(101)-rnd > 0) # original
	for (j in w){
		# print(j)
		Ti[j] = uniroot(function(x){surv(x)[j]-rnd[j]}, lower=0, upper=1e8)$root
	}

	# Ti[w] = sapply(w,function(j)uniroot(function(x){surv(x)[j]-rnd[j]},lower=0,upper=500)$root)
	Ti[Ti==0] = min(Ti[Ti>0])/2
	Ci = rexp(n, 0.01)
	Delta = as.numeric(Ti < Ci)
	Ti2 = pmin(Ti, Ci)

	list(Ti=matrix(Ti2, ncol=1), event=matrix(Delta, ncol=1), H=H, W=W, U=u)
}

# test_Ti = sim_Ti(alpha1=1, alpha2=1, delta=c(-1,1), gamma=c(-1,1))

# test_Ti[[1]]
# mean(test_Ti[[2]])
# sum(test_Ti[[1]]>=2)
# sum(test_Ti[[1]]>=4)
# sum(test_Ti[[1]]>=8)
# sum(test_Ti[[1]]>=16)
# summary(test_Ti[[1]])


#----- function to simulate longitudinal data --------------- #
sim_longitudinal_data = function(survival_data=surdata, n=500, time=c(0, 0.25, 0.5, 1, 1.5, 3), tau, sigma=1, beta, delta, normal_data=FALSE){
	# survival_data - data simulated from survival model
	# n - # of subjects
	# time - time points of observations
	# tau - quantile
	# sigma - scale parameter

	# random generation of ALD(0,sigma,p)
	rald = function(n, location=0, scale, p){
		u = rnorm(n)
		z = rexp(n)
		v = scale*z
		theta = (1-2*p)/(p*(1-p))
		tau = sqrt(2/(p*(1-p)))
		samples = theta * z + tau * sqrt(z) * u
		samples
	}

	y = matrix(NA, nrow=n, ncol=length(time)) # wide format
	Ti = survival_data$Ti
	U = survival_data$U # random effects
	H = survival_data$H
	X = cbind(1, rnorm(n))
	count = sapply(Ti, function(x) sum(x > time)) # number of observations after drop-outs

	for (i in 1:n){
		for (j in 1:count[i]){
			location = beta %*% X[i, ] + delta %*% c(H[i,1], H[i,2]*time[j]) + U[i,] %*% c(1, time[j])
            if(normal_data==FALSE) y[i,j] = location + rald(1, scale=sigma, p=tau)
            else y[i,j] = location + rnorm(1)
		}
	}
	list(y = y, X = X, J=matrix(count, ncol=1))
}


#----- function to simulate multiple joint data sets ---------#
prep_sim_data = function(N, sur_fun=sim_Ti, longi_fun=sim_longitudinal_data, alpha, tau=NULL, beta, delta, gamma, sam_size=520, sam_selected=500, normal_data=FALSE){
	# N - number of data sets to generate
	# sur_fun - function to simulate survival data
	# longi_fun - function to simulate longitudinal data
	# alpha - association mechanism for JM
	# tau - quantile
	# sam_size - sample size for each data set
	# selected_id - id for subjects that are selected in fitting model
	outdata = vector(mode='list', N)
	# Each element of outdata is a list of four elements: random sample id, full data, selected data and validation data
	for (i in 1:N){
		print(i)
		selected_id = sort(sample(c(1:sam_size), sam_selected))
		sur_data = sur_fun(n=sam_size, alpha1=alpha[1], alpha2=alpha[2], delta=delta, gamma=gamma)
		longi_data = longi_fun(n=sam_size, sur_data, tau=tau, beta=beta, delta=delta, normal_data=normal_data)
		sur_data_selected = lapply(sur_data, function(x, id){x[id,]}, selected_id)
		longi_data_selected = lapply(longi_data, function(x, id){x[id,]}, selected_id)

		full_data = list(survival_data=sur_data, longitudinal_data=longi_data)
		selected_sub = list(survival_data=sur_data_selected, longitudinal_data=longi_data_selected)
		validation_sub = list(survival_data=lapply(sur_data, function(x, id){x[id,]}, -selected_id), longitudinal_data=lapply(longi_data, function(x, id){x[id,]}, -selected_id))
		outdata[[i]] = list(in_id=selected_id, full=full_data, selected=selected_sub, valid=validation_sub)
	}
	outdata
}



## simulate multiple data sets
# sim_data = prep_sim_data(30, alpha=c(1,1),tau=0.5, beta=c(1, 1), delta=c(-1,1), gamma=c(-1,1))


#----- function to prepare every single validation data set for specific time window ---------#
select_valid_data = function(indata, time){
	t = length(time)
	y = indata[[2]][['y']]# upto time t
	count = apply(y, 1, function(x) sum(!is.na(x))) # count how many obs each subject has
	selected_id = which(count>=t)
	y_sel = y[selected_id, 1:t]
	X = indata[[2]][['X']][selected_id,]
	J = rep(t, dim(y_sel)[1]) # every one has same number of observations

	# don't need event variable anymore
	Ti = rep(time[end(time)[1]], dim(y_sel)[1])
	H = indata[[1]][['H']][selected_id,]
	W = indata[[1]][['W']][selected_id,]
	U = indata[[1]][['U']][selected_id,]

	# output - a list, same format as original input data
	list('survival_data'=list(Ti=Ti, H=H, W=W, U=U), 'longitudina_data'=list(y=y_sel, X=X, J=J, time=time, id=selected_id))
}



#### ------------------------------------------------ ####
#### Functions to fit the data using different models ####
#### ------------------------------------------------ ####

## define function to fit the joint model using quantile regression ##
QRJM_jags = function(data,tau,I=500){
	library(R2jags)
	setwd("/work/02784/myang3/QRJM/model")
	model = "/work/02784/myang3/QRJM/model/QRJM_Cholesky.txt"
	# file.show("/work/02784/myang3/QRJM/model/QRJM_Cholesky.txt")
	# load longitudinal data
	y = data[[2]][['y']]
	X = data[[2]][['X']]
	J = data[[2]][['J']]

	# load survival data
	Ti = data[[1]][['Ti']]
	H = data[[1]][['H']]
	W = data[[1]][['W']]
	event =  data[[1]][['event']]

	zeros = rep(0, I)
	zero = c(0,0)

	jags.data = list(y=y, X=X, H=H, W=W, Ti=Ti, event=event, qt=tau, t=c(0, 0.25, 0.5, 1, 1.5, 3), I=I, J=J, zeros=zeros, zero=zero)
 	jags.params = c("beta", "delta","gamma","alpha1","alpha2","sigma", "c", "w11","w21","w22")
 	jags.inits = function(){	list(beta=c(0.1, 0.2), delta=c(0.1, 0.2), gamma=c(0.1, 0.2), sigma=0.1, alpha1=0.1, alpha2=0.1, c=0.1)
  	}
  jags(data=jags.data, inits=jags.inits, jags.params, n.chain=3, n.iter=15000, n.burnin=12000, n.thin=10, model.file=model)
}


## define function to fit the joint model using mean regression ##
NormJM_jags = function(data){
	library(R2jags)
	setwd("/work/02784/myang3/QRJM/model")
	model = "/work/02784/myang3/QRJM/model/NormalJM_Cholesky.txt" #JM using LMM
	# file.show(model)
	# load longitudinal data
	y = data[[2]][['y']]
	X = data[[2]][['X']]
	J = data[[2]][['J']]
	I = dim(y)[1]
	# load survival data
	Ti = data[[1]][['Ti']]
	H = data[[1]][['H']]
	W = data[[1]][['W']]
	event =  data[[1]][['event']]

	# alpha2 = 1
	zeros = rep(0, I)
	zero = c(0,0)

	jags.data = list(y=y, X=X, H=H, W=W, Ti=Ti, event=event, t=c(0, 0.25, 0.5, 1, 1.5, 3), I=I, J=J, zeros=zeros, zero=zero)
 	jags.params = c("beta", "delta","gamma","alpha1","alpha2","sigma", "c", "w11","w21","w22")
 	jags.inits = function(){	list(beta=c(0.1,0.1), delta=c(0.1, 0.1), gamma=c(0.1, 0.1), sigma=0.1, alpha1=0.1, alpha2=0.1, c=0.1)
  	}
  jags(data=jags.data, inits=jags.inits, jags.params, n.chain=3, n.iter=15000, n.burnin=12000, n.thin=10, model.file=model)
}



#### ------------------------------------------------ ####
###### following are functions for making prediction #####
#### ------------------------------------------------ ####




##################################################################################
## gold standard part ##
#################################################################################t
# function to calculated the gold standard of predicted survival probability given t1 at time t1+dt for the simulated data
gold_std_pred_surv = function(indata, t1, dt, true_par=list(alpha=c(1,1), delta=c(-1,1), gamma=c(-1,1))){

    # indata - simulated multiple data sets, a list
    # par - true parameter list

    # number of qualified subjects in each data set
    no_of_subj = function(data, dataid){
        length(data[[dataid]][['longitudina_data']]$id)
    }
    no_data = length(indata) # find out how many simulated data sets
    no_sub = sapply(seq(1,no_data), no_of_subj, data=indata) # number of qualified subjects in each data set


    # function to calcuate gold standard survival probability at time t for a single simulated data set
	gold_std_surv = function(data, id, t, par, dataid){
		# data - a list of multiple simulated data sets
		# id - subject id
		# t - time when the survival probability is calcualted for
		# par - true parameters used in simulation
		# dataid - id of the data set in the data list
		onedata = data[[dataid]]$survival_data
		u = matrix(onedata[['U']][id,], ncol=2)
		H = matrix(onedata[['H']][id,], ncol=2)
		W = matrix(onedata[['W']][id,], ncol=2)

		alpha1 = par$alpha[1]
		alpha2 = par$alpha[2]
		gamma = par$gamma
		delta = par$delta

		if(alpha1==0 & alpha2==0) {res=t*exp(W%*%gamma)}
		else{
			res = ((exp(alpha2*u[,1]+alpha2*t*u[,2]+alpha1*delta[1]*H[,1]+alpha1*delta[2]*H[,2]*t+W%*%gamma)-exp(alpha1*H[,1]*delta[1]+W%*%gamma+alpha2*u[,1]))/(alpha2*u[,2]+alpha1*H[,2]*delta[2]))
		}
		return(as.numeric(exp(-res)))
	}

    # calculate the predicted suvival probability at t+dt given survived up to time t for all the qualified subject in indata
    pred_surv =  vector(length=no_data, mode='list')
    for (i in seq_along(pred_surv)){
        pred_surv[[i]] = sapply(seq(1,no_sub[i]), gold_std_surv, data=indata, t=t1+dt, dataid=i, par=true_par)/sapply(seq(1,no_sub[i]), gold_std_surv, data=indata, t=t1, dataid=i, par=true_par)
    }

    # return final result, should be of length = length(indata)
    return(unlist(pred_surv))
}


##################################################################################
## prediction part ##
##################################################################################

# ------- 1. Below are for QRJM ----------#

### function to make predictions of random effects based on posterior samples for single subject in certain data set###
QRJM_jags_pred = function(data, tau, par, id){
	# data here should already be "truncated" by some time point for current available historical information
	# data - a list of two lists, the fist list contains survival data and the second list contains longitudinal data
	# par - a list of posterior samples, each element of the list is in matrix format
	# id - subject id in the data set
	library(R2jags)
	setwd("/work/02784/myang3/QRJM_sim2/model")
	model = "/work/02784/myang3/QRJM_sim2/model/QRJM_Cholesky_pred.txt"

	y = as.vector(data[[2]][['y']][id, ])
	X = as.vector(data[[2]][['X']][id,])
	J = data[[2]][['J']][id]
	time = data[[2]][["time"]]

	# load survival data
	Ti = data[[1]][['Ti']][id]
	H = as.vector(data[[1]][['H']][id,])
	W = as.vector(data[[1]][['W']][id,])

	# other data info
	# I = dim(y)[1]
	# I = length(y)
	zeros = rep(0, length(id))
	zero = c(0,0)

	# get posterior samples of the parameters, should loop over variable par
	alpha1 = par$alpha1
	alpha2 = par$alpha2
	beta = c(par$'beta[1]', par$'beta[2]')
	gamma = c(par$'gamma[1]', par$'gamma[2]')
	sigma = par$sigma
	delta = c(par$'delta[1]', par$'delta[2]')
	c = par$c
	w11 = par$w11
	w21 = par$w21
	w22 = par$w22

	jags.data = list(y=y, X=X, H=H, W=W, Ti=Ti, qt=tau, t=time, J=J, zeros=zeros, zero=zero, beta=beta, delta=delta, gamma=gamma, sigma=sigma, alpha1=alpha1, alpha2=alpha2, c=c, w11=w11, w21=w21, w22=w22)
 	jags.params = c("u")
 	# jags.inits = function(){	list() }
  	jags(data=jags.data, parameters.to.save=jags.params, n.chain=1, n.iter=1000, n.burnin=500, model.file=model)
}

## function to implement the prediction of random effects for multiple simulated data sets and every subject in the data set using QRJM
QRJM_pred_re = function(indata, post_par, tau){
    # indata - validation data sets with qualified subjects only, a list, each element in the list is a list as well that contains two lists of longitudinal and survival data
    # post_par - posterior samples of the parameters, a list
    N = length(indata) # number of multiple simulated data
    K = dim(post_par[[1]])[1] # number of posterior samples
    predicted_RE = vector(length=N, mode="list")

    ### functions to collect predicted RE ###
    collect_re = function(pred_data){
        re = unlist(lapply(pred_data, colMeans))
        u1 = re[names(re)=='u[1]']
        u2 = re[names(re)=='u[2]']
        cbind(u1, u2) # output a matrix of K times 2 MC samples of REs
    }

    for (i in 1:N){
        for (j in 1:length(indata[[i]][['longitudina_data']]$id)){
            current_re = foreach(k=1:K,.packages=c("rjags","R2jags","foreach"),.verbose=T, .export=c("QRJM_jags_pred"))%dopar%{ # iterate over posterior sampels
                    set.seed(123)
                    temp <- QRJM_jags_pred(indata[[i]], tau=tau, par=post_par[[i]][k,], id=j)
                    temp$BUGSoutput$sims.matrix# collect predictions random effects
            }
        current_re = collect_re(current_re) # a matrix of K times 2
        predicted_RE[[i]][[j]]= current_re
        }
    }
    return(predicted_RE)
}
# ------- 1. Above are for QRJM ----------#


# ------- 2. Below are for LMM JM ----------#
NormJM_jags_pred = function(data, par, id){
	library(R2jags)
	setwd("/work/02784/myang3/QRJM_sim2/model")
	model = "/work/02784/myang3/QRJM_sim2/model/NormalJM_Cholesky_pred.txt" #JM using LMM

	y = as.vector(data[[2]][['y']][id,])
    X = as.vector(data[[2]][['X']][id,])
	J = data[[2]][['J']][id]
	time = data[[2]][["time"]]
	# load survival data
	Ti = data[[1]][['Ti']][id]
    H = as.vector(data[[1]][['H']][id,])
    W = as.vector(data[[1]][['W']][id,])

	# alpha2 = 1
	zeros = rep(0, length(id))
	zero = c(0,0)

    # get posterior samples of the parameters, should loop over variable par
    alpha1 = par$alpha1
    alpha2 = par$alpha2
    beta = c(par$'beta[1]', par$'beta[2]')
    gamma = c(par$'gamma[1]', par$'gamma[2]')
    sigma = par$sigma
    delta = c(par$'delta[1]', par$'delta[2]')
    c = par$c
    w11 = par$w11
    w21 = par$w21
    w22 = par$w22

    jags.data = list(y=y, X=X, H=H, W=W, Ti=Ti, t=time, J=J, zeros=zeros, zero=zero, beta=beta, delta=delta, gamma=gamma, sigma=sigma, alpha1=alpha1, alpha2=alpha2, c=c, w11=w11, w21=w21, w22=w22)
    jags.params = c("u")
    jags(data=jags.data, parameters.to.save=jags.params, n.chain=1, n.iter=1000, n.burnin=500, model.file=model)
}
# testing
# pred_re_test_lm = NormJM_jags_pred(curt_data[[1]], sim2_25_30_75fit_thetahat[[1]][1,], id=1)


## function to implement the prediction of random effects for multiple simulated data sets and every subject in the data set using LMJM
LMJM_pred_re = function(indata, post_par){
    # indata - validation data sets with qualified subjects only, a list
    # post_par - posterior samples of the parameters, a list
    N = length(indata) # number of multiple simulated data
    K = dim(post_par[[1]])[1] # number of posterior samples
    predicted_RE = vector(length=N, mode="list")

    ### functions to collect predicted RE ###
    collect_re = function(pred_data){
        re = unlist(lapply(pred_data, colMeans))
        u1 = re[names(re)=='u[1]']
        u2 = re[names(re)=='u[2]']
        cbind(u1, u2) # output a matrix of K times 2 MC samples of REs
    }

    for (i in 1:N){
        for (j in 1:length(indata[[i]][['longitudina_data']]$id)){
            current_re = foreach(k=1:K,.packages=c("rjags","R2jags","foreach"),.verbose=T, .export=c("NormJM_jags_pred"))%dopar%{ # iterate over posterior sampels
                    set.seed(123)
                    temp <- NormJM_jags_pred(indata[[i]], par=post_par[[i]][k,], id=j)
                    temp$BUGSoutput$sims.matrix# collect predictions random effects
            }
        current_re = collect_re(current_re) # a matrix of K times 2
        predicted_RE[[i]][[j]]= current_re
        }
    }
    return(predicted_RE)
}
# ------- 2. Above are for LMJM ----------#


## function to calculate the predicted survival probability based on model fitting and predicted REs for certain time window
model_pred_surv = function(indata, t1, dt, post_par, pred_re){
    # indata - simulated multiple data sets, a list
    # post_par - posterior samples of parameters, a list
    # pred_re - predicted random effects, a list of list

    # number of qualified subjects in each data set
    no_of_subj = function(data, dataid){
        length(data[[dataid]][['longitudina_data']]$id)
    }
    no_data = length(indata) # find out how many simulated data sets
    no_sub = sapply(seq(1,no_data), no_of_subj, data=indata) # number of qualified subjects in each data set


    ### define function to calculate survival probability for subject i in a single data set - indata[[dataid]]###
    ### return a vector of K MC prediction of survival probabilities. ##
	pred_surv = function(indata, dataid, subid, t, post_par, RE){
		# indata - validation data with covariates, a lit
		# dataid - id of the data set
		# post_par - values of parameters, a list
		# RE - predicted randome effects, a list of lists of matrices
		# subid - subject id in certain validation data
		# t - the survival time
		onedata = indata[[dataid]]$survival_data
		H = matrix(onedata[['H']][subid,], ncol=2)
		W = matrix(onedata[['W']][subid,], ncol=2)

		# get posterior samples of the parameters, should loop over variable par
		out = numeric(0)
		for (i in 1:dim(post_par[[dataid]])[1]){
			u = RE[[dataid]][[1]][[subid]][i,] # note here the dimension of RE! ith row of the matrix
			alpha1 = post_par[[dataid]]$alpha1[i]
			alpha2 = post_par[[dataid]]$alpha2[i]
			gamma = c(post_par[[dataid]]$'gamma[1]'[i], post_par[[dataid]]$'gamma[2]'[i])
			delta = c(post_par[[dataid]]$'delta[1]'[i], post_par[[dataid]]$'delta[2]'[i])
			if(alpha1==0 & alpha2==0) {res=t*exp(W%*%gamma)}
			else{
				res = ((exp(alpha2*u[1]+alpha2*t*u[2]+alpha1*delta[1]*H[,1]+alpha1*delta[2]*H[,2]*t+W%*%gamma)-exp(alpha1*H[,1]*delta[1]+W%*%gamma+alpha2*u[1]))/(alpha2*u[2]+alpha1*H[,2]*delta[2]))
			}
			out = c(out, exp(-res))
		}
		return(out)
	}

    # calculate the predicted suvival probability at t+dt given survived up to time t for all the qualified subject in indata
    model_pred_surv =  vector(length=no_data, mode='list')
    for (i in seq_along(model_pred_surv)){
        model_pred_surv[[i]] = sapply(seq(1,no_sub[i]), pred_surv, indata=indata, t=t1+dt, dataid=i, post_par=post_par, RE=pred_re)/sapply(seq(1,no_sub[i]), pred_surv, indata=indata, t=t1, dataid=i, post_par=post_par, RE=pred_re)
    }
    # return final result, should be of length = length(indata)
    return(model_pred_surv)
}








