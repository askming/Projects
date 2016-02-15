
# To genereate data sets for making prediction validtion
# In each data set total sample size =520, random sample 500 of them to fit the model
# use the rest 20 to make validatioin.
rm(list=ls())
setwd("/Users/askming/Documents/github/Projects/3. JM Dynamic Prediction/Simulation/Simulation2/R code")

source("sim2data.R")

selected_sim_data = function(N, sur_fun=sim_Ti, longi_fun=sim_longitudinal_data, alpha, tau, beta, delta, gamma, sam_size=520, sam_seleted=500){
	# N - number of data sets to generate
	# sur_fun - function to simulate survival data
	# longi_fun - function to sumualte longitudinal data
	# alpha - association mechanism for JM
	# tau - quantile
	# sam_size - sample size for each data set
	# selected_id - id for subjects that are selected in fitting model
	outdata = vector(mode='list', N)
	# Each element of outdata is a list of four elements: random sample id, full data, selected data and validation data
	for (i in 1:N){
		selected_id = sort(sample(c(1:sam_size), sam_seleted))
		sur_data = sur_fun(n=sam_size, alpha1=alpha[1], alpha2=alpha[2], delta=delta, gamma=gamma)
		longi_data = longi_fun(n=sam_size, sur_data, tau=tau, beta=beta, delta=delta)
		sur_data_selected = lapply(sur_data, function(x, id){x[id,]}, selected_id)
		longi_data_selected = lapply(longi_data, function(x, id){x[id,]}, selected_id)

		full_data = list(survival_data=sur_data, longitudinal_data=longi_data)
		selected_sub = list(survival_data=sur_data_selected, longitudinal_data=longi_data_selected)
		validation_sub = list(survival_data=lapply(sur_data, function(x, id){x[id,]}, -selected_id), longitudinal_data=lapply(longi_data, function(x, id){x[id,]}, -selected_id))
		outdata[[i]] = list(in_id=selected_id, full=full_data, selected=selected_sub, valid=validation_sub)
	}
	outdata
}

# step to sumulate data
t1 = Sys.time()
sel_sim_data = selected_sim_data(N=12, alpha=c(1,1),tau=0.25, beta=c(1,-1), delta=c(-2,1), gamma=c(-3,-1))
Sys.time()-t1
# Time difference of 4.806859 mins

input_data = lapply(sel_sim_data ,`[[`, 'selected')

# loglik = function(data, beta, alpha, delta, gamma, c){
#     # y = data[[2]][['y']]
#     # X = data[[2]][['X']]
#     # J = data[[2]][['J']]

#     # load survival data
#     u = data[[1]][['U']]
#     Ti = data[[1]][['Ti']]
#     H = data[[1]][['H']]
#     W = data[[1]][['W']]
#     event =  data[[1]][['event']]

#     N = 500
#     alpha1=alpha[1]; alpha2=alpha[2]
#     A = B = S = h = L = numeric(N)
#     for (i in 1:N){
#         A[i] <- alpha2*u[i,2] + alpha1*delta[2]*H[i,2]
#         B[i] <- alpha1*delta[1]*H[i,1] + alpha2*u[i,1] + gamma%*%W[i,]
#         S[i] <- exp(- c*exp(B[i])*(round(exp(A[i]^Ti[i]),3)-1)/A[i])
#         h[i] <- c*exp(gamma%*%W[i,] + alpha1*(delta[1]*H[i,1] + delta[2]*H[i,2]*Ti[i]) + alpha2*(u[i,1] + u[i,2]*Ti[i]))
#         L[i] <- (h[i]^event[i])*S[i]/1.0E+08
#     }
#     list(A, B, Ti, S, h, L)
# }

# testout = loglik(input_data[[1]], beta=c(0.1, 0.2), delta=c(0.1, 0.2), gamma=c(0.1, 0.2), alpha=c(0.1, 0.1), c=0.1)

