for (myloop in 1:nsim)
{
  set.seed(mySeeds[myloop])
  # generate random effects
  U <- mvrnorm(N, mu=rep(0,2), Sigma=Sigma.U)

  # generate latent variable theta
  theta <- NULL
  for (j in 1:T) theta <- cbind(theta, U[,1]+ (beta1[1] + beta1[2]*treat.pts + U[,2])*time[j])

  # generate longitudinal outcomes
  Y.conti <- array(NA,dim = c(N, T))
  Y.conti.2 <- array(NA,dim = c(N, T))
  Y.ordi <- array(NA,dim = c(N, T, K.ordi))
  for (i in 1:N) {
    for (j in 1:T) {
      Y.conti[i,j] <- rnorm(1, mean=a.conti+b.conti*theta[i,j],
                            sd=sd.conti)
      for (k in 1:K.ordi) {
        psi <- c(rep(0,n.ordi[k]-1),1)
        prob.y <- rep(0,n.ordi[k])
        for (l in 1:(n.ordi[k]-1)) psi[l] <- expit(a.ordi[k,l]
                                  - b.ordi[k]*theta[i,j])
        prob.y[1] <- psi[1]
        for (l in 2:n.ordi[k]) prob.y[l] <- psi[l]-psi[l-1]
        Y.ordi[i,j,k] <- sample(1:n.ordi[k],1,prob=prob.y)
} }
}
  # turn response Y and theta into long format
  Y.conti.temp <-  Y.conti.2.temp <- Y.ordi.temp <- theta.temp <- NULL
  for (i in 1:N) {
    Y.conti.temp <- c(Y.conti.temp, Y.conti[i,])
    Y.ordi.temp <- rbind(Y.ordi.temp, Y.ordi[i,,])
    theta.temp <- c(theta.temp, theta[i,])
  }
  Y.conti <- Y.conti.temp; rm(Y.conti.temp)
  Y.ordi <- Y.ordi.temp; rm(Y.ordi.temp)
  theta <- theta.temp; rm(theta.temp)

  # simulation for survival part
S <- runif(N) #simualte survival function from unif(0,1)
Ti.obs <- -log(S)/(exp(gam*treat.pts+ omega2*U[,1] + omega3*U[,2])*h0)
C <- runif(N, 10, 20)
tee <- pmin(Ti.obs, C) #observed time is min of censored and true
event <- ifelse(tee==Ti.obs, 1, 0) #0: censored ; 1: event;
# Combined longitudinal and survival data,
# measures after event/censor are deleted
subject <- rep(1:N, each = T)
treat <- rep(treat.pts, each = T)
t <- rep(time, N)
obs <- length(t)  # obs: total number of observations across subjects
zero <- rep(0, 2)
month <- rep(c(0,1,3,9,15), N)
all <- cbind(Y.conti, Y.ordi, subject, treat, t, month)
unbalance.all <- NULL
for (i in 1:obs)
{
  #keep the measures before event/censor
  if (all[i, 7] < tee[subject[i]])
  {unbalance.all <- rbind(unbalance.all, all[i,])}
  #delete measures after/at the same time of event/censor
  else {unbalance.all <- unbalance.all}
}
obs <- nrow(unbalance.all)
subject <- unbalance.all[,4]
treat <- unbalance.all[,5]
t <- unbalance.all[,6]
Y.conti <- unbalance.all[,1]
Y.ordi <- unbalance.all[,2:3]
# output BUGS data
zeros <-  numeric(N)
n <- n.ordi
data <- list("N", "n", "subject", "obs", "K.ordi",
             "Y.conti", "Y.ordi", "treat", "zero",
             "t", "treat.pts", "tee","event","zeros")
## generate data for BUGS
outname <- myloop
outputname <- paste("BUGSdata/", "data", outname, ".txt", sep="")
  bugs.data(data, data.file=outputname)
  ## transfer BUGS data to JAGS format
  bugs.data.path <- paste("./", "BUGSdata/", "data", outname, ".txt", sep="")
  jags.data.path <- paste("./", "JAGSdata/", "data", outname, ".txt", sep="")
  bugs2jags(bugs.data.path, jags.data.path)
}