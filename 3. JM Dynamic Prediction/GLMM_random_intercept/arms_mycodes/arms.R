
rm(list=ls())
library(HI)

expit <- function(x) exp(x)/(1+exp(x))
# log density of beta vector
logden.beta <- function(x, myu) {
  myu.long <- as.vector(mapply(rep, myu, nobs))
  temp <- as.vector(design.matrix %*% as.matrix(x, ncol=1)) + myu.long
  return(sum(y*temp - log(1+exp(temp))))
}

# log density of random intercept u_i
# x: the random intercept u_i
logden.u <- function(x, mybeta, myy, mycovariate, mytau) {
  temp <- mycovariate %*% as.matrix(mybeta, ncol=1) + x
  return(sum(myy * temp - log(1+exp(temp))) - 0.5*x^2*mytau)
}


dat <- read.table("../mydata.txt", header=T)
I <- length(unique(dat$ID))
myID <- unique(dat$ID)
nobs <- NULL # number of observations from each subject
subject <- NULL # subject number of each observation
for (i in 1:I) {
  temp <- subset(dat, ID == myID[i])
  temp.nrow <- nrow(temp)
  nobs <- c(nobs, temp.nrow)
  subject <- c(subject, rep(i,temp.nrow))
}
ncycles <- 5000
p <- 2
design.matrix <- cbind(rep(1,I), dat$x)
y <- dat$y
myshape <- 0.5*I-1 # the alpha parameter for the full conditional of parameter tau


##############################################################################
## MCMC cycles start from here
##############################################################################
set.seed(1353)
temp.beta <- rep(0, p); temp.u <- rep(0.1, I); temp.tau <- 1
write.table(t(c(temp.beta, temp.tau, temp.u)), file="MCMCout.txt", sep="\t",
            col.names=F, row.names=F)

time1 <- Sys.time()
for (n in 2:ncycles) {
  # update beta coefficients
  temp.beta <- arms(temp.beta, logden.beta, function(x,...) all(x>(-10))*all(x<(10)), 1, myu=temp.u)
  # update tau.u
  temp.tau <- rgamma(1, shape=myshape, rate=0.5*sum(temp.u^2))
  # update random intercept U_i
  save.u <- NULL; counter <- 0
  for (i in 1:I) {
    start <- counter + 1
    end <- counter + nobs[i]
    temp <- arms(temp.u[i], logden.u, function(x,...) (x>(-100))*(x<100), 1,
                 mybeta=temp.beta, myy=y[start:end], mycovariate=design.matrix[start:end,],
                 mytau=temp.tau)
    save.u <- c(save.u, temp)
    counter <- end
  }
  temp.u <- save.u; rm(save.u)
  write.table(t(c(temp.beta, temp.tau, temp.u)), file="MCMCout.txt", sep="\t",
              col.names=F, row.names=F, append=T)
}
Sys.time() - time1  
# Time difference of 68 mins

thinning <- 2500
mysamples <- read.table("MCMCout.txt", sep="\t") # 5000 X 1003
final.samples <- mysamples[-(1:thinning),]
str(final.samples)

#beta0.samples <- final.samples[,1]
#beta1.samples <- final.samples[,2]
#sigma.samples <- 1/final.samples[,3]

par.samples <- final.samples[,1:3]
par.samples[,3] <- 1/par.samples[,3]
colnames(par.samples) <- c("beta0", "beta1", "sigmau2")
head(par.samples)

out <- rbind(apply(par.samples, 2, mean), apply(par.samples, 2, sd),
             apply(par.samples, 2, quantile, 0.025),
             apply(par.samples, 2, quantile, 0.9755))
rownames(out) <- c("mean", "sd", "2.5", "97.5")
out <- t(out)
round(out, 3)
#          mean    sd    2.5   97.5
#beta0   -0.485 0.059 -0.600 -0.374
#beta1    1.101 0.070  0.955  1.231
#sigmau2  1.878 0.198  1.501  2.274



# prepare the dataset for WinBUGS
library(R2WinBUGS)
N <- sum(nobs); N
covariate <- dat$x
data <- list(N=N, I=I, x.long=covariate, Y.long=y, subject=subject)
#bugs.data(data, data.file="myBUGSdata.txt")

inits <- function(){
  list(myU=rep(0.0,I), tau.u=1.0, beta=rep(0,p))
}

#inits <- list(myU=rep(0.0,I), tau.u=1.0, beta=rep(0,p))
#bugs.data(inits, data.file="myBUGSinits.txt")

time1 <- Sys.time()
BUGSrun <- bugs(data, inits=inits, model.file = "GLMM.odc",
                parameters = c("beta", "sigma.u2"),
                n.chains = 1, n.iter = 5000, n.burnin=2500, n.thin=1,
                bugs.directory = "C:/WinBUGS14/",debug=F)
Sys.time() - time1 # 2 mins
print(BUGSrun,digits.summary = 3)
#             mean     sd     2.5%      25%      50%      75%    97.5%
#beta[1]    -0.475  0.056   -0.583   -0.515   -0.476   -0.438   -0.366
#beta[2]     1.106  0.066    0.975    1.061    1.105    1.149    1.235
#sigma.u2    1.855  0.193    1.503    1.729    1.844    1.972    2.284






# 1/rgamm is identical to rigamma

alpha <- 2; beta <- 10
junk1 <- 1/rgamma(100000, shape=alpha, rate=beta)
junk2 <- rigamma(100000, alpha=alpha, beta=beta)
mean(junk1); sd(junk1); summary(junk1)
mean(junk2); sd(junk2); summary(junk2)

par(mfrow=c(1,2))
plot(density(junk1))
plot(density(junk2))
