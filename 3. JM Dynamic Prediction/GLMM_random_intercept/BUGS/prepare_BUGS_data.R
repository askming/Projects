
rm(list=ls())
library(R2WinBUGS)
expit <- function(x) exp(x)/(1+exp(x))


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


# prepare the dataset for WinBUGS
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

save(BUGSrun, file="BUGSrun.rda")

