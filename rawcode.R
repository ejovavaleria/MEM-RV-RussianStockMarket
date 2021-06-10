#libraries
library(matrixcalc)
library(Matrix)
library(foreach)

#data
rv_dta <- read.csv("~/THESIS/SBER.csv", sep = ";")
rv_dta = RV_sber

train1 <- 1:400
valid1 <- 401:nrow(rv_dta)

rvall <- rv_dta$rv
rv  <- rv_dta$rv[train1]
rvv <- rv_dta$rv[valid1]

#modeling and estimating

### to modify:
# 1) mu(...): lmu eq-n, add factor
# 2) same with dmu(...)
# 3) in fits <- ... change dimensionality of paro
# 4) same change in muf(...)

mu <- function(par) {
  muo <- mean(rv)
  rvo <- muo
  
  lmu <- rep(NA,length(rv))
  lmu[1] <- par[1]+par[2]*log(rvo)+par[3]*log(muo)
  for(t in 2:length(rv)) {
    lmu[t] <- par[1]+par[2]*log(rv[t-1])+par[3]*lmu[t-1]
  }
  
  return(exp(lmu))
}
#mu(c(0.2,-0.2,+0.05))

dmu <- function(par) {
  mu1 <- mu(par)
  
  foreach(t=2:length(rv), .combine=rbind) %do% {
    mu1[t]*c(1, log(rv[t-1]), log(mu1[t-1]))
  }
}
#dmu(c(0.2,-0.2,+0.05))

nLL <- function(par) {
  mu1 <- mu(par)
  ll  <- -rv/mu1 - log(mu1)
  return(-sum(ll))
}
#nLL(c(0.2,-0.2,+0.05))

eps <- function(par) {
  rv/mu(par)
}
#eps(c(0.2,-0.2,+0.05))

s <- function(par) {
  mu1 <- mu(par)
  dmu1 <- dmu(par)
  foreach(t=2:length(rv),.combine=rbind) %do% {
    dmu1[t-1,]*(rv[t]/mu1[t]^2-1/mu1[t])
  }
}
#s(c(0.2,-0.2,+0.05))

nLLgr <- function(par) {
  -colSums(s(par))
}
nLLgr(c(0.2,-0.2,+0.05))

vcov <- function(par,Hess) {
  s1 <- s(par)
  n <- nrow(s1)
  
  bunn <- (-Hess/n)^(-1)
  
  meat <- Reduce('+',lapply(1:n, function(i) matrix((s1[i,]),ncol(s1),1)%*%matrix(s1[i,],1,ncol(s1))))/n

  forceSymmetric(bunn %*% meat %*% bunn)
}

### estimation
B <- 300
library(doParallel)

cl <- makePSOCKcluster(8)
registerDoParallel(cl)
fits <- foreach(b=1:B,
                .errorhandling = 'remove',
                .packages = c("foreach")) %dopar% { 
  paro <- runif(3,-5,5)
  fit1 <- optim(paro, nLL, nLLgr, method="BFGS",hessian=TRUE)
}
stopCluster(cl)

nLLs <- foreach(fit=fits,.combine=c) %do% {
  if(fit1$convergence == 0) {
    fit$value
  } else {
    +9999
  }
}
fit1 <- fits[[which.min(nLLs)]]

is.positive.definite(fit1$hessian)
colMeans(s(fit1$par))

#vcov0 <- (fit1$hessian/length(rv))^(-1)
vcov1 <- vcov(fit1$par,fit1$hessian)

is.positive.definite(matrix(as.numeric(vcov1),ncol(vcov1),ncol(vcov1)))
det(vcov1)

se1 <- sqrt(diag(vcov1)/length(rv))
t1 <- fit1$par / se1
pv1 <- (1-pnorm(abs(t1)))/2

eps1 <- eps(fit1$par)
acf(eps1)

### validation
muf <- function(par) {
  lmu <- log(rvall)
  for(t in valid1) {
    lmu[t] <- par[1]+par[2]*log(rvall[t-1])+par[3]*lmu[t-1]
  }
  return(exp(lmu[valid1]))
}
sqrt(mean((log(rvv) - log(muf(fit1$par)))^2)) #validation RMSE for log(RV)


