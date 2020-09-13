######### Confidence Intervals, Equivalence ########
######### and noninferiority tests ##################
######### for time-to-event outcomes ################
############## CASE STUDY veteran #################
########### Author: Kathrin Möllenhoff ##############

source("conf_weibull.R")

library(survival)

data(veteran)

delta <- 0.2 #threshold for non-inferiority/equivalence test for the difference of survival curves

N1 <- sum(veteran$trt==1)
N2 <- sum(veteran$trt==2)
N <- N1+N2

timerange <- c(0,max(veteran$time)) 
t0 <- 80 #time point of intereset (e.g. t0-days survival)

gridc <- seq(0,600,2)
veteran_r <- veteran[veteran$trt==1,] #data standard treatment
veteran_t <- veteran[veteran$trt==2,] #data test treatment
maxt <- max(veteran$time) #end of study

alpha <- 0.05 #significance level/Type I error
B <- 1000 #bootstrap repetitions

conf <- matrix(NA,nrow=4,ncol=length(gridc))
conf2 <- matrix(NA,nrow=4,ncol=length(gridc))

#########################
#########################

dat_r <- veteran_r
dat_t <- veteran_t

#building models
#Weibull models
mod1w <- survreg(Surv(time,status)~1,data=dat_r)
mod2w <- survreg(Surv(time,status)~1,data=dat_t)

#exponential models
mod1e <- survreg(Surv(time,status)~1,data=dat_r,dist="exponential")
mod2e <- survreg(Surv(time,status)~1,data=dat_t,dist="exponential")

#gaussian models
mod1g <- survreg(Surv(time,status)~1,data=dat_r,dist="gaussian")
mod2g <- survreg(Surv(time,status)~1,data=dat_t,dist="gaussian")

#logistic models
mod1l <- survreg(Surv(time,status)~1,data=dat_r,dist="logistic")
mod2l <- survreg(Surv(time,status)~1,data=dat_t,dist="logistic")

#lognormal models
mod1ln <- survreg(Surv(time,status)~1,data=dat_r,dist="lognormal")
mod2ln <- survreg(Surv(time,status)~1,data=dat_t,dist="lognormal")

#loglogistic models
mod1ll <- survreg(Surv(time,status)~1,data=dat_r,dist="loglogistic")
mod2ll <- survreg(Surv(time,status)~1,data=dat_t,dist="loglogistic")

#using AIC as model selection criterion
extractAIC(mod1w)[2]
extractAIC(mod1e)[2]
extractAIC(mod1g)[2]
extractAIC(mod1l)[2]
extractAIC(mod1ln)[2]
extractAIC(mod1ll)[2]
extractAIC(mod2w)[2]
extractAIC(mod2e)[2]
extractAIC(mod2g)[2]
extractAIC(mod2l)[2]
extractAIC(mod2ln)[2]
extractAIC(mod2ll)[2]

#estimate the underlying distributions, distributions of censor times can be estimated as well
g <- function(psi,v) {
  dexp(v,rate=psi)
}
G <- function(psi,v) {
  pexp(v,rate=psi)
}

#assume exponential distribution for censoring times
likelihood <- function(data){
  f<- function(psi){-sum(log(
    g(psi,data$time)^(1-data$status)*(1-G(psi,data$time))^(data$status)
  ))
  }
  return(f)
}

rate_r <- optimize(f=likelihood(dat_r),interval=c(0,5/max(dat_r$time)))$minimum
rate_t <- optimize(f=likelihood(dat_t),interval=c(0,5/max(dat_t$time)))$minimum


#note that survreg uses another parametrization than dweibull!!!
theta_r <- c(1/mod1w$scale,exp(mod1w$coefficients))
theta_t <- c(1/mod2w$scale,exp(mod2w$coefficients))

hazrate_r <- function(t){
  f(c(theta_r[1],theta_r[2]),t)/(1-F(c(theta_r[1],theta_r[2]),t))
}
hazrate_t <- function(t){
  f(c(theta_t[1],theta_t[2]),t)/(1-F(c(theta_t[1],theta_t[2]),t))
}

ratio <- function(t){hazrate_r(t)/hazrate_t(t)}
logratio <- function(t){log(hazrate_r(t)/hazrate_t(t))}

surv_r <- function(t){
  1-F(c(theta_r[1],theta_r[2]),t)
}
surv_t <- function(t){
  1-F(c(theta_t[1],theta_t[2]),t)
}

diff <- function(t){surv_r(t)-surv_t(t)} #S_1-S_2 

#data simulating function (Weibull distribution)
simulWeib <- function(N, lambda, gamma, rateC)
{
  # Weibull latent event times
  # survival time assuming a Weibull distribution
  Tlat <- rweibull(N, shape = gamma, scale = lambda) 
  
  # censoring times (exponential distributed)
  C <- rexp(n=N, rate=rateC)
  
  # follow-up times and event indicators
  study_end <- rep(maxt,N)
  time <- pmin(Tlat, C, study_end)
  status <- as.numeric(Tlat <= C & Tlat <= study_end) #uncensored (experiences an event)
  
  # data set
  data.frame(id=1:N,
             time=time,
             status=status
  )
}

#asymptotic variances
var_loghaz <- function(t,theta,C,N=1){
  gamma=theta[1]
  lambda=theta[2]
  N*t(c(1/gamma+log(t/lambda),-gamma/lambda))%*%C%*%c(1/gamma+log(t/lambda),-gamma/lambda)
}


var_df <- function(t,theta,C,N=1){
  gamma=theta[1]
  lambda=theta[2]
  N*t(c(-exp(-(t/lambda)^gamma)*((t/lambda)^gamma*log(t/lambda)),exp(-(t/lambda)^gamma)*(gamma*(t^gamma)*lambda^(-gamma-1))))%*%C%*%c(-exp(-(t/lambda)^gamma)*((t/lambda)^gamma*log(t/lambda)),exp(-(t/lambda)^gamma)*(gamma*(t^gamma)*lambda^(-gamma-1)))
}

#estimating model parameters by hand (in order to obtain Fisher info for asymptotic bands)
f <- function(theta,v) {
  dweibull(v,shape=theta[1],scale=theta[2])
}
F <- function(theta,v) {
  pweibull(v,shape=theta[1],scale=theta[2])
}

loglikelihood_r <- function(w) {
  theta <- w[1:2]
  -log(
    prod(
      f(theta,dat_r$time)^(dat_r$status)*(1-F(theta,dat_r$time))^(1-dat_r$status)
    )
  )
}

loglikelihood_t <- function(w) {
  theta <- w[1:2]
  -log(
    prod(
      f(theta,dat_t$time)^(dat_t$status)*(1-F(theta,dat_t$time))^(1-dat_t$status)
    )
  )
}

optim1 <- optim(par=theta_r,fn=loglikelihood_r,hessian=TRUE)
optim2 <- optim(par=theta_t,fn=loglikelihood_t,hessian=TRUE)
#calculating observed Fisher info
C1 <- solve(optim1$hessian[1:2,1:2])
C2 <- solve(optim2$hessian[1:2,1:2])

#estimating rates of censoring distribution
rate_r <- optim(par=c(theta_r,0.0001),fn=loglikelihood_r2)$par[3]
rate_t <- optim(par=c(theta_t,0.0001),start,fn=loglikelihood_t2)$par[3]

for(m in 1:length(gridc)){
  t0 <- gridc[m]
  #lower and upper confidence bounds in t0 for the log-hazard-ratio and the difference of survival functions
  conf_haz_d <- conf_weibull(method="delta",alpha=alpha,t0=t0,theta_r=theta_r[1:2],theta_t=theta_t[1:2],C1=C1[1:2,1:2],C2=C2[1:2,1:2])[1:2]
  conf_diff_d <- conf_weibull(method="delta",alpha=alpha,t0=t0,theta_r=theta_r[1:2],theta_t=theta_t[1:2],C1=C1[1:2,1:2],C2=C2[1:2,1:2])[3:4]
  
  bootconf <- conf_weibull(method="bootstrap",alpha=alpha,t0=t0,theta_r=c(theta_r,rate_r),theta_t=c(theta_t,rate_t),C1=C1,C2=C2)
  conf_haz_b <-  bootconf[1:2]
  conf_diff_b <- bootconf[3:4]
  
  #lower bounds
  conf[1,m] <- conf_haz_d[1]
  conf[2,m] <- conf_haz_b[1]
  conf[3,m] <- conf_diff_d[1]
  conf[4,m] <- conf_diff_b[1]
  #upper bounds
  conf2[1,m] <- conf_haz_d[2]
  conf2[2,m] <- conf_haz_b[2]
  conf2[3,m] <- conf_diff_d[2]
  conf2[4,m] <- conf_diff_b[2]
}

#Kaplan Meier estimation
km_fit <- survfit(Surv(time,status) ~ trt, data = veteran)
survdiff(Surv(time,status) ~ trt, data = veteran) #logrank-test


#plot the bands for the difference
curve(diff(x),xlim=c(0,600),xlab="Time (days)",ylab="Difference in Survival",ylim=c(-0.2,0.2))
points(gridc,conf[3,],type="l",lty=2)
points(gridc,conf2[3,],type="l",lty=2)
points(gridc,conf[4,],type="l",lty=3)
points(gridc,conf2[4,],type="l",lty=3)
legend("bottomright",legend=c(expression(paste("Difference in Survival ", S[1](t)-S[2](t))),"Asymptotic confidence bands","Bootstrap confidence bands"),lty=c(1,2,3),cex=0.6)

#plot the bands for the log hazard ratio
curve(logratio(x),xlim=c(0,600),xlab="Time (days)",ylab="log hazard ratio",ylim=c(-0.5,1.2))
points(gridc,conf[1,],type="l",lty=2)
points(gridc,conf2[1,],type="l",lty=2)
points(gridc,conf[2,],type="l",lty=3)
points(gridc,conf2[2,],type="l",lty=3)
legend("bottomright",legend=c(expression(paste("log hazarad ratio, log ", frac(h[1](t),h[2](t)))),"Asymptotic confidence bands","Bootstrap confidence bands"),lty=c(1,2,3),cex=0.4)


