# (pointwise) confidence bands for the (log)-hazard ratio and the difference of survival curves
#for method=delta entire confidence bands are calculated (as a function) and the value at t=t0 returned
#for method=bootstrap bounds are calculated pointwise 

conf_weibull <- function(method,alpha,t0,theta_r,theta_t,C1,C2){
  if(method=="delta"){
    #upper bound for pointwise lower and upper (1-alpha) CI (Delta Method)
    upper_haz <- function(t){
      logratio(t)+qnorm(1-alpha)*sqrt((var_loghaz(t,theta_r,C1,1)+var_loghaz(t,theta_t,C2,1)))
    }
    
    lower_haz <- function(t){
      logratio(t)-qnorm(1-alpha)*sqrt((var_loghaz(t,theta_r,C1,1)+var_loghaz(t,theta_t,C2,1)))
    }
    upper_df <- function(t){
      diff(t)+qnorm(1-alpha)*sqrt((var_df(t,theta_r,C1,1)+var_df(t,theta_t,C2,1)))
    }
    
    lower_df <- function(t){
      diff(t)-qnorm(1-alpha)*sqrt((var_df(t,theta_r,C1,1)+var_df(t,theta_t,C2,1)))
    }
    
    return(c(lower_haz(t0),upper_haz(t0),lower_df(t0),upper_df(t0)))
  }
  if(method=="bootstrap"){
    
    boothaz <- numeric()
    bootsur <- numeric()
    
    #bootstrap to obtain the standard error
    for(l in 1:B){
      dat_rb <- simulWeib(N=N1, lambda=theta_r[2], gamma=theta_r[1], rateC=theta_r[3])
      dat_tb <- simulWeib(N=N2, lambda=theta_t[2], gamma=theta_t[1], rateC=theta_t[3])
      loglikelihood_rb <- function(w) {
        theta <- w[1:2]
        -log(
          prod(
            f(theta,dat_rb$time)^(dat_rb$status)*(1-F(theta,dat_rb$time))^(1-dat_rb$status)
          )
        )
      }
      
      loglikelihood_tb <- function(w) {
        theta <- w[1:2]
        -log(
          prod(
            f(theta,dat_tb$time)^(dat_tb$status)*(1-F(theta,dat_tb$time))^(1-dat_tb$status)
          )
        )
      }
      
      theta_rb <- optim(par=theta_r,fn=loglikelihood_rb)$par
      theta_tb <- optim(par=theta_t,fn=loglikelihood_tb)$par
      
      hazrate_rb <- function(t){
        f(c(theta_rb[1],theta_rb[2]),t)/(1-F(c(theta_rb[1],theta_rb[2]),t))
      }
      hazrate_tb <- function(t){
        f(c(theta_tb[1],theta_tb[2]),t)/(1-F(c(theta_tb[1],theta_tb[2]),t))
      }
      logratiob <- function(t){log(hazrate_rb(t)/hazrate_tb(t))}
      
      surv_rb <- function(t){
        1-F(c(theta_rb[1],theta_rb[2]),t)
      }
      surv_tb <- function(t){
        1-F(c(theta_tb[1],theta_tb[2]),t)
      }
      
      diffb <- function(t){surv_tb(t)-surv_rb(t)}
      
      if(is.na(logratiob(t0))){boothaz[l] <- mean(boothaz,na.rm=TRUE);next}
      if(logratiob(t0)<Inf & -Inf<logratiob(t0)){
      boothaz[l] <- logratiob(t0)
      }else{boothaz[l] <- mean(boothaz)}
      if(is.na(diffb(t0))){bootsur[l] <- mean(bootsur,na.rm=TRUE);next}
      if(diffb(t0)<Inf & -Inf<diffb(t0)){bootsur[l] <- diffb(t0)}else{bootsur[l] <- mean(bootsur,na.rm=TRUE)}
    }
    
    up_haz <- logratio(t0)+qnorm(1-alpha)*sqrt(var(boothaz,na.rm=TRUE))
    lo_haz <- logratio(t0)-qnorm(1-alpha)*sqrt(var(boothaz,na.rm=TRUE))
    up_df <- diff(t0)+qnorm(1-alpha)*sqrt(var(bootsur,na.rm=TRUE))
    lo_df <- diff(t0)-qnorm(1-alpha)*sqrt(var(bootsur,na.rm=TRUE))
    
    return(c(lo_haz,up_haz,lo_df,up_df))
  }
}