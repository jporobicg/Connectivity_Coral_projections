## Title: 4_01_model_parameters.R
## Author: Mandy Cheung
## Date: 06-09-2024
## Description: This script interpolate the model parameters for 
##              survival and competency functions

#-------------------------------------------------------------------------------#

# Load package
library(tidyverse)
library(ggplot2)


##-----------------------------------------------------------------------------##
# ----------------------------- (1) Write data ----------------------------------
##-----------------------------------------------------------------------------##

# record model parameters

# survival model parameters - generalised weibull
# bathtub-shaped mortality curves: initial mortality is high and decrease with age,
# and increase again as larvae become very old
temp<-c(25,27,29)
lambda<-c(0.02954, 0.000138, 0.0000125)       # mortality rate (scale parameter)
lambda_ucl<-c(0.04098, 0.00234, 0.00185)      # upper cl
lambda_lcl<-c(0.02117, 0.000002, 0)           # lower cl
v<-c(0.4612,0.2069,0.1386)                    # shape parameter of weibull in bathtub (shape parameter)
v_ucl<-c(0.508,0.307,0.2019)
v_lcl<-c(0.3789,0.134,0.0892)
sigma<-c(0.0000001275,2.1545,2.3833)          # shape parameter of bathtub curve (location parameter)
sigma_ucl<-c(0.3617, 2.1634, 2.3557)
sigma_lcl<-c(0, 1.6623, 2.2419)


surv<-cbind(temp,lambda,lambda_ucl,lambda_lcl,v,v_ucl,v_lcl,sigma,sigma_ucl,sigma_lcl)
surv<-surv%>%data.frame()
remove(temp,lambda,lambda_ucl,lambda_lcl,v,v_ucl,v_lcl,sigma,sigma_ucl,sigma_lcl)


# competency model parameters - weibull 
temp<-c(25,27,29)
tc<-c(5.38,4.89,3.87)                   # the age at which larvae reaches maximum competency
tc_ucl<-c(5.74,5.01,4.01)
tc_lcl<-c(5.04,4.77,3.74)
alpha<-c(0.4497,0.4497,0.4497)          # scale parameter (competence acquisition rate)
alpha_ucl<-c(0.5497,0.5497,0.5497)
alpha_lcl<-c(0.3672,0.3672,0.3672)
beta<-c(0.01623,0.02669,0.02996)        # shape parameter (competence loss rate)
beta_ucl<-c(0.02559,0.03770,0.04094)
beta_lcl<-c(0.01047,0.01847,0.02126)
n<-c(0.3981,0.3981,0.3981)              # exponential decay parameter in early decline phase
n_ucl<-c(0.4988,0.4988,0.4988)
n_lcl<-c(0.3236,0.3236,0.3236)

compet<-cbind(temp,tc,tc_ucl,tc_lcl,alpha,alpha_ucl,alpha_lcl,beta,beta_ucl,beta_lcl,n,n_ucl,n_lcl)
compet<-compet%>%data.frame()
remove(temp,tc,tc_ucl,tc_lcl,alpha,alpha_ucl,alpha_lcl,beta,beta_ucl,beta_lcl,n,n_ucl,n_lcl)



##-----------------------------------------------------------------------------##
# --------------------------- (2) Plot to check --------------------------------
##-----------------------------------------------------------------------------##


## (A) create the generalised weibull survival function ------------------------

bathtub_curve<-function(t, lambda, v, sigma) {
  (lambda*v*(lambda*t)^(v-1))/(1-sigma*(lambda*t)^v)
}



decayfn<-function (t, lambda, v, sigma){
  # bathtub_curve<-function(t, lambda, v, sigma) {
  #  (lambda*v*(lambda*t)^(v-1))/(1-sigma*(lambda*t)^v) # attempt 1-2, 4
  # }
  bathtub_curve <- function(t, lambda, v, sigma) { # bathtub-shaped mortality rate #third attempt
    term <- 1 - sigma * (lambda * t)^v
    if (term <= 0) {               # If the denominator <= 0
      term <- .Machine$double.eps  # Change the denominator to the smallest positive number 2.220446e-16
    }
    (lambda * v * (lambda * t)^(v - 1)) / term
  }
  decay <- numeric(length(t)) # Initialize decay vector
  # Calculate decay
  for (age in seq_along(t)) {
    area <- integrate(Vectorize(bathtub_curve), lower=0, 
                      upper=t[age], lambda=lambda, 
                      v=v, sigma=sigma)$value
    decay[age] <- exp(-area)
  }
  return(decay)
}




# the age/time df
x<-seq(1,140,by=1)
# plot
plot(x, decayfn(x, surv$lambda[1], surv$v[1], surv$sigma[1]), 
     type="l", lwd=3, col="steelblue",
     xlab="Age",ylab="Decay")
lines(decayfn(x, surv$lambda[2], surv$v[2], surv$sigma[2]), lwd=3,col="pink")
lines(decayfn(x, surv$lambda[3], surv$v[3], surv$sigma[3]),lwd=3,col="orange")
legend("topright", legend = c("25C", "27C", "29C"), 
       col = c("steelblue", "pink", "orange"), lwd = 3)




## (B) create the weibull competence function ----------------------------------

competence<-function(t, tc, alpha, beta, n){
  competence <- c()
  for (age in  seq_along(t)) {
    if (t[age] < tc) {
      area <- 0
    } else {
      fxtau_early <- function(tau) { # a*exp(-a*(tau-tc))*exp(-((b1*(time-tau))^v1))
        alpha * exp(-alpha * (tau - tc)) * exp(- ((beta * (t[age] - tau))^n))
      }
      area <- integrate(fxtau_early, tc, t[age])$value
    }
    competence <- c(competence, area)
  }
  return(competence)
}

# the time/age df
x<-seq(1,100,by=1)

# plot
plot(competence(x, compet$tc[1], compet$alpha[1], compet$beta[1], compet$n[1]),
     type="l", lwd=3, col="steelblue",
     xlab="Age",ylab="Competency")
lines(competence(x, compet$tc[2], compet$alpha[2], compet$beta[2], compet$n[2]),
      lwd=3, col="pink")
lines(competence(x, compet$tc[3], compet$alpha[3], compet$beta[3], compet$n[3]),
      lwd=3, col="orange")
legend("topright", legend = c("25C", "27C", "29C"), 
       col = c("steelblue", "pink", "orange"), lwd = 3)



##-----------------------------------------------------------------------------##
# ------------------------ (3) Extrapolate parameters --------------------------
##-----------------------------------------------------------------------------##


## (A) Survival -----------------------------------------------------------------

# lambda
# First check the distribution of model parameters 
plot(x=surv$temp, y=surv$lambda)  # exponential curve
lamba.int<-splinefun(x = c(25, 27, 29), y = log(surv$lambda), method = "monoH.FC")
lamba.int <- exp(lamba.int(c(25, 26, 27, 28, 29, 30, 31, 32)))
lamba.int<-data.frame(cbind(c(25,26,27,28,29,30,31, 32),lamba.int))
plot(x=lamba.int$V1, y=lamba.int$lamba.int)  # exponential curve
#v
plot(x=surv$temp, y=surv$v) # exponential curve
v.int <- splinefun(x = c(25, 27, 29), y = log(surv$v), method = "monoH.FC")
v.int <- exp(v.int(c(25, 26, 27, 28, 29, 30, 31, 32)))
v.int<-data.frame(cbind(c(25,26,27,28,29,30,31, 32),v.int))
plot(x=v.int$V1, y=v.int$v.int)  # exponential curve
# sigma 
plot(x=surv$temp, y=surv$sigma) # logarithm curve
sigma.int <- splinefun(x = c(25, 27, 29), y = surv$sigma, method = "monoH.FC")
sigma.int <- sigma.int(c(25, 26, 27, 28, 29, 30, 31, 32))
sigma.int<-data.frame(cbind(c(25,26,27,28,29,30,31, 32),sigma.int))
plot(x=sigma.int$V1, y=sigma.int$sigma.int)  # exponential curve

extpl_spline<-data.frame(cbind(lamba.int,v.int[2],sigma.int[2]))
colnames(extpl_spline)<-c("temp","lambda","v","sigma")

# the age/time df
x<-seq(1,120,by=1)

# plot

plot(x, decayfn(x, extpl_spline$lambda[1], extpl_spline$v[1], extpl_spline$sigma[1]), 
     type="l", lwd=3, col="steelblue",
     xlab="Age",ylab="Decay", ylim=c(0,1))
lines(decayfn(x, extpl_spline$lambda[2], extpl_spline$v[2], extpl_spline$sigma[2]),lwd=3,col="lightblue")
lines(decayfn(x, extpl_spline$lambda[3], extpl_spline$v[3], extpl_spline$sigma[3]),lwd=3,col="darkolivegreen")
lines(decayfn(x, extpl_spline$lambda[4], extpl_spline$v[4], extpl_spline$sigma[4]),lwd=3,col="orange")
lines(decayfn(x, extpl_spline$lambda[5], extpl_spline$v[5], extpl_spline$sigma[5]),lwd=3,col="pink")
lines(decayfn(x, extpl_spline$lambda[6], extpl_spline$v[6], extpl_spline$sigma[6]),lwd=3,col="salmon")
lines(decayfn(x, extpl_spline$lambda[7], extpl_spline$v[7], extpl_spline$sigma[7]),lwd=3,col="tomato")
lines(decayfn(x, extpl_spline$lambda[8], extpl_spline$v[8], extpl_spline$sigma[8]),lwd=3,col="red")
legend("topright", legend = c("25C","26C","27C","28C","29C","30C","31C","32C"), 
       col = c("steelblue", "lightblue","darkolivegreen",
               "orange","pink","salmon","tomato","red"), lwd = 3)


# (B) Competency ---------------------------------------------------------------

## linear extrapolation


tc.int <- splinefun(x = c(25, 27, 29), y = compet$tc, method = "monoH.FC")
tc.int <- tc.int(c(25, 26, 27, 28, 29, 30, 31,32))
extpl_tc<-data.frame(cbind(c(25,26,27,28,29,30,31,32),tc.int))
plot(x=extpl_tc$V1, y=extpl_tc$tc.int)  # exponential curve

# temp-independent
alpha.int <- splinefun(x = c(25, 27, 29), y = compet$alpha, method = "monoH.FC")
alpha.int<-alpha.int(c(25, 26, 27, 28, 29, 30, 31,32))
extpl_alpha<-data.frame(cbind(c(25,26,27,28,29,30,31,32),alpha.int))
plot(x=extpl_alpha$V1, y=extpl_alpha$alpha.int)

beta.int <- splinefun(x = c(25, 27, 29), y = compet$beta, method = "monoH.FC")
beta.int <- beta.int(c(25, 26, 27, 28, 29, 30, 31,32))
extpl_beta<-data.frame(cbind(c(25,26,27,28,29,30,31,32),beta.int))
plot(x=extpl_beta$V1, y=extpl_beta$beta.int)  # exponential curve

# temp-independent
n.int <- splinefun(x = c(25, 27, 29), y = compet$n, method = "monoH.FC")
n.int<-n.int(c(25, 26, 27, 28, 29, 30, 31,32))
extpl_n<-data.frame(cbind(c(25,26,27,28,29,30,31,32),n.int))
plot(x=extpl_n$V1, y=extpl_n$n.int)

# combine
extpl_compet<-cbind(data.frame(extpl_tc), extpl_alpha$alpha.int, extpl_beta$beta.int, extpl_n$n.int)
colnames(extpl_compet)<-c("temp","tc","alpha","beta","n")

# the age/time df
x<-seq(1,100,by=1)
plot(competence(x, extpl_compet$tc[1], extpl_compet$alpha[1], extpl_compet$beta[1], extpl_compet$n[1]), 
     type="l", lwd=3, col="steelblue",
     xlab="Age",ylab="Competency")
lines(competence(x, extpl_compet$tc[2], extpl_compet$alpha[2], extpl_compet$beta[2], extpl_compet$n[2]),lwd=3,col="lightblue")
lines(competence(x, extpl_compet$tc[3], extpl_compet$alpha[3], extpl_compet$beta[3], extpl_compet$n[3]),lwd=3,col="darkolivegreen")
lines(competence(x, extpl_compet$tc[4], extpl_compet$alpha[4], extpl_compet$beta[4], extpl_compet$n[4]),lwd=3,col="orange")
lines(competence(x, extpl_compet$tc[5], extpl_compet$alpha[5], extpl_compet$beta[5], extpl_compet$n[5]),lwd=3,col="pink")
lines(competence(x, extpl_compet$tc[6], extpl_compet$alpha[6], extpl_compet$beta[6], extpl_compet$n[6]),lwd=3,col="salmon")
lines(competence(x, extpl_compet$tc[7], extpl_compet$alpha[7], extpl_compet$beta[7], extpl_compet$n[7]),lwd=3,col="tomato")
lines(competence(x, extpl_compet$tc[8], extpl_compet$alpha[8], extpl_compet$beta[8], extpl_compet$n[7]),lwd=3,col="red")
legend("topright", legend = c("25C","26C","27C","28C","29C","30C","31C","32C"), 
       col = c("steelblue", "lightblue","darkolivegreen",
               "orange","pink","salmon","tomato","red"), lwd = 3)



save.image(file="./memory/4_01_model_para_WS.RData")


## parameters 

extpl_compet
extpl_spline


# for use in python ------------------

decay_para={'temp':[25,26,27,28,29,30,31,32], 
            'lmbda1':[2.954000e-02, 1.677537e-03, 1.380000e-04, 3.450817e-05, 1.250000e-05, 3.762058e-06, 1.132246e-06, 3.407661e-07],
            'v1':[0.46120000, 0.30126014, 0.20690000, 0.16515002, 0.13860000, 0.11343958, 0.09284659, 0.07599191],
            'sigma1':[0.0000001275, 1.2623016324, 2.1545000000, 2.3369825925, 2.3833000000, 2.4480132916, 2.5127265832, 2.5774398748]}

compet_para={'temp':[25,26,27,28,29,30,31,32], 
              'tc':[5.380000, 5.168125, 4.890000, 4.413125, 3.870000, 3.360000, 2.850000, 2.340000],
              'alpha':[0.4497, 0.4497, 0.4497, 0.4497, 0.4497, 0.4497, 0.4497, 0.4497],
              'beta1':[0.01623000, 0.02190937, 0.02669000, 0.02877437, 0.02996000, 0.03159500, 0.03323000, 0.03486500],
              'v':[0.3981, 0.3981, 0.3981, 0.3981, 0.3981, 0.3981, 0.3981, 0.3981]}

