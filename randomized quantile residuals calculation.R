# because the output of the function residuals produces values Inf 
# we write our own code to compute the randomized quantile residuals

library("gamlss")
library("gamlss.mx")
library("gamlss.data")
library(ggplot2)
library(qqplotr)

rm(list=ls())
load("BestModel-BI-RIAS-6MP.RData")
set.seed(9)

residuals.model <- residuals(model)
n.obs <- length(residuals.model)
observed.data <- model$y[1:n.obs]

sim.res <- rep(NA, n.obs)
for (i in 1:n.obs) {
  if(residuals.model[i] == Inf) {
    prob1 <- pbinom(29, 30, fitted(model)[i])
    prob1 <- 0.5*(prob1 + 1)
    sim.res[i] <- qnorm(runif(1, prob1, 0.9999999999999999))
  } else { 
    if(observed.data[i] != 0) {
      prob1 <- pbinom(observed.data[i]-1, 30, 
                      fitted(model)[i])
      prob1 <- 0.5*(prob1 + pbinom(observed.data[i], 30, 
                                   fitted(model)[i]))
      prob2 <- pbinom(observed.data[i]+1, 30, 
                      fitted(model)[i])
      prob2 <- 0.5*(prob2 + pbinom(observed.data[i], 30, 
                                   fitted(model)[i]))
      sim.res[i] <- qnorm(runif(1, prob1, prob2))
    } else { 
      prob2 <- 0.5*(pbinom(observed.data[i], 30, 
                           fitted(model)[i]))
      sim.res[i] <- qnorm(runif(1, 0.000000000001, 
                                prob2))
    }
  }
  
}



