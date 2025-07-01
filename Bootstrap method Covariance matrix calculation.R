library("gamlss")
library("gamlss.mx")
library("gamlss.data")

rm(list=ls())

# number of the simulations for the calculation of the covariance matrix
n.sim <- 1000 

# to load the best model
load("BestModel-BI-RIAS-6MP.RData")

# to store the estimates of the regression coefficients
res <- array(NA, dim = c(length(model$mu.coefficients), 4, n.sim))
# to store the estimates of the probabilities of the random effects
res.prob.random.effects <- 
  array(NA, c(n.sim, length(model$prob)))

# the predictor variables
bootstrap.data <- model$data[,c("RID", "DX.bl", "PTGENDER", 
                                "PTEDUCAT", "APOE4", 
                                "married.binary", "TIME",
                                "AGE.centered")]

fit <- matrix(NA, nrow = nrow(bootstrap.data),
              ncol = length(model$prob))
for(i in 1:length(model$prob)) {
  fit[,i] <- fitted(model, K=i)
}

# the number of the observations per subject
num.id <- NA
k <- 0
for(i in unique(model$data$RID)) {
  k <- k+1
  num.id[k] <- sum(model$data$RID==i)
}


for(j in 1:n.sim) {

  rand.eff <- sample(1:length(model$prob), 
                     length(unique(model$data$RID)),
                     prob = model$prob, replace = T)
  
  rand.eff <- rep(rand.eff, num.id)
  
  rand.fit <- NA
  for (i in 1:nrow(bootstrap.data)) {
    rand.fit[i] <- fit[i, rand.eff[i]] 
  } 
    
  bootstrap.data$MMSE.bootstrap <- 
    rbinom(nrow(bootstrap.data), 30, rand.fit)
  
  model.bootstrap <- 
    gamlssNP(cbind(MMSE.bootstrap,30-MMSE.bootstrap)~
               TIME * DX.bl * AGE.centered + 
               PTGENDER + APOE4 + PTEDUCAT + 
               married.binary, K = 6, 
             random=~TIME|RID, mixture="np", 
             family = BI, data = na.omit(bootstrap.data),
             control = NP.control(EMn.cyc = 500))
  
  res[,,j] <- summary(model.bootstrap)
  res.prob.random.effects[j, ] <- model.bootstrap$prob
  if(j%%10==0) save(res, res.prob.random.effects,
                    file = paste("ResultsBootstrapMethodNew.",n.sim,".RData", sep = ""))
  
}

