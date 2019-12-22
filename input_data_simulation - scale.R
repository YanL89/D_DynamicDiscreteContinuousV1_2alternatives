library(dcmodels)
rm(list=ls())
setwd("C:/Users/Yan Liu/Documents/UMD research/recursive probit/2dynamic binary D-C model - series covariances/simulate_data_scale")
source("dynData.R")
source("create_X.R")
source("dynUtils.R")

#model description
#U1 = time0*beta1 + cost0*beta2
#U2 = time1*beta1 + cost1*beta3 + ASC*beta4 + income*beta5
#Y = attr1*beta6 + attr2*beta7

#data generation---------------------------------------------------------------------
nobs_sim = 1000
nalt_sim = 2
nvarD_disc = 2
nvarG_disc = 2
nvarD_reg = 2
time_sim = 20
n_sim = 1000

#simulate "global" data
ASC_sim = rep(1, nobs_sim)
income_sim = rnorm (nobs_sim, 1.5, 1)
global = matrix(data = NA, nrow = nobs_sim, ncol = nvarG_disc, byrow = FALSE)
global[, 1] = ASC_sim
global[, 2] = income_sim
global = rbind(c("ASC", "income"), global)
write.table(global, file = "global.txt", row.names=FALSE, col.names=FALSE)
#global = read.table("global.txt", sep = " ", head = T)

#simulate "dyn" data
dyn = list()
for(t in 1:time_sim){
  dyn[[t]] = list()
  dyn[[t]] = matrix(data = NA, nrow = nobs_sim, ncol = nvarD_disc * nalt_sim + nvarD_reg, byrow = FALSE)
  time0_sim = runif(nobs_sim, min = 0, max = 2 - 0.02 * t + 0.002)
  time1_sim = runif(nobs_sim, min = 2 - 0.02 * t, max = 4 - 0.02 * t + 0.002)
  cost0_sim = rnorm(nobs_sim, 2 + 0.02 * t - 0.002, 1)
  cost1_sim = rnorm(nobs_sim, 1 + 0.02 * t - 0.002, 1)
  attr1 = rnorm(nobs_sim, 3 + 0.02 * t - 0.002, 1)
  attr2 = rnorm(nobs_sim, 4 + 0.02 * t - 0.002, 1)
  dyn[[t]][, 1] = time0_sim
  dyn[[t]][, 2] = time1_sim
  dyn[[t]][, 3] = cost0_sim 
  dyn[[t]][, 4] = cost1_sim
  dyn[[t]][, 5] = attr1 
  dyn[[t]][, 6] = attr2  
  dyn[[t]] = rbind(c("time.0", "time.1", "cost.0", "cost.1", "attr1", "attr2"), dyn[[t]])
  data_name = paste("dyn", t, sep = "")
  write.table(dyn[[t]], file = paste(data_name, "txt", sep = "."), row.names=FALSE, col.names=FALSE)
  #dyn[[t]] = read.table(paste(data_name, "txt", sep = "."), sep = " ", head = T)
}


#define true coefficients - beta-------------------------------------------------------
#beta = c(time, cost0, ASC1, cost1, income1, attr1, attr2) where "time" is generic, "ASC, "cost" and "income" are specific  
beta = c(1, -0.5, -2, -3, 1, 1, 0.5)
beta_disc = beta[1:5]
beta_reg = beta[6:7]


#define time-dependent D-C covariance matrix
S_disc = diag(nalt_sim)
M = matrix(c(-1, 1), nalt_sim-1, nalt_sim)
Diff_S_disc = M %*% S_disc %*% t(M)
M1 = matrix(c(-1, -1, 1, 0, 0, 1), nalt_sim, nalt_sim+1)
S = list()
Diff_S = list()
Diff_S_1 = list()
TrueDiff_S = list()

for(t in 1:time_sim){
  S_reg = diag(1)*(1 + 0.05*t)
  S12 = as.matrix(c(0, 0.02*t), nalt_sim, 1)
  S[[t]] = diag(nalt_sim+1)
  S[[t]][1:nalt_sim, 1:nalt_sim] = S_disc
  S[[t]][1:nalt_sim, nalt_sim+1] = S12
  S[[t]][nalt_sim+1, 1:nalt_sim] = t(S12)
  S[[t]][nalt_sim+1, nalt_sim+1] = S_reg
  Diff_S[[t]] = M1 %*% S[[t]] %*% t(M1)
  Diff_S_1[[t]] = Diff_S[[t]][1, 1]
}
scale = sqrt(Diff_S_1[[1]])
#Diff_S[[t]] is the true covariance in difference at time t


#calculate utilities-------------------------------------------------------------------
#define model specification
spec = list(
  D = "dyn",
  Global = "global.txt",
  nvarD_disc = 2,
  nvarG_disc = 2,
  nvarD_reg = 2,
  
  generic = list(c("time.0", "time.1")),
  specific = list(c("cost.0"), c("ASC", "cost.1", "income")),
  reg = c("attr1", "attr2"),
  
  modifyD = function(D,Global,t){ 
    Dyan = cbind(D[[t]], Global)
    return(Dyan)
  },
  nTime = 17,  #actual time, not total time
  nLook = 3,  #look ahead period
  transition = rbind(c(1, 1), c(1, 1)),
  stopAlt = c() , #if you have more than one alternative 
  # that halts the decision process, put that in 
  # a vector (like stopAlt = c(1,2,10) if 
  # alt 1, 2 and 10 are stopping alternatives)
  method = "pmvnorm"
)

spec = checkFillDynSpec(spec)

#define function computeArgs
computeArgs = function(spec, D){
  X = list()
  RegX = list()
  
  for(t in 1:spec$nTime){
    Dreg = spec$modifyD(spec$D, spec$Global, t)
    RegX[[t]] = as.matrix(Dreg[, spec$reg])
    X[[t]] = list()
    for(l in 1:(spec$nLook+1)){
      Dt = spec$modifyD(spec$D, spec$Global, t+l-1)   #get all attributes at t
      X[[t]][[l]] = create_X(spec$generic, spec$specific, Dt)
    }
  }
  
  # check if we reached stop point already
  stopMat = matrix(0, spec$nObs, spec$nTime) 
  
  list(spec = spec, X = X, RegX = RegX, stop = stopMat, method = spec$method)    
}

args = computeArgs(spec, D)  
args = c(args, spec)

#define function computeMisc
computeMisc = function(spec, D){
  n = c()
  comIndex = 1
  for(com in spec$generic){
    n = c(n,paste("common_time",comIndex,sep="_"))
    comIndex = comIndex + 1
  }
  for(s in spec$specific)
    n = c(n,s)
  for(s in spec$reg)
    n = c(n,s)
  
  ndim = spec$nAlt - 1 + 1
  if(ndim > 1){
    for(t in 1:spec$nTime){
      for(i in 2:ndim)
        for(j in 1:i){
          n = c(n, paste("L_",i,j,"_",t,sep=""))
        }
      }
    }
  list(names = n)
}

misc = computeMisc(spec, D)


#calculate utility and choice
U = list()
choice = matrix(data = NA, nrow = nobs_sim, ncol = spec$nTime, byrow = FALSE)
Y = matrix(data = NA, nrow = nobs_sim, ncol = spec$nTime, byrow = FALSE)

for(t in 1:args$spec$nTime){
  #calculate utility for discrte part
  for(l in 1:(args$spec$nLook+1)){
    U[[l]] = matrix(args$X[[t]][[l]] %*% beta_disc, args$spec$nObs, args$spec$nAlt, byrow = T)
    U[[l]] = U[[l]]/scale
  }
  
  if(nalt_sim > 2){
    fU = getfU(U, args$spec$transition)  #logsum of recursive logit
  } else if(nalt_sim == 2) {
    fU = getProbitU(U, Diff_S_1[[t]], args$spec$transition)  #obtain instant utility plus downstream utility at  time t
  }
  
  #calculate Y for continuous part
  Y[, t] = args$RegX[[t]] %*% beta_reg
  
  #simulate error components AND obatin final utilities
  Err = matrix(data = NA, nrow = nobs_sim, ncol = nalt_sim + 1, byrow = FALSE)
  Utility = matrix(data = NA, nrow = nobs_sim, ncol = nalt_sim, byrow = FALSE)
  Prob = rep(0, nobs_sim)
  #RegErr = rep(0, nobs_sim)
  for(j in 1:n_sim){
  Err = rmvnorm(n = nobs_sim, mean = rep(0, nalt_sim + 1), sigma = S[[t]])
  Utility = fU + Err[, 1:nalt_sim]
    for(h in 1:nobs_sim){
      if(Utility[h, 1] > Utility[h, 2]){
        Prob[h] = Prob[h] + 1/n_sim
      }
    }
  }
  Y[, t] = Y[, t] + Err[, nalt_sim + 1]
  
  for(h in 1:nobs_sim){
    temp = runif(1, min=0, max=1) #general a random  probability
    if(temp <= Prob[h]){
      choice[h, t] = 0
    } else {choice[h, t] = 1
    }
  }
}

#save discrete dependent variable - "choice"
choice_name = paste("t", seq(1, spec$nTime, 1), sep="")
choice = rbind(choice_name, choice)
write.table(choice, file = "choice.txt", row.names=FALSE, col.names=FALSE)
#save continuous dependent variable - "miles"
miles_name = paste("vmt", seq(1, spec$nTime, 1), sep="")
miles = rbind(miles_name, Y)
write.table(miles, file = "miles.txt", row.names=FALSE, col.names=FALSE)
#mile = read.table(paste("miles", "txt", sep = "."), sep = " ", head = T)
