# install package "dcmodels_1.0.01" if you have not installed it before
#install.packages("dcmodels_1.0.01.tar.gz", repos=NULL, type="source")

library(dcmodels)
rm(list=ls())
setwd("C:/Users/Yan Liu/Documents/UMD research/recursive probit/2dynamic binary D-C model - series covariances")

source("dcmodels/R/probitUtils.R")
source("dcmodels/R/models.R")
source("dcmodels/R/dynUtils.R")
source("dcmodels/R/dyn.R")
source("dcmodels/R/dyn_dc.R")
source("dcmodels/R/dynData.R")

spec = list(
  D = "simulate_data_scale/dyn",
  Global = "simulate_data_scale/global.txt",
  miles = "simulate_data_scale/miles.txt",
  choices = "simulate_data_scale/choice.txt",
  
  generic = list(c("time.0", "time.1")),
  specific = list(c(c("cost.0")), c("ASC", "cost.1", "income")),
  reg = c("attr1", "attr2"),
  
  nvar_altspec = 2,
  modifyD = function(D,Global,t){ 
    Dyan = cbind(D[[t]], Global)
    return(Dyan)
  },
  SD = "hessian",
  ASC = FALSE, # specify this now
  nTime = 17,  #actual time, not total time
  nLook = 3,  #look ahead period
  transition = rbind(c(1, 1), c(1, 1)),
  stopAlt = c(), # HERE, if you have more than one alternative 
  # that halts the decision process, put that in 
  # a vector (like stopAlt = c(1,2,10) if 
  # alt 1, 2 and 10 are stopping alternatives)
  reltol = 1e-5,
  delta = 1e-1,
  method = "pmvnorm",
  SD = "hessian",    #bootstrap OR none
  nboot = 100,
  verbose = TRUE
)

spec = checkFillDynSpec(spec)
modelFns = dyn_dc
D = NULL

Sys.time()
Mdyn = model_Yan(modelFns, spec, D)
Sys.time()
Mdyn

#calculate LL at 0
spec = checkFillDynSpec(spec)
Args = dyn_dc$computeArgs(spec,NULL)
#zero = rep(0, length(dyn_dc$computeStart(spec,NULL)))
b = rep(0, 7)
for(t in 1:Args$spec$nTime){
  b = c(b, seq(0, 1, 1))
}
cat("LL at 0: " ,sum(log(dyn_dc$LLVec(b,Args))),"\n")

