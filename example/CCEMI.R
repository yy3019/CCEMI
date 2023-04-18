library(bootImpute)
library(tidyverse)
library(rstan)
library(mice)
library(investr)
library(MASS)
library(dplyr)
library(knitr)
library(matrixStats)
library(foreach)
library(AER)
library(MCMCglmm)


#### M = Imputation numbers, normally we set it as 2.
#### maxit = 20
#### method, details see 'mice' package
#### predM, imputation matrix for missing data

impOnce <- function(inputData,M) {
  miceImps <- mice::mice(inputData, m=M, maxit = 20, print = F, method = "norm", predictorMatrix = predM) 
  for (i in 1:M) {
    imps[[i]] <- mice::complete(miceImps,i)
  }
  imps
}

#### imp, the missing data need to be imupted
#### M = bn, imputation numbers, normally we set it as 2.
#### ss, set seed for reproducibility

impst <- bootImpute(imp, impOnce, nBoot=200, nImp = bn,seed=ss, M = bn) 


analyseImp <- function(inputData) {
  coef(lm(y~ x1 + x2 + x3 + z,data=inputData)) ### set model for analysis
}


resultt = bootImputeAnalyse(impst, analyseImp)


a_bmt =  cbind(resultt$ests, resultt$ci)