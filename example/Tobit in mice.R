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


mice.impute.RUTR <- function(y, ry, type, x, wy = NULL, upper_bound, type_c,...){ 
  
  ### Choose upper bound for y
  '%ni%' <- Negate('%in%')
  upper_bound_y <- upper_bound[which(type_c %ni% attributes(type)$name)[1]]
  
  ### Data for tobit
  y1 = y
  y1[wy] = upper_bound_y - 1 
  dt = data.frame(y1, x)
  
  ### Tobit regression, Mu & Sigma
  utobit <- tobit(y1 ~ ., data = dt, left = upper_bound_y, dist = "gaussian")
  mu <- predict(utobit) 
  mu_u <- mu[wy] 
  sigma_u <- utobit$scale 
  
  ### Impute LOD
  lod = length(which(wy*1 == 1))
  imp_u = c()
  for(i in 1:lod){
    imp_u[i] = rtnorm(1, mean = mu_u[i],sd = sigma_u, upper = upper_bound_y)}
  
  return(imp_u)
}

