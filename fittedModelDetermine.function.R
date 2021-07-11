
library(tidyverse)
library(data.table)
library(MASS)

fittedModelDetermine <- function(fittedModel){
  
  fit.model1 <- function(data){
    tryCatch(glmmPQL(data = data, ratio ~ trt_ind, random = ~ 1 | cluster, 
                     family = binomial, weights = n), 
             warning = function(w) { "warning" }, 
             error = function(e) { "error" })
  }
  
  fit.model2 <- function(data){
    tryCatch(glmmPQL(data = data, ratio ~ trt_ind + time, random = ~ 1 | cluster, 
                     family = binomial, weights = n), 
             warning = function(w) { "warning" }, 
             error = function(e) { "error" })
  }
  
  fit.model3 <- function(data){
    tryCatch(glmmPQL(data = data, ratio ~ trt_ind + factor(time), random = ~ 1 | cluster, 
                     family = binomial, weights = n), 
             warning = function(w) { "warning" }, 
             error = function(e) { "error" })
  }
  
  fit.model4 <- function(data){
    tryCatch(glmmPQL(data = data, ratio ~ trt_ind, random = list(cluster = ~1, trt_ind = ~1), 
                     family = binomial, weights = n), 
             warning = function(w) { "warning" }, 
             error = function(e) { "error" })
  }
  
  fit.model5 <- function(data){
    tryCatch(glmmPQL(data = data, ratio ~ trt_ind + time, random = list(cluster = ~1, trt_ind = ~1), 
                     family = binomial, weights = n), 
             warning = function(w) { "warning" }, 
             error = function(e) { "error" })
  }
  
  fit.model6 <- function(data){
    tryCatch(glmmPQL(data = data, ratio ~ trt_ind + factor(time), random = list(cluster = ~1, trt_ind = ~1), 
                     family = binomial, weights = n), 
             warning = function(w) { "warning" }, 
             error = function(e) { "error" })
  }
  
  
  if (fittedModel == 1){
    fit.model = fit.model1
  } else if (fittedModel == 2){
    fit.model = fit.model2
  } else if (fittedModel == 3){
    fit.model = fit.model3
  } else if (fittedModel == 4){
    fit.model = fit.model4
  } else if (fittedModel == 5){
    fit.model = fit.model5
  } else {fit.model = fit.model6}
  
  return(fit.model)
  
}

