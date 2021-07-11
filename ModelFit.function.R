
library(tidyverse)
library(data.table)
library(MASS)

ModelFit <- function(glmFun, 
                     nCluster, nGroup, nTransit, m, 
                     beta0, beta4, rho34, rho24, rho13, t, 
                     sigma, sigma2=NA, 
                     secularTrend=0, randomTrt=0,
                     x1=NULL, x2=NULL, lincom=c(1, -1)) {
  
  
  result_vec = c()
  
  data_small = gData(nCluster, nGroup, nTransit, m, 
                     beta0, beta4, rho34, rho24, rho13, t, 
                     sigma, sigma2, 
                     secularTrend, randomTrt)
  
  glmfit <- glmFun(data = data_small)
  
  if (! is.character(glmfit)) {
    
    coef_trt1 <- coef(summary(glmfit))["trt_indtrt1", "Value"]
    lci_trt1 <- coef(summary(glmfit))["trt_indtrt1", "Value"] - 1.96 * coef(summary(glmfit))["trt_indtrt1", "Std.Error"]
    uci_trt1 <- coef(summary(glmfit))["trt_indtrt1", "Value"] + 1.96 * coef(summary(glmfit))["trt_indtrt1", "Std.Error"]
    pvalue_trt1 <- coef(summary(glmfit))["trt_indtrt1", "p-value"]
    
    coef_trt2 <- coef(summary(glmfit))["trt_indtrt2", "Value"]
    lci_trt2 <- coef(summary(glmfit))["trt_indtrt2", "Value"] - 1.96 * coef(summary(glmfit))["trt_indtrt2", "Std.Error"]
    uci_trt2 <- coef(summary(glmfit))["trt_indtrt2", "Value"] + 1.96 * coef(summary(glmfit))["trt_indtrt2", "Std.Error"]
    pvalue_trt2 <- coef(summary(glmfit))["trt_indtrt2", "p-value"]
    
    coef_trt3 <- coef(summary(glmfit))["trt_indtrt3", "Value"]
    lci_trt3 <- coef(summary(glmfit))["trt_indtrt3", "Value"] - 1.96 * coef(summary(glmfit))["trt_indtrt3", "Std.Error"]
    uci_trt3 <- coef(summary(glmfit))["trt_indtrt3", "Value"] + 1.96 * coef(summary(glmfit))["trt_indtrt3", "Std.Error"]
    pvalue_trt3 <- coef(summary(glmfit))["trt_indtrt3", "p-value"]
    
    coef_trt4 <- coef(summary(glmfit))["trt_indtrt4", "Value"]
    lci_trt4 <- coef(summary(glmfit))["trt_indtrt4", "Value"] - 1.96 * coef(summary(glmfit))["trt_indtrt4", "Std.Error"]
    uci_trt4 <- coef(summary(glmfit))["trt_indtrt4", "Value"] + 1.96 * coef(summary(glmfit))["trt_indtrt4", "Std.Error"]
    pvalue_trt4 <- coef(summary(glmfit))["trt_indtrt4", "p-value"]
    
    if ( "time" %in% rownames(coef(summary(glmfit))) ){
      
      coef_time <- coef(summary(glmfit))["time", "Value"]
      lci_time <- coef(summary(glmfit))["time", "Value"] - 1.96 * coef(summary(glmfit))["time", "Std.Error"]
      uci_time <- coef(summary(glmfit))["time", "Value"] + 1.96 * coef(summary(glmfit))["time", "Std.Error"]
      pvalue_time <- coef(summary(glmfit))["time", "p-value"]
      
    } else {
      
      coef_time <- NA
      lci_time <- NA
      uci_time <- NA
      pvalue_time <- NA
      
    }
    
    
    test1Result <- diff.waldtest(x1=x1, x2=x2, model=glmfit)
    test1_lci_x1 <- test1Result$sim.lwr.x1
    test1_uci_x1 <- test1Result$sim.upr.x1
    test1_lci_x2 <- test1Result$sim.lwr.x2
    test1_uci_x2 <- test1Result$sim.upr.x2
    test1_pvalue <- test1Result$pchi
    
    test2Result <- diff.ttest(x1=x1, x2=x2, lincom=lincom, model=glmfit)
    test2_lci <- test2Result$lwr
    test2_uci <- test2Result$upr
    test2_pvalue <- test2Result$p
    
    
    
    nPeriod = nCluster + 1 + nTransit
    N = m * (nGroup*nCluster*nPeriod)
    
    beta3 = beta4 * rho34
    beta1 = beta3 * rho13
    beta2 = beta4 * rho24
    
    result_vec <- c(N,
                    m,
                    nCluster,
                    nGroup,
                    nTransit,
                    beta0,
                    beta1,
                    beta2,
                    beta3,
                    beta4,
                    rho13,
                    rho24,
                    rho34,
                    sigma,
                    sigma2,
                    coef_trt1,
                    lci_trt1,
                    uci_trt1,
                    pvalue_trt1,
                    coef_trt2,
                    lci_trt2,
                    uci_trt2,
                    pvalue_trt2,
                    coef_trt3,
                    lci_trt3,
                    uci_trt3,
                    pvalue_trt3,
                    coef_trt4,
                    lci_trt4,
                    uci_trt4,
                    pvalue_trt4,
                    coef_time,
                    lci_time,
                    uci_time,
                    pvalue_time,
                    test1_lci_x1,
                    test1_uci_x1,
                    test1_lci_x2,
                    test1_uci_x2,
                    test1_pvalue,
                    test2_lci,
                    test2_uci,
                    test2_pvalue)
    
  }
    return(result_vec)
}