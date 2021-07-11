
library(tidyverse)
library(data.table)


gData = function(nCluster, nGroup, nTransit, m, beta0, beta4, rho34, rho24, rho13, t=0, sigma, sigma2=NA, secularTrend=0, randomTrt=0){
  
  nPeriod = (nCluster + 1 + nTransit)
  
  if(secularTrend == 2 & length(t) !=  nPeriod) stop('the length of the time effect vector should match the number of total time points')
  if(secularTrend != 2 & length(t) > 1) stop('the length of the time effect vector should be 1')
  if(secularTrend == 0 & length(t) == 1  & t != 0) stop('the time effect can only be zero if no secular trend is selected')
  if(randomTrt == 1 & is.na(sigma2) == T) stop('please specify the variance of the random treatment effect')
  if(randomTrt == 0 & is.na(sigma2) == F) {
    sigma2 = NA
    warning('When randomTrt option is FALSE (default), the input of sigma2 is ignored.')
  }
  
  beta3 = beta4 * rho34
  beta2 = beta4 * rho24
  beta1 = beta3 * rho13

  cluster = rep(1:nCluster, each = nPeriod*m*nGroup)
  group = rep(1:(nCluster*nGroup), each=nPeriod*m) # cluster id
  period = rep(1:nPeriod, time=nCluster*nGroup, each = m) # period
  time = period - 1
  subject = seq(1:(nCluster*nPeriod*m*nGroup)) # subject id
  
  if((nCluster %% 2) == 0) {
    seq = rep(c(0,1), time=nCluster/2, each=nPeriod*m*nGroup)  
  } else {
    seq = c(rep(c(0,1), time=floor(nCluster/2), each=nPeriod*m*nGroup), rep(0, nPeriod*m*nGroup))
  }
  
  trt = vector()
  for (j in 1:nCluster){
    trt_j = rep(rep(c(rep(0,j), rep((2-(j%%2)), nTransit), rep(3, nPeriod-j-nTransit)), each=m), time=nGroup)
    # C:trt=0; C+P:trt=1; C+IPE:trt=2; C+P+IPE:trt=3
    # if j is odd, then the "transition" is C+P and x=1; otherwise it's C+IPE and x=2
    trt = append(trt, trt_j)
  }
  
  inttime = vector() # the starting time for each intervention
  for (j in 1:nCluster){
    inttime_j = rep(rep(c(rep(1,j), 1+j, rep(2+j, nPeriod-j-1)), each=m), time=nGroup)
    inttime = append(inttime, inttime_j)
  }
  
  intduration = period-inttime
  
  # generate 
  
  data = data.frame(cluster, group, period, time, trt, subject, seq, inttime, intduration) %>%
    mutate(trt_ind = ifelse(trt==1, "trt1", 
                            ifelse(trt==2, "trt2", 
                                   ifelse(trt==3 & seq==0, "trt3", 
                                          ifelse(trt==3 & seq==1, "trt4", "control")))),
           trt_ind = factor(trt_ind, levels = c("control", "trt1", "trt2", "trt3", "trt4")))
  
  
  data$u1 = rep(rnorm(nCluster, 0, sigma), each=nPeriod*m*nGroup) # for different icc
  
  
  if (randomTrt==1){
    
    u2 = rep(c(rep(rnorm(1, 0, sigma2), 1*m), rep(rnorm(1, 0, sigma2), nTransit*m), rep(rnorm(1, 0, sigma2), (nPeriod-1-nTransit)*m)), time=nGroup)
    for (j in 2:nCluster){
      u2_new = rep(c(rep(rnorm(1, 0, sigma2), j*m), rep(rnorm(1, 0, sigma2), nTransit*m), rep(rnorm(1, 0, sigma2), (nPeriod-j-nTransit)*m)), time=nGroup)
      u2 = append(u2, u2_new)
    }
    data$u2 = u2
    
  }
  
  
  if (secularTrend==2){
    
    # nonlinear time effect
    data$time_effect = NULL # time_effect is the categorical time effect
    for (j in 1:nPeriod){
      data$time_effect[data$period == j] <- t[j]
    }

  }
  
  
  # linear time + nonrandom trt effect
  # linear time + random trt effect
  # nonlinear time + nonrandom trt effect
  # nonlinear time + random trt effect
  
  
  if (secularTrend!=2 & randomTrt==0){
    log_odds = with(data, beta0 + beta1*(trt_ind == "trt1") + beta2*(trt_ind == "trt2") + beta3*(trt_ind == "trt3") +
                      beta4*(trt_ind == "trt4") + u1 + t*time)
  } else if (secularTrend!=2 & randomTrt==1){
    log_odds = with(data, beta0 + beta1*(trt_ind == "trt1") + beta2*(trt_ind == "trt2") + beta3*(trt_ind == "trt3") +
                      beta4*(trt_ind == "trt4") + u1 + u2 + t*time)
  } else if (secularTrend==2 & randomTrt==0){
    log_odds = with(data, beta0 + beta1*(trt_ind == "trt1") + beta2*(trt_ind == "trt2") + beta3*(trt_ind == "trt3") +
                      beta4*(trt_ind == "trt4") + u1 + time_effect)
  } else if (secularTrend==2 & randomTrt==1) {
    log_odds = with(data, beta0 + beta1*(trt_ind == "trt1") + beta2*(trt_ind == "trt2") + beta3*(trt_ind == "trt3") +
                      beta4*(trt_ind == "trt4") + u1 + u2 + time_effect)
  }

  
  prop = plogis(log_odds)
  data$y = rbinom(n = nCluster*nPeriod*m*nGroup, size = 1, prob = prop)
  
  data_small = data.table(data)
  data_small[, ratio := mean(y), by = c("cluster", "time", "trt_ind")]
  data_small[, n := length(y), by = c("cluster", "time",  "trt_ind")]
  data_small <- unique(data_small[, c("cluster", "time", "trt_ind", "ratio", "n")])

  return(data_small)
  


}

# # # example
# data = gData(nCluster=1, nGroup=1, nTransit=1, m=1,
#              beta0=0, beta4=0.8, rho34=0.5, rho24=0.5, rho13=0.5,
#              t=0, sigma=0.05, sigma2=0.05, secularTrend = 0, randomTrt = 0)





