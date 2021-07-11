

# testing whether x1 = x2
diff.ttest <- function(x1, x2, lincom=c(1,-1), model){
  diffest <- lincom[1] * coef(summary(model))[x1, "Value"] + lincom[2] * coef(summary(model))[x2, "Value"]
  vardiff <- ( (lincom[1] * coef(summary(model))[x1, "Std.Error"])^2 + 
                 (lincom[2] * coef(summary(model))[x2, "Std.Error"])^2 ) + (2 * lincom[1] * lincom[2] * vcov(model)[x1, x2]) 
  # variance of x1 + variance of x2 - 2*covariance of x1 and x2
  diffse <- sqrt(vardiff)
  tdiff <- (diffest)/(diffse)
  DF <- coef(summary(model))[x1, "DF"]
  ptdiff <- 2*(1- pt(abs(tdiff), DF, lower.tail=T) )
  upr <- diffest + qt(.975, df = DF)*diffse
  lwr <- diffest + qt(.025, df = DF)*diffse
  return(list(
    p=round(ptdiff, digits = 4), 
    lwr=round(lwr, digits = 4), 
    upr=round(upr, digits = 4)
  ))
}


# testing whether x1=x2=0
diff.waldtest <- function(x1, x2, model){
  
  var1 <- (coef(summary(model))[x1, "Std.Error"])^2
  var2 <- (coef(summary(model))[x2, "Std.Error"])^2
  cov12 <- vcov(model)[x1, x2]
  est1 <- coef(summary(model))[x1, "Value"]
  est2 <- coef(summary(model))[x2, "Value"]
  delta <- var1 * var2 - (cov12)^2
  wald.stat <- (1/delta) * ((est1^2 * var2) - (2 * est1 * est2 * cov12) + (est2^2 * var1) )
  DF <- coef(summary(model))[x1, "DF"]
  
  pchi <- 1- pchisq(wald.stat, 2, lower.tail=T) 
  sim.upr.x1 <- est1 + qt(.9875, df = DF) * sqrt(var1)
  sim.lwr.x1 <- est1 + qt(.0125, df = DF) * sqrt(var1)
  sim.upr.x2 <- est2 + qt(.9875, df = DF) * sqrt(var2)
  sim.lwr.x2 <- est2 + qt(.0125, df = DF) * sqrt(var2)
  
  return(list(pchi=round(pchi, digits = 4), 
              sim.upr.x1=round(sim.upr.x1, digits = 4), 
              sim.lwr.x1=round(sim.lwr.x1, digits = 4),
              sim.upr.x2=round(sim.upr.x2, digits = 4), 
              sim.lwr.x2=round(sim.lwr.x2, digits = 4)
  ))
}