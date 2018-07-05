# This file contains tools for working with Type S/M errors under simple
# hypothosis test scenarios.
# In future, as more complex designs such as ANOVA are added, they should each
# be given their own .R file.

# Function originally provided in Gelman and Carlin (2014), reused with
# permission under the MIT license
retrodesign <- function(A, s, alpha=.05, df=Inf, n.sims=1000000){
  z <- qt(1-alpha/2, df)
  p.hi <- 1 - pt(z-A/s, df)
  p.lo <- pt(-z-A/s, df)
  power <- p.hi + p.lo
  typeS <- p.lo/power
  estimate <- A + s*rt(n.sims,df)
  significant <- abs(estimate) > s*z
  exaggeration <- mean(abs(estimate)[significant])/A
return(list(power=power, typeS=typeS, exaggeration=exaggeration))
}

