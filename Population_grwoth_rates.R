#########################################################################
# Calculate population grwoth rates
# Author: Shengman Lyu (shengman.lyu@usys.ethz.ch)
# Date: 07.01.2022
#########################################################################
rm(list=ls())

# load IPM functions
source("Functions_for_IPMs.R")

# read bootstrapped vital rates 
p <- read.csv("vr_bootstrap_20210315.csv")
p$pgr <- as.numeric(NA)

for(a in 1:nrow(p)) {
  # a = 1
  da <- p[a,]
  
  # skip
  if(is.na(da$pair)) next
  
  # calculate population growth rates
  ipm <- mk.kernel(n=3000, L=as.numeric(da$L), U=as.numeric(da$U), par=da) 
  if(sum(is.na(ipm$K))>0) next
  p[a,"pgr"] <- lambda.k(ipm$K)
  print(a)
  
  # plot
  # plot.kernel(x=ipm$meshpts, y=ipm$meshpts, k=ipm$K)
}

# Population grwoth rates
p$pgr