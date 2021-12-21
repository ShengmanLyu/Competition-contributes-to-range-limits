#########################################################################
# IPM parameters and functions
# Spring 2020
# update 25.11.2020
#-----------------------------

#------------------------------------------
# IPM functions

# vital rate parameters
par <- c("surv.int", "surv.slope", 
         "growth.int", "growth.slope", "growth.sd", 
         "flowering.int", "flowering.slope", 
         "fecundity.int.linear", "fecundity.slope.linear",
         "germination.prob",
         "establishment.funiche", 
         "establishment.comp",
         "seedling.size.mean",
         "seedling.size.sd",
         "L", "U")

length(par)

# 1. probability of surviving (logistic)
survival.z=function(z,params) {
  # linear predictor
  u <- params$surv.int+params$surv.slope*z
  # inverse logit transformation
  survival.p <- exp(u)/(1+exp(u))
  survival.p[u>700] <- 1 # exp(710) gives Inf values
  # return the probability of survival
  return(survival.p)
}

# 2. growth function (gaussian)
growth.z1z <- function(z1,z, params) {
  # mean size time t+1, z1
  mu <- params$growth.int+params$growth.slope*z
  # sd about mean
  sig <- params$growth.sd
  # probability density of new size z1 for current size = z
  pd <- dnorm(z1,mean=mu,sd=sig)
  # return probability density
  return(pd)
}

# 2. growth function (gaussian)
# to predict size
growth.z <- function(z, params) {
  z1 <- params$growth.int+params$growth.slope*z
  return(z1)
}


# 3 probaility of flowering (logistic)
flowering.z <- function(z, params) {
  # linear predictor
  u <- params$flowering.int+params$flowering.slope*z
  # probability of flowering, inverse logit
  flowering.p <- exp(u)/(1+exp(u))
  flowering.p[u>700] <- 1 # exp(710) gives Inf values
  # return probability of flowering
  return(flowering.p) 
}

# 4 seed production (poisson)
seed.z.poisson <- function(z,params) {
  # linear predictor
  u = params$fecundity.int+params$fecundity.slope*z
  # seed prodution, exponential
  seeds <- exp(u) 
  # nagetative seeds to 0
  #u[u<0] = 0 
  #seeds <- u 
  # return seed production
  return(seeds)
}

# 4.1 seed production (gaussian)
seed.z.gaussian <- function(z,params) {
  # linear predictor
  u = params$fecundity.int.linear+params$fecundity.slope.linear*z
  # nagetative seeds to 0
  u[u<0] = 0 
  seeds <- u 
  # return seed production
  return(seeds)
}

# 5 probability of seed germination (constant)
germination <- function(params) {
  # get the probability of seedling establishment
  germination.p <- params$germination.prob
  # return
  return(germination.p)
}

# 6.1 probability of seedling establishment in the absence of competition (logistic)
establishment.fun <- function(params) {
  # establishment using FuNiche data
  u.funiche <- params$establishment.funiche
  establishment.prob.funiche <- exp(u.funiche)/(1+exp(u.funiche))
  # return
  return(establishment.prob.funiche)
}

# 6.2 probability of seedling establishment in the presence of competition (logistic)
establishment.com <- function(params) {
  # establishment with background species
  u.comp <- params$establishment.comp
  establishment.prob.comp <- exp(u.comp)/(1+exp(u.comp))
  
  # return
  return(establishment.prob.comp)
}

# 6.3 probability of seedling establishment in the presence of competition (logistic)
establishment.com0.5 <- function(params) {
  # establishment with background species
  u.comp <- params$establishment.comp
  establishment.prob.comp <- exp(u.comp)/(1+exp(u.comp))
  
  # return
  return(establishment.prob.comp*0.5)
}

# 7 seedling size (constant)
seedling.z1 <- function(z1, params) {
  # probability density of recruits
  rpd <- dnorm(z1, mean=params$seedling.size.mean, sd=params$seedling.size.sd)
  # return probability density
  return(rpd)
}

# G kernel: survival-growth
G.z1z <- function(z1,z,params) {
  # combine survival and growth
  sg <- survival.z(z, params) * growth.z1z(z1,z,params)
  # return
  return(sg)
}

# R kernel: reproduction
R.z1z=function(z1,z,params) {
  # calculate fecundity kernel
  f = flowering.z(z, params) * 
    seed.z.gaussian(z, params) * 
    germination(params) *
    establishment.fun(params) * establishment.com(params) *
    seedling.z1(z1, params)
  # return
  return(f)
}

# full kernel
full.z1z = function(z1,z,params) {
  # combine growth and reproduction kernel
  G.z1z(z1,z,params) + R.z1z(z1,z,params)
}

#------------------------------
# To make an IPM kernel
mk.kernel <- function(n, L, U, par) {
# n: number of meshpoints
# par: vital rates
# L, U: lower and upper limit of size
# d: density of competitor
# fun: full kernel 
# mesh points 
  h <- (U - L)/n
  meshpts <- L + ((1:n) - 1/2) * h
  G <- h * (outer(meshpts, meshpts, G.z1z, params = par))
  R <- h * (outer(meshpts, meshpts, R.z1z, params = par))
  K <- G+R
  return(list(meshpts = meshpts, G=G, R=R, K = K))
}

#------------------------------
# To plot kernels
library(plot.matrix)
plot.kernel <- function(x,y,k,... ) {
  # x,y: meshpoints
  # k: IPM kernel
  image(x, y, t(k), ...)
}

#-------------------------------
# To project an IPM
project.ipm <- function(n0, k, nstep){
  # k: value of make.kernel, with meshpoints and kernel
  # nstep: number of step to project

  # project only one step
  if(nstep==1) { 
    nt <- k %*% n0
    re <- nt
  }
  # project n steps
  else if(nstep >=2) {
    ntmat <- matrix(0, nrow=length(n0), ncol = nstep)
    ntmat[,1] <- n0; # each column is nt in one year
    for(k in 2:nstep) { ntmat[,k] <- k %*% ntmat[,k-1] }
    re <- ntmat
  }
  return(re)
}

#--------------------------------------------------
# calculate population growth rate (lambda) using kernel
lambda.k <- function(k) {
# k: IPM full kerbel
  l=Re(eigen(k, only.values = TRUE)$values[1]) 
  return(l)
}

m <- matrix(1:16, nrow=4)
eigen.m <- eigen(m, only.values = TRUE)
Re(eigen.m$values)

#--------------------------------------------------
# To perform matrix-level sensitivity, elasticity analysis on an IPM
elasticity.k <- function(k) {
  # k: full kernel of an IPM
  # right eigenvalue
  w.eigen=Re(eigen(K)$vectors[,1])
  # stable size distribution
  stable.dist=w.eigen/sum(w.eigen) 
  # left eigenvalue
  v.eigen=Re(eigen(t(K))$vectors[,1])
  # reproductive value when population is stable?
  repro.val=v.eigen/v.eigen[1]  
  
  # elasticity and sensitivity matrices
  v.dot.w=sum(stable.dist*repro.val)*h
  # sensitivity matrix
  sens=outer(repro.val,stable.dist)/v.dot.w
  # elasticity matrix
  elas=matrix(as.vector(sens)*as.vector(K)/lam,nrow=n)
  
  return(list(sensitivty=sens, elasticity = elas))
}

#--------------------------------------------------
# To perform vital rate-level sensitivity, elasticity analysis on an IPM
elasticity.vr <- function(params, lambda.true, vital.rate, delta=0.01, ...) {
  # params: paramters of vital rates in an IPM
  # lambda.true: actural population growth rate of the IPM
  # vital.rate: vital rates of interest (of which caclcate sensitivity and elastivity)
  # delta: the perturbation of vital rates
  
  nvr = length(vital.rate) # number of vital rates
  sens = numeric(nvr); names(sens) = vital.rate # vector to hold parameter sensitivities
  
  for(i in 1:nvr){ 
    m.par = params
    vr.i <- vital.rate[i]
    print(vr.i)
    
    # IPM down
    #m.par[vr.i] = m.par[vr.i] - delta
    #IPM.down = mk.kernel(par=m.par, L = as.numeric(m.par["L"]), U=as.numeric(m.par["U"]), ...)
    # IPM.down = mk_K_ceiling(m.par=m.par, L = m.par["L"], U=m.par["U"], ...) # IPM course
    #lambda.down = lambda.k(IPM.down$K)
    
    # IPM up
    m.par[vr.i]=m.par[vr.i] + 1*delta
    IPM.up = mk.kernel(par=m.par, L = as.numeric(m.par["L"]), U=as.numeric(m.par["U"]), ...)
    #IPM.up = mk_K_ceiling(m.par=m.par, L = m.par["L"], U=m.par["U"], ...) # IPM course
    if(sum(is.na(IPM.up$K))>0 | sum(is.infinite(IPM.up$K))>0) {warning("failed to calculare lambda"); print(vr.i); next}
    lambda.up = lambda.k(IPM.up$K)
    
    # calculate sensitivity
    si = (lambda.up-lambda.true)/(1*delta) 
    sens[i]=si 
  }   
  # calculate elasticity: do I have to ca
  elas = sens*abs(params[vital.rate])/lambda.true; names(elas) = vital.rate
  
  return(list(sensitivty = sens, elastivity = elas))
}

#--------------------------------------------------
# To perform parameter-level sensitivity, elasticity analysis on an IPM
elasticity.par <- function(params, lambda.true, vital.rate, delta=0.01, ...) {
  # params: paramters of vital rates in an IPM
  # lambda.true: actural population growth rate of the IPM
  # vital.rate: vital rates of interest (of which caclcate sensitivity and elastivity)
  # delta: the perturbation of vital rates
  
  nvr = length(vital.rate) # number of vital rates
  sens = numeric(nvr); names(sens) = vital.rate # vector to hold parameter sensitivities
  
  for(i in 1:nvr){ 
    m.par = params
    vr.i <- vital.rate[i]
    print(vr.i)
    
    # IPM down
    #m.par[vr.i] = m.par[vr.i] - delta
    #IPM.down = mk.kernel(par=m.par, L = as.numeric(m.par["L"]), U=as.numeric(m.par["U"]), ...)
    # IPM.down = mk_K_ceiling(m.par=m.par, L = m.par["L"], U=m.par["U"], ...) # IPM course
    #lambda.down = lambda.k(IPM.down$K)
    
    # IPM up
    m.par[vr.i]=m.par[vr.i] + 1*delta
    IPM.up = mk.kernel(n=3000, par=m.par, L = as.numeric(m.par["L"]), U=as.numeric(m.par["U"]), ...)
    #IPM.up = mk_K_ceiling(m.par=m.par, L = m.par["L"], U=m.par["U"], ...) # IPM course
    if(sum(is.na(IPM.up$K))>0 | sum(is.infinite(IPM.up$K))>0) {warning("failed to calculare lambda"); print(vr.i); next}
    lambda.up = lambda.k(IPM.up$K)
    
    # calculate sensitivity
    si = (lambda.up-lambda.true)/(1*delta) 
    sens[i]=si 
  }   
  # calculate elasticity: do I have to ca
  elas = sens*abs(params[vital.rate])/lambda.true; names(elas) = vital.rate
  
  return(list(sensitivty = sens, elastivity = elas))
}

elas = 1:10*1:10
names(elas) <- letters[1:10]

#--------------------------------------------------
# To perform LTRE analysis on IPMs
ltre <- function(vital.rate, params.ctl, lambda.ctl, params.trt, lambda.trt, params.ref, lambda.ref, sensitivity.ref=NA, ...) {
  # params.ctl: vital rates of control IPM (only one as reference level)
  # params.trt: vital rates of treatment IPMs (can be more than 1, then each row is a unique treatment)
  # params.ref: vital rates of reference IPM (can be calcuate as the mean of control and treatment IPMs)
  # lambda.ctl, lambda.trt, lambda.mean: lambdas of control, treatment and mean-IPM
  # vital.rates: vital rates of interest (of which calculate sensitivity and elastivity)
  
  # sensitivity of mean-kernel
  if(sum(is.na(sensitivity.ref[vital.rate])) != 0) {
    print("calculate sensitivity of the reference IPM")
    elas.vr <- elasticity.par(params = params.ref, lambda.true = lambda.ref, vital.rate = vital.rate, ...)
    sens.mean <- elas.vr$sensitivity
  }
  else sens.mean <- sensitivity.ref
  
  # number of vital rate parameters
  npar <- length(vital.rate)
  
  # number of treatments
  ntrt <- nrow(params.ctl)
  
  # matrix to hold LTRE for each treatment IPM
  LTRE <- matrix(NA, nrow=ntrt, ncol=npar + 2)
  for(i in 1:ntrt) {
    # parameter differences between control and treatment
    par.dif <- params.trt[vital.rate] - params.ctl[vital.rate]
    # LTRE of parameters
    ltre.i <- sens.mean * par.dif
    # combine LTRE, lambda difference, and LTRE sum
    LTRE[i,] <- c(as.numeric(ltre.i[1,]), lambda.trt[i] - lambda.ctl, sum(ltre.i))
  }
  colnames(LTRE) <- c(vital.rate, "dif.lambda", "sum.ltre")
  
  return(list(LTRE=LTRE, sensitivity.reference=sens.mean))
}

#--------------------------------------------------
# To standardise relaitve contribution of vital rates
# This standardisation makes relative contribution compaable across species
ltre.standardise <- function(x) {
  # x: a vector of relative contribution, including positive and negtaive values, of vital rates to be standardised
  
  # Check if NA in x
  if(sum(is.na(x)) >0) warning("NAs in the relative contribution!")
  
  # standardisation
  x.positive <- abs(x)
  x.sum <- sum(x.positive)
  x.divided <- x.positive/x.sum
  x.standard <- ifelse(x< 0, (x.divided*-1), x.divided)
  
  return(x.standard)
}

# test
ltre.standardise(c(-1,1,2,-2))

#--------------------------------------------------
# calculate sensitivity of species using intrinsic and invasion growth rates
sensitivity.igr <- function(pgr.intrinsic, pgr.invasion) {
  # pgr.intrinsic: intrinsic growth rate of the focal species
  # pgr.invasion: invasion growth rate of the focal species
  
  # log
  pgr.intrinsic.log <- log(pgr.intrinsic)
  pgr.invasion.log <- log(pgr.invasion)
  
  # calculate sensitivity
  if(pgr.intrinsic.log > 0) {
    sens <- 1- pgr.invasion.log/pgr.intrinsic.log
  }
  else if(pgr.intrinsic.log < 0) {
    sens <- pgr.invasion.log/pgr.intrinsic.log - 1
  }
  else if(pgr.intrinsic.log == 0) {
    sens <- 1- pgr.invasion.log/pgr.intrinsic.log
  }
  # output
  return(sens)
}

# when intrinic > 0
intrinsic = exp(0.5)
invasion = exp(0.5 + seq(-0.8, 0.8, length.out = 10))
sensitivity.igr(pgr.intrinsic = intrinsic, pgr.invasion = invasion)

# when intrinic < 0
intrinsic = exp(-0.5)
invasion = exp(-0.5 + seq(-0.8, 0.8, length.out = 10))
sensitivity.igr(pgr.intrinsic = intrinsic, pgr.invasion = invasion)

#--------------------------------------------------
# calculate sensitivity and coexistence outcome of 2 species using invasion growth rate
coex.igr <- function(sps1,sps2, igr10, igr20, igr12, igr21) {
  # igr10, igr20: intrinsic growth rate in the absence of competitors
  # igr11: invasion growth rate of species 1 invading species 2
  # igr12: invasion growth rates of species 2 invaing species 1
  
  # calculate sensitivity
  s1 <- sensitivity.igr(igr10, igr12); names(s1) <- sps1
  s2 <- sensitivity.igr(igr20, igr21); names(s2) <- sps2
  ss <- c(s1, s2)
  
  # log
  igr10 <- log(igr10)
  igr20 <- log(igr20)
  igr12 <- log(igr12)
  igr21 <- log(igr21)
  
  # coexsitence outcome independent of the intrinsic growth rates
  if(igr12 > 0 & igr21 >0) {
    outcome <- "Coexistence"
    winner <- "sps1_sps2"
  }
  else if(igr12 > 0 & igr21 < 0) { 
    outcome <- "Competitive exclusion"
    winner <- "sps1"
  }
  else if(igr12 < 0 & igr21 > 0) { 
    outcome <- "Competitive exclusion"
    winner <- "sps2"
  }
  else if(igr12 < 0 & igr21 < 0) { 
    outcome <- "Priority effect"
    winner <- "none"
  }
  
  # superior
  if(s1 > s2) superior <- "sps2"
  else superior <- "sps1"
  
  return(list(sensitivity = ss, outcome = outcome, winner = winner, superior = superior)) 
}

#--------------------------------------------------
# calculate ND and FD and coexistence outcome of 2 species
coex.ndfd <- function(s12, s21) {
  # s12: sensitivity of species 1 invading species 2
  s <- c(s12, s21)
  
  # ND
  nd <- 1 - sqrt(s12*s21)
  names(nd) <- "ND"
  
  # FD
  # afd <- exp( sqrt( (log(s12/gm) + log(s21/gm))/2 ) )
  # afd <- exp(sqrt(mean(log(s)**2) - (mean(log(s)))**2))
  rfd <- sqrt(s21/s12)
  names(rfd)  <- "RFD"
  
  # superior
  if(rfd > 1) superior <- "sps1"
  else superior <- "sps2"
  
  # Coexistence outcomes
  # make sure FD always > 1
  if(rfd < 1) rfd2 <- 1/rfd
  else rfd2 <- rfd
  
  # coexistence metrics
  coex.metric <- 1/ (rfd2*(1-nd))
  names(coex.metric) <- "Coexistence metric"
  
  if(rfd2 <= 1/(1-nd)) { 
    outcome <- "Coexistence"
    winner <- "sp1_sps2"
    }
  else if(rfd > 1) {
    outcome <- "Competitive exclusion"
    winner <- "sps1"
    }
  else if(rfd < 1) {
    outcome <- "Competitive exclusion"
    winner <- "sps2"
  }
  return(list(ndfd = c(nd, rfd, coex.metric), outcome = outcome, winner = winner, superior = superior) )
}

# test
coex.ndfd(0.6,2)
coex.ndfd(2, 0.6)

#-----------------------------------------------------------
# Plot ND and RFD (log-scale) frame indicating coexistence outcomes
coex.frame <- function(x1=-0.5,x2=0.9) {
  # x1, x2 lower and upper limits of ND, in which x1 < -0.001
  x <- c(seq(x1,-0.001, length.out = 100), seq(0.001,x2,length.out = 100))
  y <- log(1/(1-x))
  p <- ggplot() + 
    geom_line(data=data.frame(x=x, y=y), aes(x=x, y=y), inherit.aes = FALSE) +
    geom_line(data=data.frame(x=x, y=-y), aes(x=x, y=y), inherit.aes = FALSE) +
    geom_line(data=data.frame(x=x, y=rep(0,length(x))), aes(x=x, y=y), inherit.aes = FALSE, linetype="dashed")
  return(p)
}
coex.frame()

#-----------------------------------------------------------
# Plot ND and RFD (log-scale) frame indicating coexistence outcomes
# add frame
# x1, x2 lower and upper limits of ND, in which x1 < -0.001
frame.coex = function(x1 = -0.01, x2= 0.91) {
  x <- c(seq(x1,-0.001, length.out = 100), seq(0.001,x2,length.out = 100))
  y <- log(1/(1-x))
  return(data.frame(x=x,y=y))
}
frame.coex()

# add polygan
# polygan of coexistence area
polygan.coex = function(x1 = 0.0001, x2= 0.91) {
  xx1 = seq(x1,x2, length.out = 100)
  xx2 = rep(x2,100)
  xx3 = seq(x2,x1,length.out = 100)
  xx <- c(xx1,xx2,xx3)
  yy <- c(log(1/(1-xx1)), log(1/(1-xx2)),-log(1/(1-xx3)))
  return(data.frame(x=xx,y=yy))
}
ggplot(polygan.coex(), aes(x=x,y=y)) +
  geom_polygon()

polygan.prio = function(x1 = -3, x2= -0.001) {
  xx1 = seq(x1,x2, length.out = 100)
  xx2 = rep(-0.001,100)
  xx3 = seq(x2,x1,length.out = 100)
  xx <- c(xx1,xx2,xx3)
  yy <- c(log(1/(1-xx1)), log(1/(1-xx2)),-log(1/(1-xx3)))
  return(data.frame(x=xx,y=yy))
}
ggplot(polygan.prio(), aes(x=x,y=y)) +
  geom_polygon()


frame.coex = function(x1 = -0.01, x2= 0.91) {
  # x1, x2 lower and upper limits of ND, in which x1 < -0.001
  x <- c(seq(x1,-0.001, length.out = 100), seq(0.001,x2,length.out = 100))
  y <- log(1/(1-x))
  return(data.frame(x=x,y=y))
}
frame.coex()


#-----------------------------------------------------------
# To extrac statics from lmer
library(car)
coef.lmer <- function(model, ...) {
  fit <- Anova(model, ...)
  return(as.data.frame(fit))
}

#-----------------------------------------------------------
#mean_sd <- function(x) {c(mean(x), sd(x))}

#-----------------------------------------------------------
#mean_se <- function(x) {c(mean(x), sqrt(var(x)/length(x)))}

#-----------------------------------------------------------
# confidence interval using linear model
# same as: Hmisc::mean_cl_normal
mean_ci <- function(x) {
  fit <- lm(x~1)
  c(mean(x), confint(fit)[1, ])
  }

#-----------------------------------------------------------
# confidence interval using SE
mean_ci_se <- function(x) {
  u = mean(x)
  se = sqrt(var(x)/length(x))
  data.frame(y=u, ymin=u-1.96*se, ymax=u+1.96*se)
}

#-----------------------------------------------------------
# confidence interval using SD
mean_ci_sd <- function(x, ...) {
  u = mean(x, ...)
  sd = sqrt(var(x, ...))
  data.frame(y=u, ymin=u-1.96*sd, ymax=u+1.96*sd)
}

#-----------------------------------------------------------
# confidence interval (Hanna's)
mean_ci_hanna <- function(x, ...) {
  u = mean(x, ...)
  x.min = quantile(x, probs = 0.025, ...)
  x.max = quantile(x, probs = 0.975, ...)
  data.frame(y=u, ymin=x.min, ymax=x.max) 
}