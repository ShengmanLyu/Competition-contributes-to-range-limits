#*****************************************************************************
# Functions ----
# Shengman Lyu
# Update 2022.01.13
# See more detail in the manuscript "Competition contributes to both warm and cold range edges"
# E-mail:shengman.lyu@usys.ethz.ch
#*****************************************************************************

# The codes provided here define functions used for 
#  - vital rate model disagnostics
#  - vital rate model selection
#  - builing intergral projection models (IPM)
#  - IPM diagnostics: size eviction, number of bins, uncertainty (non-parametric bootstrap)
#  - basic analyses of IPMs: projection, population growth rates
#  - perturbation analyses of IPM: sensitivity, elasticity, life-table response experiment, vital rate replacement analyses
#  - coexistence analyses: coexistence outcomes, niche differences, relative fitness differences
#  - statistical analyses: different types of confidence intervals


#************************************************************
# ****** Load R packages ******----
#************************************************************
library(MuMIn)
library(car)
library(ggplot2)
library(Matrix)

#************************************************************
# ****** Functions for vital rate model diagnostics ******----
# codes are adapted from Ellner et al. (2016). Data-Driven Modelling of Structured Populations : A Practical Guide to the Integral Projection Model. Cham : Springer.
#************************************************************

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Diagnostics plots of linear models  ----
# this applies to growth and fecundity
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot.diag.lm <- function(x, y, data, fc) {
  # x: size t as character
  # y: vital rate to test as character
  # data: data including x and y variables
  # fc: name of species
  
  # lm model 
  m.lm = lm(as.formula(paste(y, "~", x)) , data=data)
  
  par(mfrow=c(2,2))
  # Plot residuals versus fitted for growth model
  # "Residuals should have constant variance and no trend in mean"
  plot(as.formula(paste(y, "~", x)) , data=data, xlab = x, ylab = y, main = fc)
  
  # Absolute residuals versus fitted
  # Residuals should have no trend in mean
  zhat <- fitted(m.lm)
  resid <- residuals(m.lm)
  sresid <- rstandard(m.lm)
  
  plot(zhat, sqrt(abs(sresid)), xlab = "Fitted values", ylab = "sqrt(|Std Residuals|)")
  gam.sresid <- gam(sqrt(abs(sresid)) ~ s(zhat), method = "REML")
  ss <- seq(min(zhat), max(zhat), length.out = 100)
  points(ss, predict(gam.sresid, newdata =data.frame(zhat = ss), type = "response"), type = "l", lwd=2)
  
  # Normal qq-plot
  # "Residuals should be Gaussian (< 5% out of the band)."
  qqPlot(sresid, main = "", xlab = "Normal quantiles", ylab = "Standardized residual quantiles", 
         col.lines = "black", lwd = 1)
  
  # compare to a gam fit
  # Linear or non-linear?
  m.gam <- gam(as.formula(paste(y, "~ s(", x, ")")), data = data, method = "REML")
  ss2 <- seq(min((data$size.autumn0), max((data$size.autumn0))), length.out = 100)
  plot(log(ss2), predict(m.lm, newdata = data.frame(size.autumn0 = ss2), type = "response"), type="l",xlab = x, ylab = paste(y, "_fitted"), lwd=2)
  points(log(ss2), predict(m.gam, newdata = data.frame(size.autumn0 = ss2), type = "response"), col="red", type="l", lwd=2, lty=2)
  
  # unlog
  #plot(ss2, predict(m.lm, newdata = data.frame(size.autumn0 = ss2), type = "response"), type="l",xlab = "Size t", ylab = "Fitted size t+1", lwd=2)
  #points(ss2, predict(m.gam, newdata = data.frame(size.autumn0 = ss2), type = "response"), col="red", type="l", lwd=2, lty=2)
  
  par(mfrow=c(1,1))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Diagnostics plot of generalised linear models  ----
# this applies to survival
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot.diag.sur <- function(data, fc) {
  # data: data including size (size.autumn0) and vital rate (survival) variables
  # fc: name of species
  
  # fit models to the reduced data set
  data$log.size = log(data$size.autumn0)
  m.glm = glm(survival ~ log.size, data = data, family = binomial)
  m.gam = gam(survival ~ s(log.size), data = data, family = binomial, method = "REML")
  
  # compare fitted glm, grouped data, and fitted gam
  plot(data$log.size,data$survival, xlab = "log(Size z)", ylab = "Survival probability", main = fc)
  ss <- seq(min(data$log.size), max(data$log.size), length = 200)
  points(ss, predict(m.glm, newdata = data.frame(log.size = ss), type="response"), type="l")
  ss2 <- seq(min(data$log.size), max(data$log.size), length = 20)
  points(ss2, predict(m.gam, newdata = data.frame(log.size = ss2), type = "response"), type = "p", col = "red", pch = 16)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Diagnostics plot of generalised linear models  ----
# this applies to flowering
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot.diag.flo <- function(data, fc) {
  # data: data including size (size.autumn0) and vital rate (flowering.autumn0) variables
  # fc: name of species
  
  # fit models to the reduced data set
  data$log.size = log(data$size.autumn0)
  m.glm = glm(flowering.autumn0 ~ log.size, data = data, family = binomial)
  m.gam = gam(flowering.autumn0 ~ s(log.size), data = data, family = binomial, method = "REML")
  
  # compare fitted glm, grouped data, and fitted gam
  plot(data$log.size,data$flowering.autumn0, xlab = "log(Size z)", ylab = "Flowering probability", main = fc)
  ss <- seq(min(data$log.size), max(data$log.size), length = 200)
  points(ss, predict(m.glm, newdata = data.frame(log.size = ss), type="response"), type="l")
  ss2 <- seq(min(data$log.size), max(data$log.size), length = 20)
  points(ss2, predict(m.gam, newdata = data.frame(log.size = ss2), type = "response"), type = "p", col = "red", pch = 16)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Linearity----
# To calculate the AICc of lm, lm with a quadratic term and gam models
# Similar AICc indicate linear relationships 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(MuMIn)
linearity.lm <- function(x, y, data) {
  # x: size t as character
  # y: vital rate to test as character
  # data: data including x and y variables
  
  # linear or quadratic
  m.linear = lm(as.formula(paste(y, "~", x)) , data=data); AIC.linear = AICc(m.linear)
  m.quadratic =  lm(as.formula(paste(y, "~", x, "+ I(", x,"^2)")), data=data); AIC.quadratic = AICc(m.quadratic)
  m.gam = gam(as.formula(paste(y, "~ s(", x, ")")), data = data, method = "REML"); AIC.gam = AICc(m.gam)
  
  data.frame(n = nrow(m.linear$model), AIC.linear = AIC.linear, AIC.quadratic = AIC.quadratic, AIC.gam = AIC.gam, AIC.delata = AIC.linear - AIC.quadratic)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Linearity----
# To compare the AICc of glm, glm with a quadratic term and gam models
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
linearity.glm <- function(x, y, data) {
  # x: size t as character
  # y: vital rate to test as character
  # data: data including x and y variables
  
  # linear or quadratic
  m.linear = glm(as.formula(paste(y, "~", x)) , data=data, family = binomial); AIC.linear = AICc(m.linear)
  m.quadratic =  glm(as.formula(paste(y, "~", x, "+ I(", x,"^2)")), data=data, family = binomial); AIC.quadratic = AICc(m.quadratic)
  m.gam = gam(as.formula(paste(y, "~ s(", x, ")")), data = data, family = binomial, method = "REML"); AIC.gam = AICc(m.gam)
  
  data.frame(n = nrow(m.linear$model), AIC.linear = AIC.linear, AIC.quadratic = AIC.quadratic, AIC.gam = AIC.gam, AIC.delata = AIC.linear - AIC.quadratic)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Normality----
# Not applicable to glm
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
normality.lm <- function(x, y, data) {
  # x: size t as character
  # y: vital rate to test as character
  # data: data including x and y variables
  
  m.linear = lm(as.formula(paste(y, "~", x)) , data=data)
  sresid = rstandard(m.linear)
  test = shapiro.test(sresid)
  data.frame(test$statistic, test$p.value)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Constant variance ----
# Not applicable to glm
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
constant.var.lm <- function(x, y, data) {
  # x: size t as character
  # y: vital rate to test as character
  # data: data including x and y variables
  
  m.linear = lm(as.formula(paste(y, "~", x)) , data=data)
  zhat <- fitted(m.linear)
  sresid = rstandard(m.linear)
  test = cor.test(zhat,sqrt(abs(sresid)),method="k")
  data.frame(test$statistic, test$p.value)
}

#************************************************************
# ****** Functions for vital rate model selection ******----
#************************************************************
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# function to compare candidate models of lm (family = gaussian) and glm
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

modelcmp.glm <- function(data, response, candidate, ...) {
  # data
  # response: response variable as a character, e.g. "survival"
  # candidate: candidate models as a string vector, e.g. "size.autumn0 * site * background.species"
  # ...: family = gaussian for linear models
  
  model.candidate <- paste(response, candidate, sep="~")
  aic <- NULL
  for(i in 1:length(model.candidate)) {
    #i =17
    fit <- glm(formula = as.formula(model.candidate[i]), data = data,...)
    aic[i] <- AICc(fit)
  }
  # delta AIC
  delta.aic <- aic - min(aic)
  
  # AIC weight
  weight.aic <- (exp(-0.5*delta.aic))/sum(exp(-0.5*delta.aic))
  
  out <- data.frame(y = response, n=nrow(data), models = candidate, aic = aic, aic.delta =delta.aic, aic.weight = weight.aic)
  out <- out[order(out$aic),]
  return(out)
}

#************************************************************
# ****** Functions for making intergral projection model (IPMs) ****** ----
#************************************************************

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Vital rate parameters----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par <- c("surv.int", "surv.slope", 
         "growth.int", "growth.slope", "growth.sd", 
         "flowering.int", "flowering.slope", 
         "fecundity.int", "fecundity.slope",
         "germination.prob",
         "establishment.funiche.int", 
         "establishment.comp.int",
         "seedling.size.mean",
         "seedling.size.sd",
         "L", "U")
length(par)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Vital rate functions ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

# 2. growth function (ceiling)
growth.z1z_ceiling <- function(z1,z, params) {
  
  # mean size time t+1, z1
  mu <- params$growth.int+params$growth.slope*z
  
  # ceiling
  mu[mu > params$U] = params$U
  print(sum(mu > params$U))
  
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

# 4. seed production (gaussian and log)
seed.z.gaussian <- function(z,params) {
  # linear predictor
  u = params$fecundity.int+params$fecundity.slope*z
  seeds <- exp(u) 
  # nagetative seeds to 0
  seeds[seeds < 0] = 0 
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
  u.funiche <- params$establishment.funiche.int
  establishment.prob.funiche <- exp(u.funiche)/(1+exp(u.funiche))
  # return
  return(establishment.prob.funiche)
}

# 6.2 probability of seedling establishment in the presence of competition (logistic)
establishment.com <- function(params) {
  # establishment with background species
  u.comp <- params$establishment.comp.int
  establishment.prob.comp <- exp(u.comp)/(1+exp(u.comp))
  
  # return
  return(establishment.prob.comp)
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

# G kernel: survival-growth
G.z1z_ceiling <- function(z1,z,params) {
  # combine survival and growth
  sg <- survival.z(z, params) * growth.z1z_ceiling(z1,z,params)
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to make an IPM kernel ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# To make an IPM kernel
mk.kernel <- function(L, U, par, n, ceiling = FALSE) {
# n: number of meshpoints
# par: vital rates
# L, U: lower and upper limit of size
# d: density of competitor
# fun: full kernel 
# mesh points 
  h <- (U - L)/n
  meshpts <- L + ((1:n) - 1/2) * h
  
  # no ceiling
  if(ceiling) {
    # ceiling
    G <- h * (outer(meshpts, meshpts, G.z1z_ceiling, params = par))
    R <- h * (outer(meshpts, meshpts, R.z1z, params = par))
    K <- G+R
  }
  else {
    G <- h * (outer(meshpts, meshpts, G.z1z, params = par))
    R <- h * (outer(meshpts, meshpts, R.z1z, params = par))
    K <- G+R
  }
  
  return(list(meshpts = meshpts, G=G, R=R, K = K))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to plot an IPM kernel ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# To plot kernels
plot.kernel <- function(x,y,k,... ) {
  # x,y: meshpoints
  # k: IPM kernel
  image(x, y, t(k), ...)
}

#************************************************************
# ****** Functions for basic analyses of IPMs ******----
#************************************************************

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to project an IPM kernel ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to calculate population growth rate (lambda) by iteration----
# This is faster for 1000 x 1000 or bigger matrices
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
lambda.iter = function(k, tol=1e-8) {
  # k: IPM full kerbel
  if(sum(k == 0) == length(k)) { return(list(lambda=NA, w = NA)) }
  else {
    qmax=10*tol
    lam=1
    x=rep(1,nrow(k))   
    k=Matrix(k)
    while(qmax>tol) {
      x1 = k%*%x
      qmax=sum(abs(x1-lam*x))  
      lam=sum(x1)
      # lambda is 0 when population is 0
      if(lam == 0) { x = 0; qmax = 0 }
      else x=x1/lam # one individual in total
    } 
    # lambda: converged population grwoth rate
    # w: stablized size distribution
    return(list(lambda=lam,w=x/sum(x)))
  }
} 

# test
#k = ipm.i$K

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to calculate population growth rate (lambda) using eigenvalue----
# This is faster than using iteration for 250 x 250 or smaller matrices
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
lambda.k <- function(k, only.lambda = TRUE) {
  # k: IPM full kerbel
  # only.lambda: whether or not calcaulte lamdab only or also stable size distribution
  
  # lambda is NA when the kernel is null
  if(sum(k == 0) == length(k)) { return(list(lambda = NA, w=NA)) }
  else {
    eigen.k = eigen(k, only.values = only.lambda)
    lam = abs(eigen.k$values[1]) 
    if(only.lambda) return(list(lambda=lam, w=NA))
    else { 
      eigen.vec = eigen.k$vectors[,1]
      w = abs(eigen.vec)/sum(abs(eigen.vec))
      return(list(lambda=lam, w=w)) 
    }
  }
}

#************************************************************
# ****** Functions for IPM diagnostics ******----
#************************************************************
# Checklist of IPM modelling----
#  - Vital rate models: linearity, normality, constant variance,
#    overfitting (when sample is small), whether covariates other than size improve the models (e.g. age, year, blocks)
#  - IPMs: model structure, size eviction, number of bins, uncertainty (bootstrap)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to determine the number of bins ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
number.bins = function(params, start = 100, end = 3000, step.n = 100, tol=1e-8) {
  # k: IPM full kerbel
  # This is faster for 1000 x 1000 or bigger matrices
  
  # seq of n to try
  nn <- seq(start, end, by = step.n)
  ll <- rep(NA, length(nn))
  k1 <- mk.kernel(L = as.numeric(params$L), U = as.numeric(params$U), par=params, n=nn[1])
  ll[1] <- lambda.iter(k1$K)$lambda
  
  for(i in 2:length(nn)) {
    ki = mk.kernel(L = as.numeric(params$L), U = as.numeric(params$U), par=params, n=nn[i])
    ll[i] = lambda.iter(ki$K)$lambda
    dif.lam = ll[i] - ll[i-1]
    # check
    if(dif.lam < tol) {
      print("IPM converged")
      out <- list(n.converged = nn[i], 
                  data.frame(number.bins = nn, lambdas = ll))
      break
    }
  }
  
  # check after loop
  if(nn[i] == end) {
    print("IPM failed to converge")
    out <- list(n.converged = NA, 
                data.frame(number.bins = nn, lambdas = ll))
  }
  return(out)
} 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to calculate the probability of size eviction ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
size.eviction <- function(params = NA, L = NA, U = NA, n, growth = FALSE, ...) {
  # params: vital rates parameters
  # L, U: the lower and upper bounds implemented in the IPM
  # n: number of bins implemented in the IPM
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Probability of recruit eviction
  p_eviction_recruit_lower <- pnorm(L,mean=params$seedling.size.mean, sd=params$seedling.size.sd)
  p_eviction_recruit_upper <- pnorm(U,mean=params$seedling.size.mean, sd=params$seedling.size.sd,lower.tail = FALSE)
  p_eviction_recruit <- data.frame(p_eviction_recruit_lower = p_eviction_recruit_lower, p_eviction_recruit_upper = p_eviction_recruit_upper)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Probability of growth eviction
  if(growth) {
    # mesh points
    h <- (U - L)/n
    meshpts <- L + ((1:n) - 1/2) * h
    # using sapply
    # calculate the probability of eviction for each initial size z
    p_eviction_growth.full = 1 - sapply(meshpts, function(z) integrate(function(u) G.z1z(u,z,params), L, U, ...)$value); 
    #p_eviction_growth_upper = sapply(meshpts,function(z) integrate(function(u) G.z1z(u,z,params), U, U*10, ...)$value); 
    #p_eviction_growth_lower = sapply(meshpts,function(z) integrate(function(u) G.z1z(u,z,params), -L*10, L, ...)$value); 
    p_eviction_growth_upper = sapply(meshpts,function(z) integrate(function(u) G.z1z(u,z,params), U, Inf, ...)$value); 
    p_eviction_growth_lower = sapply(meshpts,function(z) integrate(function(u) G.z1z(u,z,params), -Inf, L, ...)$value); 
    p_eviction_growth <- data.frame(z = meshpts,
                                    p_eviction_growth.full = p_eviction_growth.full,
                                    p_eviction_growth_lower = p_eviction_growth_lower,
                                    p_eviction_growth_upper = p_eviction_growth_upper)
  }
  else {
    p_eviction_growth <- data.frame(z = NA,
                                    p_eviction_growth.full = NA,
                                    p_eviction_growth_lower = NA,
                                    p_eviction_growth_upper = NA)
  }
  # output
  list(p_eviction_recruit = p_eviction_recruit,
       p_eviction_growth = p_eviction_growth)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to estimate population growth rate after correcting size eviction by iteraction ----
# Eviction measures for "floor-ceiling" solution to eviction, 
# using iteration to compute dominant eigenvalue/vectors. 
# Computes approximation to the effect of expanding the size range from
# (minsize, maxsize) by applying a demographic floor at minsize and
# demographic ceiling at maxsize. Computes everything internally, given the 
# kernels & size range, using midpoint rule. 
# codes are adaped from: Williams, J.L., Miller, T.E.X. & Ellner, S.P. (2012). Avoiding unintentional eviction from integral projection models. Ecology, 93, 2008-2014.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
eviction_delta.lambda_iter=function(growthKernel = NA, kernel = NA, survivalFunction = NA,
                                 minsize = NA, maxsize = NA, params = NA, n.big.matrix = NA) {
  # REQUIRED ARGUMENTS 
  #  growthKernel = growth kernel function g(new.size, old.size, params)
  #  kernel = complete kernel function K(new.size,old.size,params) 
  #  minsize,maxsize = size limits in the model, aka [L,U] or [xmin,xmax] 
  #  params = parameter vector passed to the kernels. params must be defined
  #      and used in the call even if it is not used by kernel or growthKernel.
  #  survivalFunction = size-dependent survival function s(x).
  # OPTIONAL ARGUMENTS
  #  n.big.matrix = linear dimension of approximating matrix for the kernel 
  
  growthKernel=match.fun(growthKernel); 
  survivalFunction=match.fun(survivalFunction);
  kernel=match.fun(kernel); 
  
  h = (maxsize-minsize)/n.big.matrix;
  y = minsize + h*c(1:n.big.matrix)-(h/2);
  
  # kernel,v,w for uncorrected model 
  Kmat=matrix(0,n.big.matrix,n.big.matrix);
  for(j in 1:n.big.matrix) {
    Kmat[,j]=h*kernel(y,y[j],params); 
  }
  domeigK=lambda.iter(Kmat); 
  lambda=domeigK$lambda; w=domeigK$w; 
  v=lambda.iter(t(Kmat))$w; 
  
  # Calculate the probability of size eviction
  eps = 1-sapply(y,function(z) integrate(function(u) growthKernel(u,z,params), minsize, maxsize)$value); 
  eps.U = sapply(y,function(z) integrate(function(u) growthKernel(u,z,params), maxsize, Inf)$value); 
  eps.L = sapply(y,function(z) integrate(function(u) growthKernel(u,z,params), -Inf, minsize)$value); 
  
  sx=survivalFunction(y,params);
  rho=eps*sx; rho.U=eps.U*sx; rho.L=eps.L*sx; 
  
  # Construct the iteration matrix for the corrected model 
  # add big and small classes as 2 columns at the right 
  Kmat2=cbind(Kmat, kernel(y,maxsize,params),kernel(y,minsize,params));
  
  # add the bottom rows: evictees are sent to the small or large class
  Kmat2=rbind(Kmat2, c(h*rho.U,rho.U[n.big.matrix],rho.U[1])); 
  Kmat2=rbind(Kmat2, c(h*rho.L,rho.L[n.big.matrix],rho.L[1])); 
  
  lambda2=lambda.iter(Kmat2)$lambda; 
  
  vnew=v+(1/lambda)*(v[n.big.matrix]*rho.U+v[1]*rho.L);  
  dlambdaU=vnew[n.big.matrix]*sum(rho.U*w)/sum(vnew*w); 
  dlambdaL=vnew[1]*sum(rho.L*w)/sum(vnew*w); 
  
  return(list(evict=eps,evict.U=eps.U,evict.L=eps.L,lambda=lambda,
              lambda2=lambda2,dlambda=lambda2-lambda,
              dlambdaU=dlambdaU,dlambdaL=dlambdaL))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to estimate population growth rate after correcting size eviction using eigenvalue ----
# Eviction measures for "floor-ceiling" solution to eviction. 
# Computes approximation to the effect of expanding the size range from (minsize, maxsize) by applying a demographic floor at minsize and
# demographic ceiling at maxsize. Computes everything internally, given the  kernels & size range 
# codes are adaped from: Williams, J.L., Miller, T.E.X. & Ellner, S.P. (2012). Avoiding unintentional eviction from integral projection models. Ecology, 93, 2008-2014.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
eviction_delta.lambda = function(growthKernel = NA, kernel = NA, survivalFunction = NA,
                            minsize = NA, maxsize = NA, params = NA, n.big.matrix=NA) {
  # REQUIRED ARGUMENTS 
  #  growthKernel = growth kernel function g(new.size, old.size, params)
  #  kernel = complete kernel function K(new.size,old.size,params) 
  #  minsize, maxsize = size limits in the model, aka [L,U] or [xmin,xmax] 
  #  params = parameter vector passed to the kernels. params must be defined
  #      and used in the call even if it is not used by kernel or growthKernel.
  #  survivalFunction = size-dependent survival function s(x).
  # OPTIONAL ARGUMENTS
  #  n.big.matrix = linear dimension of approximating matrix for the kernel 
  # see Williams, J.L., Miller, T.E.X. & Ellner, S.P. (2012). Avoiding unintentional eviction from integral projection models. Ecology, 93, 2008-2014.

  growthKernel=match.fun(growthKernel); 
  survivalFunction=match.fun(survivalFunction);
  kernel=match.fun(kernel); 
  
  # define mesh for midpoing rule 
  h = (maxsize-minsize)/n.big.matrix;
  b = minsize+c(0:n.big.matrix)*h;
  y = 0.5*(b[1:n.big.matrix]+b[2:(n.big.matrix+1)]);
  
  # kernel,v,w for uncorrected model 
  Kmat = h*outer(y,y,kernel,params=params); 
  w=eigen(Kmat)$vectors[,1]; w=abs(w)/sum(abs(w)); 
  v=eigen(t(Kmat))$vectors[,1]; v=abs(v)/max(abs(v));
  
  eps = 1-sapply(y,function(z) integrate(function(u) growthKernel(u,z,params), minsize, maxsize)$value); 
  eps.U = sapply(y,function(z) integrate(function(u) growthKernel(u,z,params), maxsize, Inf)$value); 
  eps.L = sapply(y,function(z) integrate(function(u) growthKernel(u,z,params), -Inf, minsize)$value); 
  
  sx=survivalFunction(y,params);
  rho=eps*sx; rho.U=eps.U*sx; rho.L=eps.L*sx; 
  
  # Construct the iteration matrix for the corrected model 
  # add small and big classes as 2 columns at the right 
  Kmat2=cbind(Kmat, kernel(y,minsize,params),kernel(y,maxsize,params));
  
  # add the bottom rows: evictees are sent to the small or large class
  Kmat2=rbind(Kmat2, c(h*rho.L,rho.L[1],rho.L[n.big.matrix])); 
  Kmat2=rbind(Kmat2, c(h*rho.U,rho.U[1],rho.U[n.big.matrix])); 
  lambda=abs(eigen(Kmat)$values[1]);
  lambda2=abs(eigen(Kmat2)$values[1]);
  
  vnew=v+(1/lambda)*(v[n.big.matrix]*rho.U+v[1]*rho.L);  
  dlambdaU=vnew[n.big.matrix]*sum(rho.U*w)/sum(vnew*w); 
  dlambdaL=vnew[1]*sum(rho.L*w)/sum(vnew*w); 
  
  return(list(evict=eps,evict.U=eps.U,evict.L=eps.L,lambda=lambda,
              lambda2=lambda2,dlambda=lambda2-lambda,
              dlambdaU=dlambdaU,dlambdaL=dlambdaL))
}

#************************************************************
# ****** Functions for IPM perturbation analyses ******----
#************************************************************

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to perform matrix-level sensitivity, elasticity analyses ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to perform vital rate-level sensitivity, elasticity analyses ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to perform parameter-level sensitivity, elasticity analyses ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to perform life-table response experiment, LTRE ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to standardise contributions of vital rates ----
# This standardisation makes relative contribution compaable across species
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

#************************************************************
# ****** Functions for coexistence analyses ******----
#************************************************************

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to calculate sensitivity of species using intrinsic and invasion growth rates ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to calculate sensitivity and coexistence outcome using invasion growth rate ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to calculate ND and FD and coexistence outcome ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to plot ND, FD and coexistence outcome ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# add frame
frame.coex = function(x1 = -0.01, x2= 0.91) {
  # x1, x2 lower and upper limits of ND, in which x1 < -0.001
  x <- c(seq(x1,-0.001, length.out = 100), seq(0.001,x2,length.out = 100))
  y <- log(1/(1-x))
  return(data.frame(x=x,y=y))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# add polygan
# polygan of coexistence
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

 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# polygan of priority effects
polygan.prio = function(x1 = -3, x2= -0.001) {
  xx1 = seq(x1,x2, length.out = 100)
  xx2 = rep(-0.001,100)
  xx3 = seq(x2,x1,length.out = 100)
  xx <- c(xx1,xx2,xx3)
  yy <- c(log(1/(1-xx1)), log(1/(1-xx2)),-log(1/(1-xx3)))
  return(data.frame(x=xx,y=yy))
}

#************************************************************
# ****** Functions for statistical analyses ******----
#************************************************************

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to extrac statics from lmer ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
coef.lmer <- function(model, ...) {
  fit <- Anova(model, ...)
  return(as.data.frame(fit))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to compuate confidence interval using linear model ----
# same as: Hmisc::mean_cl_normal
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mean_ci <- function(x) {
  fit <- lm(x~1)
  c(mean(x), confint(fit)[1, ])
  }

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to compuate SE-based confidence interval ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mean_ci_se <- function(x) {
  u = mean(x)
  se = sqrt(var(x)/length(x))
  data.frame(y=u, ymin=u-1.96*se, ymax=u+1.96*se)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to compuate SD-based confidence interval ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mean_ci_sd <- function(x, ...) {
  u = mean(x, ...)
  sd = sqrt(var(x, ...))
  data.frame(y=u, ymin=u-1.96*sd, ymax=u+1.96*sd)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to compuate percentile-based confidence interval ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mean_ci_hanna <- function(x, ...) {
  u = median(x, ...)
  x.min = quantile(x, probs = 0.025, ...)
  x.max = quantile(x, probs = 0.975, ...)
  data.frame(y=u, ymin=x.min, ymax=x.max) 
}
