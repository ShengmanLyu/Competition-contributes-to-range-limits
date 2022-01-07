#########################################################################
# Functions for IPMs and coexistence analyses
# Author: Shengman Lyu (shengman.lyu@usys.ethz.ch)
# Date: 07.01.2022
# See Ellner, S.P., Childs, D.Z. & Rees, M. (2016). Data-Driven Modelling of Structured Populations : A Practical Guide to the Integral Projection Model. Cham : Springer.
#########################################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 0. vital rate parameters----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. probability of surviving (logistic)----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
survival.z=function(z,params) {
  # linear predictor
  u <- params$surv.int+params$surv.slope*z
  # inverse logit transformation
  survival.p <- exp(u)/(1+exp(u))
  survival.p[u>700] <- 1 # exp(710) gives Inf values
  # return the probability of survival
  return(survival.p)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. growth function (gaussian)----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. growth function (gaussian)----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# to predict size z1
growth.z <- function(z, params) {
  # z: initial size
  # z1: size at time t+1
  z1 <- params$growth.int+params$growth.slope*z
  return(z1)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3 probaility of flowering (logistic)----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
flowering.z <- function(z, params) {
  # linear predictor
  u <- params$flowering.int+params$flowering.slope*z
  # probability of flowering, inverse logit
  flowering.p <- exp(u)/(1+exp(u))
  flowering.p[u>700] <- 1 # exp(710) gives Inf values
  # return probability of flowering
  return(flowering.p) 
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4 seed production (poisson)----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4.1 seed production (gaussian)----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
seed.z.gaussian <- function(z,params) {
  # linear predictor
  u = params$fecundity.int.linear+params$fecundity.slope.linear*z
  # nagetative seeds to 0
  u[u<0] = 0 
  seeds <- u 
  # return seed production
  return(seeds)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 5 probability of seed germination (constant)----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
germination <- function(params) {
  # get the probability of seedling establishment
  germination.p <- params$germination.prob
  # return
  return(germination.p)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 6.1 probability of seedling establishment in the absence of competition (logistic)----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
establishment.fun <- function(params) {
  # establishment using FuNiche data
  u.funiche <- params$establishment.funiche
  establishment.prob.funiche <- exp(u.funiche)/(1+exp(u.funiche))
  # return
  return(establishment.prob.funiche)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 6.2 probability of seedling establishment in the presence of competition (logistic)----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
establishment.com <- function(params) {
  # establishment with background species
  u.comp <- params$establishment.comp
  establishment.prob.comp <- exp(u.comp)/(1+exp(u.comp))
  
  # return
  return(establishment.prob.comp)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 7 seedling size (constant)-----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
seedling.z1 <- function(z1, params) {
  # probability density of recruits
  rpd <- dnorm(z1, mean=params$seedling.size.mean, sd=params$seedling.size.sd)
  # return probability density
  return(rpd)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# G kernel: survival-growth----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
G.z1z <- function(z1,z,params) {
  # combine survival and growth
  sg <- survival.z(z, params) * growth.z1z(z1,z,params)
  # return
  return(sg)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# R kernel: reproduction-----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# full kernel----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
full.z1z = function(z1,z,params) {
  # combine growth and reproduction kernel
  G.z1z(z1,z,params) + R.z1z(z1,z,params)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# To make an IPM kernel----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
  return(list(meshpts = meshpts, G=G, R=R, K=K))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# To plot a kernel----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot.kernel <- function(x,y,k,... ) {
  # x,y: meshpoints
  # k: IPM kernel
  image(x, y, t(k), ...)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# To project an IPM-----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# To calculate population growth rate (lambda) using kernel----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
lambda.k <- function(k) {
# k: IPM full kerbel
  l=Re(eigen(k, only.values = TRUE)$values[1]) 
  return(l)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# To detect size eviction-----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
size.eviction <- function(params = NA, L = NA, U = NA, n, growth = FALSE, ...) {
  # params: vital rates parameters
  # L, U: the lower and upper bounds implemented in the IPM
  # n: number of bins implemented in the IPM
  # growth: whether perform size eviction on growth
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Probability of recruit eviction
  p_eviction_recruit_lower <- pnorm(L,mean=params$seedling.size.mean, sd=params$seedling.size.sd)
  p_eviction_recruit_upper <- pnorm(U,mean=params$seedling.size.mean, sd=params$seedling.size.sd,lower.tail = FALSE)
  p_eviction_recruit <- data.frame(p_eviction_recruit_lower = p_eviction_recruit_lower, p_eviction_recruit_upper = p_eviction_recruit_upper)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Probability of growth eviction
  p_eviction_growth <- data.frame(z = NA,
                                  p_eviction_growth_lower = NA,
                                  p_eviction_growth_upper = NA)
  if(growth) {
    # mesh points
    h <- (U - L)/n
    meshpts <- L + ((1:n) - 1/2) * h
    # data frame for output
    p_eviction_growth <- data.frame(z = meshpts,
                                    p_eviction_growth_lower = NA,
                                    p_eviction_growth_upper = NA)
    # calculate the probability of eviction for each initial size z
    for(i in 1:length(p_eviction_growth$z)) {
      # initial size z
      z.i <- p_eviction_growth$z[i]
      # make a new growth function for integral
      G.z1 <- function(z1) G.z1z(z1, z=z.i, params)
      # integrate
      p_eviction_grwoth_lower.i <- integrate(G.z1, lower = -L*10, upper = L, ...)$value
      p_eviction_grwoth_upper.i <- integrate(G.z1, lower = U, upper = U*10, ...)$value
      #p_eviction_grwoth_lower.i <- integrate(G.z1, lower = -Inf, upper = L, ...)$value
      #p_eviction_grwoth_upper.i <- integrate(G.z1, lower = U, upper = Inf, ...)$value
      p_eviction_growth[i,"p_eviction_growth_lower"] = p_eviction_grwoth_lower.i
      p_eviction_growth[i,"p_eviction_growth_upper"] = p_eviction_grwoth_upper.i
    }
  }
  
  # output
  list(p_eviction_recruit = p_eviction_recruit,
       p_eviction_growth = p_eviction_growth)
}
