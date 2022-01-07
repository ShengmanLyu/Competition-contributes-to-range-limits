#########################################################################
# Functions for IPMs and coexistence analyses
# Author: Shengman Lyu (shengman.lyu@usys.ethz.ch)
# Date: 07.01.2022
# See Carroll, I.T., Cardinale, B.J. & Nisbet, R.M. (2011). Niche and fitness differences relate the maintenance of diversity to ecosystem function. Ecology, 92, 1157-1165.

#########################################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# To calculate sensitivity of species using intrinsic and invasion growth rates----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# To calculate sensitivity and coexistence outcome of 2 species using invasion growth rate----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# To calculate ND and FD and coexistence outcome of 2 species----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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