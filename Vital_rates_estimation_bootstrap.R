#****************************************************************************
# Estimating vital rates and bootstrap
# Author: Shengman Lyu (shengman.lyu@usys.ethz.ch)
# Date: 07.01.2022
# See paper Competition contributes to both warm and cold range edges
#****************************************************************************
library(readxl)
library(tidyverse)
library(lme4)
library(lmerTest)
library(MASS)
library(Hmisc)
library(ggpubr)

rm(list=ls())

#****************************************************************************
# ** 1. Estimate vital rates using the best model and bootstrap ** ------
#****************************************************************************

# 0 read data---------------------------------------
# NICH data
d <- read_excel("data_181920_v6.xlsx", na="NA",col_names = TRUE)
d
# remove seedling planted in autumn 2019 because mostly of them died: 16755
d <- d %>% filter(!(!is.na(note.size.autumn0) & note.size.autumn0 == "planting.autumn0" & year == 2020))
d

# convert to factor in order to relevel them
d$site <- factor(d$site, levels=c("Les Posses", "Solalex", "Anzeindaz"))
d$focal.species <- factor(d$focal.species)
d$focal.before <- factor(d$focal.before)
d$background.species <- factor(d$background.species)
d

# biomass of greenhouse seedling
b <- read_excel("/Users/slyu/LVSM/NICH/Data/Biomass/focal biomass/biomass data for correlation_2020.xlsx",na="NA", col_names = TRUE)
b$greenhouse <- ifelse(b$site == "greenhouse" | b$site == "chngreenhouse", "yes", "no")
b$flowering <- ifelse(b$stem.number == 0, "no", "yes")
b

# FuNiche data
funiche <- read_csv("/Users/slyu/LVSM/NICH/Data/germination/FunNiche/FunNiche_data_2019_2020_dec2020_20210114.csv")
f <- filter(funiche, site %in% c(4,13,19))
f$site <- as.character(f$site)
f$site <- dplyr::recode(f$site, "4" = "Les Posses", "13" = "Solalex", "19"="Anzeindaz")
f$species <- factor(f$species)
f$site <- factor(f$site)
f

# vital rates
vr <- read_excel("/Users/slyu/LVSM/NICH/Results/IPM/vital rates.xlsx", col_names = TRUE, na="NA")
vr <- vr[-1,]
vr

# bootstrapped vital rates
nboot=99
vr.boot <- array(NA, dim = c(nrow(vr), ncol(vr), nboot))

# standard error of vital rates
vr.se <- read_excel("/Users/slyu/LVSM/NICH/Results/IPM/vital rates_se.xlsx", col_names = TRUE, na="NA")
vr.se <- vr.se[-1,]
vr.se

# best vital rates models
vr.models <- read_excel("/Users/slyu/LVSM/NICH/Results/IPM/vital rates_model comparison.xlsx", na="NA") # best.model_20210105
table(vr.models$response)

# species
sps <- c("Anal","Armo","Asal","Plal", "Poal", "Seca", "Trba",        # 7 alpine species 
         "Brer", "Crbi", "Daca","Melu" , "Plla", "Potr", "Sapr", # 7 lowland species
         "none")
sps

# 0 functions --------------------------------------
mean_sd <- function(x) {c(mean(x), sd(x))}

for(n in 1:nrow(vr)) {
  # n = 1
  # fc = "Brer"
  # bg = "none"
  # st = "Les Posses"
  vi <- vr[n,]
  fc <- as.character(vi$focal.species)
  bg <- as.character(vi$background.species)
  st <- as.character(vi$site)
  
  # skip when no pair
  if(is.na(vi$pair) | bg == "site" | bg == "site_none" | bg == "species" | bg == "species_none") next
  
  # print
  print(n)
  print(paste(fc,st,bg, sep="_"))
  
  # relevel focal species, background species and site
  d$focal.species <- relevel(d$focal.species, ref=fc)
  d$focal.before <- relevel(d$focal.before, ref=fc)
  d$background.species <- relevel(d$background.species, ref=bg)
  d$site <- relevel(d$site, ref=st)
  
  # relevel FuNiche data where Daca is absent
  if(fc != "Daca") {
    f$species <- relevel(f$species, ref=fc)
    f$site <- relevel(f$site, ref=st)
  }
  
  ## 1 survival-------------------------------------------------
  # size model produced very similar results with probability
  d.sur <- filter(d, focal.before == fc &
                    background.species %in% sps &
                    !is.na(size.autumn0) & !is.na(survival) &
                    (!is.na(if.same) & if.same==TRUE) )
  
  sample.size.pair <- table(d.sur[,c("background.species", "site")])[1,1]
  vr[n, "n.survival"] <- sample.size.pair
  
  y.sur <- as.character(filter(vr.models, species == fc & response == "survival")[3])
  m.sur <- as.character(filter(vr.models, species == fc & response == "survival")[5])
  
  # glm
  glm.sur <- glm(as.formula(paste(y.sur, m.sur, sep="~")), family="binomial", data=d.sur)
  
  # intercept and sigma
  vr[n, "surv.int"] <- coef(glm.sur)[1]
  vr[n, "surv.sd.global"] <- sigma(glm.sur)
  
  # se of intercept
  vr.se[n, "surv.int"] <- summary(glm.sur)$coefficients[1,2]
  
  # slope if model include size.autumn0
  if("0" %in% strsplit(m.sur, split = "")[[1]]) { 
    # bootstrap
    boot.sur <- mvrnorm(nboot,mu=coef(glm.sur)[1:2], Sigma=vcov(glm.sur)[1:2,1:2])
    
    vr[n, "surv.slope"] <- coef(glm.sur)[2]
    vr.boot[n, match("surv.int", colnames(vr)), 1:nboot] = boot.sur[,1]
    vr.boot[n, match("surv.slope", colnames(vr)), 1:nboot] = boot.sur[,2]
  } else { 
    # bootstrap
    boot.sur <- mvrnorm(nboot,mu=coef(glm.sur)[1], Sigma=vcov(glm.sur)[1])
    
    vr[n, "surv.slope"] <- 0 
    vr.boot[n, match("surv.int", colnames(vr)), 1:nboot] = boot.sur[,1]
    vr.boot[n, match("surv.slope", colnames(vr)), 1:nboot] = 0
  }
  
  rm(m.sur)
  rm(glm.sur)
  rm(sur.boot)
  print("survival")
  
  ## 2 growth-----------------------------------------------------
  d.gro <- filter(d, focal.species == fc &
                    background.species %in% sps &
                    !is.na(survival) & survival == 1 &
                    !is.na(size.autumn0) & !is.na(size.autumn1) &
                    (!is.na(if.same) & if.same==TRUE) )
  # sample size
  sample.size.pair <- table(d.gro[,c("background.species", "site")])[1,1]
  vr[n, "n.growth"] <- sample.size.pair
  
  # model
  y.gro <- as.character(filter(vr.models, species == fc & response == "growth")[3])
  m.gro <- as.character(filter(vr.models, species == fc & response == "growth")[5])
  
  # full model skip when sample < 3
  #if(sample.size.pair > 2) {
  
  # glm
  lm.gro <- lm(as.formula(paste(y.gro, m.gro, sep="~")), data=d.gro)
  
  # bootstrap
  boot.gro <- mvrnorm(nboot,mu=coef(lm.gro)[1:2], Sigma=vcov(lm.gro)[1:2,1:2])
  
  # intrcept and slope
  vr[n, "growth.int"] <- coef(lm.gro)[1]
  vr[n, "growth.slope"] <- coef(lm.gro)[2]
  vr.boot[n, match("growth.int", colnames(vr)), 1:nboot] = boot.gro[,1]
  vr.boot[n, match("growth.slope", colnames(vr)), 1:nboot] = boot.gro[,2]
  
  # SD of residual at pair level (no bootstapping)
  if(sample.size.pair > 1) {
    d.gro.sd <- filter(d.gro, site == st & background.species == bg)
    resi.gro <- d.gro.sd$size.autumn1 - (as.numeric(vr[n, "growth.int"]) + d.gro.sd$size.autumn0 * as.numeric(vr[n, "growth.slope"]))
    vr[n,"growth.sd"] <- sd(resi.gro)
    vr.boot[n, match("growth.sd", colnames(vr)), 1:nboot] = sd(resi.gro)
  }
  else{
    vr[n,"growth.sd"] <- 0.01
    vr.boot[n, match("growth.sd", colnames(vr)), 1:nboot] = 0.01
  }
  
  rm(m.gro)
  rm(lm.gro)
  rm(boot.gro)
  print("growth") 
  
  # 3 flowering-----------------------------------------------------
  d.flo <- filter(d, focal.before == fc &
                    background.species %in% sps &
                    !is.na(size.autumn0) & !is.na(flowering.autumn0) &
                    (!is.na(if.same) & if.same==TRUE) ) 
  
  sample.size.pair <- table(d.flo[,c("background.species", "site")])[1,1]
  vr[n, "n.flowering"] <- sample.size.pair
  
  y.flo <- as.character(filter(vr.models, species == fc & response == "flowering")[3])
  m.flo <- as.character(filter(vr.models, species == fc & response == "flowering")[5])
  
  # full model skip when sample < 2
  if(sample.size.pair > 1 | m.flo != "size.autumn0 * site * background.species") {
    # glm
    glm.flo <- glm(as.formula(paste(y.flo, m.flo, sep="~")), family="binomial", data=d.flo)
    
    # slope if model include size.autumn0
    if("0" %in% strsplit(m.flo, split = "")[[1]]) { 
      # bootstrap
      boot.flo <- mvrnorm(nboot,mu=coef(glm.flo)[1:2], Sigma=vcov(glm.flo)[1:2,1:2])
      
      vr[n, "flowering.int"] <- coef(glm.flo)[1]  
      vr[n, "flowering.slope"] <- coef(glm.flo)[2] 
      vr.boot[n, match("flowering.int", colnames(vr)), 1:nboot] = boot.flo[,1]
      vr.boot[n, match("flowering.slope", colnames(vr)), 1:nboot] = boot.flo[,2]
    } else { 
      # bootstrap
      boot.flo <- mvrnorm(nboot,mu=coef(glm.flo)[1], Sigma=vcov(glm.flo)[1])
      
      vr[n, "flowering.int"] <- coef(glm.flo)[1]
      vr[n, "flowering.slope"] <- 0
      vr.boot[n, match("flowering.int", colnames(vr)), 1:nboot] = boot.flo[,1]
      vr.boot[n, match("flowering.slope", colnames(vr)), 1:nboot] = 0
    }
    
    rm(m.flo)
    rm(glm.flo)
    rm(boot.flo)
    print("flowering")
  }
  
  #**************************************************
  # 4.1 fecundity_poisson ----
  d.fec <- filter(d, focal.before == fc &
                    background.species %in% sps &
                    !is.na(size.autumn0) & !is.na(seeds.autumn0) &
                    (!is.na(if.same) &if.same==TRUE) )
  
  # sample size
  sample.size.pair <- table(d.fec[,c("background.species", "site")])[1,1]
  sample.size.site <- table(d.fec[,"site"])[1]
  sample.size.bg <- table(d.fec[,"background.species"])[1]
  vr[n, "n.fecundity"] <- sample.size.pair
  
  y.fec <- as.character(filter(vr.models, species == fc & comment1 == "linear model" & response == "fecundity")[3])
  m.fec <- as.character(filter(vr.models, species == fc & comment1 == "linear model" & response == "fecundity")[5])
  
  # skip when model is NA or full model and sample < 3
  if(!is.na(m.fec)) {
    # glm
    glm.fec <- glm(as.formula(paste(y.fec, m.fec, sep="~")), family=poisson(link="log"), data=d.fec)
    
    # bootstrap
    boot.fec <- mvrnorm(nboot,mu=coef(glm.fec)[1:2], Sigma=vcov(glm.fec)[1:2,1:2])
    
    # coefs
    vr[n, "fecundity.int"] <- coef(glm.fec)[1]
    vr[n, "fecundity.slope"] <- coef(glm.fec)[2]
    vr.boot[n, match("fecundity.int", colnames(vr)), 1:nboot] = boot.fec[,1]
    vr.boot[n, match("fecundity.slope", colnames(vr)), 1:nboot] = boot.fec[,2]
    
    rm(glm.fec)
    rm(boot.fec)
    print("fecunidty")
  }
  
  #**************************************************
  # 4.2 fecundity_linear ----
  d.fec <- filter(d, focal.before == fc &
                    background.species %in% sps &
                    !is.na(size.autumn0) & !is.na(seeds.autumn0) &
                    (!is.na(if.same) &if.same==TRUE) )
  
  # sample size
  sample.size.pair <- table(d.fec[,c("background.species", "site")])[1,1]
  sample.size.site <- table(d.fec[,"site"])[1]
  sample.size.bg <- table(d.fec[,"background.species"])[1]
  
  y.fec <- as.character(filter(vr.models, species == fc & comment1 == "linear model" & response == "fecundity")[3])
  m.fec <- as.character(filter(vr.models, species == fc & comment1 == "linear model" & response == "fecundity")[5])
  
  # skip when model is NA or full model and sample < 3
  if(!is.na(m.fec)) {
    # glm
    lm.fec <- glm(as.formula(paste(y.fec, m.fec, sep="~")), family=gaussian, data=d.fec)
    
    # bootstrap
    boot.fec.lm <- mvrnorm(nboot,mu=coef(lm.fec)[1:2], Sigma=vcov(lm.fec)[1:2,1:2])
    
    # coefs
    vr[n, "fecundity.int.linear"] <- coef(lm.fec)[1]
    vr[n, "fecundity.slope.linear"] <- coef(lm.fec)[2]
    vr.boot[n, match("fecundity.int.linear", colnames(vr)), 1:nboot] = boot.fec.lm[,1]
    vr.boot[n, match("fecundity.slope.linear", colnames(vr)), 1:nboot] = boot.fec.lm[,2]
    
    rm(m.fec)
    rm(lm.fec)
    rm(boot.fec.lm)
    print("fecunidty")
  }
  
  # 5 germination (no bootstrap) ----------------------------------------------
  #if(fc == "Seca") f.ger <- filter(f, species == fc & site %in% c("Solalex", "Anzeindaz"))
  #else f.ger <- filter(f, species == fc)
  f.ger <- filter(f, species == fc)
  
  sample.size.pair <- table(f.ger[,"site"])[1]
  vr[n,"n.germination"] <- sample.size.pair
  
  y.ger <- as.character(filter(vr.models, species == fc & response == "germination")[3])
  m.ger <- as.character(filter(vr.models, species == fc & response == "germination")[5])
  
  if(!is.na(m.ger)) {
    # lm
    lm.ger <- lm(as.formula(paste(y.ger, m.ger, sep="~")), data=f.ger)
    
    # intercept
    vr[n,"germination.prob"] <- coef(lm.ger)[1]
    vr.boot[n, match("germination.prob", colnames(vr)), 1:nboot] = coef(lm.ger)[1]
    
    rm(m.ger)
    rm(lm.ger)
    print("germination")
  }
  
  #**************************************************
  # 6.1 seedling establishment of FuNiche (no bootstrap) ----
  f.est <- filter(f, species == fc & !is.na(establishement_summer_19))
  sample.size.site <- table(f.est[,"site"])[1]
  sample.size.sp <- nrow(f.est)
  vr[n,"n.establishment.funiche"] <- sample.size.site
  
  y.est <- as.character(filter(vr.models, species == fc & response == "establishment_FuNiche")[3])
  m.est <- as.character(filter(vr.models, species == fc & response == "establishment_FuNiche")[5])
  
  #if(m.est == "site" & fc != "Daca" & sample.size.site >3) {
  if(!is.na(m.est)) {
    # lm
    glm.est <- glm(as.formula(paste(y.est, m.est, sep="~")), data=f.est, family = "binomial")
    
    # intercept
    vr[n,"establishment.funiche"] <- coef(glm.est)[1]
    vr.boot[n, match("establishment.funiche", colnames(vr)), 1:nboot] = coef(glm.est)[1]
    
    rm(m.est)
    rm(glm.est)
    print("establishment.funiche")
  } else{
    vr[n,"establishment.funiche"] <- 20
    vr.boot[n, match("establishment.funiche", colnames(vr)), 1:nboot] = 20
    print("establishment.funiche")
  }
  
  #**************************************************
  # 6.2 seedling establishment with competition_best model (no bootstrap) ----
  d.est <- filter(d, focal.species == fc & 
                    (!is.na(planting.spring1) & planting.spring1=="yes") &
                    !is.na(survival.summer1.seedling) &
                    background.species %in% sps)
  
  # remove Poal and Potr because of different bias in site
  if(fc == "Armo") d.est <- filter(d.est, !background.species %in% c("Poal", "Potr"))
  
  sample.size.pair <- table(d.est[,c("background.species", "site")])[1,1]
  vr[n,"n.establishment.comp"] <- sample.size.pair
  vr[n,"n.establishment.comp1"] <- sample.size.pair
  
  y.est <- as.character(filter(vr.models, species == fc & comment1 == "logistic model" & response == "establishment")[3])
  m.est <- as.character(filter(vr.models, species == fc & comment1 == "logistic model" & response == "establishment")[5])
  
  #if(sample.size.pair > 2 & !is.na(m.est)) {
  
  # glm
  glm.est2 <- glm(as.formula(paste(y.est, m.est, sep="~")), data=d.est, family="binomial")
  
  # intercept
  vr[n,"establishment.comp"] <- coef(glm.est2)[1]
  vr.boot[n, match("establishment.comp", colnames(vr)), 1:nboot] = coef(glm.est2)[1]
  
  #} 
  
  rm(m.est)
  rm(glm.est2)
  print("establishment12") 
  
  #**************************************************
  # 7 recruit size (no bootstrap) ----
  b.rec <- filter(b, species == fc & greenhouse == "yes" & flowering == "no" & !is.na(drymass.g))
  sample.size.sp <- nrow(b.rec)
  vr[n, "n.seedling.size"] <- sample.size.sp
  
  lm.rec <- b.rec %>% 
    .$drymass.g %>%
    mean_sd
  vr[n, "seedling.size.mean"] = lm.rec[1] 
  vr[n, "seedling.size.sd"] = lm.rec[2]
  vr.boot[n, match("seedling.size.mean", colnames(vr)), 1:nboot] = lm.rec[1] 
  vr.boot[n, match("seedling.size.sd", colnames(vr)), 1:nboot] = lm.rec[2] 
  
  rm(lm.rec)
  print("recruit size")
  
  #**************************************************
  # 8 size range: L, U (no bootstrap) ----
  d.siz <- filter(d, focal.species == fc &
                    !is.na(size.autumn0) &
                    background.species %in% sps)
  
  b.siz <- filter(b, species == fc & greenhouse == "yes" & flowering == "no" & !is.na(drymass.g))
  lm.siz.gr <- range(b.siz$drymass.g)
  
  sample.size.pair <- table(d.siz[,c("background.species", "site")])[1,1]
  sample.size.site <- table(d.siz[,"site"])[1]
  sample.size.bg <- table(d.siz[,"background.species"])[1]
  vr[n, "n.size.range"] <- sample.size.pair
  
  y.siz <- as.character(filter(vr.models, species == fc & response == "size range")[3])
  m.siz <- as.character(filter(vr.models, species == fc & response == "size range")[5])
  
  if(m.siz == "site*background.species" | m.siz == "site + background.species") {
    lm.siz <- d.siz %>% 
      filter(background.species == bg & site == st) %>%
      .$size.autumn0 %>%
      range
    vr[n, "L"] = min(lm.siz, lm.siz.gr, na.rm = TRUE)
    vr[n, "U"] = max(lm.siz, lm.siz.gr, na.rm = TRUE)
    vr.boot[n, match("L", colnames(vr)), 1:nboot] = as.numeric(vr[n, "L"])
    vr.boot[n, match("U", colnames(vr)), 1:nboot] = as.numeric(vr[n, "U"])
    
    rm(lm.siz.gr); rm(lm.siz)
    print("size range")
  } else if(m.siz == "background.species") {
    lm.siz <- d.siz %>% 
      filter(background.species == bg) %>%
      .$size.autumn0 %>%
      range
    vr[n, "L"] = min(lm.siz, lm.siz.gr, na.rm = TRUE)
    vr[n, "U"] = max(lm.siz, lm.siz.gr, na.rm = TRUE)
    vr.boot[n, match("L", colnames(vr)), 1:nboot] = as.numeric(vr[n, "L"])
    vr.boot[n, match("U", colnames(vr)), 1:nboot] = as.numeric(vr[n, "U"])
    rm(lm.siz.gr); rm(lm.siz)
    print("size range")
  }
  
} # for loop

#****************************************************************************
# ** 2 Species and site level vital rates ** ------
#****************************************************************************
for(n in 1:nrow(vr)) {
  # n = 1
  vi <- vr[n,]
  fc <- as.character(vi$focal.species)
  bg <- as.character(vi$background.species)
  st <- as.character(vi$site)
  
  # skip
  if(bg != "site" & bg != "site_none" & bg != "species" & bg != "species_none") next
  
  # print
  print(n)
  print(paste(fc,st,bg, sep="_"))
  
  # relevel focal species, background species and site
  d$focal.species <- relevel(d$focal.species, ref=fc)
  d$focal.before <- relevel(d$focal.before, ref=fc)
  # d$background.species <- relevel(d$background.species, ref=bg)
  d$site <- relevel(d$site, ref=st)
  
  # site level with/without non-competition focals
  if(bg=="site" | bg == "species") { 
    sps <- c("Anal","Armo","Asal","Plal", "Poal", "Seca", "Trba",        # 7 alpine species 
             "Brer", "Crbi", "Daca","Melu" , "Plla", "Potr", "Sapr") # 7 lowland species 
  } else if(bg == "site_none" | bg == "species_none") { 
    sps <- c("Anal","Armo","Asal","Plal", "Poal", "Seca", "Trba",        # 7 alpine species 
             "Brer", "Crbi", "Daca","Melu" , "Plla", "Potr", "Sapr", # 7 lowland species
             "none")
  }
  
  ## 1 survival-------------------------------------------------
  # size is not importamt here!
  d.sur <- filter(d, focal.before == fc &
                    background.species %in% sps &
                    !is.na(size.autumn0) & !is.na(survival) &
                    (!is.na(if.same) & if.same==TRUE) )
  
  sample.size.site <- table(d.sur[,c("site")])[1]
  sample.size.sps <- nrow(d.sur)
  
  if(bg == "site" | bg == "site_none") {
    vr[n, "n.survival"] <- sample.size.site
    glm.sur <- glm(survival ~ size.autumn0 + site, family="binomial", data=d.sur)
  } else if(bg == "species" | bg == "species_none") {
    vr[n, "n.survival"] <- sample.size.sps
    glm.sur <- glm(survival ~ size.autumn0, family="binomial", data=d.sur)
  }
  
  # bootstrap
  boot.sur <- mvrnorm(nboot,mu=coef(glm.sur)[1:2], Sigma=vcov(glm.sur)[1:2,1:2])
  
  # coefs
  vr[n, "surv.int"] <- coef(glm.sur)[1]
  vr[n, "surv.slope"] <- coef(glm.sur)[2]
  vr.boot[n, match("surv.int", colnames(vr)), 1:nboot] = boot.sur[,1]
  vr.boot[n, match("surv.slope", colnames(vr)), 1:nboot] = boot.sur[,2]
  
  
  rm(glm.sur)
  rm(boot.sur)
  print("survival")
  
  ## 2 growth-----------------------------------------------------
  d.gro <- filter(d, focal.species == fc &
                    background.species %in% sps &
                    !is.na(survival) & survival == 1 &
                    !is.na(size.autumn0) & !is.na(size.autumn1) &
                    (!is.na(if.same) & if.same==TRUE) )
  
  # sample size
  sample.size.site <- table(d.gro[,c("site")])[1]
  sample.size.sps <- nrow(d.gro)
  
  if(bg == "site" | bg == "site_none") {
    vr[n, "n.growth"] <- sample.size.site
    lm.gro <- lm(size.autumn1~size.autumn0 + site, data=d.gro)
  } else if(bg == "species" | bg == "species_none") {
    vr[n, "n.growth"] <- sample.size.sps
    lm.gro <- lm(size.autumn1~size.autumn0, data=d.gro)
  }
  
  # bootstrap
  boot.gro <- mvrnorm(nboot,mu=coef(lm.gro)[1:2], Sigma=vcov(lm.gro)[1:2,1:2])
  
  # coefs and sigma
  vr[n, "growth.int"] <- coef(lm.gro)[1]
  vr[n, "growth.slope"] <- coef(lm.gro)[2]
  vr[n,"growth.sd"] <- sd(resid(lm.gro))
  vr.boot[n, match("growth.int", colnames(vr)), 1:nboot] = boot.gro[,1]
  vr.boot[n, match("growth.slope", colnames(vr)), 1:nboot] = boot.gro[,2]
  vr.boot[n, match("growth.sd", colnames(vr)), 1:nboot] = sd(resid(lm.gro))
  
  rm(lm.gro)
  rm(boot.gro)
  print("growth")
  
  # 3 flowering-----------------------------------------------------
  d.flo <- filter(d, focal.before == fc &
                    background.species %in% sps &
                    !is.na(size.autumn0) & !is.na(flowering.autumn0) &
                    (!is.na(if.same) & if.same==TRUE) ) 
  
  # sample size
  sample.size.site <- table(d.flo[,c("site")])[1]
  sample.size.sps <- nrow(d.flo)
  
  if(bg == "site" | bg == "site_none") {
    vr[n, "n.flowering"] <- sample.size.site
    glm.flo <- glm(flowering.autumn0 ~ size.autumn0 + site, family="binomial", data=d.flo)
  } else if(bg == "species" | bg == "species_none") {
    vr[n, "n.flowering"] <- sample.size.sps
    glm.flo <- glm(flowering.autumn0 ~ size.autumn0, family="binomial", data=d.flo)
  }
  
  # bootstrap
  boot.flo <- mvrnorm(nboot,mu=coef(glm.flo)[1:2], Sigma=vcov(glm.flo)[1:2,1:2])
  
  # intercept and sigma
  vr[n, "flowering.int"] <- coef(glm.flo)[1]
  vr[n, "flowering.slope"] <- coef(glm.flo)[2]
  vr.boot[n, match("flowering.int", colnames(vr)), 1:nboot] = boot.flo[,1]
  vr.boot[n, match("flowering.slope", colnames(vr)), 1:nboot] = boot.flo[,2]
  
  rm(glm.flo)
  rm(boot.flo)
  print("flowering")
  
  # 4 fecundity_poisson -------------------------------------------
  d.fec <- filter(d, focal.before == fc &
                    background.species %in% sps &
                    !is.na(size.autumn0) & !is.na(seeds.autumn0) &
                    (!is.na(if.same) &if.same==TRUE) )
  
  # sample size
  sample.size.site <- table(d.fec[,c("site")])[1]
  sample.size.sps <- nrow(d.fec)
  
  if(fc != "Armo" & (bg == "site" | bg == "site_none")) {
    vr[n, "n.fecundity"] <- sample.size.site
    glm.fec <- glm(seeds.autumn0 ~ size.autumn0 + site, family="poisson", data=d.fec)
    
    # bootstrap
    boot.fec <- mvrnorm(nboot,mu=coef(glm.fec)[1:2], Sigma=vcov(glm.fec)[1:2,1:2])
    
    # coefs and sigma
    vr[n, "fecundity.int"] <- coef(glm.fec)[1]
    vr[n, "fecundity.slope"] <- coef(glm.fec)[2]
    vr.boot[n, match("fecundity.int", colnames(vr)), 1:nboot] = boot.fec[,1]
    vr.boot[n, match("fecundity.slope", colnames(vr)), 1:nboot] = boot.fec[,2]
    
  } else if(fc != "Armo" & (bg == "species" | bg == "species_none")) {
    vr[n, "n.fecundity"] <- sample.size.sps
    glm.fec <- glm(seeds.autumn0 ~ size.autumn0, family="poisson", data=d.fec)
    
    # bootstrap
    boot.fec <- mvrnorm(nboot,mu=coef(glm.fec)[1:2], Sigma=vcov(glm.fec)[1:2,1:2])
    
    # coefs and sigma
    vr[n, "fecundity.int"] <- coef(glm.fec)[1]
    vr[n, "fecundity.slope"] <- coef(glm.fec)[2]
    vr.boot[n, match("fecundity.int", colnames(vr)), 1:nboot] = boot.fec[,1]
    vr.boot[n, match("fecundity.slope", colnames(vr)), 1:nboot] = boot.fec[,2]
    
  }
  rm(glm.fec)
  rm(boot.fec)
  print("fecunidty")
  
  
  # 4 fecundity_linear -----------------------------------------------
  d.fec <- filter(d, focal.before == fc &
                    background.species %in% sps &
                    !is.na(size.autumn0) & !is.na(seeds.autumn0) &
                    (!is.na(if.same) &if.same==TRUE) )
  
  # sample size
  sample.size.site <- table(d.fec[,c("site")])[1]
  sample.size.sps <- nrow(d.fec)
  
  if(fc != "Armo" & (bg == "site" | bg == "site_none")) {
    vr[n, "n.fecundity"] <- sample.size.site
    lm.fec <- glm(seeds.autumn0 ~ size.autumn0 + site, data=d.fec)
    
    # bootstrap
    boot.fec.lm <- mvrnorm(nboot,mu=coef(lm.fec)[1:2], Sigma=vcov(lm.fec)[1:2,1:2])
    
    # coefs and sigma
    vr[n, "fecundity.int.linear"] <- coef(lm.fec)[1]
    vr[n, "fecundity.slope.linear"] <- coef(lm.fec)[2]
    vr.boot[n, match("fecundity.int.linear", colnames(vr)), 1:nboot] = boot.fec.lm[,1]
    vr.boot[n, match("fecundity.slope.linear", colnames(vr)), 1:nboot] = boot.fec.lm[,2]
    
  } else if(fc != "Armo" & (bg == "species" | bg == "species_none")) {
    vr[n, "n.fecundity"] <- sample.size.sps
    lm.fec <- glm(seeds.autumn0 ~ size.autumn0, data=d.fec)
    
    # bootstrap
    boot.fec.lm <- mvrnorm(nboot,mu=coef(lm.fec)[1:2], Sigma=vcov(lm.fec)[1:2,1:2])
    
    # coefs and sigma
    vr[n, "fecundity.int.linear"] <- coef(lm.fec)[1]
    vr[n, "fecundity.slope.linear"] <- coef(lm.fec)[2]
    vr.boot[n, match("fecundity.int.linear", colnames(vr)), 1:nboot] = boot.fec.lm[,1]
    vr.boot[n, match("fecundity.slope.linear", colnames(vr)), 1:nboot] = boot.fec.lm[,2]
  }
  
  rm(lm.fec)
  rm(boot.fec.lm)
  print("fecunidty2")
  
  # 5 germination ----------------------------------------------
  #if(fc == "Seca") f.ger <- filter(f, species == fc & site %in% c("Solalex", "Anzeindaz"))
  #else f.ger <- filter(f, species == fc)
  f.ger <- filter(f, species == fc)
  
  sample.size.pair <- table(f.ger[,"site"])[1]
  vr[n,"n.germination"] <- sample.size.pair
  
  y.ger <- as.character(filter(vr.models, species == fc & response == "germination")[3])
  m.ger <- as.character(filter(vr.models, species == fc & response == "germination")[5])
  
  if(!is.na(m.ger)) {
    # lm
    lm.ger <- lm(as.formula(paste(y.ger, m.ger, sep="~")), data=f.ger)
    # intercept and sigma
    vr[n,"germination.prob"] <- coef(lm.ger)[1]
    vr.boot[n, match("germination.prob", colnames(vr)), 1:nboot] = coef(lm.ger)[1]
    
    rm(m.ger)
    rm(lm.ger)
    print("germination")
  }
  
  # 6.1 seedling establishment of FuNiche-------------------------------------
  f.est <- filter(f, species == fc & !is.na(establishement_summer_19))
  sample.size.site <- table(f.est[,"site"])[1]
  sample.size.sp <- nrow(f.est)
  vr[n,"n.establishment.funiche"] <- sample.size.site
  
  y.est <- as.character(filter(vr.models, species == fc & response == "establishment_FuNiche")[3])
  m.est <- as.character(filter(vr.models, species == fc & response == "establishment_FuNiche")[5])
  
  #if(m.est == "site" & fc != "Daca" & sample.size.site >3) {
  if(!is.na(m.est)) {
    # lm
    glm.est <- glm(as.formula(paste(y.est, m.est, sep="~")), data=f.est, family = "binomial")
    
    # intercept and sigma
    vr[n,"establishment.funiche"] <- coef(glm.est)[1]
    vr.boot[n, match("establishment.funiche", colnames(vr)), 1:nboot] = coef(glm.est)[1]
    
    rm(m.est)
    rm(glm.est)
    print("establishment.funiche")
  } else {
    vr[n,"establishment.funiche"] <- 20
    vr.boot[n, match("establishment.funiche", colnames(vr)), 1:nboot] = 20
    print("establishment.funiche")
  }
  
  # 6.2 seedling establishment with competition----------------------------------
  d.est <- filter(d, focal.species == fc & 
                    (!is.na(planting.spring1) & planting.spring1=="yes") &
                    !is.na(survival.summer1.seedling) &
                    background.species %in% sps)
  
  # sample size
  sample.size.site <- table(d.est[,c("site")])[1]
  sample.size.sps <- nrow(d.est)
  
  if((bg == "site" | bg == "site_none") & sample.size.site > 2) {
    vr[n, "n.establishment.comp"] <- sample.size.site
    lm.est <- glm(survival.summer1.seedling ~ site, data=d.est, family="binomial")
  } else if(bg == "species" | bg == "species_none" | sample.size.site < 3) {
    vr[n, "n.establishment.comp"] <- sample.size.sps
    lm.est <- glm(survival.summer1.seedling ~ 1, data=d.est, family="binomial")
  }
  
  # coefs
  vr[n,"establishment.comp"] <- coef(lm.est)[1]
  vr.boot[n, match("establishment.comp", colnames(vr)), 1:nboot] = coef(lm.est)[1]
  
  rm(lm.est)
  print("establishment2") 
  
  # 7 recruit size--------------------------------------------
  b.rec <- filter(b, species == fc & greenhouse == "yes" & flowering == "no" & !is.na(drymass.g))
  sample.size.sp <- nrow(b.rec)
  vr[n, "n.seedling.size"] <- sample.size.sp
  
  lm.rec <- b.rec %>% 
    .$drymass.g %>%
    mean_sd
  vr[n, "seedling.size.mean"] = lm.rec[1] 
  vr[n, "seedling.size.sd"] = lm.rec[2]
  vr.boot[n, match("seedling.size.mean", colnames(vr)), 1:nboot] = lm.rec[1] 
  vr.boot[n, match("seedling.size.sd", colnames(vr)), 1:nboot] = lm.rec[2] 
  
  rm(lm.rec)
  print("recruit size")
  
  # 8 size range: L, U----------------------------------------
  d.siz <- filter(d, focal.species == fc &
                    !is.na(size.autumn0) &
                    background.species %in% sps)
  
  # sample size
  sample.size.site <- table(d.siz[,c("site")])[1]
  sample.size.sps <- nrow(d.siz)
  
  if(bg == "site" | bg == "site_none") {
    vr[n, "n.size.range"] <- sample.size.site
    lm.siz <- d.siz %>% 
      filter(site == st) %>%
      .$size.autumn0 %>%
      range
  } else if(bg == "species" | bg == "species_none") {
    vr[n, "n.size.range"] <- sample.size.sps
    lm.siz <- d.siz %>% 
      .$size.autumn0 %>%
      range
  }
  
  # greenhouse plants
  b.siz <- filter(b, species == fc & greenhouse == "yes" & flowering == "no" & !is.na(drymass.g))
  lm.siz.gr <- range(b.siz$drymass.g)
  
  vr[n, "L"] = min(lm.siz, lm.siz.gr, na.rm = TRUE)
  vr[n, "U"] = max(lm.siz, lm.siz.gr, na.rm = TRUE)
  vr.boot[n, match("L", colnames(vr)), 1:nboot] = as.numeric(vr[n, "L"])
  vr.boot[n, match("U", colnames(vr)), 1:nboot] = as.numeric(vr[n, "U"])
  
  rm(lm.siz)
  rm(lm.siz.gr)
  print("size range")
  
} # for loop

#****************************************************************************
# ** 3 Fill up vital rates gaps **------
#****************************************************************************
for(n in 1:nrow(vr)) {
  # n = 1
  vi <- vr[n,]
  fc <- as.character(vi$focal.species)
  bg <- as.character(vi$background.species)
  st <- as.character(vi$site)
  
  # print
  print(n)
  print(paste(fc,st,bg, sep="_"))
  
  # relevel focal species, background species and site
  d$focal.species <- relevel(d$focal.species, ref=fc)
  #d$background.species <- relevel(d$background.species, ref=bg)
  d$site <- relevel(d$site, ref=st)
  
  ## 1 survival -----------------------------------------------------
  
  ## 2 growth-----------------------------------------------------
  # SD of growth
  if(!is.na(vi$pair) & is.na(vi$growth.sd) ) {
    vr[n, "growth.sd"] <- 0.01
    vr.boot[n, match("growth.sd", colnames(vr)), 1:nboot] = 0.01
  }
  
  ## 3 flowering -----------------------------------------------------
  
  ## 4 fecundity -----------------------------------------------------
  # Armo no flowered plants
  if(!is.na(vi$pair) & is.na(vi$fecundity.int) & fc == "Armo") {
    # poisson
    vr[n, "fecundity.int"] <- -20 # this makes exp(-20) --> 0
    vr[n, "fecundity.slope"] <- 0
    vr.boot[n, match("fecundity.int", colnames(vr)), 1:nboot] = -20
    vr.boot[n, match("fecundity.slope", colnames(vr)), 1:nboot] = 0
    
    # linear
    vr[n, "fecundity.int.linear"] <- 0
    vr[n, "fecundity.slope.linear"] <- 0
    vr.boot[n, match("fecundity.int.linear", colnames(vr)), 1:nboot] = 0
    vr.boot[n, match("fecundity.slope.linear", colnames(vr)), 1:nboot] = 0
  }
  
  ## 5 germination ------
  # Daca germination establishment and estiablishment-FuNiche
  if(!is.na(fc) & fc == "Daca" & !is.na(vi$pair)) {
    if(st == "Les Posses") {
      vr[n, "germination.prob"] <- 0.055
      vr[n, "establishment.funiche"] <- 20 # logit(20) == 1
      vr.boot[n, match("germination.prob", colnames(vr)), 1:nboot] = 0.055
      vr.boot[n, match("establishment.funiche", colnames(vr)), 1:nboot] = 20
    }
    else if(st == "Solalex") {
      vr[n, "germination.prob"] <- 0.12
      vr[n, "establishment.funiche"] <- 20 # logit(20) == 1
      vr.boot[n, match("germination.prob", colnames(vr)), 1:nboot] = 0.12
      vr.boot[n, match("establishment.funiche", colnames(vr)), 1:nboot] = 20
    }
    else if(st == "Anzeindaz") {
      vr[n, "germination.prob"] <- 0.04
      vr[n, "establishment.funiche"] <- 20 # logit(20) == 1
      vr.boot[n, match("germination.prob", colnames(vr)), 1:nboot] = 0.04
      vr.boot[n, match("establishment.funiche", colnames(vr)), 1:nboot] = 20
    }
  }
} 

