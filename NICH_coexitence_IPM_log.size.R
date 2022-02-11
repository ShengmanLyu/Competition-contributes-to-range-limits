#*******************************************************************************
# IPM using log(size) -----
# 2022.01.13
# update: 
#*******************************************************************************
rm(list=ls())

library(readxl)
library(tidyverse)
library(lme4)
library(lmerTest)
library(MASS)
library(Hmisc)
library(ggpubr)
library(viridis)
library(MuMIn)
library(mgcv)

#****************************************************************************
# 1. Vital rates models diagnostics at the species level ------
#****************************************************************************
rm(list=ls())

# ** - 1.0 read data---------------------------------------
# load functions
source("NICH_coexistence_IPM_functions.R")

# NICH data
d <- read_excel("data_181920_v6.xlsx", na="NA",col_names = TRUE)
d
# remove seedling planted in autumn 2019 because mostly of them died: 16755
d <- d %>% filter(!(!is.na(note.size.autumn0) & note.size.autumn0 == "planting.autumn0" & year == 2020))
d

# ** - 1.1 size at time 0 -------------------------------------
fig.size0 <- d %>%
  filter(!is.na(size.autumn0) & !is.na(survival)) %>%
  filter((!is.na(if.same) & if.same==TRUE)) %>%
  filter(!is.na(site) & !is.na(focal.species) & !is.na(background.species)) %>%
  filter(!(focal.species %in% c("Anvu", "Trca"))) %>%
  filter(!(background.species %in% c("Anvu", "Trca"))) %>%
  mutate(focal.species = factor(focal.species, levels = c("Anal", "Armo","Asal", "Plal", "Poal","Seca", "Trba", "Brer", "Crbi", "Daca", "Melu", "Plla", "Potr", "Sapr"))) %>%  
  mutate(background.species = factor(background.species, levels = c("Anal","Asal", "Plal", "Poal","Seca", "Trba", "Brer", "Crbi", "Melu", "Plla", "Potr", "Sapr", "none"))) %>%  
  ggplot(aes(x=log(size.autumn0))) +
  geom_histogram() +
  facet_wrap(~focal.species, nrow=2, scales = "free") +
  scale_x_continuous(name = "log(Size at year t)") + 
  scale_y_continuous(name="Count")
fig.size0

# ** - 1.2 size at time 1 -------------------------------------
fig.size1 <- d %>%
  filter(!is.na(size.autumn0) & !is.na(survival)) %>%
  filter((!is.na(if.same) & if.same==TRUE)) %>%
  filter(!is.na(site) & !is.na(focal.species) & !is.na(background.species)) %>%
  filter(!(focal.species %in% c("Anvu", "Trca"))) %>%
  filter(!(background.species %in% c("Anvu", "Trca"))) %>%
  mutate(focal.species = factor(focal.species, levels = c("Anal", "Armo","Asal", "Plal", "Poal","Seca", "Trba", "Brer", "Crbi", "Daca", "Melu", "Plla", "Potr", "Sapr"))) %>%  
  mutate(background.species = factor(background.species, levels = c("Anal","Asal", "Plal", "Poal","Seca", "Trba", "Brer", "Crbi", "Melu", "Plla", "Potr", "Sapr", "none"))) %>%  
  ggplot(aes(x=log(size.autumn1))) +
  geom_histogram() +
  facet_wrap(~focal.species, nrow=2, scales = "free") +
  scale_x_continuous(name = "log(Size at year t + 1)") + 
  scale_y_continuous(name="Count")
fig.size1

# ** - 1.3 growth -------------------------------------
fig.growth <- d %>%
  filter(!is.na(size.autumn0) & !is.na(size.autumn1)) %>%
  filter((!is.na(if.same) & if.same==TRUE)) %>%
  filter(!is.na(site) & !is.na(focal.species) & !is.na(background.species)) %>%
  filter(!(focal.species %in% c("Anvu", "Trca"))) %>%
  filter(!(background.species %in% c("Anvu", "Trca"))) %>%
  mutate(focal.species = factor(focal.species, levels = c("Anal", "Armo","Asal", "Plal", "Poal","Seca", "Trba", "Brer", "Crbi", "Daca", "Melu", "Plla", "Potr", "Sapr"))) %>%  
  ggplot(aes(x=log(size.autumn0), y = log(size.autumn1), col=site)) +
  geom_point() +
  geom_rug(sides = "b", col="grey50") +
  facet_wrap(~focal.species, nrow=2, scales = "free") +
  scale_x_continuous(name = "log(Size at year t)") +
  scale_y_continuous(name="log(Size at year t+1)")
fig.growth

# plot
for(fc in c("Anal", "Armo","Asal", "Plal", "Poal","Seca", "Trba", "Brer", "Crbi", "Daca", "Melu", "Plla", "Potr", "Sapr")) {
  # fc = "Sapr"
  d %>%
    filter(!is.na(size.autumn0) & !is.na(size.autumn1)) %>%
    filter((!is.na(if.same) & if.same==TRUE)) %>%
    filter(!is.na(site) & !is.na(focal.species) & !is.na(background.species)) %>%
    filter(!(background.species %in% c("Anvu", "Trca"))) %>%
    #filter(background.species != "none") %>%
    filter(focal.species == fc) %>%
    plot.diag.lm("log(size.autumn0)"," log(size.autumn1)",data=., fc)
}

# linearity
linearity.gro <- d %>%
  filter(!is.na(size.autumn0) & !is.na(size.autumn1)) %>%
  filter((!is.na(if.same) & if.same==TRUE)) %>%
  filter(!is.na(site) & !is.na(focal.species) & !is.na(background.species)) %>%
  filter(!(background.species %in% c("Anvu", "Trca"))) %>%
  mutate(focal.species = factor(focal.species, levels = c("Anal", "Armo","Asal", "Plal", "Poal","Seca", "Trba", "Brer", "Crbi", "Daca", "Melu", "Plla", "Potr", "Sapr"))) %>%  
  #filter(background.species != "none") %>%
  #filter(focal.species == fc) %>%
  split(.$focal.species) %>%
  purrr::map_dfr(~linearity.lm("log(size.autumn0)"," log(size.autumn1)", data=.))
linearity.gro

# normality
normality.gro <- d %>%
  filter(!is.na(size.autumn0) & !is.na(size.autumn1)) %>%
  filter((!is.na(if.same) & if.same==TRUE)) %>%
  filter(!is.na(site) & !is.na(focal.species) & !is.na(background.species)) %>%
  filter(!(background.species %in% c("Anvu", "Trca"))) %>%
  mutate(focal.species = factor(focal.species, levels = c("Anal", "Armo","Asal", "Plal", "Poal","Seca", "Trba", "Brer", "Crbi", "Daca", "Melu", "Plla", "Potr", "Sapr"))) %>%  
  #filter(background.species != "none") %>%
  #filter(focal.species == fc) %>%
  split(.$focal.species) %>%
  purrr::map_dfr(~normality.lm("log(size.autumn0)"," log(size.autumn1)", data=.))
normality.gro

# constant variacne
cv.gro <- d %>%
  filter(!is.na(size.autumn0) & !is.na(size.autumn1)) %>%
  filter((!is.na(if.same) & if.same==TRUE)) %>%
  filter(!is.na(site) & !is.na(focal.species) & !is.na(background.species)) %>%
  filter(!(background.species %in% c("Anvu", "Trca"))) %>%
  mutate(focal.species = factor(focal.species, levels = c("Anal", "Armo","Asal", "Plal", "Poal","Seca", "Trba", "Brer", "Crbi", "Daca", "Melu", "Plla", "Potr", "Sapr"))) %>%  
  #filter(background.species != "none") %>%
  #filter(focal.species == fc) %>%
  split(.$focal.species) %>%
  purrr::map_dfr(~constant.var.lm("log(size.autumn0)"," log(size.autumn1)", data=.))
cv.gro

diag.gro <- cbind(linearity.gro, normality.gro, cv.gro)
diag.gro$species = c("Anal", "Armo","Asal", "Plal", "Poal","Seca", "Trba", "Brer", "Crbi", "Daca", "Melu", "Plla", "Potr", "Sapr")
diag.gro$vital.rate = "growth"
diag.gro

# ** - 1.4 fecundity -------------------------------------
# distributions
d %>%
  filter(!is.na(size.autumn0) & !is.na(size.autumn1)) %>%
  filter((!is.na(if.same) & if.same==TRUE)) %>%
  filter(!is.na(site) & !is.na(focal.species) & !is.na(background.species)) %>%
  filter(!(focal.species %in% c("Anvu", "Trca"))) %>%
  filter(!(background.species %in% c("Anvu", "Trca"))) %>%
  .$seeds.autumn0 %>% log() %>% hist()

fig.fecundity <- d %>%
  filter(!is.na(size.autumn0) & !is.na(size.autumn1)) %>%
  filter((!is.na(if.same) & if.same==TRUE)) %>%
  filter(!is.na(site) & !is.na(focal.species) & !is.na(background.species)) %>%
  filter(!(focal.species %in% c("Anvu", "Trca"))) %>%
  filter(!(background.species %in% c("Anvu", "Trca"))) %>%
  mutate(focal.species = factor(focal.species, levels = c("Anal", "Armo","Asal", "Plal", "Poal","Seca", "Trba", "Brer", "Crbi", "Daca", "Melu", "Plla", "Potr", "Sapr"))) %>%  
  ggplot(aes(x=log(size.autumn0), y = log(seeds.autumn0), col=site)) +
  geom_point() +
  geom_rug(sides = "b", col="grey50") +
  facet_wrap(~focal.species, nrow=2, scales = "free") +
  scale_x_continuous(name = "log(Size at year t)") +
  scale_y_continuous(name="log(Seed number)")
fig.fecundity

# plot
for(fc in c("Anal","Asal", "Plal", "Poal","Seca", "Trba", "Brer", "Crbi", "Daca", "Melu", "Plla", "Potr", "Sapr")) {
  #fc = "Anal"
  d %>%
    filter(!is.na(size.autumn0) & !is.na(seeds.autumn0) & seeds.autumn0 != 0) %>%
    filter((!is.na(if.same) & if.same==TRUE)) %>%
    filter(!is.na(site) & !is.na(focal.species) & !is.na(background.species)) %>%
    filter(!(background.species %in% c("Anvu", "Trca"))) %>%
    #filter(background.species != "none") %>%
    filter(focal.species == fc) %>%
    plot.diag.lm("log(size.autumn0)","log(seeds.autumn0)",data=., fc)
}
#dev.off()

# linearity
linearity.fec <- d %>%
  filter(!is.na(size.autumn0) & !is.na(seeds.autumn0) & seeds.autumn0 != 0) %>%  
  filter((!is.na(if.same) & if.same==TRUE)) %>%
  filter(!is.na(site) & !is.na(focal.species) & !is.na(background.species)) %>%
  filter(!(background.species %in% c("Anvu", "Trca"))) %>%
  filter(focal.species != "Armo") %>%
  mutate(focal.species = factor(focal.species, levels = c("Anal", "Asal", "Plal", "Poal","Seca", "Trba", "Brer", "Crbi", "Daca", "Melu", "Plla", "Potr", "Sapr"))) %>%  
  #filter(background.species != "none") %>%
  split(.$focal.species) %>%
  purrr::map_dfr(~linearity.lm("log(size.autumn0)"," log(seeds.autumn0)", data=.))
linearity.fec

# normality
normality.fec <- d %>%
  filter(!is.na(size.autumn0) & !is.na(seeds.autumn0) & seeds.autumn0 != 0) %>%
  filter((!is.na(if.same) & if.same==TRUE)) %>%
  filter(!is.na(site) & !is.na(focal.species) & !is.na(background.species)) %>%
  filter(!(background.species %in% c("Anvu", "Trca"))) %>%
  filter(focal.species != "Armo") %>%
  mutate(focal.species = factor(focal.species, levels = c("Anal","Asal", "Plal", "Poal","Seca", "Trba", "Brer", "Crbi", "Daca", "Melu", "Plla", "Potr", "Sapr"))) %>%  
  #filter(background.species != "none") %>%
  #filter(focal.species == fc) %>%
  split(.$focal.species) %>%
  purrr::map_dfr(~normality.lm("log(size.autumn0)", "log(seeds.autumn0)", data=.))
normality.fec

# constant variacne
cv.fec<- d %>%
  filter(!is.na(size.autumn0) & !is.na(seeds.autumn0) & seeds.autumn0 != 0) %>%
  filter((!is.na(if.same) & if.same==TRUE)) %>%
  filter(!is.na(site) & !is.na(focal.species) & !is.na(background.species)) %>%
  filter(!(background.species %in% c("Anvu", "Trca"))) %>%
  filter(focal.species != "Armo") %>%
  mutate(focal.species = factor(focal.species, levels = c("Anal","Asal", "Plal", "Poal","Seca", "Trba", "Brer", "Crbi", "Daca", "Melu", "Plla", "Potr", "Sapr"))) %>%  
  #filter(background.species != "none") %>%
  #filter(focal.species == fc) %>%
  split(.$focal.species) %>%
  purrr::map_dfr(~constant.var.lm("log(size.autumn0)","log(seeds.autumn0)", data=.))
cv.fec

diag.fec <- cbind(linearity.fec, normality.fec, cv.fec)
diag.fec$species = c("Anal", "Asal", "Plal", "Poal","Seca", "Trba", "Brer", "Crbi", "Daca", "Melu", "Plla", "Potr", "Sapr")
diag.fec$vital.rate = "fecundity"

# ** - 1.5 survival -------------------------------------
fig.survival <- d %>%
  filter(!is.na(size.autumn0) & !is.na(survival)) %>%
  filter((!is.na(if.same) & if.same==TRUE)) %>%
  filter(!is.na(site) & !is.na(focal.species) & !is.na(background.species)) %>%
  filter(!(focal.species %in% c("Anvu", "Trca"))) %>%
  filter(!(background.species %in% c("Anvu", "Trca"))) %>%
  mutate(focal.species = factor(focal.species, levels = c("Anal", "Armo","Asal", "Plal", "Poal","Seca", "Trba", "Brer", "Crbi", "Daca", "Melu", "Plla", "Potr", "Sapr"))) %>%  
  mutate(background.species = factor(background.species, levels = c("Anal","Asal", "Plal", "Poal","Seca", "Trba", "Brer", "Crbi", "Melu", "Plla", "Potr", "Sapr", "none"))) %>%  
  ggplot(aes(x=log(size.autumn0), y = survival, col = site)) +
  geom_jitter(height=0.1) +
  geom_rug(sides = "b", col="grey50") +
  facet_wrap(~focal.species, nrow=2, scales = "free") +
  scale_x_continuous(name = "log(Size at year t)") +
  scale_y_continuous(name="Prob (survival)")
fig.survival

# plot
for(fc in c("Anal","Armo","Asal", "Plal", "Poal","Seca", "Trba", "Brer", "Crbi", "Daca", "Melu", "Plla", "Potr", "Sapr")) {
  #fc = "Armo"
  d %>%
    filter(!is.na(size.autumn0) & !is.na(survival)) %>%
    filter((!is.na(if.same) & if.same==TRUE)) %>%
    filter(!is.na(site) & !is.na(focal.species) & !is.na(background.species)) %>%
    filter(!(background.species %in% c("Anvu", "Trca"))) %>%
    #filter(background.species != "none") %>%
    filter(focal.species == fc) %>%
    plot.diag.sur(data=., fc)
}
#dev.off()

# linearity
linearity.sur <- d %>%
  filter(!is.na(size.autumn0) & !is.na(survival)) %>%
  filter((!is.na(if.same) & if.same==TRUE)) %>%
  filter(!is.na(site) & !is.na(focal.species) & !is.na(background.species)) %>%
  filter(!(focal.species %in% c("Anvu", "Trca"))) %>%
  filter(!(background.species %in% c("Anvu", "Trca"))) %>%
  mutate(focal.species = factor(focal.species, levels = c("Anal", "Armo","Asal", "Plal", "Poal","Seca", "Trba", "Brer", "Crbi", "Daca", "Melu", "Plla", "Potr", "Sapr"))) %>%  
  mutate(background.species = factor(background.species, levels = c("Anal","Asal", "Plal", "Poal","Seca", "Trba", "Brer", "Crbi", "Melu", "Plla", "Potr", "Sapr", "none"))) %>%  
  #filter(background.species != "none") %>%
  split(.$focal.species) %>%
  purrr::map_dfr(~linearity.glm("log(size.autumn0)","survival", data=.))
linearity.sur$species = c("Anal", "Armo","Asal", "Plal", "Poal","Seca", "Trba", "Brer", "Crbi", "Daca", "Melu", "Plla", "Potr", "Sapr")
linearity.sur$vital.rate = "survival"
linearity.sur

# ** - 1.6 flowering -------------------------------------
fig.flowering <- d %>%
  filter(!is.na(size.autumn0) & !is.na(flowering.autumn0)) %>%
  filter((!is.na(if.same) & if.same==TRUE)) %>%
  filter(!is.na(site) & !is.na(focal.species) & !is.na(background.species)) %>%
  filter(!(focal.species %in% c("Anvu", "Trca"))) %>%
  filter(!(background.species %in% c("Anvu", "Trca"))) %>%
  mutate(focal.species = factor(focal.species, levels = c("Anal", "Armo","Asal", "Plal", "Poal","Seca", "Trba", "Brer", "Crbi", "Daca", "Melu", "Plla", "Potr", "Sapr"))) %>%  
  ggplot(aes(x=log(size.autumn0), y = flowering.autumn0, col = site)) +
  geom_jitter(height=0.1) +
  geom_rug(sides = "b", col="grey50") +
  facet_wrap(~focal.species, nrow=2, scales = "free") +
  scale_x_continuous(name = "log(Size at year t)") +
  scale_y_continuous(name="Prob (flowering)")
fig.flowering

# plot
for(fc in c("Anal", "Armo","Asal", "Plal", "Poal","Seca", "Trba", "Brer", "Crbi", "Daca", "Melu", "Plla", "Potr", "Sapr")) {
  # fc= "Armo"
  d %>%
    filter(!is.na(size.autumn0) & !is.na(flowering.autumn0)) %>%
    filter((!is.na(if.same) & if.same==TRUE)) %>%
    filter(!is.na(site) & !is.na(focal.species) & !is.na(background.species)) %>%
    filter(!(focal.species %in% c("Anvu", "Trca"))) %>%
    filter(!(background.species %in% c("Anvu", "Trca"))) %>%
    filter(focal.species == fc) %>%
    plot.diag.flo(data=., fc)
}
#dev.off()
  
data = d %>%
  filter(!is.na(size.autumn0) & !is.na(flowering.autumn0)) %>%
  filter((!is.na(if.same) & if.same==TRUE)) %>%
  filter(!is.na(site) & !is.na(focal.species) & !is.na(background.species)) %>%
  filter(!(focal.species %in% c("Anvu", "Trca"))) %>%
  filter(!(background.species %in% c("Anvu", "Trca"))) %>%
  filter(focal.species != "Armo") %>%
  mutate(focal.species = factor(focal.species, levels = c("Anal", "Asal", "Plal", "Poal","Seca", "Trba", "Brer", "Crbi", "Daca", "Melu", "Plla", "Potr", "Sapr"))) %>%  
  filter(focal.species == "Poal")
  
# linearity
linearity.flo <- d %>%
  filter(!is.na(size.autumn0) & !is.na(flowering.autumn0)) %>%
  filter((!is.na(if.same) & if.same==TRUE)) %>%
  filter(!is.na(site) & !is.na(focal.species) & !is.na(background.species)) %>%
  filter(!(focal.species %in% c("Anvu", "Trca"))) %>%
  filter(!(background.species %in% c("Anvu", "Trca"))) %>%
  mutate(focal.species = factor(focal.species, levels = c("Anal","Armo", "Asal", "Plal", "Poal","Seca", "Trba", "Brer", "Crbi", "Daca", "Melu", "Plla", "Potr", "Sapr"))) %>%  
  split(.$focal.species) %>%
  purrr::map_dfr(~linearity.glm("log(size.autumn0)", "flowering.autumn0", data=.))

linearity.flo$species = c("Anal","Armo", "Asal", "Plal", "Poal","Seca", "Trba", "Brer", "Crbi", "Daca", "Melu", "Plla", "Potr", "Sapr")
linearity.flo$vital.rate = "flowering"
linearity.flo

#****************************************************************************
# 2. Vital rates model selection ------
#****************************************************************************
# rm(list=ls())

# load functions
source("NICH_coexistence_IPM_functions.R")

# raw data: 17679
d <- read_excel("data_181920_v6.xlsx", na="NA",col_names = TRUE)
d 
# remove seedling planted in autumn 2019 because mostly of them died: 16755
d <- d %>% filter(!(!is.na(note.size.autumn0) & note.size.autumn0 == "planting.autumn0" & year == 2020))
d$site <- factor(d$site)
d

# log size
d$size.autumn0 = log(d$size.autumn0); hist(d$size.autumn0)
d$size.autumn1 = log(d$size.autumn1); hist(d$size.autumn1)

# focal species
sps <- c("Anal","Armo","Asal","Plal", "Poal", "Seca", "Trba",      # 7 alpine species 
         "Brer", "Crbi", "Daca","Melu" , "Plla", "Potr", "Sapr", # 7 lowland species
         "none")   

# candidate models without size
models.1 <- c("1", 
              "site",
              "background.species",
              "site + background.species",
              "site*background.species")

# candidate models with size
models.size <- c("1", 
                 "size.autumn0",
                 "site",
                 "background.species",
                 "size.autumn0 + site",
                 "size.autumn0 + background.species",
                 "site + background.species",
                 "size.autumn0 + site + background.species",
                 "size.autumn0*site + background.species",
                 "size.autumn0*background.species + site",
                 "size.autumn0 + site*background.species",
                 "size.autumn0 * site * background.species")

aic <- aic.best <- NULL
for(fc in sps[1:14]) {
  #fc <- "Sapr"
  
  #**************************************************
  # ** - 2.1 survival----
  d.sur <- d %>%
    filter(focal.before == fc &
             background.species %in% sps &
             !is.na(size.autumn0) & !is.na(survival) &
             !(!is.na(note.size.autumn0) & note.size.autumn0 == "planting.autumn0" & year == 2020) &
             (!is.na(if.same) & if.same==TRUE) )
  
  if(nrow(sample.pair(d.sur)) > 5) {
    # AIC comparison
    aic.sur <- modelcmp.glm(response = "survival", candidate = models.size, d.sur, family="binomial")
    aic.sur$response <- "survival"
    print("survival")
  }
  
  #**************************************************
  # ** - 2.2 growth----
  d.gro <- d %>%
    filter(focal.species == fc &
             !is.na(survival) & survival == 1 &
             background.species %in% sps &
             !is.na(size.autumn0) & !is.na(size.autumn1) &
             !(!is.na(note.size.autumn0) & note.size.autumn0 == "planting.autumn0" & year == 2020) &
             (!is.na(if.same) & if.same==TRUE) )
  
  if(nrow(sample.pair(d.gro)) > 5) {
    # AIC comparison
    aic.gro <- modelcmp.glm(response = "size.autumn1", candidate = models.size, data=d.gro, family= "gaussian")
    aic.gro$response <- "growth"
    print("growth")
  }
  
  #**************************************************
  # ** - 2.3 flowering----
  d.flo <- d %>%
    filter(focal.before == fc &
             background.species %in% sps &
             !is.na(size.autumn0) & !is.na(flowering.autumn0) &
             !(!is.na(note.size.autumn0) & note.size.autumn0 == "planting.autumn0" & year == 2020) &           
             (!is.na(if.same) & if.same==TRUE) ) 
  
  if(nrow(sample.pair(d.flo)) > 5) {
    # AIC comparison
    aic.flo <- modelcmp.glm(response = "flowering.autumn0", candidate = models.size, data=d.flo, family = "binomial")
    aic.flo$response <- "flowering"
    print("flowering")
  }
  
  #**************************************************
  # ** - 2.4 fecundity----
  d.fec <- d %>%
    filter(focal.before == fc &
             background.species %in% sps &
             !is.na(size.autumn0) & !is.na(seeds.autumn0) & seeds.autumn0 != 0 &
             !(!is.na(note.size.autumn0) & note.size.autumn0 == "planting.autumn0" & year == 2020) &           
             (!is.na(if.same) &if.same==TRUE) )
  
  if(nrow(d.fec) == 0) aic.fec <- data.frame(y = "seeds.autumn0", n = 0, models=NA, aic = NA, response = "fecundity")
  else if(nrow(sample.pair(d.fec)) > 5) {
    # AIC comparison
    aic.fec <- modelcmp.glm(response = "log(seeds.autumn0)", candidate = models.size, data=d.fec, family="gaussian")
    aic.fec$response <- "fecundity"
    print("fecundity")
  }
  
  #**************************************************
  # ** - 2.5 germination----
  # Using FuNiche data
  
  #**************************************************
  # ** - 2.6 seedling establishment----
  d.est <- d %>%
    filter(focal.species == fc & 
             (!is.na(planting.spring1) & planting.spring1=="yes") &
             !(!is.na(note.size.autumn0) & note.size.autumn0 == "planting.autumn0" & year == 2020) &   
             !is.na(survival.summer1.seedling) &
             background.species %in% sps)
  
  #if(nrow(sample.pair(d.est)) > 5 & length(unique(sample.pair(d.est)$site)) >1) {
  # AIC comparison
  aic.est <- modelcmp.glm(response = "survival.summer1.seedling", candidate = models.1, data=d.est, family="binomial")
  aic.est$response <- "establishment"
  print("establishment")
  #}
  
  #**************************************************
  # ** - 2.7 recruit size----
  # species level using greenhouse plants
  
  #**************************************************
  # ** - 2.8 size range: L, U----
  # only compare the biggest individuals
  d.siz <- d %>% 
    filter(focal.species == fc &
             !is.na(size.autumn1) &
             !(!is.na(note.size.autumn0) & note.size.autumn0 == "planting.autumn0" & year == 2020) &   
             background.species %in% sps)
  
  if(nrow(sample.pair(d.siz)) > 5) {
    # AIC comparison
    aic.siz <- modelcmp.glm(response = "size.autumn1", candidate = models.1, data=d.siz, family = "gaussian")
    aic.siz$response <- "size range"
    print("size.range")
  }
  
  # AIC comparison
  # combine aic results
  aic.fc <- bind_rows(aic.sur, aic.gro, aic.flo, aic.fec, aic.est, aic.siz)
  aic.fc$species <- fc
  
  aic.best.fc <- bind_rows(aic.sur[1,], aic.gro[1,], aic.flo[1,], aic.fec[1,], aic.est[1,], aic.siz[1,])
  aic.best.fc$species <- fc
  
  aic.sur <- aic.gro <- aic.flo <- aic.fec <- aic.est <- aic.siz <- NULL
  
  aic <- rbind(aic, aic.fc)
  aic.best <- rbind(aic.best, aic.best.fc)
}

aic
aic.best

#****************************************************************************
# ** - 2.9. Germination using FuNiche data ------
# Les Posses: 4; Solalex: 13; Anzeindaz: 19
#****************************************************************************
seeds <- read_excel("FunNiche Seed Counts.xlsx")
seeds
colnames(seeds) <- c("rep","Plal", "Poal", "Anal", "Trba", "Seca", "Asal", "Armo", "Plla", "Brer", "Melu", "Trca", "Crbi", "Potr", "Sapr")

funiche <- read_excel("FunNiche_data_2019_2020_dec2020.xlsx", na="NA")
funiche

# make sps names consistent
funiche$species <- dplyr::recode(funiche$species, 
                                 "crebie" = "Crbi",
                                 "tribad" = "Trba",
                                 "salpra" = "Sapr",
                                 "astalp" = "Asal",
                                 "tricam" = "Trca",
                                 "plalan" = "Plla",
                                 "sescae" = "Seca",
                                 "broere" = "Brer",
                                 "medlup" = "Melu",
                                 "arnmon" = "Armo",
                                 "antalp" = "Anal",
                                 "plaalp" = "Plal",
                                 "poaalp" = "Poal",
                                 "poatri" = "Potr")

unique(funiche$species)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# seeds per stick
seeds %>% 
  pivot_longer(cols = 2:15, names_to = "species", values_to = "seeds.per.stick") %>%
  ggplot(aes(x=species, y=seeds.per.stick)) + geom_boxplot() + geom_point()

sn <- seeds %>% 
  pivot_longer(cols = 2:15, names_to = "species", values_to = "seeds.per.stick") %>%
  split(.$species) %>%
  purrr::map_dbl(~mean(.$seeds.per.stick, na.rm=TRUE))
sn

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# germination and seedling establishment
funiche$seeds_autumn_18 <- as.numeric(NA)
funiche$germination_spring_19 <- as.numeric(NA)
funiche$establishement_summer_19 <- as.numeric(NA)

# calculate germination rates
for(a in 1:nrow(funiche)) {
  # a = 1
  fi <- funiche[a, ]
  sp <- fi$species
  
  # seeds per stick
  funiche[a, "seeds_autumn_18"] = as.numeric(sn[sp])
  
  # germination
  funiche[a, "germination_spring_19"] = as.numeric(funiche[a, "seedlings_total_summer_19"]/sn[sp])
  
  # establishment
  if(funiche[a, "presence_spring_19"] == 1) {
    if(funiche[a, "presence_summer_19"]==1) funiche[a, "establishement_summer_19"] = 1
    else if(funiche[a, "presence_summer_19"]==0) funiche[a, "establishement_summer_19"] = 0
  }
  print(a)
}

# plot
funiche %>% 
  filter(site %in% c(4,13,19)) %>%
  ggplot(aes(y=log(germination_spring_19), x=factor(elevation), col=factor(elevation))) + geom_boxplot()+ geom_jitter(height=0.1) + facet_wrap(~species)

# germination against elevation for all sites
funiche %>% 
  ggplot(aes(y=log(germination_spring_19), x=elevation, col=elevation)) + geom_jitter(height=0.1) + facet_wrap(~species, scales="free")

# plot
funiche %>% 
  filter(site %in% c(4,13,19)) %>%
  ggplot(aes(y=log(establishement_summer_19), x=factor(elevation), col=factor(elevation))) + geom_boxplot()+ geom_jitter(height=0.1) + facet_wrap(~species)

# check
funiche %>% 
  filter(site %in% c(13,19)) %>%
  filter(species == "Seca") %>%
  .$germination_spring_19 %>%
  mean()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# species or site level?
# germination
funiche$site <- as.character(funiche$site)
sps <- c("Anal","Armo","Asal","Plal", "Poal", "Seca", "Trba",        # 7 alpine species 
         "Brer", "Crbi", "Trca","Melu" , "Plla", "Potr", "Sapr") # 7 lowland species)
model.site <-  c("1", "site")

# distribution of germination rates
funiche %>%
  filter(species %in% sps) %>%
  #mutate(germination_spring_19 = ifelse(germination_spring_19 == 0, 0.001, ) 
  filter(site %in% c("4","13","19")) %>%
  ggplot(aes(x=germination_spring_19)) +
  geom_histogram() +
  facet_wrap(~species)

aic.ger <- aic.ger.best <- NULL
for(sp in sps) {
  #sp = "Anal"
  aic.germination <- funiche %>%
    filter(species == sp &
             site %in% c("4","13","19")) %>%
    modelcmp.glm(data=., response = "germination_spring_19", candidate = model.site,family="gaussian")
  aic.germination$response = "germination"
  aic.germination$species = sp
  aic.ger <- rbind(aic.ger, aic.germination)
  aic.ger.best <- rbind(aic.ger.best, aic.germination[1,])
}
aic.ger
aic.ger.best

# establishment
sps <- c("Anal","Asal","Plal", "Poal", "Trba","Armo", "Seca",        # 5 alpine species, no Armo and Seca
         "Brer", "Crbi", "Trca","Melu" , "Plla", "Potr", "Sapr") # 7 lowland species)
model.site <-  c("1", "site")

aic.est <- aic.est.best <- NULL
for(sp in sps) {
  # sp = "Seca"
  if(sp %in% c("Armo", "Seca")) aic.establish <- data.frame(y="establishement_summer_19", n =3, models= model.site, aic=NA, aic.delta = NA, aic.weight= NA, response= "establishment_FuNiche", species=sp)
  else {
    aic.establish <- funiche %>%
      filter(species == sp &
               !is.na(establishement_summer_19) &
               site %in% c("4","13","19")) %>%
      modelcmp.glm(response = "establishement_summer_19", candidate = model.site,  family="binomial")
    aic.establish$response = "establishment_FuNiche"
    aic.establish$species = sp
  }
  aic.est <- rbind(aic.est, aic.establish)
  aic.est.best <- rbind(aic.est.best, aic.establish[1,])
}
aic.est
aic.est.best

aic.funiche <- bind_rows(aic.ger, aic.est)
aic.best.funiche <- bind_rows(aic.ger.best, aic.est.best)

#****************************************************************************
# 3. Vital rates estimation ------
#****************************************************************************

#****************************************************************************
# ** - 3.1 Vital rates models: fitting and bootstrap ** ------
#****************************************************************************

rm(list=ls())
# NICH data
# load functions
source("NICH_coexistence_IPM_functions.R")

# raw data: 17679
d <- read_excel("data_181920_v6.xlsx", na="NA",col_names = TRUE)
d 
# remove seedling planted in autumn 2019 because mostly of them died
# 16755
d <- d %>% filter(!(!is.na(note.size.autumn0) & note.size.autumn0 == "planting.autumn0" & year == 2020))
d$site <- factor(d$site)
d

# log size
d$size.autumn0 = log(d$size.autumn0); hist(d$size.autumn0)
d$size.autumn1 = log(d$size.autumn1); hist(d$size.autumn1)

# convert to factor in order to relevel them
d$site <- factor(d$site, levels=c("Les Posses", "Solalex", "Anzeindaz"))
d$focal.species <- factor(d$focal.species)
d$focal.before <- factor(d$focal.before)
d$background.species <- factor(d$background.species)
d

# FuNiche data for germination and establishment
funiche <- read_csv("FunNiche_data_2019_2020_dec2020_20210114.csv")
f <- filter(funiche, site %in% c(4,13,19))
f$site <- as.character(f$site)
f$site <- dplyr::recode(f$site, "4" = "Les Posses", "13" = "Solalex", "19"="Anzeindaz")
f$species <- factor(f$species)
f$site <- factor(f$site)
f

# biomass of greenhouse seedling
g <- read_excel("biomass data for correlation_2020.xlsx",na="NA", col_names = TRUE)
g$greenhouse <- ifelse(g$site == "greenhouse" | g$site == "chngreenhouse", "yes", "no")
g$flowering <- ifelse(g$stem.number == 0, "no", "yes")
g

# vital rates
vr <- read_excel("vital.rates.xlsx", col_names = TRUE, na="NA")
vr <- vr[-1,]
vr

# new data frame to hold bootstrapped vital rates
nboot=999
vr.boot <- array(NA, dim = c(nrow(vr), ncol(vr), nboot))

# best vital rate models
vr.models <- read_excel("model selection.xlsx", na="NA") # best.model_20210105
vr.models
table(vr.models$response)

# species
sps <- c("Anal","Armo","Asal","Plal", "Poal", "Seca", "Trba",        # 7 alpine species 
         "Brer", "Crbi", "Daca","Melu" , "Plla", "Potr", "Sapr", # 7 lowland species
         "none")
sps

for(n in 1:nrow(vr)) {
  # n = 1
  vi <- vr[n,]
  fc <- as.character(vi$focal.species)
  bg <- as.character(vi$background.species)
  st <- as.character(vi$site)
  
  # skip when no pair
  if(is.na(vi$pair) | vi$pair == "no" | bg == "site" | bg == "site_none" | bg == "species" | bg == "species_none") next
  
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
  
  #**************************************************
  #**** - 3.1.1 survival----
  d.sur <- filter(d, focal.before == fc &
                    background.species %in% sps &
                    !is.na(size.autumn0) & !is.na(survival) &
                    (!is.na(if.same) & if.same==TRUE) )
  
  sample.size.site <- sum(table(d.sur[,c("background.species", "site")])[,1])
  sample.size.pair <- table(d.sur[,c("background.species", "site")])[1,1]
  
  y.sur <- as.character(filter(vr.models, species == fc & response == "survival")[3])
  m.sur <- as.character(filter(vr.models, species == fc & response == "survival")[5])
  
  # glm
  glm.sur <- glm(as.formula(paste(y.sur, m.sur, sep="~")), family="binomial", data=d.sur)
  
  if(m.sur == "site") {
    # bootstrap
    boot.sur <- mvrnorm(nboot,mu=coef(glm.sur)[1], Sigma=vcov(glm.sur)[1])
    
    vr[n, "n.survival"] <- sample.size.site
    vr[n, "surv.int"] <- coef(glm.sur)[1]
    vr[n, "surv.slope"] <- 0 
    vr[n, "surv.sd.global"] <- sigma(glm.sur)
    
    vr.boot[n, match("surv.int", colnames(vr)), 1:nboot] = boot.sur[,1]
    vr.boot[n, match("surv.slope", colnames(vr)), 1:nboot] = 0
  }
  else {
    # intercept, slope and sigma
    vr[n, "n.survival"] <- sample.size.pair
    vr[n, "surv.int"] <- coef(glm.sur)[1]
    vr[n, "surv.slope"] <- coef(glm.sur)[2]
    vr[n, "surv.sd.global"] <- sigma(glm.sur)
    
    # bootstrap
    sur.boot <- mvrnorm(nboot,mu=coef(glm.sur)[1:2], Sigma=vcov(glm.sur)[1:2,1:2])
    vr.boot[n, match("surv.int", colnames(vr)), 1:nboot] = sur.boot[,1]
    vr.boot[n, match("surv.slope", colnames(vr)), 1:nboot] = sur.boot[,2]
    
  }
  
  rm(m.sur)
  rm(glm.sur)
  rm(sur.boot)
  print("survival")
  
  #**************************************************
  #**** - 3.1.2 growth----
  # also try non-linear growth
  d.gro <- filter(d, focal.species == fc &
                    background.species %in% sps &
                    !is.na(survival) & survival == 1 &
                    !is.na(size.autumn0) & !is.na(size.autumn1) &
                    (!is.na(if.same) & if.same==TRUE) )
  # sample size
  sample.size.bg <- sum(table(d.gro[,c("background.species", "site")])[1,])
  sample.size.pair <- table(d.gro[,c("background.species", "site")])[1,1]
  vr[n, "n.growth"] <- sample.size.pair
  
  # model
  y.gro <- as.character(filter(vr.models, species == fc & response == "growth")[3])
  m.gro <- as.character(filter(vr.models, species == fc & response == "growth")[5])
  
  # growth model
  lm.gro <- lm(as.formula(paste(y.gro, m.gro, sep="~")), data=d.gro)
  
  # intrcept, slope and sigma
  vr[n, "growth.int"] <- coef(lm.gro)[1]
  vr[n, "growth.slope"] <- coef(lm.gro)[2]
  vr[n, "growth.sd.global"] <- sigma(lm.gro)
  
  # bootstrap
  boot.gro <- mvrnorm(nboot,mu=coef(lm.gro)[1:2], Sigma=vcov(lm.gro)[1:2,1:2])
  vr.boot[n, match("growth.int", colnames(vr)), 1:nboot] = boot.gro[,1]
  vr.boot[n, match("growth.slope", colnames(vr)), 1:nboot] = boot.gro[,2]
  
  # SD of residual at BG level
  if(m.gro == "size.autumn0 + background.species")  {
    if(sample.size.bg > 2) {
      d.gro.sd <- filter(d.gro, background.species == bg)
      resi.gro <- d.gro.sd$size.autumn1 - (as.numeric(vr[n, "growth.int"]) + d.gro.sd$size.autumn0 * as.numeric(vr[n, "growth.slope"]))
      vr[n,"growth.sd"] <- sd(resi.gro)
      for(b in 1:nboot) {
        # b = 1
        resi.gro.b <- d.gro.sd$size.autumn1 - (as.numeric(vr.boot[n,match("growth.int", colnames(vr)), b]) + d.gro.sd$size.autumn0 * as.numeric(vr.boot[n, match("growth.slope", colnames(vr)), b]))
        vr.boot[n, match("growth.sd", colnames(vr)), b] = sd(resi.gro.b)
      }
    }
    else{
      d.gro.sd <- filter(d.gro, site == st)
      lm.gro.sd <- lm(size.autumn1 ~ size.autumn0, data=d.gro.sd)
      resid.gro <- resid(lm.gro.sd)
      
      vr[n,"n.growth.sd"] <- nrow(d.gro.sd)
      vr[n,"growth.sd"] <- sd(resid.gro)
    }
  }
  # sd of growth residual at pair level
  else {
    if(sample.size.pair > 2) {
      d.gro.sd <- filter(d.gro, site == st & background.species == bg)
      resi.gro <- d.gro.sd$size.autumn1 - (as.numeric(vr[n, "growth.int"]) + d.gro.sd$size.autumn0 * as.numeric(vr[n, "growth.slope"]))
      vr[n,"growth.sd"] <- sd(resi.gro)
      for(b in 1:nboot) {
        # b = 1
        resi.gro.b <- d.gro.sd$size.autumn1 - (as.numeric(vr.boot[n,match("growth.int", colnames(vr)), b]) + d.gro.sd$size.autumn0 * as.numeric(vr.boot[n, match("growth.slope", colnames(vr)), b]))
        vr.boot[n, match("growth.sd", colnames(vr)), b] = sd(resi.gro.b)
      }
    }
    else{
      d.gro.sd <- filter(d.gro, site == st)
      lm.gro.sd <- lm(size.autumn1 ~ size.autumn0, data=d.gro.sd)
      resid.gro <- resid(lm.gro.sd)
      
      vr[n,"n.growth.sd"] <- nrow(d.gro.sd)
      vr[n,"growth.sd"] <- sd(resid.gro)
    }
  }
  
  # SD of residual at site level independent of size
  # also use size dependent and species level
  
  rm(m.gro)
  rm(lm.gro)
  rm(boot.gro)
  print("growth") 
  
  #**************************************************
  # **** - 3.1.3 flowering----
  d.flo <- filter(d, focal.before == fc &
                    background.species %in% sps &
                    !is.na(size.autumn0) & !is.na(flowering.autumn0) &
                    (!is.na(if.same) & if.same==TRUE) ) 
  
  sample.size.pair <- table(d.flo[,c("background.species", "site")])[1,1]
  vr[n, "n.flowering"] <- sample.size.pair
  
  y.flo <- as.character(filter(vr.models, species == fc & response == "flowering")[3])
  m.flo <- as.character(filter(vr.models, species == fc & response == "flowering")[5])
  
  # when the model is not NA
  if(!is.na(m.flo)) {
    # glm
    glm.flo <- glm(as.formula(paste(y.flo, m.flo, sep="~")), family="binomial", data=d.flo)
    
    # intecept and slope
    vr[n, "flowering.int"] <- coef(glm.flo)[1]  
    vr[n, "flowering.slope"] <- coef(glm.flo)[2] 
    
    # bootstrap
    boot.flo <- mvrnorm(nboot,mu=coef(glm.flo)[1:2], Sigma=vcov(glm.flo)[1:2,1:2])
    vr.boot[n, match("flowering.int", colnames(vr)), 1:nboot] = boot.flo[,1]
    vr.boot[n, match("flowering.slope", colnames(vr)), 1:nboot] = boot.flo[,2]
    
    rm(m.flo)
    rm(glm.flo)
    rm(boot.flo)
    print("flowering")
  }
  # Armo don't flower
  else { 
    vr[n, "flowering.int"] <- -400 # logit(-400) == 0 in R
    vr[n, "flowering.slope"] <- 0
    vr.boot[n, match("flowering.int", colnames(vr)), 1:nboot] = -400
    vr.boot[n, match("flowering.slope", colnames(vr)), 1:nboot] = 0
    print("Armo")
  }
  
  #**************************************************
  # **** - 3.1.4. fecundity ----
  d.fec <- filter(d, focal.before == fc &
                    background.species %in% sps &
                    !is.na(size.autumn0) & !is.na(seeds.autumn0) & seeds.autumn0 != 0 &
                    (!is.na(if.same) &if.same==TRUE) )
  
  # sample size
  sample.size.pair <- table(d.fec[,c("background.species", "site")])[1,1]
  sample.size.site <- table(d.fec[,"site"])[1]
  sample.size.bg <- table(d.fec[,"background.species"])[1]
  
  vr[n, "n.fecundity"] <- sample.size.pair
  
  y.fec <- as.character(filter(vr.models, species == fc  & response == "fecundity")[3])
  m.fec <- as.character(filter(vr.models, species == fc &  response == "fecundity")[5])
  
  # skip when model is NA
  if(!is.na(m.fec)) {
    # glm
    lm.fec <- glm(as.formula(paste(y.fec, m.fec, sep="~")), family="gaussian", data=d.fec)
    
    # coefs
    vr[n, "fecundity.int"] <- coef(lm.fec)[1]
    vr[n, "fecundity.slope"] <- coef(lm.fec)[2]
    
    # bootstrap
    boot.fec.lm <- mvrnorm(nboot,mu=coef(lm.fec)[1:2], Sigma=vcov(lm.fec)[1:2,1:2])
    vr.boot[n, match("fecundity.int", colnames(vr)), 1:nboot] = boot.fec.lm[,1]
    vr.boot[n, match("fecundity.slope", colnames(vr)), 1:nboot] = boot.fec.lm[,2]
    
    rm(m.fec)
    rm(lm.fec)
    rm(boot.fec.lm)
    print("fecunidty")
  }
  else { # Armo produce no seeds
    vr[n, "fecundity.int"] <- -800 # exp(-800) == 0 in R
    vr[n, "fecundity.slope"] <- 0
    vr.boot[n, match("fecundity.int", colnames(vr)), 1:nboot] = -800
    vr.boot[n, match("fecundity.slope", colnames(vr)), 1:nboot] = 0
  }
  
  #**************************************************
  # **** - 3.1.5 germination (no bootstrap) ----
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
  # **** - 3.1.6.1 seedling establishment of FuNiche (no bootstrap) ----
  f.est <- filter(f, species == fc & !is.na(establishement_summer_19))
  
  sample.size.site <- table(f.est[,"site"])[1]
  sample.size.sp <- nrow(f.est)
  vr[n,"n.establishment.funiche"] <- sample.size.site
  
  y.est <- as.character(filter(vr.models, species == fc & response == "establishment_FuNiche")[3])
  m.est <- as.character(filter(vr.models, species == fc & response == "establishment_FuNiche")[5])
  
  if(!is.na(m.est)) {
    # lm
    glm.est <- glm(as.formula(paste(y.est, m.est, sep="~")), data=f.est, family = "binomial")
    
    # intercept
    vr[n,"establishment.funiche.int"] = coef(glm.est)[1]
    vr.boot[n, match("establishment.funiche.int", colnames(vr)), 1:nboot] = coef(glm.est)[1]
    
    rm(m.est)
    rm(glm.est)
    print("establishment.funiche")
  } else{
    vr[n,"establishment.funiche.int"] = 20 # exp(20) == 1 in R
    vr.boot[n, match("establishment.funiche.int", colnames(vr)), 1:nboot] = 20
    print("establishment.funiche")
  }
  
  #**************************************************
  # **** - 3.1.6.2 seedling establishment with competition_best model (no bootstrap) ----
  d.est <- filter(d, focal.species == fc & 
                    (!is.na(planting.spring1) & planting.spring1=="yes") &
                    !is.na(survival.summer1.seedling) &
                    background.species %in% sps)
  
  # remove Poal and Potr because of different bias in site
  if(fc == "Armo") d.est <- filter(d.est, !background.species %in% c("Poal", "Potr"))
  
  sample.size.pair <- table(d.est[,c("background.species", "site")])[1,1]
  vr[n,"n.establishment.comp"] <- sample.size.pair
  
  y.est <- as.character(filter(vr.models, species == fc & response == "establishment")[3])
  m.est <- as.character(filter(vr.models, species == fc & response == "establishment")[5])

  # glm
  glm.est2 <- glm(as.formula(paste(y.est, m.est, sep="~")), data=d.est, family="binomial")
  
  # intercept
  vr[n,"establishment.comp.int"] <- coef(glm.est2)[1]
  vr.boot[n, match("establishment.comp.int", colnames(vr)), 1:nboot] = coef(glm.est2)[1]
  
  rm(m.est)
  rm(glm.est2)
  print("establishment12") 
  
  #**************************************************
  # **** - 3.1.7 recruit size (no bootstrap) ----
  mean_sd <- function(x) {c(mean(x), sd(x))}

  # model
  g.rec <- filter(g, species == fc & greenhouse == "yes" & flowering == "no" & !is.na(drymass.g))
  sample.size.sp <- nrow(g.rec)
  vr[n, "n.seedling.size"] <- sample.size.sp
  
  # mean and sd
  lm.rec <- g.rec %>% 
    .$drymass.g %>%
    log() %>%
    mean_sd()
  
  vr[n, "seedling.size.mean"] = lm.rec[1] 
  vr[n, "seedling.size.sd"] = lm.rec[2]
  vr.boot[n, match("seedling.size.mean", colnames(vr)), 1:nboot] = lm.rec[1] 
  vr.boot[n, match("seedling.size.sd", colnames(vr)), 1:nboot] = lm.rec[2] 
  rm(lm.rec)
  print("recruit size")
  
  #**************************************************
  # **** - 3.1.8 size range: L, U (no bootstrap) ----
  # only compare the biggest individuals
  # size range for species with limited size?
  d.siz <- filter(d, focal.species == fc &
                    !is.na(size.autumn0) &
                    background.species %in% sps)
  
  g.siz <- filter(g, species == fc & greenhouse == "yes" & flowering == "no" & !is.na(drymass.g))
  lm.siz.gr <- range(log(g.siz$drymass.g))
  
  sample.size.pair <- table(d.siz[,c("background.species", "site")])[1,1]
  sample.size.site <- table(d.siz[,"site"])[1]
  sample.size.bg <- table(d.siz[,"background.species"])[1]
  
  # model
  m.siz <- as.character(filter(vr.models, species == fc & response == "size range")[5])
  
  if(sample.size.pair > 2) {
    lm.siz <- d.siz %>% 
      filter(background.species == bg & site == st) %>%
      .$size.autumn0 %>%
      range()
    
    vr[n, "n.size.range"] <- sample.size.pair
    vr[n, "L"] = min(c(lm.siz, lm.siz.gr), na.rm = TRUE)
    if(is.infinite(as.numeric(vr[n, "L"]))) stop() 
    vr[n, "U"] = max(c(lm.siz, lm.siz.gr), na.rm = TRUE)
    if(is.infinite(as.numeric(vr[n, "U"]))) stop()
    vr.boot[n, match("L", colnames(vr)), 1:nboot] = as.numeric(vr[n, "L"])
    vr.boot[n, match("U", colnames(vr)), 1:nboot] = as.numeric(vr[n, "U"])
    
    rm(lm.siz.gr); rm(lm.siz)
  }
  # size range for each background species
  else {
    lm.siz <- d.siz %>% 
      filter(site == st) %>%
      .$size.autumn0 %>%
      range()
    
    vr[n, "n.size.range"] <- sample.size.site
    vr[n, "L"] = min(c(lm.siz, lm.siz.gr), na.rm = TRUE)
    if(is.infinite(as.numeric(vr[n, "L"]))) stop() 
    vr[n, "U"] = max(c(lm.siz, lm.siz.gr), na.rm = TRUE)
    if(is.infinite(as.numeric(vr[n, "U"]))) stop()
    vr.boot[n, match("L", colnames(vr)), 1:nboot] = as.numeric(vr[n, "L"])
    vr.boot[n, match("U", colnames(vr)), 1:nboot] = as.numeric(vr[n, "U"])
    
    rm(lm.siz.gr); rm(lm.siz)
  }
  print("size range")
} # for loop

#****************************************************************************
# ** - 3.2 Species and site level vital rates ------
#****************************************************************************
for(n in 1:nrow(vr)) {
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
  d$site <- relevel(d$site, ref=st)
  
  # relevel FuNiche data where Daca is absent
  if(fc != "Daca") {
    f$species <- relevel(f$species, ref=fc)
    f$site <- relevel(f$site, ref=st)
  }
  
  # site level with/without non-competition focals
  if(bg=="site" | bg == "species") { 
    sps <- c("Anal","Armo","Asal","Plal", "Poal", "Seca", "Trba",        # 7 alpine species 
             "Brer", "Crbi", "Daca","Melu" , "Plla", "Potr", "Sapr") # 7 lowland species 
  } else if(bg == "site_none" | bg == "species_none") { 
    sps <- c("Anal","Armo","Asal","Plal", "Poal", "Seca", "Trba",        # 7 alpine species 
             "Brer", "Crbi", "Daca","Melu" , "Plla", "Potr", "Sapr", # 7 lowland species
             "none")
  }
  
  ## **** - 3.2.1 survival-------------------------------------------------
  d.sur <- filter(d, focal.before == fc &
                    background.species %in% sps &
                    !is.na(size.autumn0) & !is.na(survival) &
                    (!is.na(if.same) & if.same==TRUE) )
  
  sample.size.site <- table(d.sur[,c("site")])[1]
  sample.size.sps <- nrow(d.sur)
  
  if(bg == "site" | bg == "site_none") {
    vr[n, "n.survival"] <- sample.size.site
    glm.sur <- glm(survival ~ size.autumn0 * site, family="binomial", data=d.sur)
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
  
  ## **** - 3.2.2 growth-----------------------------------------------------
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
    lm.gro <- lm(size.autumn1~size.autumn0 * site, data=d.gro)
    
    # bootstrap
    boot.gro <- mvrnorm(nboot,mu=coef(lm.gro)[1:2], Sigma=vcov(lm.gro)[1:2,1:2])
    
    # coefs
    vr[n, "growth.int"] <- coef(lm.gro)[1]
    vr[n, "growth.slope"] <- coef(lm.gro)[2]
    vr.boot[n, match("growth.int", colnames(vr)), 1:nboot] = boot.gro[,1]
    vr.boot[n, match("growth.slope", colnames(vr)), 1:nboot] = boot.gro[,2]
    
    # sigma of mean model
    d.gro.sd <- filter(d.gro, site == st)
    resi.gro <- d.gro.sd$size.autumn1 - (as.numeric(vr[n, "growth.int"]) + d.gro.sd$size.autumn0 * as.numeric(vr[n, "growth.slope"]))
    vr[n,"growth.sd"] <- sd(resi.gro)
    
    # sigma of bootstraps
    for(b in 1:nboot) {
      # b = 1
      resi.gro.b <- d.gro.sd$size.autumn1 - (as.numeric(vr.boot[n,match("growth.int", colnames(vr)), b]) + d.gro.sd$size.autumn0 * as.numeric(vr.boot[n, match("growth.slope", colnames(vr)), b]))
      vr.boot[n, match("growth.sd", colnames(vr)), b] = sd(resi.gro.b)
    }
  } else if(bg == "species" | bg == "species_none") {
    vr[n, "n.growth"] <- sample.size.sps
    lm.gro <- lm(size.autumn1~size.autumn0, data=d.gro)
    
    # bootstrap
    boot.gro <- mvrnorm(nboot,mu=coef(lm.gro)[1:2], Sigma=vcov(lm.gro)[1:2,1:2])
    
    # coefs and sigma
    vr[n, "growth.int"] <- coef(lm.gro)[1]
    vr[n, "growth.slope"] <- coef(lm.gro)[2]
    vr[n,"growth.sd"] <- sd(resid(lm.gro))
    vr.boot[n, match("growth.int", colnames(vr)), 1:nboot] = boot.gro[,1]
    vr.boot[n, match("growth.slope", colnames(vr)), 1:nboot] = boot.gro[,2]
    vr.boot[n, match("growth.sd", colnames(vr)), 1:nboot] = sd(resid(lm.gro))
  }
  
  rm(lm.gro)
  rm(boot.gro)
  print("growth")
  
  # **** - 3.2.3 flowering-----------------------------------------------------
  d.flo <- filter(d, focal.before == fc &
                    background.species %in% sps &
                    !is.na(size.autumn0) & !is.na(flowering.autumn0) &
                    (!is.na(if.same) & if.same==TRUE) ) 
  
  # sample size
  sample.size.site <- table(d.flo[,c("site")])[1]
  sample.size.sps <- nrow(d.flo)
  
  if(fc != "Armo" & (bg == "site" | bg == "site_none")) {
    vr[n, "n.flowering"] <- sample.size.site
    glm.flo <- glm(flowering.autumn0 ~ size.autumn0 * site, family="binomial", data=d.flo)
    
    # bootstrap
    boot.flo <- mvrnorm(nboot,mu=coef(glm.flo)[1:2], Sigma=vcov(glm.flo)[1:2,1:2])
    
    # intercept and sigma
    vr[n, "flowering.int"] <- coef(glm.flo)[1]
    vr[n, "flowering.slope"] <- coef(glm.flo)[2]
    vr.boot[n, match("flowering.int", colnames(vr)), 1:nboot] = boot.flo[,1]
    vr.boot[n, match("flowering.slope", colnames(vr)), 1:nboot] = boot.flo[,2]
    
    rm(glm.flo)
    rm(boot.flo)
    
  } else if(fc != "Armo" & (bg == "species" | bg == "species_none")) {
    vr[n, "n.flowering"] <- sample.size.sps
    glm.flo <- glm(flowering.autumn0 ~ size.autumn0, family="binomial", data=d.flo)
    
    # bootstrap
    boot.flo <- mvrnorm(nboot,mu=coef(glm.flo)[1:2], Sigma=vcov(glm.flo)[1:2,1:2])
    
    # intercept and sigma
    vr[n, "flowering.int"] <- coef(glm.flo)[1]
    vr[n, "flowering.slope"] <- coef(glm.flo)[2]
    vr.boot[n, match("flowering.int", colnames(vr)), 1:nboot] = boot.flo[,1]
    vr.boot[n, match("flowering.slope", colnames(vr)), 1:nboot] = boot.flo[,2]
    
    rm(glm.flo)
    rm(boot.flo)
    
  }
  else if(fc == "Armo") {
    # intercept and sigma
    vr[n, "flowering.int"] <- -400
    vr[n, "flowering.slope"] <- 0
    vr.boot[n, match("flowering.int", colnames(vr)), 1:nboot] = -400
    vr.boot[n, match("flowering.slope", colnames(vr)), 1:nboot] = 0
  }
  
  print("flowering")
  
  # **** - 3.2.4 fecundity -------------------------------------------
  d.fec <- filter(d, focal.before == fc &
                    background.species %in% sps &
                    !is.na(size.autumn0) & !is.na(seeds.autumn0) & seeds.autumn0 != 0 &
                    (!is.na(if.same) &if.same==TRUE) )
  
  # sample size
  sample.size.site <- table(d.fec[,c("site")])[1]
  sample.size.sps <- nrow(d.fec)
  
  if(fc != "Armo" & (bg == "site" | bg == "site_none")) {
    vr[n, "n.fecundity"] <- sample.size.site
    glm.fec <- glm(log(seeds.autumn0) ~ size.autumn0 * site, family="gaussian", data=d.fec)
    
    # bootstrap
    boot.fec <- mvrnorm(nboot,mu=coef(glm.fec)[1:2], Sigma=vcov(glm.fec)[1:2,1:2])
    
    # coefs and sigma
    vr[n, "fecundity.int"] <- coef(glm.fec)[1]
    vr[n, "fecundity.slope"] <- coef(glm.fec)[2]
    vr.boot[n, match("fecundity.int", colnames(vr)), 1:nboot] = boot.fec[,1]
    vr.boot[n, match("fecundity.slope", colnames(vr)), 1:nboot] = boot.fec[,2]
    rm(glm.fec)
    rm(boot.fec)
    
  } else if(fc != "Armo" & (bg == "species" | bg == "species_none")) {
    vr[n, "n.fecundity"] <- sample.size.sps
    glm.fec <- glm(log(seeds.autumn0) ~ size.autumn0, family="gaussian", data=d.fec)
    
    # bootstrap
    boot.fec <- mvrnorm(nboot,mu=coef(glm.fec)[1:2], Sigma=vcov(glm.fec)[1:2,1:2])
    
    # coefs and sigma
    vr[n, "fecundity.int"] <- coef(glm.fec)[1]
    vr[n, "fecundity.slope"] <- coef(glm.fec)[2]
    vr.boot[n, match("fecundity.int", colnames(vr)), 1:nboot] = boot.fec[,1]
    vr.boot[n, match("fecundity.slope", colnames(vr)), 1:nboot] = boot.fec[,2]
    rm(glm.fec)
    rm(boot.fec)
  }
  else if(fc == "Armo") {
    vr[n, "fecundity.int"] <- -800 # exp(-800) == 0 in R
    vr[n, "fecundity.slope"] <- 0
    vr.boot[n, match("fecundity.int", colnames(vr)), 1:nboot] = -800
    vr.boot[n, match("fecundity.slope", colnames(vr)), 1:nboot] = 0
  }
  print("fecunidty")
  
  # **** - 3.2.5 germination ----------------------------------------------
  f.ger <- filter(f, species == fc)
  
  sample.size.pair <- table(f.ger[,"site"])[1]
  vr[n,"n.germination"] <- sample.size.pair
  
  y.ger <- as.character(filter(vr.models, species == fc & response == "germination")[3])
  m.ger <- as.character(filter(vr.models, species == fc & response == "germination")[5])
  
  if(!is.na(m.ger) & (bg == "site" | bg == "site_none")) {
    # lm
    lm.ger <- lm(as.formula(paste(y.ger, m.ger, sep="~")), data=f.ger)
    # intercept and sigma
    vr[n,"germination.prob"] <- coef(lm.ger)[1]
    vr.boot[n, match("germination.prob", colnames(vr)), 1:nboot] = coef(lm.ger)[1]
    rm(m.ger)
    rm(lm.ger)
  }
  else if(!is.na(m.ger) & (bg == "species" | bg == "species_none")) {
    lm.ger <- lm(as.formula(paste(y.ger, "1", sep="~")), data=f.ger)
    # intercept and sigma
    vr[n,"germination.prob"] <- coef(lm.ger)[1]
    vr.boot[n, match("germination.prob", colnames(vr)), 1:nboot] = coef(lm.ger)[1]
    rm(m.ger)
    rm(lm.ger)
  }
  print("germination")
  
  # **** - 3.2.6.1 seedling establishment of FuNiche-------------------------------------
  f.est <- filter(f, species == fc & !is.na(establishement_summer_19))
  
  sample.size.site <- table(f.est[,"site"])[1]
  sample.size.sp <- nrow(f.est)
  vr[n,"n.establishment.funiche"] <- sample.size.site
  
  y.est <- as.character(filter(vr.models, species == fc & response == "establishment_FuNiche")[3])
  m.est <- as.character(filter(vr.models, species == fc & response == "establishment_FuNiche")[5])
  
  if(!is.na(m.est) & (bg == "site" | bg == "site_none")) {
    glm.est <- glm(as.formula(paste(y.est, m.est, sep="~")), data=f.est, family="binomial")
    vr[n,"establishment.funiche.int"] <- coef(glm.est)[1]
    vr.boot[n, match("establishment.funiche.int", colnames(vr)), 1:nboot] = coef(glm.est)[1]
    rm(m.est)
    rm(glm.est)
  }
  else if(!is.na(m.est) & (bg == "species" | bg == "species_none")) {
    glm.est <- glm(as.formula(paste(y.est, "1", sep="~")), data=f.est, family="binomial")
    vr[n,"establishment.funiche.int"] <- coef(glm.est)[1]
    vr.boot[n, match("establishment.funiche.int", colnames(vr)), 1:nboot] = coef(glm.est)[1]
    rm(m.est)
    rm(glm.est)
  }
  else if(is.na(m.est)) {
    vr[n,"establishment.funiche.int"] <- 20
    vr.boot[n, match("establishment.funiche.int", colnames(vr)), 1:nboot] = 20
  }
  print("establishment.funiche")

  # **** - 3.2.6.2 seedling establishment with competition----------------------------------
  d.est <- filter(d, focal.species == fc & 
                    (!is.na(planting.spring1) & planting.spring1=="yes") &
                    !is.na(survival.summer1.seedling) &
                    background.species %in% sps)
  
  # sample size
  sample.size.site <- table(d.est[,c("site")])[1]
  sample.size.sps <- nrow(d.est)
  
  y.est <- as.character(filter(vr.models, species == fc & response == "establishment")[3])
  m.est <- as.character(filter(vr.models, species == fc & response == "establishment")[5])
  
  if((bg == "site" | bg == "site_none") & sample.size.site > 2) {
    vr[n, "n.establishment.comp"] <- sample.size.site
    glm.est2 <- glm(as.formula(paste(y.est, m.est, sep="~")), data=d.est, family="binomial")
    
  } else if(bg == "species" | bg == "species_none" | sample.size.site < 3) {
    vr[n, "n.establishment.comp"] <- sample.size.sps
    glm.est2 <- glm(survival.summer1.seedling ~ 1, data=d.est, family="binomial")
  }
  
  # coefs
  vr[n,"establishment.comp.int"] <- coef(glm.est2)[1]
  vr.boot[n, match("establishment.comp.int", colnames(vr)), 1:nboot] = coef(glm.est2)[1]
  
  rm(glm.est2)
  print("establishment2") 
  
  # **** - 3.2.7 recruit size--------------------------------------------
  mean_sd <- function(x) {c(mean(x), sd(x))}
  
  # model
  g.rec <- filter(g, species == fc & greenhouse == "yes" & flowering == "no" & !is.na(drymass.g))
  sample.size.sp <- nrow(g.rec)
  vr[n, "n.seedling.size"] <- sample.size.sp
  
  # mean and sd
  lm.rec <- g.rec %>% 
    .$drymass.g %>%
    log() %>%
    mean_sd()
  
  vr[n, "seedling.size.mean"] = lm.rec[1] 
  vr[n, "seedling.size.sd"] = lm.rec[2]
  vr.boot[n, match("seedling.size.mean", colnames(vr)), 1:nboot] = lm.rec[1] 
  vr.boot[n, match("seedling.size.sd", colnames(vr)), 1:nboot] = lm.rec[2] 
  rm(lm.rec)
  print("recruit size")
  
  # **** - 3.2.8 size range: L, U----------------------------------------
  # sample size
  d.siz <- filter(d, focal.species == fc &
                    !is.na(size.autumn0) &
                    background.species %in% sps)
  
  sample.size.site <- table(d.siz[,c("site")])[1]
  sample.size.sps <- nrow(d.siz)
  
  if((bg == "site" | bg == "site_none")) {
    vr[n, "n.size.range"] <- sample.size.site
    d.siz <- filter(d, focal.species == fc &
                      site == st &
                      !is.na(size.autumn0) &
                      background.species %in% sps)
    g.siz <- filter(g, species == fc & greenhouse == "yes" & flowering == "no" & !is.na(drymass.g))
    
    vr[n, "L"] = min(c(range(log(g.siz$drymass.g)), range(d.siz$size.autumn0)), na.rm = TRUE)
    vr[n, "U"] = max(c(range(log(g.siz$drymass.g)), range(d.siz$size.autumn0)), na.rm = TRUE)
    vr.boot[n, match("L", colnames(vr)), 1:nboot] = as.numeric(vr[n, "L"])
    vr.boot[n, match("U", colnames(vr)), 1:nboot] = as.numeric(vr[n, "U"])
  } 
  else if(bg == "species" | bg == "species_none") {
    vr[n, "n.size.range"] <- sample.size.sps
    
    d.siz <- filter(d, focal.species == fc &
                      !is.na(size.autumn0) &
                      background.species %in% sps)
    g.siz <- filter(g, species == fc & greenhouse == "yes" & flowering == "no" & !is.na(drymass.g))
    
    vr[n, "L"] = min(range(log(g.siz$drymass.g)), range(d.siz$size.autumn0), na.rm = TRUE)
    vr[n, "U"] = max(range(log(g.siz$drymass.g)), range(d.siz$size.autumn0), na.rm = TRUE)
    vr.boot[n, match("L", colnames(vr)), 1:nboot] = as.numeric(vr[n, "L"])
    vr.boot[n, match("U", colnames(vr)), 1:nboot] = as.numeric(vr[n, "U"])
  }
  print("size range")
  
} # for loop

#****************************************************************************
# **  - 3.3 Fill up vital rates gaps **------
#****************************************************************************
for(n in 1:nrow(vr)) {
  # n = 1
  vi <- vr[n,]
  fc <- as.character(vi$focal.species)
  bg <- as.character(vi$background.species)
  st <- as.character(vi$site)
  
  # relevel focal species, background species and site
  d$focal.species <- relevel(d$focal.species, ref=fc)
  #d$background.species <- relevel(d$background.species, ref=bg)
  d$site <- relevel(d$site, ref=st)
  
  ## **** - 3.3.1 survival -----------------------------------------------------
  
  ## **** - 3.3.2 growth-----------------------------------------------------
  # SD of growth
  #if(!is.na(vi$pair) & !is.na(vi$growth.int) & is.na(vi$growth.sd) ) {
  #  vr[n, "growth.sd"] <- 0.01
  #  vr.boot[n, match("growth.sd", colnames(vr)), 1:nboot] = 0.01
  #}
  
  ## **** - 3.3.3 flowering -----------------------------------------------------
  
  ## **** - 3.3.4 fecundity -----------------------------------------------------
  
  ## **** - 3.3.5 germination ------
  # Daca germination establishment and estiablishment-FuNiche
  if(!is.na(fc) & fc == "Daca" & !is.na(vi$pair) & vi$pair != "no") {
    if(st == "Les Posses") {
      vr[n, "germination.prob"] <- 0.055
      vr[n, "establishment.funiche.int"] <- 20 # logit(20) == 1
      vr.boot[n, match("germination.prob", colnames(vr)), 1:nboot] = 0.055
      vr.boot[n, match("establishment.funiche.int", colnames(vr)), 1:nboot] = 20
      
      # print
      print(n)
      print(paste(fc,st,bg, sep="_"))
      
    }
    else if(st == "Solalex") {
      vr[n, "germination.prob"] <- 0.12
      vr[n, "establishment.funiche.int"] <- 20 # logit(20) == 1
      vr.boot[n, match("germination.prob", colnames(vr)), 1:nboot] = 0.12
      vr.boot[n, match("establishment.funiche.int", colnames(vr)), 1:nboot] = 20
      
      # print
      print(n)
      print(paste(fc,st,bg, sep="_"))
    }
    else if(st == "Anzeindaz") {
      vr[n, "germination.prob"] <- 0.04
      vr[n, "establishment.funiche.int"] <- 20 # logit(20) == 1
      vr.boot[n, match("germination.prob", colnames(vr)), 1:nboot] = 0.04
      vr.boot[n, match("establishment.funiche.int", colnames(vr)), 1:nboot] = 20
      
      # print
      print(n)
      print(paste(fc,st,bg, sep="_"))
    }
  }
  ## **** - 3.3.6 size range: plots transplanted in 2018!----
} 

#****************************************************************************
# combine into one data frame ----
vr$bootstrap = 0
vr.bootstrap <- vr
for(i in 1:499) {
  # i = 1
  vr.i <- as.data.frame(vr.boot[,,i])
  vr.i$bootstrap = i
  colnames(vr.i) = colnames(vr)
  vr.i[1:9] = vr[1:9]
  vr.bootstrap = rbind(vr.bootstrap, vr.i)
  print(i)
}
vr.bootstrap

#****************************************************************************
# ** - 3.4 Plot vital rates against raw data ** ------
#****************************************************************************
# **** - 3.4.0 make predictions -------------------------------------
# vital rate paramters of all pairs
vr.par <- vr %>% 
  filter(!(background.species %in% c("site", "site_none", "species", "species_none"))) %>%
  filter(!is.na(pair) & pair != "no") 
vr.par

# make predictions based on estimated vital rates and size ranges
vr.prediction <- NULL

for(i in 1:nrow(vr.par)) {
  # i = 1
  par.i <- vr.par[i,]
  
  # get pair information
  fc.i <- par.i$focal.species
  bg.i <- par.i$background.species
  st.i <- par.i$site
  
  # make size range based on lower and upper size bounds
  size.i <- seq(as.numeric(par.i$L), as.numeric(par.i$U), length.out = 100)
  
  # survival
  survival.i <- survival.z(z=size.i, params = par.i)
  
  # growth
  growth.i <- growth.z(z=size.i, params = par.i)
  
  # flowering
  flowering.i <- flowering.z(z = size.i, params = par.i)
  
  # fecundity
  fecundity.i <- seed.z.gaussian(z=size.i, params = par.i)
  
  # combine all vital rates into one data frame
  df.i <- data.frame(focal.species = fc.i,
                     background.species = bg.i,
                     site = st.i,
                     xx=rep(size.i,4),
                     yy = c(survival.i, growth.i, flowering.i, fecundity.i),
                     vital.rate = rep(c("survival", "growth", "flowering", "fecundity"), each=100))
  
  # combine all pairs into one data
  vr.prediction <- rbind(vr.prediction, df.i)
}
vr.prediction
vr.prediction$focal.species <- factor(vr.prediction$focal.species, levels = c("Anal", "Armo","Asal", "Plal", "Poal","Seca", "Trba", "Brer", "Crbi", "Daca", "Melu", "Plla", "Potr", "Sapr"))
vr.prediction$background.species <- factor(vr.prediction$background.species, levels = c("Anal","Asal", "Plal", "Poal","Seca", "Trba", "Brer", "Crbi", "Melu", "Plla", "Potr", "Sapr", "none"))

# chose focal species
sps = c("Anal", "Armo","Asal", "Plal", "Poal","Seca", "Trba", "Brer", "Crbi", "Daca", "Melu", "Plla", "Potr", "Sapr")
fc = "Daca"

# **** - 3.4.1 survival -------------------------------------
fig.survival <- d %>%
  filter(!is.na(size.autumn0) & !is.na(survival)) %>%
  filter((!is.na(if.same) & if.same==TRUE)) %>%
  filter(!is.na(site) & !is.na(focal.species) & !is.na(background.species)) %>%
  filter(!(focal.species %in% c("Anvu", "Trca"))) %>%
  filter(!(background.species %in% c("Anvu", "Trca"))) %>%
  mutate(background.species = factor(background.species, levels = c("Anal","Asal", "Plal", "Poal","Seca", "Trba", "Brer", "Crbi", "Melu", "Plla", "Potr", "Sapr", "none"))) %>%  
  # for each focal species
  #filter(focal.species %in% fc) %>%
  #ggplot(aes(x=size.autumn0, y = survival, col=site)) +
  #geom_jitter(height=0.1, alpha = 0.5) +
  #geom_line(data=filter(vr.prediction, vital.rate == "survival" & focal.species %in% fc), aes(x=xx, y =yy)) +
  #scale_color_manual(values=c("Les Posses" = "orange", "Solalex"="green", "Anzeindaz"="blue")) + 
  #facet_wrap(~background.species, nrow=2) +
  # all focal species
  mutate(focal.species = factor(focal.species, levels = c("Anal", "Armo","Asal", "Plal", "Poal","Seca", "Trba", "Brer", "Crbi", "Daca", "Melu", "Plla", "Potr", "Sapr"))) %>%  
  ggplot(aes(x=size.autumn0, y = survival, col = background.species, shape = site, linetype=site)) +
  geom_jitter(height=0.1, width=0, col="grey") +
  geom_line(data=filter(vr.prediction, vital.rate == "survival"), aes(x=xx, y =yy)) +
  scale_color_manual(values=c("Anal" = turbo(12)[1], "Asal"=turbo(12)[2], "Plal"=turbo(12)[3], "Poal"=turbo(12)[4], "Seca"=turbo(12)[5], "Trba"=turbo(12)[6], 
                              "Brer" = turbo(12)[7], "Crbi"=turbo(12)[8], "Melu"=turbo(12)[9], "Plla"=turbo(12)[10], "Potr"=turbo(12)[11], "Sapr"=turbo(12)[12],
                              "none" = "black")) + 
  facet_wrap(~focal.species, nrow=2, scales = "free") +
  geom_rug(sides = "b", col="grey50") +
  scale_x_continuous(name = "log(Size at year t)") +
  scale_y_continuous(name="Prob (Survival)")
fig.survival

# **** - 3.4.2 growth -------------------------------------
fig.growth <- d %>%
  filter(!is.na(size.autumn0) & !is.na(survival)) %>%
  filter((!is.na(if.same) & if.same==TRUE)) %>%
  filter(!is.na(site) & !is.na(focal.species) & !is.na(background.species)) %>%
  filter(!(focal.species %in% c("Anvu", "Trca"))) %>%
  filter(!(background.species %in% c("Anvu", "Trca"))) %>%
  mutate(background.species = factor(background.species, levels = c("Anal","Asal", "Plal", "Poal","Seca", "Trba", "Brer", "Crbi", "Melu", "Plla", "Potr", "Sapr", "none"))) %>%  
  # for each focal species
  #filter(focal.species %in% fc) %>%
  #ggplot(aes(x=size.autumn0, y = size.autumn1, col=site)) +
  #geom_abline(intercept = 0, slope=1, linetype="dashed", col="grey") +
  #geom_jitter(height=0.1, alpha = 0.5) +
  #geom_line(data=filter(vr.prediction, vital.rate == "growth" & focal.species %in% fc), aes(x=xx, y =yy)) +
  #scale_color_manual(values=c("Les Posses" = "orange", "Solalex"="green", "Anzeindaz"="blue")) + 
  #facet_wrap(~background.species, nrow=2) +
  # all focal species
  mutate(focal.species = factor(focal.species, levels = c("Anal", "Armo","Asal", "Plal", "Poal","Seca", "Trba", "Brer", "Crbi", "Daca", "Melu", "Plla", "Potr", "Sapr"))) %>%  
  ggplot(aes(x=size.autumn0, y = size.autumn1, col = background.species, shape = site, linetype=site)) +
  geom_point(col="grey") +
  geom_line(data=filter(vr.prediction, vital.rate == "growth"), aes(x=xx, y =yy)) +
  scale_color_manual(values=c("Anal" = turbo(12)[1], "Asal"=turbo(12)[2], "Plal"=turbo(12)[3], "Poal"=turbo(12)[4], "Seca"=turbo(12)[5], "Trba"=turbo(12)[6], 
                              "Brer" = turbo(12)[7], "Crbi"=turbo(12)[8], "Melu"=turbo(12)[9], "Plla"=turbo(12)[10], "Potr"=turbo(12)[11], "Sapr"=turbo(12)[12],
                              "none" = "black")) + 
  facet_wrap(~focal.species, nrow=2, scales = "free") +
  geom_rug(sides = "b", col="grey50") +
  scale_x_continuous(name = "log(Size at year t)") +
  scale_y_continuous(name="log(Size at year t+1)")
fig.growth

# **** - 3.4.3 flowering -------------------------------------
fig.flowering <- d %>%
  filter(!is.na(size.autumn0) & !is.na(survival)) %>%
  filter((!is.na(if.same) & if.same==TRUE)) %>%
  filter(!is.na(site) & !is.na(focal.species) & !is.na(background.species)) %>%
  filter(!(focal.species %in% c("Anvu", "Trca"))) %>%
  filter(!(background.species %in% c("Anvu", "Trca"))) %>%
  mutate(background.species = factor(background.species, levels = c("Anal","Asal", "Plal", "Poal","Seca", "Trba", "Brer", "Crbi", "Melu", "Plla", "Potr", "Sapr", "none"))) %>%  
  # for each focal species
  #filter(focal.species %in% fc) %>%
  #ggplot(aes(x=size.autumn0, y = flowering.autumn0, col=site)) +
  #geom_jitter(height=0.1, alpha = 0.5) +
  #geom_line(data=filter(vr.prediction, vital.rate == "flowering" & focal.species %in% fc), aes(x=xx, y =yy)) +
  #scale_color_manual(values=c("Les Posses" = "orange", "Solalex"="green", "Anzeindaz"="blue")) + 
  #facet_wrap(~background.species, nrow=2) +
  # all focal species
  mutate(focal.species = factor(focal.species, levels = c("Anal", "Armo","Asal", "Plal", "Poal","Seca", "Trba", "Brer", "Crbi", "Daca", "Melu", "Plla", "Potr", "Sapr"))) %>%  
  ggplot(aes(x=size.autumn0, y = flowering.autumn0, col = background.species, shape = site, linetype=site)) +
  geom_jitter(height=0.1, width=0,col="grey") +
  geom_line(data=filter(vr.prediction, vital.rate == "flowering"), aes(x=xx, y =yy)) +
  scale_color_manual(values=c("Anal" = turbo(12)[1], "Asal"=turbo(12)[2], "Plal"=turbo(12)[3], "Poal"=turbo(12)[4], "Seca"=turbo(12)[5], "Trba"=turbo(12)[6], 
                              "Brer" = turbo(12)[7], "Crbi"=turbo(12)[8], "Melu"=turbo(12)[9], "Plla"=turbo(12)[10], "Potr"=turbo(12)[11], "Sapr"=turbo(12)[12],
                              "none" = "black")) + 
  facet_wrap(~focal.species, nrow=2, scales = "free") +
  geom_rug(sides = "b", col="grey50") +
  scale_x_continuous(name = "log(Size at year t)") +
  scale_y_continuous(name="Prob (Flowering)")
fig.flowering

# **** - 3.4.4 fecundity -------------------------------------
fig.fecundity <- d %>%
  filter(!is.na(size.autumn0) & !is.na(survival)) %>%
  filter((!is.na(if.same) & if.same==TRUE)) %>%
  filter(!is.na(site) & !is.na(focal.species) & !is.na(background.species)) %>%
  filter(!(focal.species %in% c("Anvu", "Trca"))) %>%
  filter(!(background.species %in% c("Anvu", "Trca"))) %>%
  mutate(background.species = factor(background.species, levels = c("Anal","Asal", "Plal", "Poal","Seca", "Trba", "Brer", "Crbi", "Melu", "Plla", "Potr", "Sapr", "none"))) %>%  
  # for each focal species
  #filter(focal.species %in% fc) %>%
  #ggplot(aes(x=size.autumn0, y = log(seeds.autumn0), col=site)) +
  #geom_jitter(height=0.1, alpha = 0.5) +
  #geom_line(data=filter(vr.prediction, vital.rate == "fecundity" & focal.species %in% fc), aes(x=xx, y =log(yy))) +
  #scale_color_manual(values=c("Les Posses" = "orange", "Solalex"="green", "Anzeindaz"="blue")) + 
  #facet_wrap(~background.species, nrow=2) +
  # all focal species
  mutate(focal.species = factor(focal.species, levels = c("Anal", "Armo","Asal", "Plal", "Poal","Seca", "Trba", "Brer", "Crbi", "Daca", "Melu", "Plla", "Potr", "Sapr"))) %>%  
  ggplot(aes(x=size.autumn0, y = log(seeds.autumn0), col = background.species, shape = site, linetype=site)) +
  geom_point(col="grey") +
  geom_line(data=filter(vr.prediction, vital.rate == "fecundity"), aes(x=xx, y =log(yy))) +
  scale_color_manual(values=c("Anal" = turbo(12)[1], "Asal"=turbo(12)[2], "Plal"=turbo(12)[3], "Poal"=turbo(12)[4], "Seca"=turbo(12)[5], "Trba"=turbo(12)[6], 
                              "Brer" = turbo(12)[7], "Crbi"=turbo(12)[8], "Melu"=turbo(12)[9], "Plla"=turbo(12)[10], "Potr"=turbo(12)[11], "Sapr"=turbo(12)[12],
                              "none" = "black")) + 
  facet_wrap(~focal.species, nrow=2, scales = "free") +
  geom_rug(sides = "b", col="grey50") +
  scale_x_continuous(name = "log(Size at year t)") +
  scale_y_continuous(name="log(Numer of seeds)")
fig.fecundity

# **** - 3.4.5 germination ----------------------------------------------
vr.par$site2 <- dplyr::recode(vr.par$site, "Les Posses" = "Low", "Solalex" = "Middle", "Anzeindaz" = "High")
vr.par$focal.species <- factor(vr.par$focal.species, levels = c("Anal", "Armo","Asal", "Plal", "Poal","Seca", "Trba", "Brer", "Crbi", "Daca", "Melu", "Plla", "Potr", "Sapr"))
fig.germination  <- f %>%
  mutate(focal.species = species) %>%
  filter(focal.species != "Trca") %>%
  mutate(site = dplyr::recode(site, "Les Posses" = "Low", "Solalex" = "Middle", "Anzeindaz" = "High")) %>%
  mutate(site = factor(site, levels= c("Low", "Middle", "High"))) %>%
  mutate(focal.species = factor(focal.species, levels = c("Anal", "Armo","Asal", "Plal", "Poal","Seca", "Trba", "Brer", "Crbi", "Daca", "Melu", "Plla", "Potr", "Sapr"))) %>%
  ggplot(aes(x=site, y = germination_spring_19)) +
  facet_wrap(~focal.species,scales="free", nrow=2) +
  geom_point(col="grey") +
  geom_pointrange(data=vr.par,aes(x=site2,y=germination.prob, ymin=germination.prob-germination.sd.global, ymax=germination.prob+germination.sd.global)) +
  coord_cartesian(ylim=c(0,0.8)) +
  scale_x_discrete(name = "Site") +
  scale_y_continuous(name="Prob (Germination)")
fig.germination

# **** - 3.4.6.1 establishement without competition (FuNiche) --------------------------------
vr.par$establish.fun <- as.numeric(NA)
vr.par$establish.fun.min <- as.numeric(NA)
vr.par$establish.fun.max <- as.numeric(NA)

for(i in 1:nrow(vr.par)) {
  # i = 1
  par.i <- vr.par[i,]
  vr.par[i, "establish.fun"] = establishment.fun(params = par.i)
  vr.par[i, "establish.fun.max"] = exp(par.i$establishment.funiche.int + par.i$establishment.funiche.sd.global)/(1+exp(par.i$establishment.funiche.int + par.i$establishment.funiche.sd.global))
  vr.par[i, "establish.fun.min"] = exp(par.i$establishment.funiche.int - par.i$establishment.funiche.sd.global)/(1+exp(par.i$establishment.funiche.int - par.i$establishment.funiche.sd.global))
}

fig.establish.fun  <- f %>%
  mutate(focal.species = species) %>%
  filter(focal.species != "Trca") %>%
  mutate(focal.species = factor(focal.species, levels = c("Anal", "Armo","Asal", "Plal", "Poal","Seca", "Trba", "Brer", "Crbi", "Daca", "Melu", "Plla", "Potr", "Sapr"))) %>%
  mutate(site = dplyr::recode(site, "Les Posses" = "Low", "Solalex" = "Middle", "Anzeindaz" = "High")) %>%
  mutate(site = factor(site, levels= c("Low", "Middle", "High"))) %>%
  ggplot(aes(x=site, y = establishement_summer_19)) +
  geom_jitter(col="grey", height = 0.05, width=0.1) +
  geom_pointrange(data=vr.par,aes(x=site2,y=establish.fun, ymin=establish.fun.min, ymax=establish.fun.max)) +
  facet_wrap(~focal.species, scales="free",nrow=2) +
  scale_x_discrete(name = "Site") +
  scale_y_continuous(name="Prob (Establishment no-comp)")
fig.establish.fun

# **** - 3.4.6.2 seedling establishment under competition----------------
vr.par$establish.com <- as.numeric(NA)
vr.par$establish.com.min <- as.numeric(NA)
vr.par$establish.com.max <- as.numeric(NA)

for(i in 1:nrow(vr.par)) {
  # i = 1
  par.i <- vr.par[i,]
  vr.par[i, "establish.com"] = establishment.com(params = par.i)
  vr.par[i, "establish.com.max"] = exp(par.i$establishment.comp.int + par.i$establishment.comp.sd.global)/(1+exp(par.i$establishment.comp.int + par.i$establishment.comp.sd.global))
  vr.par[i, "establish.com.min"] = exp(par.i$establishment.comp.int - par.i$establishment.comp.sd.global)/(1+exp(par.i$establishment.comp.int - par.i$establishment.comp.sd.global))
}
vr.par$site = dplyr::recode(vr.par$site, "Les Posses" = "Low", "Solalex" = "Middle", "Anzeindaz" = "High")
vr.par$site = factor(vr.par$site, levels= c("Low", "Middle", "High"))

fig.establish.competition <- d %>%
  filter(!is.na(planting.spring1)) %>%
  filter(planting.spring1=="yes") %>%
  filter(!is.na(survival.summer1.seedling)) %>%
  filter(!(focal.species %in% c("Anvu", "Trca"))) %>%
  filter(!(background.species %in% c("Anvu", "Trca"))) %>%
  mutate(site = dplyr::recode(site, "Les Posses" = "Low", "Solalex" = "Middle", "Anzeindaz" = "High")) %>%
  mutate(site = factor(site, levels= c("Low", "Middle", "High"))) %>%
  # one focal species
  #filter(focal.species == fc) %>%
  #ggplot(aes(x=site, y=survival.summer1.seedling, col=site)) +
  #geom_jitter(col="grey", height=0.05, width=0.1) +
  #geom_pointrange(data=filter(vr.par, focal.species == fc),aes(x=site,y=establish.com, ymin=establish.com.min, ymax=establish.com.max), alpha=0.6) +
  #scale_color_manual(values=c("Low" = "orange", "Middle" = "green", "High" = "blue")) + 
  #facet_wrap(~background.species) +
  
  # all species
  mutate(focal.species = factor(focal.species, levels = c("Anal", "Armo","Asal", "Plal", "Poal","Seca", "Trba", "Brer", "Crbi", "Daca", "Melu", "Plla", "Potr", "Sapr"))) %>%
  ggplot(aes(x=site, y=survival.summer1.seedling, col=background.species)) +
  geom_jitter(col="grey", height=0.05, width=0.1) +
  geom_pointrange(data=vr.par,aes(x=site,y=establish.com, ymin=establish.com.min, ymax=establish.com.max), alpha=0.6) +
  scale_color_manual(values=c("Anal" = turbo(12)[1], "Asal"=turbo(12)[2], "Plal"=turbo(12)[3], "Poal"=turbo(12)[4], "Seca"=turbo(12)[5], "Trba"=turbo(12)[6], 
                              "Brer" = turbo(12)[7], "Crbi"=turbo(12)[8], "Melu"=turbo(12)[9], "Plla"=turbo(12)[10], "Potr"=turbo(12)[11], "Sapr"=turbo(12)[12],
                              "none" = "black")) + 
  facet_wrap(~focal.species, scales="free", nrow=2) +
  
  scale_x_discrete(name = "Site") +
  scale_y_continuous(name="Prob (Establishment under comp)")
fig.establish.competition
 
#**************************************************
# **** - 3.4.7 recruit size----
vr.recruit <- filter(vr, background.species == "species")
vr.recruit$species = vr.recruit$focal.species
rec.prediciton <- NULL
for(i in 1:nrow(vr.recruit)) {
  # i = 1
  vr.recruit.i <- vr.recruit[i,]
  xx <- seq(vr.recruit.i$L,vr.recruit.i$U, length.out = 100)
  yy <- seedling.z1(xx, params = vr.recruit.i)
  d.i <- data.frame(xx=xx, yy=yy, focal.species = vr.recruit.i$focal.species)
  rec.prediciton <- rbind(rec.prediciton, d.i)
}
rec.prediciton$focal.species <- factor(rec.prediciton$focal.species,levels = c("Anal", "Armo","Asal", "Plal", "Poal","Seca", "Trba", "Brer", "Crbi", "Daca", "Melu", "Plla", "Potr", "Sapr") )

fig.recruit <- g %>%
  filter(greenhouse == "yes") %>%
  filter(flowering == "no") %>%
  filter(!is.na(drymass.g)) %>%
  mutate(focal.species = factor(species, levels = c("Anal", "Armo","Asal", "Plal", "Poal","Seca", "Trba", "Brer", "Crbi", "Daca", "Melu", "Plla", "Potr", "Sapr"))) %>%
  ggplot(aes(x=log(drymass.g))) +
  geom_histogram(aes(y = ..density..), fill="grey") +
  geom_line(data=rec.prediciton, aes(x=xx, y=yy))  +
  facet_wrap(~focal.species, scales="free", nrow=2)  +
  scale_x_continuous(name = "log(Seedling size)") +
  scale_y_continuous(name="Density")
fig.recruit

fig.germination
fig.establish.fun
fig.establish.competition
fig.recruit
fig.survival
fig.growth
fig.flowering
fig.fecundity
fig.establish.competition

#****************************************************************************
# ** - 3.5 Check strange vital rates parameters ** ------
#****************************************************************************
# if there are strange intercept or slopes
vr %>%
  pivot_longer(cols = par, names_to="vital.rate", values_to="estimate") %>%
  #filter(estimate < 100) %>%
  ggplot(aes(x=estimate)) +
    geom_histogram() +
  facet_wrap(~vital.rate, scales="free") +
  stat_bin(aes(y=..count.., label=..count..), geom="text", vjust=-.5) 

#****************************************************************************
# 4. IPM implementation ------
#****************************************************************************

#****************************************************************************
# ** - 4.1 Size eviction ------
#****************************************************************************
# vital rates data
vr.eviction <- vr
vr.eviction$p_eviction_recruit_lower = as.numeric(NA)
vr.eviction$p_eviction_recruit_upper = as.numeric(NA)
vr.eviction$p_eviction_growth_lower = as.numeric(NA)
vr.eviction$p_eviction_growth_upper = as.numeric(NA)

for(i in 2:nrow(vr.eviction)) {
  # i = 2
  vr.i <- vr.eviction[i,]
  if(is.na(vr.i$pgr) | vr.i$background.species %in% c("site", "site_none", "species", "species_none")) {
    next
  }
  else {
    eviction.test.i <- size.eviction(params = vr.i, L = as.numeric(vr.i$L), U= as.numeric(vr.i$U), n = 1000, growth=TRUE, rel.tol = 1e-10)

    # recruit
    vr.eviction[i,"p_eviction_recruit_lower"] <- eviction.test.i$p_eviction_recruit$p_eviction_recruit_lower
    vr.eviction[i,"p_eviction_recruit_upper"] <- eviction.test.i$p_eviction_recruit$p_eviction_recruit_upper

    # growth
    # plot
    #ggplot(eviction.test.i$p_eviction_growth, aes(x=z, y = p_eviction_growth_lower)) + geom_point()
    #ggplot(eviction.test.i$p_eviction_growth, aes(x=z, y = p_eviction_growth_upper)) + geom_point()
    vr.eviction[i,"p_eviction_growth_lower"] <- mean(eviction.test.i$p_eviction_growth$p_eviction_growth_lower)
    vr.eviction[i,"p_eviction_growth_upper"] <- mean(eviction.test.i$p_eviction_growth$p_eviction_growth_upper)
  }
}

# recruit
# upper bound: is fine
r.upper <- vr.eviction %>% 
  filter(pair == "yes") %>%
  filter(!(background.species %in% c("site","site_none","species", "species_none"))) %>%
  .$p_eviction_recruit_upper
hist(r.upper)

# lower bound
r.lower <- vr.eviction %>% 
  filter(pair == "yes") %>%
  filter(!(background.species %in% c("site","site_none","species", "species_none"))) %>%
  .$p_eviction_recruit_lower
hist(r.lower)

# growth
# lower bound
g.lower <- vr.eviction %>% 
  filter(pair == "yes") %>%
  filter(!(background.species %in% c("site","site_none","species", "species_none"))) %>%
  .$p_eviction_growth_lower
hist(g.lower)

# upper bound
g.upper <- vr.eviction %>% 
  filter(pair == "yes") %>%
  filter(!(background.species %in% c("site","site_none","species", "species_none"))) %>%
  .$p_eviction_growth_upper
hist(g.upper)

# mean
mean(c(r.upper, r.lower, g.lower, g.upper))

#****************************************************************************
# ** - 4.2 Compare population grwoth rates with vs without ceiling ------
#****************************************************************************
vr.eviction$pgr
vr.eviction$lambda <- as.numeric(NA)
vr.eviction$lambda.ceiling <- as.numeric(NA)
vr.eviction$lambda.delta <- as.numeric(NA)
vr.eviction

for(a in 1:nrow(vr.eviction)) {
  # a = 1
  da <- vr.eviction[a,]
  da = vr.boot[a,]
  
  # skip species not needed
  if(is.na(da$pgr) | da$background.species %in% c("site", "site_none", "species", "species_none")) next
  
  # skip when kernel is null
  ipm.i <- mk.kernel(n = 1000, L = as.numeric(da$L), U = da$U, par = da)
  if(sum(ipm.i$K != 0) == 0) next
  
  # calculate population grwoth with and without ceiling
  lambda.delta.i <- eviction_delta.lambda_iter(growthKernel = growth.z1z, kernel = full.z1z, survivalFunction = survival.z,
                                               minsize = as.numeric(da$L), maxsize = da$U, params = da, n.big.matrix = 3000)
    
  vr.eviction[a,"lambda"] <- lambda.delta.i$lambda
  vr.eviction[a,"lambda.ceiling"] <- lambda.delta.i$lambda2
  vr.eviction[a,"lambda.delta"] <- lambda.delta.i$dlambda
  print(a)
}

#************************************************
# compare pgr with and without ceiling
vr.eviction %>%
  mutate(lambda = ifelse(lambda > 0.001, lambda, 0.001)) %>%
  mutate(lambda.ceiling = ifelse(lambda.ceiling > 0.001, lambda.ceiling, 0.001)) %>%
  ggplot(aes(x=lambda, y=lambda.ceiling)) +
  #ggplot(aes(x=log(lambda), y=log(lambda.ceiling))) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_abline(intercept = 0, slope=1) +
  geom_point()

#****************************************************************************
# ** - 4.3 Determine number of bins ------
#****************************************************************************
vr.test <-  vr %>% 
  filter(pair == "yes") %>%
  filter(!is.na(growth.int))
vr.test$n.converge <- as.numeric(NA)

# select 100 IPMs randomly
set.seed(1)
ipms <- sample(1:nrow(vr.test), 100)

# compute number of bins when lambdas converge
test <- sapply(ipms, function(i) {number.bins(params = vr.test[i,], tol=1e-8)})

#****************************************************************************
# ** - 4.4 calculate lambda using cluster ------
#****************************************************************************
library(parallel)

# on Euler
# source("/cluster/home/slyu/ipm_bootstrap_20220121/NICH_coexistence_IPM_functions.R")
# vr.boot <- read.csv("/cluster/home/slyu/ipm_bootstrap_20220121/vr.bootstrap_20220121.csv")
# n.cores = 20
# n.seq <- 1:nrow(vr)

# on Mac
# source("/Users/slyu/LVSM/R codes/NICH_coexistence_IPM_functions.R")
# vr.boot <- read.csv("/Users/slyu/LVSM/NICH/Results/IPM_log.size/vital.rate/vr.bootstrap_20220121.csv")
# n.cores = 2
# n.seq <- 5509

lambda.boot <- function(a) {
  # a = 276
  da <- vr.boot[a,]
  
  if(is.na(da$pair) | 
     da$pair == "no" | 
     da$background.species %in% c("site","site_none","species", "species_none")) {
    lambda.a <- list(lambda = "skip")
    print(paste(a, "_skip"))
  }
  else {
    ipm <- mk.kernel(n=3000, L=as.numeric(da$L), U=as.numeric(da$U), par=da)
    if(sum(is.na(ipm$K))>0) { 
      lambda.a <- list(lambda = "null")
      print(paste(a, "_null")) 
    } 
    else {
      cal.a <- lambda.eigen(ipm$K)
      lam.a <- cal.a$lambda
      lambda.a <- list(lambda = lam.a)
      print(a)
    }
  }
  return(lambda.a)
}

pgr.bootstrap <- mclapply(n.seq, lambda.boot, mc.cores = n.cores)
pgr.bootstrap

# output
# save(pgr.bootstrap, file="/cluster/home/slyu/ipm_bootstrap_20220121/lambda_bootstrap.RData")

#save.image("/Users/slyu/LVSM/NICH/Results/IPM_bootstrap/ipm_bootstrap_image_20210315.RData")
#save(p, file="pgr_boot_20210315.R")

#****************************************************************************
# ** - 4.5 combine bootstrapped vital rates and lambdas ------
#****************************************************************************
rm(list=ls())

# bootstraped vital rates 
vr.boot <- read_csv("/Users/slyu/LVSM/NICH/Results/IPM_log.size/vital.rate/vr.bootstrap_20220121.csv")
vr.boot$pgr <- as.numeric(NA)
vr.boot

# bootstraped lambdas
load("/Users/slyu/LVSM/NICH/Results/IPM_log.size/vital.rate/lambda_bootstrap.RData")
length(pgr.bootstrap)

# combine bootstraped vital rates and lambdas
for(i in 1:nrow(vr.boot)) {
  # i = 1
  pgr.boot.i <- pgr.bootstrap[[i]]
  if(pgr.boot.i$lambda != "skip" & pgr.boot.i$lambda != "null") {
    vr.boot[i,"pgr"] = pgr.boot.i$lambda
    print(i)
  }
}

#****************************************************************************
# 5. Bootstraped lambdas ------
#****************************************************************************

#****************************************************************************
# ** - 5.0 Explore bootstraped lambdas ------
#****************************************************************************
rm(list=ls())

# functions
source("/Users/slyu/LVSM/R codes/NICH_coexistence_IPM_functions.R")

# bootstraped lambdas 
vr.boot <- read_csv("SD2_lambda.only_bootstrap.csv", na="NA", col_names = TRUE)
vr.boot

#************************************************
# is there any extremely big lambdas?
# 198 PGR > 20: Brer (1), Crbi (193)
# 198 PGR > 50: Crbi (154)
# 132 PGR  > 100: Crbi (122)
# 125 PGR > 150: Crbi: 117

vr.boot %>%
  filter(!(background.species %in% c("site", "site_none", "species", "species_none"))) %>%
  filter(pgr > 20 & !is.na(pgr)) %>%
  .$background.species %>% table()
  .[c("focal.species", "background.species", "site", "bootstrap", "pgr")]

#************************************************
# is there any extremely small lambdas?
# floor: 0.01
# PGR  == 0 (< 1e-324): Crbi (9)
# 1885 PGR < exp(-3): Brer (2201), Trba (1093), Crbi (118), Asal (33), Potr (13)
# 1584 PGR < 0.01: Brer (2199), Crbi (101), Potr (4), Trba (13)
# 1584 PGR < 0.001: Brer (2197), Crbi (75)
# 1584 PGR < 0.0001: Brer (2192), Crbi (58)
# 1460 PGR < 0.00001: Brer (2043), Crbi (49)

vr.boot %>%
  filter(!(background.species %in% c("site", "site_none", "species", "species_none"))) %>%
  filter(pgr< exp(-4) & !is.na(pgr)) %>%
  .$focal.species %>% table()
.[c("focal.species", "background.species", "site", "bootstrap", "pgr")]
  
#****************************************************************************
# ** - 5.1 Calculate mean and standard error across all bootstraps ------
#****************************************************************************
vr.mean <- filter(vr.boot, bootstrap == 0)

vr.mean$pgr.mean <- as.numeric(NA)
vr.mean$pgr.median <- as.numeric(NA)
vr.mean$pgr.ci.min<- as.numeric(NA)
vr.mean$pgr.ci.max <- as.numeric(NA)
vr.mean$pgr.se <- as.numeric(NA)
vr.mean$pgr.cv <- as.numeric(NA)

vr.mean$pgr.log.mean <- as.numeric(NA)
vr.mean$pgr.log.median <- as.numeric(NA)
vr.mean$pgr.log.ci.min <- as.numeric(NA)
vr.mean$pgr.log.ci.max <- as.numeric(NA)
vr.mean$pgr.log.se <- as.numeric(NA)
vr.mean$pgr.log.cv <- as.numeric(NA)

for(i in 1:nrow(vr.mean)) {
  # i = 1
  di = vr.mean[i,]
  fc.i = di$focal.species
  bg.i = di$background.species
  st.i = di$site
  
  # get bootstrapped PGR for each pair
  pgr.boot.i = filter(vr.boot,focal.species == fc.i & background.species == bg.i & site == st.i)$pgr
    
  # drop small or huge PGRs
  pgr.boot.i = pgr.boot.i[pgr.boot.i < 20]
  pgr.boot.i = pgr.boot.i[pgr.boot.i > 0.01]
  #pgr.boot.i = pgr.boot.i[pgr.boot.i > 0.1]
  
  # log
  pgr.boot.i.log = log(pgr.boot.i)
  print(sum(!is.na(pgr.boot.i.log)))
  
  if(sum(!is.na(pgr.boot.i.log)) > 1) {
    # raw PGR
    mean.se.i <- ggplot2::mean_se(pgr.boot.i)
    mean.ci.i <- median_ci_quantile(pgr.boot.i, na.rm=TRUE)
    cv.i <- sd((pgr.boot.i))/mean.se.i[1]
    # store calculated mean, se and cv
    vr.mean[i,"pgr.mean"] = mean.se.i[1]
    vr.mean[i,"pgr.median"] = mean.ci.i[1]
    vr.mean[i,"pgr.ci.min"] = mean.ci.i[2]
    vr.mean[i,"pgr.ci.max"] = mean.ci.i[3]
    vr.mean[i,"pgr.se"] = mean.se.i[3] - mean.se.i[1]
    vr.mean[i,"pgr.cv"] = cv.i
    
    # log PGR
    mean.log.se.i <- ggplot2::mean_se(pgr.boot.i.log)
    mean.log.ci.i <- median_ci_quantile(pgr.boot.i.log, na.rm=TRUE)
    cv.log.i <- sd((pgr.boot.i.log))/mean.log.se.i[1]
    
    # store calculated mean, se and cv
      vr.mean[i,"pgr.log.mean"] = mean.log.se.i[1]
    vr.mean[i,"pgr.log.median"] = mean.log.ci.i[1]
    vr.mean[i,"pgr.log.ci.min"] = mean.log.ci.i[2]
    vr.mean[i,"pgr.log.ci.max"] = mean.log.ci.i[3]
    vr.mean[i,"pgr.log.se"] = mean.log.se.i[3] - mean.log.se.i[1]
    vr.mean[i,"pgr.log.cv"] = cv.log.i
    
    mean.se.i <- cv.i <- NA
  }
}

#************************************************
# competitor type: intra, inter, none
# vr mean
vr.mean$competitor <- as.character(NA)
vr.mean$neighbor <- as.character(NA)

for(i in 1:nrow(vr.mean)) {
  di <- vr.mean[i,]
  # Intrinsic,intra, inter
  if(di$background.species == "none") {
    vr.mean[i,"competitor"] = "None" 
    vr.mean[i,"neighbor"] = "Neighbors absent" 
  }
  else if(di$focal.species == di$background.species) {
    vr.mean[i,"competitor"] = "Intraspecific" 
    vr.mean[i,"neighbor"] = "Neighbors present" 
  }
  else if(di$focal.species != di$background.species) {
    vr.mean[i,"competitor"] = "Interspecific"
    vr.mean[i,"neighbor"] = "Neighbors present" 
  }
}

# vr boot
vr.boot$competitor <- as.character(NA)
vr.boot$neighbor <- as.character(NA)

# intrinsic
vr.boot[vr.boot$background.species == "none", "competitor"] = "None"
vr.boot[vr.boot$background.species == "none", "neighbor"] = "Neighbors absent"
vr.boot[vr.boot$background.species != "none", "neighbor"] = "Neighbors present"

# intra vs inter
vr.boot[vr.boot$background.species != "none" & vr.boot$background.species == vr.boot$focal.species, "competitor"] = "Intraspecific"
vr.boot[vr.boot$background.species != "none" & vr.boot$background.species != vr.boot$focal.species, "competitor"] = "Interspecific"

#************************************************
# plot PGR by focal species
fig.pgr_focal <- vr.mean %>%
  filter(!background.species %in% c("site","site_none","species", "species_none")) %>%
  filter(!is.na(pgr.mean)) %>%
  mutate(competitor = factor(competitor, levels=c("None", "Interspecific", "Intraspecific"))) %>%
  mutate(focal.species = factor(focal.species, levels=c("Anal","Armo","Asal","Plal", "Poal", "Seca", "Trba",    # 7 alpine species
                                                        "Brer", "Crbi", "Daca","Melu" , "Plla", "Potr", "Sapr"))) %>% # 7 lowland species
  ggplot(aes(x=elevation, y=pgr.log.median, ymin=pgr.log.ci.min, ymax=pgr.log.ci.max, col=competitor)) +
  geom_hline(yintercept = 0, linetype="dotted") +
  geom_pointrange(alpha=0.8, size=0.3) +
  geom_line(aes(group=background.species),size=0.3) +
  scale_color_manual(values=c("None" = "#6DCD59FF",  "Intraspecific" = "#FB9A06FF", "Interspecific" = "grey40")) +
  scale_y_continuous(name="ln(Population growth rate)") +
  scale_x_continuous(breaks=c(890,1400, 1900), name="Elevation (m)") +
  facet_wrap(~focal.species, scales = "free", ncol=4)
fig.pgr_focal

#****************************************************************************
# ** - 5.2 Intrinsic growth rate ------
#****************************************************************************

#***********************************************
# intrinsic PGR based on mean and SE
# test: only intrinsic
vr.mean %>%
  filter(background.species == "none") %>%
  lmer(pgr.log.median ~ elevation * origin.focal+ (1|focal.species), data=.) %>%
  #summary
  Anova(type=2)

# coefficients of mean PGR
coef.intrinsic.mean <- 
  vr.mean %>%
  filter(background.species == "none") %>%
  lmer(pgr.log.median ~ elevation * origin.focal+ (1|focal.species), data=.) %>%
  fixef()

#***********************************************
# intrinsic PGR for each bootstrap
# test for each bootstrap: only intrinsic
test.intrinsic.boot <- vr.boot %>%
  filter(background.species == "none") %>%
  mutate(pgr = ifelse(pgr < 0.01, 0.01, pgr)) %>% mutate(pgr = ifelse(pgr > 20, 20, pgr)) %>%
  split(.$bootstrap) %>%
  purrr::map(~lmer(log(pgr) ~ elevation * origin.focal+ (1|focal.species) , data=.)) 

coef.intrisinc <- test.intrinsic.boot %>%
  purrr::map_dfr(fixef, .id="bootstrap")

# get coefficient for each bootstrap
coef.intrinsic.boot <- test.intrinsic.boot %>% purrr::map(fixef)
coef.intrinsic.boot

# make fitted line
pre.intrinsic.bootstrap <- NULL
coef.intrinsic.boot.mean <- data.frame(bootstrap = 1:500,
                                       int.highland = NA,
                                       slo.highland = NA,
                                       int.lowland = NA,
                                       slo.lowland = NA)

for(i in 1:length(coef.intrinsic.boot)) {
  coef.i <- coef.intrinsic.boot[[i]]
  # elevation.sequence
  ee = seq(890,1900, length.out = 500)
  
  # extract bootstrapped linear coefficients
  coef.intrinsic.boot.mean[i,"int.highland"] = coef.i["(Intercept)"]
  coef.intrinsic.boot.mean[i,"slo.highland"] = coef.i["elevation"]
  coef.intrinsic.boot.mean[i,"int.lowland"] = coef.i["(Intercept)"] + coef.i["origin.focalLowland"]
  coef.intrinsic.boot.mean[i,"slo.lowland"] = coef.i["elevation"] + coef.i["elevation:origin.focalLowland"]
  
  # highland
  pre.highland = coef.intrinsic.boot.mean[i,"int.highland"] + coef.intrinsic.boot.mean[i,"slo.highland"]*ee
  # lowland
  pre.lowland = coef.intrinsic.boot.mean[i,"int.lowland"] + coef.intrinsic.boot.mean[i,"slo.lowland"]*ee
  # combine
  pre.i <- data.frame(bootstrap = i,
                      elevation = c(ee,ee),
                      pgr = c(pre.highland, pre.lowland),
                      origin.focal =rep(c("Highland", "Lowland"), each=500),
                      neighbor = "Neighbors absent"
  )
  pre.intrinsic.bootstrap <- rbind(pre.intrinsic.bootstrap, pre.i)
}

coef.intrinsic.mean
coef.intrinsic.mean.global = apply(coef.intrinsic.boot.mean,2,mean)[2:5]
coef.intrinsic.mean.global

# CI of coefficients
coef.invasion.intrisinc <- test.intrinsic.boot %>%
  purrr::map_dfr(fixef) 
# highland: 0.0007047861 0.0004836991 0.0008923608
median_ci_quantile(coef.invasion.intrisinc$elevation)
# lowland: -0.00112471 -0.00137776 -0.0008852111
median_ci_quantile(coef.invasion.intrisinc$`elevation:origin.focalLowland`)

# add ribbon
# SD CI
d.ribbon_pgr.intrinsic <- pre.intrinsic.bootstrap %>%
  group_by(elevation, origin.focal) %>%
  dplyr::summarise(mean_ci_sd(pgr)[2], mean_ci_sd(pgr)[3])
d.ribbon_pgr.intrinsic$neighbor = "Neighbors absent"

# add mean line
d.seg.intrinsic <- data.frame(x=890,
                              xend=1900,
                              y = c(coef.intrinsic.mean.global["int.highland"] + coef.intrinsic.mean.global["slo.highland"]*890,
                                    coef.intrinsic.mean.global["int.lowland"] + coef.intrinsic.mean.global["slo.lowland"]*890),
                              yend = c(coef.intrinsic.mean.global["int.highland"] + coef.intrinsic.mean.global["slo.highland"]*1900,
                                       coef.intrinsic.mean.global["int.lowland"] + coef.intrinsic.mean.global["slo.lowland"]*1900),
                              origin.focal = c("Highland", "Lowland"),
                              neighbor = "Neighbors absent",
                              significant = "marginal")

# fitted line for each bootstrap
fig.pgr.intrinsic_boot <- pre.intrinsic.bootstrap %>%
  ggplot(aes(col=origin.focal, fill=origin.focal)) +
  geom_hline(yintercept = 0, linetype="dashed") +
  
  # add CI
  geom_ribbon(data=d.ribbon_pgr.intrinsic, aes(x=elevation, ymin=ymin, ymax=ymax), col="white",alpha=0.3,show.legend = FALSE) +
  
  # add
  stat_summary(data=filter(vr.mean, background.species == "none"), aes(x=elevation, y=pgr.log.median), 
               fun.data= mean_se, alpha=0.8, show.legend = FALSE, position = position_dodge(width=20)) +
  
  # general trend of highland species
  geom_segment(data=d.seg.intrinsic, aes(x=x, xend=xend, y= y, yend = yend), size=1, show.legend = FALSE) +
  coord_cartesian(ylim = c(-3,3)) +
  scale_color_manual(values = c("Highland" = "#3E4A89FF",  "Lowland" = "#FB9A06FF")) +
  scale_fill_manual(values = c("Highland" = "#3E4A89FF",  "Lowland" = "#FB9A06FF")) +
  scale_x_continuous(name = "Elevation (m)", breaks = c(890,1400, 1900)) +
  scale_y_continuous(name = "ln(Population growth rate)")
fig.pgr.intrinsic_boot

#****************************************************************************
# ** - 5.2 Intraspecific invasion growth rate ------
#****************************************************************************

#***********************************************
# plot: intra-specific
fig.pgr.intra_ci_focal <- vr.mean %>%
  filter(pair=="yes") %>%
  filter(background.species == focal.species) %>% 
  mutate(elevation = factor(elevation, levels=c("890", "1400", "1900"))) %>%
  mutate(focal.species = factor(focal.species, levels=c("Anal","Armo","Asal","Plal", "Poal", "Seca", "Trba",    # 7 alpine species
                                                        "Brer", "Crbi", "Daca","Melu" , "Plla", "Potr", "Sapr"))) %>% # 7 lowland species
  ggplot(aes(x=focal.species, y=pgr.log.median, col=elevation, ymin=pgr.log.ci.min, ymax = pgr.log.ci.max)) +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_pointrange(position = position_dodge2(width=0.5)) +
  scale_color_manual(values=c("890" = "#FB9A06FF", "1400" = "#6DCD59FF", "1900" = "#3E4A89FF")) +
  facet_wrap(~origin.focal, scale="free") +
  coord_cartesian(ylim=c(-4,2))  +
  scale_x_discrete(name = "Focal species") +
  scale_y_continuous(name = "ln(Population growth rate)")
fig.pgr.intra_ci_focal

#***********************************************
# 11/34 (32%)intrinsic PGR have 95 CIs including zero, i.e., equilibrium
# 18/34 (53%) intrinsic PGR have 95 CIs including (-0.1, 0.1)
# 26/34 (76%) intrinsic PGR have 95 CIs of lambda including (0.5, 1.5)
d.temp <- vr.mean %>%
  filter(background.species == focal.species) %>% 
  mutate(focal.species = factor(focal.species, levels=c("Anal","Armo","Asal","Plal", "Poal", "Seca", "Trba",    # 7 alpine species
                                                        "Brer", "Crbi", "Daca","Melu" , "Plla", "Potr", "Sapr"))) %>% # 7 lowland species
  .$pgr.log.median

d.temp1 <- vr.mean %>%
  filter(background.species == focal.species) %>% 
  mutate(focal.species = factor(focal.species, levels=c("Anal","Armo","Asal","Plal", "Poal", "Seca", "Trba",    # 7 alpine species
                                                        "Brer", "Crbi", "Daca","Melu" , "Plla", "Potr", "Sapr"))) %>% # 7 lowland species
  .$pgr.log.ci.min

d.temp2 <- vr.mean %>%
  filter(background.species == focal.species) %>% 
  mutate(focal.species = factor(focal.species, levels=c("Anal","Armo","Asal","Plal", "Poal", "Seca", "Trba",    # 7 alpine species
                                                        "Brer", "Crbi", "Daca","Melu" , "Plla", "Potr", "Sapr"))) %>% # 7 lowland species
  .$pgr.log.ci.max

sum(!is.na(d.temp))
sum(d.temp1 < 0 & d.temp2 > 0, na.rm = TRUE)

#***********************************************
# 10 intra-PGRs were predicted to be over equilibrium
# 13 intra-PGRs were predicted to be under equilibrium

#****************************************************************************
# ** - 5.3 Invasion growth rate ------
#****************************************************************************
#***********************************************
# invasion PGR based on mean and SE
# distribution of invasion grwoth rates
vr.mean %>%
  filter(!(background.species %in% c("none","site","site_none","species", "species_none"))) %>%
  .$pgr.log.median %>%
  hist()

# test
test.pgr.invasion <- vr.mean %>%
  filter(!(background.species %in% c("none","site","site_none","species", "species_none"))) %>%
  filter(pair == "yes") %>%
  filter(!is.na(pgr.log.mean)) %>%
  filter(!is.infinite(pgr.log.mean)) %>%
  lmer(pgr.log.mean ~ elevation * origin.focal+ (1|focal.species) + (1|background.species), data=.) 

test.pgr.invasion %>%
    #summary() # 312 pairs
  Anova(type=2)

#****************************************
# invasion PGR test for each bootstrap
# test
test.invasion.boot <- vr.boot %>%
  filter(!(background.species %in% c("none","site","site_none","species", "species_none"))) %>%
  filter(pair == "yes") %>%
  filter(!is.na(pgr)) %>%
  
  # floor small PGRs
  mutate(pgr = ifelse(pgr < 0.01, 0.01, pgr)) %>%
  mutate(pgr = ifelse(pgr > 20, 20, pgr)) %>%
  
  split(.$bootstrap) %>%
  purrr::map(~lmer(log(pgr) ~ elevation * origin.focal+ (1|focal.species) + (1|background.species), data=.) )

# get ANOVA
p.value_pgr.invasion <- 
  test.invasion.boot %>% 
  purrr::map_dfr(coef.lmer)

# get coefficient for each bootstrap
coef.invasion.boot <- 
  test.invasion.boot %>% 
  purrr::map(fixef)
coef.invasion.boot

# make prediction for each bootstrap
pre.invasion.bootstrap <- NULL
coef.invasion.boot.mean <- data.frame(bootstrap = 1:500,
                                      int.highland = NA,
                                      slo.highland = NA,
                                      int.lowland = NA,
                                      slo.lowland = NA)

for(i in 1:length(coef.invasion.boot)) {
  #  i = 6
  coef.i <- coef.invasion.boot[[i]]
  # elevation.sequence
  ee = seq(890,1900, length.out = 500)
  
  # extract bootstrapped linear coefficients
  coef.invasion.boot.mean[i,"int.highland"] = coef.i["(Intercept)"]
  coef.invasion.boot.mean[i,"slo.highland"] = coef.i["elevation"]
  coef.invasion.boot.mean[i,"int.lowland"] = coef.i["(Intercept)"] + coef.i["origin.focalLowland"]
  coef.invasion.boot.mean[i,"slo.lowland"] = coef.i["elevation"] + coef.i["elevation:origin.focalLowland"]
  
  # highland
  pre.highland = coef.invasion.boot.mean[i,"int.highland"] + coef.invasion.boot.mean[i,"slo.highland"]*ee
  # lowland
  pre.lowland = coef.invasion.boot.mean[i,"int.lowland"] + coef.invasion.boot.mean[i,"slo.lowland"]*ee
  # combine
  pre.i <- data.frame(bootstrap = i,
                      elevation = c(ee,ee),
                      pgr = c(pre.highland, pre.lowland),
                      origin.focal =rep(c("Highland", "Lowland"), each=500),
                      neighbor = "Neighbors present"
  )
  pre.invasion.bootstrap <- rbind(pre.invasion.bootstrap, pre.i)
}
pre.invasion.bootstrap
coef.invasion.boot.mean
coef.invasion.mean.global <- apply(coef.invasion.boot.mean,2,mean)[2:5]
which(coef.invasion.boot.mean$slo.lowland > 0.0001)

# CI of coefficients
coef.invasion.interaction <- test.invasion.boot %>%
  purrr::map_dfr(fixef) 
# highland: 0.0005369556 0.0004219872 0.0006567842
hist(coef.invasion.interaction$elevation)
median_ci_quantile(coef.invasion.interaction$elevation)
# lowland: -0.0008081642 -0.001127529 -0.0005109716
hist(coef.invasion.interaction$`elevation:origin.focalLowland`)
median_ci_quantile(coef.invasion.interaction$`elevation:origin.focalLowland`)

# add ribbon
# SD CI
d.ribbon_pgr.invasion <- pre.invasion.bootstrap %>%
  group_by(elevation, origin.focal) %>%
  dplyr::summarise(mean_ci_sd(pgr)[2], mean_ci_sd(pgr)[3])
d.ribbon_pgr.invasion$neighbor = "Neighbors present"

# add mean trend
d.seg.invasion <- data.frame(x=890,
                             xend=1900,
                             y = c(coef.invasion.mean.global["int.highland"] + coef.invasion.mean.global["slo.highland"]*890,
                                   coef.invasion.mean.global["int.lowland"] + coef.invasion.mean.global["slo.lowland"]*890),
                             yend = c(coef.invasion.mean.global["int.highland"] + coef.invasion.mean.global["slo.highland"]*1900,
                                      coef.invasion.mean.global["int.lowland"] + coef.invasion.mean.global["slo.lowland"]*1900),
                             origin.focal = c("Highland", "Lowland"),
                             neighbor = "Neighbors present",
                             significant = "significant")

# fitted line for each bootstrap
fig.pgr.invasion_boot <-  pre.invasion.bootstrap %>%
  ggplot(aes( col=origin.focal, fill=origin.focal)) +
  geom_hline(yintercept = 0, linetype="dashed") +

  # add ribbon
  geom_ribbon(data=d.ribbon_pgr.invasion, aes(x=elevation, ymin=ymin, ymax=ymax), col="white", alpha=0.4, show.legend = FALSE) +
  
  # add mean points
  stat_summary(data = filter(vr.mean, !(background.species %in% c("none","site","site_none","species", "species_none"))),
               aes(x=elevation, y=pgr.log.median), fun.data = mean_se, alpha=0.8, show.legend = FALSE, position=position_dodge(width=20)) +
  
  # general trend of highland species
  geom_segment(data=d.seg.invasion, aes(x=x, y=y, xend=xend, yend = yend), size=1, show.legend = FALSE) +
  
  coord_cartesian(ylim = c(-3,3)) +
  scale_color_manual(values = c("Highland" = "#3E4A89FF",  "Lowland" = "#FB9A06FF")) +
  scale_fill_manual(values = c("Highland" = "#3E4A89FF",  "Lowland" = "#FB9A06FF")) +
  
  scale_x_continuous(name = "Elevation (m)", breaks = c(890,1400, 1900)) +
  scale_y_continuous(name = "ln(Population growth rate)")
fig.pgr.invasion_boot
fig.pgr.intrinsic_boot

#****************************************************************************
# ** - 5.4 Combine intrinsic and invasion PGRs ------
#****************************************************************************
fig.pgr_boot <-  bind_rows(pre.intrinsic.bootstrap, pre.invasion.bootstrap)  %>%
  ggplot(aes( col=origin.focal, fill=origin.focal)) +
  geom_hline(yintercept = 0, linetype="dashed") +

  # add ribbon
  geom_ribbon(data=bind_rows(d.ribbon_pgr.intrinsic, d.ribbon_pgr.invasion), 
              aes(x=elevation, ymin=ymin, ymax=ymax), 
              col="white", alpha=0.4, show.legend = FALSE) +
  
  # add mean points
  stat_summary(data = filter(vr.mean, !(background.species %in% c("site","site_none","species", "species_none"))),
               aes(x=elevation, y=pgr.log.median), 
               fun.data = mean_se, show.legend = FALSE) +
  
  # general trend of highland species
  geom_segment(data = bind_rows(d.seg.intrinsic, d.seg.invasion), 
               aes(x=x, y=y, xend=xend, yend = yend), 
               size=1, show.legend = FALSE) +
  coord_cartesian(ylim = c(-1,2.4)) +
  scale_color_manual(values = c("Highland" = "#3E4A89FF",  "Lowland" = "#FB9A06FF")) +
  scale_fill_manual(values = c("Highland" = "#3E4A89FF",  "Lowland" = "#FB9A06FF")) +
  scale_x_continuous(name = "Elevation (m)", breaks = c(890,1400, 1900)) +
  scale_y_continuous(name = "ln(Population growth rate)") +
  facet_wrap(~neighbor, scales = "free")
fig.pgr_boot

#****************************************************************************
# 6. Coexistence outcomes using bootstrapped lambdas ------
#****************************************************************************

#*************************************************
# ** - 6.1 Calculate sensitivity -------
#*************************************************
# make columns for sensitivity
# Intrinsic lambdas
vr.intrinsic <- filter(vr.boot, background.species == "none")
vr.intrinsic$ID_focal.species <- paste(vr.intrinsic$focal.species, vr.intrinsic$site, vr.intrinsic$bootstrap)

# match intrinsic lambdas
vr.boot$ID_focal.species <- paste(vr.boot$focal.species, vr.boot$site, vr.boot$bootstrap)
vr.boot$pgr.intrinsic.focal <- vr.intrinsic[match(vr.boot$ID_focal.species, vr.intrinsic$ID_focal.species), ]$pgr
vr.boot

# set ceiling and floor for PGRs
vr.boot[vr.boot$pgr > 20 & !is.na(vr.boot$pgr), "pgr"] = 20
vr.boot[vr.boot$pgr < 0.01 & !is.na(vr.boot$pgr), "pgr"] = 0.01

# calculate sensitivty
vr.boot <- as.data.frame(vr.boot)
vr.boot$sensitivity.raw2 <- as.numeric(NA)
vr.boot[vr.boot$pgr.intrinsic.focal > 1 & !is.na(vr.boot$pgr.intrinsic.focal), "sensitivity.raw2"] = 1- log(vr.boot[vr.boot$pgr.intrinsic.focal > 1 & !is.na(vr.boot$pgr.intrinsic.focal), ]$pgr) / log(vr.boot[vr.boot$pgr.intrinsic.focal > 1 & !is.na(vr.boot$pgr.intrinsic.focal), ]$pgr.intrinsic.focal)
vr.boot[vr.boot$pgr.intrinsic.focal < 1 & !is.na(vr.boot$pgr.intrinsic.focal), "sensitivity.raw2"] = log(vr.boot[vr.boot$pgr.intrinsic.focal < 1 & !is.na(vr.boot$pgr.intrinsic.focal), ]$pgr) / log(vr.boot[vr.boot$pgr.intrinsic.focal < 1 & !is.na(vr.boot$pgr.intrinsic.focal), ]$pgr.intrinsic.focal) - 1

#*************************************************
# ** - 6.2 Coexistence outcome -----
#*************************************************
# outcomes
out <- read_excel("coexistence_bootstrap_null.xlsx")
out <- out[-1,]
out$site <- factor(out$site, levels = c("Les Posses","Solalex", "Anzeindaz"))
out$origin.pair = dplyr::recode(out$origin.pair, "Highland_Highland" = "Highland-highland", "Lowland_Lowland" = "Lowland-lowland", "Lowland_Highland" = "Lowland-highland")

# make a new dataframe for bootstrapped outcomes
out.boot <- NULL
for(i in 0:499) {
  out.i = out
  out.i$bootstrap = i
  out.boot <- rbind(out.boot, out.i)
}

# predict outcomes
for(a in 1:nrow(out.boot)) {
  # a = 1
  di <- out.boot[a, ]
  bt = di$bootstrap
  sps1 <- di$sps1
  sps2 <- di$sps2
  st <- as.character(di$site)
  type <- di$type
  origin1 <- di$origin.sps1
  origin2 <- di$origin.sps2
  
  # get IGR and sensitivity
  igr1.sps <- mean(filter(vr.boot, focal.species == sps1 & background.species == "none" & bootstrap == bt)$pgr, na.rm=TRUE)
  igr2.sps <- mean(filter(vr.boot, focal.species == sps2 & background.species == "none" & bootstrap == bt)$pgr, na.rm=TRUE)
  
  igr1 <- filter(vr.boot, focal.species == sps1 & background.species == "none" & site == st & bootstrap == bt)$pgr
  igr2 <- filter(vr.boot, focal.species == sps2 & background.species == "none" & site == st & bootstrap == bt)$pgr
  
  igr11 <- filter(vr.boot, focal.species == sps1 & background.species == sps1 & site == st & bootstrap == bt)$pgr
  igr22 <- filter(vr.boot, focal.species == sps2 & background.species == sps2 & site == st & bootstrap == bt)$pgr
  
  igr12 <- filter(vr.boot, focal.species == sps1 & background.species == sps2 & site == st & bootstrap == bt)$pgr
  igr21 <- filter(vr.boot, focal.species == sps2 & background.species == sps1 & site == st & bootstrap == bt)$pgr
  
  # get sensitivity
  # facilitation scored
  s12 <- filter(vr.boot, focal.species == sps1 & background.species == sps2 & site == st & bootstrap == bt)$sensitivity.raw2
  s21 <- filter(vr.boot, focal.species == sps2 & background.species == sps1 & site == st & bootstrap == bt)$sensitivity.raw2
  
  out.boot[a,"lambda.sps.1"] <- igr1.sps
  out.boot[a,"lambda.sps.2"] <- igr2.sps
  
  out.boot[a,"lambda1"] <- igr1
  out.boot[a,"lambda2"] <- igr2
  
  out.boot[a,"igr11"] <- igr11
  out.boot[a,"igr22"] <- igr22
  
  out.boot[a,"igr12"] <- igr12
  out.boot[a,"igr21"] <- igr21
  
  # skip when IGR is NA
  if(is.na(igr1) | is.na(igr2) | is.na(igr11) | is.na(igr22) | is.na(igr12) | is.na(igr21) | is.na(s12) | is.na(s21)) next
  
  # skip estreme values
  if(is.na(igr1) | is.na(igr2) | is.na(igr11) | is.na(igr22) | is.na(igr12) | is.na(igr21) | is.na(s12) | is.na(s21)) next
  
  # skip when sensitivity is NA
  if(is.na(s12) | is.na(s21)) next
  
  # species grows faster across three sites
  if(igr1.sps > igr2.sps) out.boot[a, "superior.lambda.sps"] <- "sps1"
  else out.boot[a, "superior.lambda.sps"] <- "sps2"
  
  # species grows faster in each site
  if(igr1 > igr2) out.boot[a, "superior.lambda"] <- "sps1"
  else out.boot[a, "superior.lambda"] <- "sps2"
  
  # score negative sensitivty
  if(s12 < 0) { s12 = 0.1; print("facilitation scored") }
  if(s21 < 0) { s21 = 0.1; print("facilitation scored") }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # coexistence by IGR
  si <- coex.igr(sps1, sps2, igr1, igr2, igr12, igr21)
  #s12 <- si$sensitivity[1]
  #s21 <- si$sensitivity[2]
  
  out.boot[a,"sensitivity.12"] <- s12
  out.boot[a,"sensitivity.21"] <- s21
  
  out.boot[a, "outcome.igr"] <- si$outcome[1]
  out.boot[a, "winner.igr"] <- si$winner[1]
  out.boot[a, "superior.igr"] <- si$superior[1]
  si <- NULL
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # coexistence by NDFD
  ni <- coex.ndfd(s12, s21)
  out.boot[a,"nd"] <- ni$ndfd[1]
  out.boot[a,"fd"] <- ni$ndfd[2]
  out.boot[a, "coexistence.metric"] <- ni$ndfd[3] # coexistence metric
  
  out.boot[a, "outcome.ndfd"] <- ni$outcome[1]
  out.boot[a, "winner.ndfd"] <- ni$winner[1]
  out.boot[a, "superior.ndfd"] <- ni$superior[1]
  
  # FD highland/lowland
  if(origin1 == "Lowland" & origin2 == "Highland") out.boot[a, "fd.highlow"] <- 1/ni$ndfd[2]
  else out.boot[a, "fd.highlow"] <- ni$ndfd[2]
  
  # superior origin
  if(ni$superior[1] == "sps1") out.boot[a, "origin.superior"] <- origin1
  else out.boot[a, "origin.superior"] <- origin2
  
  # whether faster-grower is superior? site level
  if(out.boot[a, "superior.ndfd"] == out.boot[a,"superior.lambda"]) out.boot[a, "superior.fast"] = "yes"
  else out.boot[a, "superior.fast"] = "no"
  
  # whether faster-grower is superior? species level
  if(out.boot[a, "superior.ndfd"] == out.boot[a,"superior.lambda.sps"]) out.boot[a, "superior.fast.sps"] = "yes"
  else out.boot[a, "superior.fast.sps"] = "no"
  
  # FD fast-growing/slow-growing
  if(out.boot[a,"superior.lambda.sps"] == "sps1") out.boot[a,"fd.fastslow"] = ni$ndfd[2]
  else out.boot[a,"fd.fastslow"] = 1/ni$ndfd[2]
  
  ni <- NULL
  print(a)
}

#*************************************************
# ** - 6.3 IGR outcomes across sites ------
#*************************************************
source("/Users/slyu/LVSM/R codes/NICH_coexistence_IPM_functions.R")

out.boot <- read_csv("SD3_outcome_bootstrap.csv")
out.boot$site <- factor(out.boot$site, levels=c("Les Posses", "Solalex", "Anzeindaz"))

#*************************************************
# how many pairs including facilitations
# 6282/52003 = 12%
out.boot %>%
  filter(!is.na(outcome.ndfd)) %>%
  filter(sensitivity.12 == 0.1 | sensitivity.21 == 0.1)

#*************************************************
# outcomes across the three sites
# hig-high
out.hh <- out.boot %>%
  filter(sensitivity.12 != 0.1 & sensitivity.21 != 0.1)  %>%
  filter(origin.pair == "Highland-highland") %>%
  split(.$bootstrap) %>%
  purrr::map(~.[c("outcome.igr","site")]) %>%
  purrr::map(table)

# low-low
out.ll <- out.boot %>%
  filter(sensitivity.12 != 0.1 & sensitivity.21 != 0.1)  %>%
  filter(origin.pair == "Lowland-lowland") %>%
  split(.$bootstrap) %>%
  purrr::map(~.[c("outcome.igr","site")]) %>%
  purrr::map(table)

# low-high
out.lh <- out.boot %>%
  filter(sensitivity.12 != 0.1 & sensitivity.21 != 0.1)  %>%
  filter(origin.pair == "Lowland-highland") %>%
  split(.$bootstrap) %>%
  purrr::map(~.[c("outcome.igr","site")]) %>%
  purrr::map(table)

#*************************************************
# Summary competitive outcomes
out.igr.boot <- NULL
for(i in 1:500) {
  # i = 1
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # high-high
  if("Coexistence" %in% rownames(out.hh[[i]])) coex.hhi <- out.hh[[i]]["Coexistence",] 
  else coex.hhi <- c(NA, NA, NA)
  
  if("Competitive exclusion" %in% rownames(out.hh[[i]])) exc.hhi <- out.hh[[i]]["Competitive exclusion",]
  else exc.hhi <- c(NA, NA, NA)
  
  if("Priority effect" %in% rownames(out.hh[[i]])) pri.hhi <- out.hh[[i]]["Priority effect",]
  else pri.hhi <- c(NA, NA, NA)
  out.hhi <- data.frame(bootstrap = i,
                        origin.pair = "Highland-highland",
                        elevation = rep(c(890, 1400, 1900),3),
                        outcome.igr = rep(c("Coexistence", "Competitive exclusion", "Priority effect"), each=3),
                        number.pair = c(coex.hhi, exc.hhi, pri.hhi)
  )
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # low-high
  if("Coexistence" %in% rownames(out.lh[[i]])) coex.lhi <- out.lh[[i]]["Coexistence",] 
  else coex.lhi <- c(NA, NA, NA)
  
  if("Competitive exclusion" %in% rownames(out.lh[[i]])) exc.lhi <- out.lh[[i]]["Competitive exclusion",]
  else exc.lhi <- c(NA, NA, NA)
  
  if("Priority effect" %in% rownames(out.lh[[i]])) pri.lhi <- out.lh[[i]]["Priority effect",]
  else pri.lhi <- c(NA, NA, NA)
  out.lhi <- data.frame(bootstrap = i,
                        origin.pair = "Lowland-highland",
                        elevation = rep(c(890, 1400, 1900),3),
                        outcome.igr = rep(c("Coexistence", "Competitive exclusion", "Priority effect"), each=3),
                        number.pair = c(coex.lhi, exc.lhi, pri.lhi)
  )
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # low-low
  if("Coexistence" %in% rownames(out.ll[[i]])) coex.lli <- out.ll[[i]]["Coexistence",] 
  else coex.lli <- c(NA, NA, NA)
  
  if("Competitive exclusion" %in% rownames(out.ll[[i]])) exc.lli <- out.ll[[i]]["Competitive exclusion",]
  else exc.lli <- c(NA, NA, NA)
  
  if("Priority effect" %in% rownames(out.ll[[i]])) pri.lli <- out.ll[[i]]["Priority effect",]
  else pri.lli <- c(NA, NA, NA)
  out.lli <- data.frame(bootstrap = i,
                        origin.pair = "Lowland-lowland",
                        elevation = rep(c(890, 1400, 1900),3),
                        outcome.igr = rep(c("Coexistence", "Competitive exclusion", "Priority effect"), each=3),
                        number.pair = c(coex.lli, exc.lli, pri.lli)
  )
  
  out.i <- rbind(out.lli, out.lhi, out.hhi)
  out.igr.boot <- rbind(out.igr.boot, out.i)
  
  out.i <- NA
}

# summarise
# show 2 outcomes
outcome.igr.summary <- out.igr.boot %>%
  mutate(outcome.igr = dplyr::recode(outcome.igr, "Priority effect" = "Competitive exclusion", "Competitive exclusion" = "Competitive exclusion", "Coexistence" = "Coexistence")) %>%
  group_by(origin.pair, elevation, outcome.igr) %>%
  dplyr::summarise(mean_ci_sd(number.pair, na.rm=TRUE)[1], mean_ci_sd(number.pair, na.rm=TRUE)[2], mean_ci_sd(number.pair, na.rm=TRUE)[3]) 

outcome.igr.summary$ymin[outcome.igr.summary$outcome.igr == "Competitive exclusion"] = with(outcome.igr.summary, ymin[outcome.igr == "Competitive exclusion"] + y[outcome.igr=="Coexistence"])
outcome.igr.summary$ymax[outcome.igr.summary$outcome.igr == "Competitive exclusion"] = with(outcome.igr.summary, ymax[outcome.igr == "Competitive exclusion"] + y[outcome.igr=="Coexistence"])

fig.outcome.igr_2outcomes <- outcome.igr.summary %>%
  mutate(elevation = factor(elevation, levels = c("890","1400","1900"))) %>%
  mutate(origin.pair = factor(origin.pair, levels = c("Lowland-lowland", "Lowland-highland", "Highland-highland"))) %>%
  mutate(outcome.igr = factor(outcome.igr, levels=c("Competitive exclusion", "Coexistence"))) %>%
  ggplot(aes(x=elevation, y = y, ymin=ymin, ymax=ymax, fill=outcome.igr)) + 
  scale_fill_manual(values=c("Coexistence" = "grey82", "Competitive exclusion"= "grey100" )) +
  geom_bar(stat = "identity",  position="stack", col="grey35", width = 0.7) +
  geom_errorbar(width=0.2, alpha=0.6, col="grey35") +
  facet_wrap(~origin.pair, scales = "free", nrow=1, strip.position = NULL) +
  scale_y_continuous(name = "Number of pairs", labels = scales::number_format(accuracy = 1)) +
  scale_x_discrete(name = "Elevation (m)")
fig.outcome.igr_2outcomes

#****************************************************
# ** - 6.4 NDFD outcomes across sites -----
#****************************************************

#************************************
# NDFD mean and SE----
out.mean <- filter(out.boot, bootstrap == 0)

out.mean$cm.mean <- 
  out.mean$cm.median <- 
  out.mean$cm.min <-
  out.mean$cm.max <-
  out.mean$cm.sample <- 
  
  out.mean$nd.mean <-
  out.mean$nd.median <-
  out.mean$nd.mean2 <-
  out.mean$nd.min <-
  out.mean$nd.max <-
  out.mean$nd.sample <-
  
  out.mean$fd.mean <-
  out.mean$fd.median <-
  out.mean$fd.mean2 <-
  out.mean$fd.min <- 
  out.mean$fd.max <- 
  out.mean$fd.sample  <- 
  
  out.mean$n.facilitation <- as.numeric(NA)

for(i in 1:nrow(out.mean)) {
  # i = 1
  di <- out.mean[i,]
  
  # a pair
  spi.1 <- di$sps1
  spi.2 <- di$sps2
  st.i <- as.character(di$site)
  
  # get the bootstrapped data
  out.i <- out.boot %>%
    filter(sps1 == spi.1 & sps2 == spi.2 & site == st.i) %>%
    filter(!is.infinite(coexistence.metric)) %>%
    filter(!is.infinite(nd)) %>%
    filter(!is.infinite(fd.highlow))
  
  # number of non-persisting and facilitations
  out.mean[i,"n.facilitation"] <- nrow(filter(out.i, !is.na(coexistence.metric) & (sensitivity.12 == 0.1 | sensitivity.21 == 0.1)))
  
  out.i <- out.boot %>%
    filter(sps1 == spi.1 & sps2 == spi.2 & site == st.i) %>%
    filter(!is.infinite(coexistence.metric)) %>%
    filter(!is.infinite(nd)) %>%
    filter(!is.infinite(fd.highlow)) %>%
    filter(sensitivity.12 != 0.1 & sensitivity.21 != 0.1)
  
  # skip when N < 3
  if(sum(!is.na(out.i$coexistence.metric)) < 3) next
  
  # CM
  cm.i <- mean_ci_quantile(log(out.i$coexistence.metric))
  out.mean[i, "cm.mean"] = cm.i$y
  out.mean[i, "cm.median"] = median_ci_quantile(log(out.i$coexistence.metric))$y
  out.mean[i, "cm.min"] = cm.i$ymin
  out.mean[i, "cm.max"] = cm.i$ymax
  out.mean[i, "cm.sample"] = sum(!is.na(out.i$coexistence.metric))
  rm(cm.i)
  
  # ND
  nd.i <- mean_ci_quantile(out.i$nd)
  out.mean[i, "nd.mean"] = nd.i$y
  out.mean[i, "nd.median"] = median_ci_quantile(out.i$nd)$y
  out.mean[i, "nd.min"] = nd.i$ymin
  out.mean[i, "nd.max"] = nd.i$ymax
  out.mean[i, "nd.sample"] = sum(!is.na(out.i$nd))
  out.mean[i, "nd.mean2"] = mean_ci_quantile(log(1-out.i$nd))$y
  rm(nd.i)
  
  # FD
  fd.i <- mean_ci_quantile(log(out.i$fd.highlow))
  out.mean[i, "fd.mean"] = fd.i$y
  out.mean[i, "fd.median"] = median_ci_quantile(log(out.i$fd.highlow))$y
  out.mean[i, "fd.min"] = fd.i$ymin
  out.mean[i, "fd.max"] = fd.i$ymax
  out.mean[i, "fd.sample"] = sum(!is.na(log(out.i$fd.highlow)))
  
  if(as.character(di$origin.pair) %in% c("Highland-highland", "Lowland-lowland")) {
    out.mean[i, "fd.mean2"] = mean_ci_quantile(abs(log(out.i$fd.highlow)))$y
  }
  else {
    out.mean[i, "fd.mean2"] = fd.i$y
  }
  
  rm(fd.i)
  
  print(i)
}
out.mean$ID.pair <- paste(out.mean$pair, out.mean$site)
out.boot$ID.pair <- paste(out.boot$pair, out.boot$site)
out.mean$origin.pair <- factor(out.mean$origin.pair, levels = c("Lowland-lowland", "Lowland-highland", "Highland-highland"))

hist(out.mean$cm.sample)
hist(out.mean$nd.sample)
hist(out.mean$fd.sample)

#****************************
# CM-----

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Distribtuion
hist(out.boot$coexistence.metric, breaks=100,labels = TRUE)
hist(log(out.boot$coexistence.metric), breaks=100,labels = TRUE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# test based on mean
out.mean %>%
  filter(sensitivity.12 != 0.1 & sensitivity.21 != 0.1) %>%
  lmer(cm.mean ~ origin.pair * elevation + (1|pair), data=.) %>%
  Anova()
  summary()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# test for each bootstrap
test.cm.boot <- out.boot %>% 
  filter(sensitivity.12 != 0.1 & sensitivity.21 != 0.1) %>%
  filter(!is.infinite(coexistence.metric)) %>%
  mutate(cm.log = log(coexistence.metric)) %>% filter(!is.infinite(cm.log)) %>%
  split(.$bootstrap) %>%
  purrr::map(~lmer(cm.log ~ elevation*origin.pair + (1|pair), data=.)) 

# coefficient for each bootstrap
coef.cm.boot <- 
  test.cm.boot %>%
  purrr::map(fixef)
coef.cm.boot

# predict using coefficients
pred.cm <- NULL
coef.cm.global <- data.frame(bootstrap = 1:500,
                             int.highhigh = NA,
                             slo.highhigh = NA,
                             int.lowhigh = NA,
                             slo.lowhigh = NA,
                             int.lowlow = NA,
                             slo.lowlow = NA)

for(i in 1:length(coef.cm.boot)) {
  # i = 1
  coef.i <- coef.cm.boot[[i]]
  # elevation.sequence
  ee = seq(890,1900, length.out = 100)
  
  # extract bootstrapped linear coefficients
  coef.cm.global[i,"int.highhigh"] = coef.i["(Intercept)"]
  coef.cm.global[i,"slo.highhigh"] = coef.i["elevation"]
  coef.cm.global[i,"int.lowhigh"] = coef.i["(Intercept)"] +  coef.i["origin.pairLowland-highland"] 
  coef.cm.global[i,"slo.lowhigh"] = coef.i["elevation"] +  coef.i["elevation:origin.pairLowland-highland"] 
  coef.cm.global[i,"int.lowlow"] = coef.i["(Intercept)"] +  coef.i["origin.pairLowland-lowland"] 
  coef.cm.global[i,"slo.lowlow"] = coef.i["elevation"] +  coef.i["elevation:origin.pairLowland-lowland"] 
  
  # highland-highland
  pre.hh = coef.cm.global[i,"int.highhigh"] + coef.cm.global[i,"slo.highhigh"]*ee
  # lowland-highland
  pre.lh = coef.cm.global[i,"int.lowhigh"] + coef.cm.global[i,"slo.lowhigh"]*ee
  # lowland-lowland
  pre.ll = coef.cm.global[i,"int.lowlow"] + coef.cm.global[i,"slo.lowlow"]*ee
  
  # combine
  pred.i <- data.frame(bootstrap = i,
                       elevation = c(ee,ee,ee),
                       cm = c(pre.hh, pre.lh, pre.ll),
                       origin.pair =rep(c("Highland-highland", "Lowland-highland", "Lowland-lowland"), each=100) 
  )
  pred.cm <- rbind(pred.cm, pred.i)
}
pred.cm$origin.pair = factor(pred.cm$origin.pair, levels = c("Lowland-lowland", "Lowland-highland", "Highland-highland") )
coef.cm.global
pred.cm

# CI of coefficients
coef.cm.interaction <- test.cm.boot %>%
  purrr::map_dfr(fixef) 

# high-high: not significant
#  0.0005666086 -3.884874e-05 0.001517295
hist(coef.cm.interaction$elevation)
median_ci_quantile(coef.cm.interaction$elevation)
# low-low: not significant
#  -0.0006775102 -0.001718274 4.30188e-05
hist(coef.cm.interaction$`elevation:origin.pairLowland-lowland`)
median_ci_quantile(coef.cm.interaction$`elevation:origin.pairLowland-lowland`)
# low-high: not significant
# -0.0004085 -0.0009613982 8.342681e-05
hist(coef.cm.interaction$`elevation:origin.pairLowland-highland`)
median_ci_quantile(coef.cm.interaction$`elevation:origin.pairLowland-highland`)

# add confidecne interval
# CI of SD
d.ribbon.cm <- pred.cm %>%
  group_by(elevation, origin.pair) %>%
  dplyr::summarise(mean_ci_sd(cm)[2], mean_ci_sd(cm)[3])
d.ribbon.cm$origin.pair = factor(d.ribbon.cm$origin.pair, levels = c("Lowland-lowland", "Lowland-highland", "Highland-highland") )

# add mean trend
coef.cm.mean = apply(coef.cm.global,2,mean)[2:7]
d.seg.cm <- data.frame(x=890,
                       xend=1900,
                       y=c(coef.cm.mean["int.highhigh"] + coef.cm.mean["slo.highhigh"]*890, 
                           coef.cm.mean["int.lowhigh"] + coef.cm.mean["slo.lowhigh"]*890,
                           coef.cm.mean["int.lowlow"] + coef.cm.mean["slo.lowlow"]*890),
                       yend=c(coef.cm.mean["int.highhigh"] + coef.cm.mean["slo.highhigh"]*1900,
                              coef.cm.mean["int.lowhigh"] + coef.cm.mean["slo.lowhigh"]*1900,
                              coef.cm.mean["int.lowlow"] + coef.cm.mean["slo.lowlow"]*1900),
                       origin.pair = c("Highland-highland", "Lowland-highland", "Lowland-lowland"),
                       significant = "significant")
d.seg.cm$origin.pair = factor(d.seg.cm$origin.pair, levels = c("Lowland-lowland", "Lowland-highland", "Highland-highland") )

fig.cm_boot <- pred.cm %>%
  ggplot(aes(fill = origin.pair, col=origin.pair)) +
  geom_hline(yintercept = 0, linetype="dashed") +
  #geom_path(aes(x=elevation, y=cm, group=bootstrap), col="grey", alpha=0.5) + 
  
  # CI
  geom_ribbon(data=d.ribbon.cm, aes(x=elevation, ymin=ymin, ymax=ymax), alpha=0.4, col="white", show.legend = FALSE) +
  
  # mean points
  stat_summary(data = filter(out.mean, (sensitivity.12 != 0.1 & sensitivity.21 != 0.1)),
               aes(y=cm.median, x=elevation), fun.data=mean_se, show.legend = FALSE) +
  
  # general trend
  geom_segment(data=d.seg.cm, aes(x=x, xend=xend, y= y, yend = yend), size=1, show.legend = FALSE) +
  facet_wrap(~origin.pair, scales="free") +
  #coord_cartesian(ylim=c(-1.4,0.6)) +
  scale_fill_manual(values=c("Lowland-lowland" = "#FB9A06FF", "Lowland-highland" = "#6DCD59FF", "Highland-highland" = "#3E4A89FF")) +
  scale_color_manual(values=c("Lowland-lowland" = "#FB9A06FF", "Lowland-highland" = "#6DCD59FF", "Highland-highland" = "#3E4A89FF")) +
  scale_x_continuous(name = "Elevation (m)", breaks = c(890,1400, 1900)) +
  scale_y_continuous(name = "ln(Coexistence metric)",labels = scales::number_format(accuracy = 0.1))
fig.cm_boot

#****************************
# ND----

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# test based on mean
out.mean %>%
  filter(sensitivity.12 != 0.1 & sensitivity.21 != 0.1) %>%
  lmer(nd.mean2 ~ origin.pair * elevation + (1|pair), data=.) %>%
  Anova()
  summary()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# test for each bootstrap
test.nd.boot <- out.boot %>% 
  filter(sensitivity.12 != 0.1 & sensitivity.21 != 0.1) %>%
  mutate(no = 1-nd) %>% mutate(nd.log = log(no)) %>% filter(!is.infinite(nd.log)) %>%
  split(.$bootstrap) %>%
  purrr::map(~lmer(nd.log ~ elevation*origin.pair + (1|pair), data=.)) 

# ANOVA
anova.nd.boot <- test.nd.boot %>% purrr::map_dfr(coef.lmer)
p.value_nd <- anova.nd.boot[(1:500)*3, 3]

# coefficient for each bootstrap
coef.nd.boot <- test.nd.boot %>% purrr::map(fixef)
coef.nd.boot

# predict using coefficients
pred.nd <- NULL
coef.nd.global <- data.frame(bootstrap = 1:500,
                             int.highhigh = NA,
                             slo.highhigh = NA,
                             int.lowhigh = NA,
                             slo.lowhigh = NA,
                             int.lowlow = NA,
                             slo.lowlow = NA)

for(i in 1:length(coef.nd.boot)) {
  # i = 1
  coef.i <- coef.nd.boot[[i]]
  # elevation.sequence
  ee = seq(890,1900, length.out = 100)
  
  # extract bootstrapped linear coefficients
  coef.nd.global[i,"int.highhigh"] = coef.i["(Intercept)"]
  coef.nd.global[i,"slo.highhigh"] = coef.i["elevation"]
  coef.nd.global[i,"int.lowhigh"] = coef.i["(Intercept)"] +  coef.i["origin.pairLowland-highland"] 
  coef.nd.global[i,"slo.lowhigh"] = coef.i["elevation"] +  coef.i["elevation:origin.pairLowland-highland"] 
  coef.nd.global[i,"int.lowlow"] = coef.i["(Intercept)"] +  coef.i["origin.pairLowland-lowland"] 
  coef.nd.global[i,"slo.lowlow"] = coef.i["elevation"] +  coef.i["elevation:origin.pairLowland-lowland"] 
  
  # highland-highland
  pre.hh = coef.nd.global[i,"int.highhigh"] + coef.nd.global[i,"slo.highhigh"]*ee
  # lowland-highland
  pre.lh = coef.nd.global[i,"int.lowhigh"] + coef.nd.global[i,"slo.lowhigh"]*ee
  # lowland-lowland
  pre.ll = coef.nd.global[i,"int.lowlow"] + coef.nd.global[i,"slo.lowlow"]*ee
  
  # combine
  pred.i <- data.frame(bootstrap = i,
                       elevation = c(ee,ee,ee),
                       nd = c(pre.hh, pre.lh, pre.ll),
                       origin.pair =rep(c("Highland-highland", "Lowland-highland", "Lowland-lowland"), each=100) 
  )
  pred.nd <- rbind(pred.nd, pred.i)
}
pred.nd$origin.pair = factor(pred.nd$origin.pair, levels = c("Lowland-lowland", "Lowland-highland", "Highland-highland"))

# CI coefficients
coef.nd.interaction <- test.nd.boot %>%
  purrr::map_dfr(fixef) 
# high-high: not significant
# -0.0003718221 -0.001065897 0.0001266079
hist(coef.nd.interaction$elevation)
median_ci_quantile(coef.nd.interaction$elevation)
# low-low: not significant
# 0.0005340055 -5.87114e-05 0.001244548
hist(coef.nd.interaction$`elevation:origin.pairLowland-lowland`)
median_ci_quantile(coef.nd.interaction$`elevation:origin.pairLowland-lowland`)
# low-high: not significant
# 0.0002187646 -0.0002948585 0.0006438424
hist(coef.nd.interaction$`elevation:origin.pairLowland-highland`)
median_ci_quantile(coef.nd.interaction$`elevation:origin.pairLowland-highland`)

# plot with CI
# add mean trend
coef.nd.mean = apply(coef.nd.global,2,mean)[2:7]
d.seg.nd <- data.frame(x=890,
                       xend=1900,
                       y=c(coef.nd.mean["int.highhigh"] + coef.nd.mean["slo.highhigh"]*890, 
                           coef.nd.mean["int.lowhigh"] + coef.nd.mean["slo.lowhigh"]*890,
                           coef.nd.mean["int.lowlow"] + coef.nd.mean["slo.lowlow"]*890),
                       yend=c(coef.nd.mean["int.highhigh"] + coef.nd.mean["slo.highhigh"]*1900,
                              coef.nd.mean["int.lowhigh"] + coef.nd.mean["slo.lowhigh"]*1900,
                              coef.nd.mean["int.lowlow"] + coef.nd.mean["slo.lowlow"]*1900),
                       origin.pair = c("Highland-highland", "Lowland-highland", "Lowland-lowland"),
                       significant = "significant")
d.seg.nd$origin.pair <- factor(d.seg.nd$origin.pair, levels = c("Lowland-lowland", "Lowland-highland", "Highland-highland"))

# add ribbon
# SD CI
d.ribbon.nd <- pred.nd %>%
  group_by(elevation, origin.pair) %>%
  dplyr::summarise(mean_ci_sd(nd)[2], mean_ci_sd(nd)[3])
d.ribbon.nd$origin.pair <- factor(d.ribbon.nd$origin.pair, levels = c("Lowland-lowland", "Lowland-highland", "Highland-highland"))

fig.nd_boot <- pred.nd %>%
  ggplot(aes(fill=origin.pair, color=origin.pair)) +
  #geom_path(aes(x=elevation, y=nd, group=bootstrap), col="grey", alpha=0.5) + 
  
  # CI
  geom_ribbon(data=d.ribbon.nd, aes(x=elevation, ymin=ymin, ymax=ymax), alpha=0.4, col="white", show.legend = FALSE) +
  
  # mean points
  stat_summary(data = filter(out.mean,  (sensitivity.12 != 0.1 & sensitivity.21 != 0.1)),
               aes(y=nd.mean2, x=elevation), fun.data = mean_se, show.legend = FALSE) +
  
  # general trend
  geom_segment(data=d.seg.nd, aes(x=x, xend=xend, y= y, yend = yend), size=1, show.legend = FALSE) +
  facet_wrap(~origin.pair, scales="free") +
  #coord_cartesian(ylim=c(-1.3,0.6)) +
  scale_fill_manual(values=c("Lowland-lowland" = "#FB9A06FF", "Lowland-highland" = "#6DCD59FF", "Highland-highland" = "#3E4A89FF")) +
  scale_color_manual(values=c("Lowland-lowland" = "#FB9A06FF", "Lowland-highland" = "#6DCD59FF", "Highland-highland" = "#3E4A89FF")) +
  scale_x_continuous(name = "Elevation (m)", breaks = c(890,1400, 1900)) +
  scale_y_continuous(name = "ln(Niche overlap)",labels = scales::number_format(accuracy = 0.1))
fig.nd_boot

#****************************
# FD----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# test based on mean

# overall
out.mean %>%
  filter(sensitivity.12 != 0.1 & sensitivity.21 != 0.1) %>%
  lmer(abs(fd.mean2) ~ origin.pair * elevation + (1|pair), data=.) %>%
  Anova()
  summary()

# high-high
out.mean %>%
  filter(sensitivity.12 != 0.1 & sensitivity.21 != 0.1) %>%
  filter(origin.pair %in% c("Highland-highland")) %>%
  lmer(abs(fd.mean2) ~ elevation + (1|pair), data=.) %>%
  Anova()
  summary()

# low-low
out.mean %>%
  filter(sensitivity.12 != 0.1 & sensitivity.21 != 0.1) %>%
  filter(origin.pair %in% c("Lowland-lowland")) %>%
  lmer(abs(fd.mean) ~ elevation + (1|pair), data=.) %>%
  Anova()
  summary

# low-high
out.mean %>%
  filter(sensitivity.12 != 0.1 & sensitivity.21 != 0.1) %>%
  filter(origin.pair %in% c("Lowland-highland")) %>%
  lmer(fd.mean ~ elevation + (1|pair), data=.) %>%
  Anova()
  summary()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# lmer test for each bootstrap
#~~~~~~~~~~~~~~~~~~~~~~~~
# high-high
test.fd.boot_hh <- out.boot %>% 
  filter(sensitivity.12 != 0.1 & sensitivity.21 != 0.1) %>%
  filter(origin.pair == "Highland-highland") %>%
  mutate(fd.log = abs(log(fd.highlow))) %>%
  filter(!is.infinite(fd.log)) %>%
  split(.$bootstrap) %>%
  purrr::map(~lmer(fd.log ~ elevation + (1|pair), data=.)) 

# coefficients
# -0.0002271186 -0.0006686035 0.000119544
coef.fd.boot_hh <- test.fd.boot_hh %>% purrr::map_dfr(fixef)
median_ci_quantile(coef.fd.boot_hh$elevation)

coef.fd.boot_hh <- test.fd.boot_hh %>%
  purrr::map(fixef)

#~~~~~~~~~~~~~~~~~~~~~~~~
# low-high
test.fd.boot_lh <- out.boot %>%
  filter(sensitivity.12 != 0.1 & sensitivity.21 != 0.1) %>%
  filter(origin.pair =="Lowland-highland") %>%
  filter(!is.na(fd.highlow)) %>%
  mutate(fd.log = log(fd.highlow)) %>%
  filter(!is.infinite(fd.log)) %>%
  split(.$bootstrap) %>%
  purrr::map(~lmer(fd.log ~ elevation + (1|pair), data=.)) 

# coefficients
# 0.0001346951 -0.0002440668 0.0005411766
coef.fd.boot_lh <- test.fd.boot_lh %>% purrr::map_dfr(fixef)
median_ci_quantile(coef.fd.boot_lh$elevation)

coef.fd.boot_lh <- test.fd.boot_lh %>%
  purrr::map(fixef)

#~~~~~~~~~~~~~~~~~~~~~~~~
# low-low
test.fd.boot_ll <- out.boot %>% 
  filter(sensitivity.12 != 0.1 & sensitivity.21 != 0.1) %>%
  filter(origin.pair == "Lowland-lowland") %>%
  mutate(fd.log = abs(log(fd.highlow))) %>%
  filter(!is.infinite(fd.log)) %>%
  split(.$bootstrap) %>%
  purrr::map(~lmer(fd.log ~ elevation + (1|pair), data=.)) 

# coefficients
# -6.382185e-05 -0.0003379845 0.0001890299
coef.fd.boot_ll <- test.fd.boot_ll %>% purrr::map_dfr(fixef)
median_ci_quantile(coef.fd.boot_ll$elevation)

coef.fd.boot_ll <- test.fd.boot_ll %>%
  purrr::map(fixef)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# predict using coefficients
pred.fd <- NULL
coef.fd.global <- data.frame(bootstrap = 1:500,
                             int.highhigh = NA,
                             slo.highhigh = NA,
                             int.lowhigh = NA,
                             slo.lowhigh = NA,
                             int.lowlow = NA,
                             slo.lowlow = NA)

for(i in 1:length(coef.fd.boot_hh)) {
  # i = 1
  
  # elevation.sequence
  ee = seq(890,1900, length.out = 100)
  
  # extract bootstrapped linear coefficients
  coef.fd.global[i,"int.highhigh"] = coef.fd.boot_hh[[i]]["(Intercept)"]
  coef.fd.global[i,"slo.highhigh"] = coef.fd.boot_hh[[i]]["elevation"]
  coef.fd.global[i,"int.lowhigh"] = coef.fd.boot_lh[[i]]["(Intercept)"]
  coef.fd.global[i,"slo.lowhigh"] = coef.fd.boot_lh[[i]]["elevation"]
  coef.fd.global[i,"int.lowlow"] = coef.fd.boot_ll[[i]]["(Intercept)"]
  coef.fd.global[i,"slo.lowlow"] = coef.fd.boot_ll[[i]]["elevation"]
  
  # highland-highland
  pre.hh = coef.fd.global[i,"int.highhigh"] + coef.fd.global[i,"slo.highhigh"]*ee
  # lowland-highland
  pre.lh = coef.fd.global[i,"int.lowhigh"] + coef.fd.global[i,"slo.lowhigh"]*ee
  # lowland-lowland
  pre.ll = coef.fd.global[i,"int.lowlow"] + coef.fd.global[i,"slo.lowlow"]*ee
  
  # combine
  pred.i <- data.frame(bootstrap = i,
                       elevation = c(ee,ee,ee),
                       fd = c(pre.hh, pre.lh, pre.ll),
                       origin.pair =rep(c("Highland-highland", "Lowland-highland", "Lowland-lowland"), each=100),
                       significant = rep(c("no", "significant", "no"), each = 100))
  pred.fd <- rbind(pred.fd, pred.i)
}
pred.fd
pred.fd$origin.pair =  factor(pred.fd$origin.pair, levels = c("Lowland-lowland", "Lowland-highland", "Highland-highland") )

# CI of slopes
hist(coef.fd.global$slo.highhigh)
hist(coef.fd.global$slo.lowhigh)
hist(coef.fd.global$slo.lowlow)

fd.slope.ci.hh <-  mean_ci_quantile(coef.fd.global$slo.highhigh)
fd.slope.ci.lh <- mean_ci_quantile(coef.fd.global$slo.lowhigh)
fd.slope.ci.ll <- mean_ci_quantile(coef.fd.global$slo.lowlow)

# combine three types
fd.slope <- as.data.frame(rbind(fd.slope.ci.hh, fd.slope.ci.lh, fd.slope.ci.ll))
colnames(fd.slope) <- c("y", "ymin", "ymax")
fd.slope$origin.pair = c("Highland-highland", "Lowland-highland", "Lowland-lowland")
fd.slope$origin.pair = factor(fd.slope$origin.pair, levels = c("Lowland-lowland", "Lowland-highland", "Highland-highland")) 

fig.fd_slope <- fd.slope %>%
  ggplot(aes(x=origin.pair, y=y, ymin=ymin, ymax=ymax)) +
  geom_hline(yintercept = 0,linetype="dashed") +
  geom_pointrange() +
  #geom_errorbar() +
  scale_x_discrete(name="Type of interactions") +
  scale_y_continuous(name="Slope")
fig.fd_slope

# plot with CI
# add mean trend
coef.fd.mean = apply(coef.fd.global,2,mean)[2:7]
d.seg.fd <- data.frame(x=890,
                       xend=1900,
                       y=c(coef.fd.mean["int.highhigh"] + coef.fd.mean["slo.highhigh"]*890, 
                           coef.fd.mean["int.lowhigh"] + coef.fd.mean["slo.lowhigh"]*890,
                           coef.fd.mean["int.lowlow"] + coef.fd.mean["slo.lowlow"]*890),
                       yend=c(coef.fd.mean["int.highhigh"] + coef.fd.mean["slo.highhigh"]*1900,
                              coef.fd.mean["int.lowhigh"] + coef.fd.mean["slo.lowhigh"]*1900,
                              coef.fd.mean["int.lowlow"] + coef.fd.mean["slo.lowlow"]*1900),
                       origin.pair = c("Highland-highland", "Lowland-highland", "Lowland-lowland"),
                       significant = c("significant", "significant", "no"))
d.seg.fd$origin.pair = factor(d.seg.fd$origin.pair, levels = c("Lowland-lowland", "Lowland-highland", "Highland-highland") )

# add ribbon
d.ribbon.fd <- pred.fd %>%
  group_by(elevation, origin.pair) %>%
  dplyr::summarise(mean_ci_sd(fd)[2], mean_ci_sd(fd)[3])
d.ribbon.fd$origin.pair = factor(d.ribbon.fd$origin.pair, levels = c("Lowland-lowland", "Lowland-highland", "Highland-highland") )

fig.fd_boot <- pred.fd %>%
  ggplot(aes(fill=origin.pair, col=origin.pair)) +
  geom_hline(yintercept = 0, linetype="dashed") +
  #geom_path(aes(x=elevation, y=fd, group=bootstrap), col="grey", alpha=0.5) + 
  
  # CI
  geom_ribbon(data=d.ribbon.fd, aes(x=elevation, ymin=ymin, ymax=ymax), alpha=0.4, col="white", show.legend = FALSE) +
  
  # mean points
  stat_summary(data = filter(out.mean,(sensitivity.12 != 0.1 & sensitivity.21 != 0.1)),
               aes(y=fd.mean2, x=elevation),fun.data=mean_se, show.legend = FALSE) +
  
  # general trend
  geom_segment(data=d.seg.fd, aes(x=x, xend=xend, y= y, yend = yend, linetype = significant), size=1, show.legend = FALSE) +
  
  # facet
  facet_wrap(~origin.pair, scale="free") +
  scale_fill_manual(values=c("Lowland-lowland" = "#FB9A06FF", "Lowland-highland" = "#6DCD59FF", "Highland-highland" = "#3E4A89FF")) +
  scale_color_manual(values=c("Lowland-lowland" = "#FB9A06FF", "Lowland-highland" = "#6DCD59FF", "Highland-highland" = "#3E4A89FF")) +
  scale_linetype_manual(values = c("significant" = "solid", "no" = "dotted")) +
  #coord_cartesian(ylim=c(-0.5,1)) +
  scale_x_continuous(name = "Elevation (m)", breaks = c(890,1400, 1900)) +
  scale_y_continuous(name = "ln(Fitness difference)", labels = scales::number_format(accuracy = 0.1))
fig.fd_boot

#****************************************************
# 7. figures -----
#****************************************************

#****************************************************
# ** - Figure 1: Population growth rate ####
# Combined
p <- fig.pgr_boot +
  theme_classic2(base_family = "Arial") +
  theme(axis.text=element_text(size=11),
        axis.title=element_text(size=12),
        strip.background = element_blank(),
        strip.text = element_text(size=12),
        strip.placement = "outside",
        panel.grid.major = element_line(colour="grey90"),
        legend.position = "none")

#****************************************************
# ** - Figure 2: coexistence outcomes across sites ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# show 2 outcomes
p <- fig.outcome.igr_2outcomes +
  theme_classic2(base_family = "Arial") +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=11),
        strip.background = element_blank(),
        strip.text = element_text(size=11),
        strip.placement = "outside",
        panel.grid.major = element_line(colour="grey90"),
        legend.position = "none")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CM
p <- fig.cm_boot +
  theme_classic2(base_family = "Arial") +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=11),
        strip.background = element_blank(),
        strip.text = element_text(size=11),
        strip.placement = "outside",
        panel.grid.major = element_line(colour="grey90"),
        legend.position = "none")

#****************************************************
# ** - Figure 3: ND and RFD across sites ####

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ND
p <- fig.nd_boot +
  theme_classic2(base_family = "Arial") +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=12),
        strip.background = element_blank(),
        strip.text = element_text(size=11),
        #strip.text.x = element_text(hjust = 0),
        strip.placement = "outside",
        panel.grid.major = element_line(colour="grey90"),
        legend.position = "none")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# fd
p <- fig.fd_boot +
  theme_classic2(base_family = "Arial") +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=12),
        strip.background = element_blank(),
        strip.text = element_text(size=11),
        #strip.text.x = element_text(hjust = 0),
        strip.placement = "outside",
        panel.grid.major = element_line(colour="grey90"),
        legend.position = "none")

