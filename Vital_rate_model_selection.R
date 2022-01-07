#########################################################################
# Vital rate model selection
# Author: Shengman Lyu (shengman.lyu@usys.ethz.ch)
# Date: 07.01.2022
#########################################################################
library(readxl)
library(tidyverse)
library(lme4)
library(lmerTest)

#*****************************************************************************
# ****** Model comparison: interactions of site and background species******----
#*****************************************************************************
rm(list=ls())

# data
d <- read_excel("/Users/slyu/LVSM/NICH/Data/Biomass/data_181920_v6.xlsx", na="NA",col_names = TRUE)
d
# remove seedling planted in autumn 2019 because mostly of them died: 16755
d <- d %>% filter(!(!is.na(note.size.autumn0) & note.size.autumn0 == "planting.autumn0" & year == 2020))
d$site <- factor(d$site)
d

# focal species
sps <- c("Anal","Armo","Asal","Plal", "Poal", "Seca", "Trba",      # 7 alpine species 
         "Brer", "Crbi", "Daca","Melu" , "Plla", "Potr", "Sapr", # 7 lowland species
         "none")   

# functional to compare candidate models
modelcmp.glm <- function(data, response, candidate, ...) {
  model.candidate <- paste(response, candidate, sep="~")
  aic <- NULL
  for(i in 1:length(model.candidate)) {
    #i =17
    fit <- glm(formula = as.formula(model.candidate[i]), data = data,...)
    aic[i] <- AIC(fit)
  }
  # delta AIC
  delta.aic <- aic - min(aic)
  
  # AIC weight
  weight.aic <- (exp(-0.5*delta.aic))/sum(exp(-0.5*delta.aic))
  
  out <- data.frame(y = response, n=nrow(data), models = candidate, aic = aic, aic.delta =delta.aic, aic.weight = weight.aic)
  out <- out[order(out$aic),]
  return(out)
}

modelcmp.lm <- function(data, response, candidate, ...) {
  model.candidate <- paste(response, candidate,sep="~")
  aic <- NULL
  for(i in 1:length(model.candidate)) {
    #i =17
    fit <- lm(formula = as.formula(model.candidate[i]), data = data, ...)
    aic[i] <- AIC(fit)
  }
  # delta AIC
  delta.aic <- aic - min(aic)
  
  # AIC weight
  weight.aic <- (exp(-0.5*delta.aic))/sum(exp(-0.5*delta.aic))
  
  out <- data.frame(y = response, n=nrow(data), models = candidate, aic = aic, aic.delta =delta.aic, aic.weight = weight.aic)
  out <- out[order(out$aic),]
  return(out)
}

# candidate models without size
models.1 <- c("1", 
              "site",
              "background.species",
              "site + background.species",
              "site*background.species")

# candidate models with size
models.size <- c("1", 
                 "size.autumn0",
                 #"I(size.autumn0^2)",
                 #"size.autumn0 + I(size.autumn0^2)",
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

models.size.log <- c("1", 
                     "log(size.autumn0)",
                     #"I(size.autumn0^2)",
                     #"size.autumn0 + I(size.autumn0^2)",
                     "site",
                     "background.species",
                     "log(size.autumn0) + site",
                     "log(size.autumn0) + background.species",
                     "site + background.species",
                     "log(size.autumn0) + site + background.species",
                     "log(size.autumn0)*site + background.species",
                     "log(size.autumn0)*background.species + site",
                     "log(size.autumn0) + site*background.species",
                     "log(size.autumn0) * site * background.species")

# function to select pairs 
sample.pair <- function(data, n=1) {
  bg.lp <- names(table(data[,c("background.species", "site")])[,1])[table(data[,c("background.species", "site")])[,1]>n]
  bg.sl <- names(table(data[,c("background.species", "site")])[,1])[table(data[,c("background.species", "site")])[,2]>n]
  bg.az <- names(table(data[,c("background.species", "site")])[,1])[table(data[,c("background.species", "site")])[,3]>n]
  
  out <- filter(data, 
                site == "Les Posses" & background.species %in% bg.lp | 
                  site == "Solalex" & background.species %in% bg.sl |
                  site == "Anzeindaz" & background.species %in% bg.az)
  return(out)
}

## Include competition-only focals
# d <- filter(d, background.species != "none")
d

aic <- aic.best <- NULL

for(fc in sps[1:14]) {
  #fc <- "Sapr"
  
  #**************************************************
  ## 1 survival----
  d.sur <- d %>%
    filter(focal.before == fc &
             background.species %in% sps &
             !is.na(size.autumn0) & !is.na(survival) &
             !(!is.na(note.size.autumn0) & note.size.autumn0 == "planting.autumn0" & year == 2020) &
             (!is.na(if.same) & if.same==TRUE) )
  sample.pair(d.sur, n=0)
  
  if(nrow(sample.pair(d.sur)) > 5) {
    # AIC comparison
    # log size
    #aic.sur <- modelcmp.glm(response = "survival", candidate = models.size.log, d.sur, family="binomial")
    aic.sur <- modelcmp.glm(response = "survival", candidate = models.size, data=sample.pair(d.sur, n=3), family="binomial")
    
    aic.sur <- modelcmp.glm(response = "survival", candidate = models.size, d.sur, family="binomial")
    aic.sur$response <- "survival"
    print("survival")
  }
  
  #**************************************************
  # 2 growth----
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
  # 3 flowering----
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
  # 4 fecundity----
  d.fec <- d %>%
    filter(focal.before == fc &
             background.species %in% sps &
             !is.na(size.autumn0) & !is.na(seeds.autumn0) &
             !(!is.na(note.size.autumn0) & note.size.autumn0 == "planting.autumn0" & year == 2020) &           
             (!is.na(if.same) &if.same==TRUE) )
  
  if(nrow(d.fec) == 0) aic.fec <- data.frame(y = "seeds.autumn0", n = 0, models=NA, aic = NA, response = "fecundity")
  else if(nrow(sample.pair(d.fec)) > 5) {
    # AIC comparison
    aic.fec <- modelcmp.glm(response = "seeds.autumn0", candidate = models.size, data=d.fec, family="gaussian")
    aic.fec$response <- "fecundity"
    print("fecundity")
  }
  
  
  #**************************************************
  # 5 germination----
  # Using FuNiche data
  
  #**************************************************
  # 6 seedling establishment----
  d.est <- d %>%
    filter(focal.species == fc & 
             (!is.na(planting.spring1) & planting.spring1=="yes") &
             !(!is.na(note.size.autumn0) & note.size.autumn0 == "planting.autumn0" & year == 2020) &   
             !is.na(survival.summer1.seedling) &
             background.species %in% sps)
  
  #if(nrow(sample.pair(d.est)) > 5 & length(unique(sample.pair(d.est)$site)) >1) {
  # AIC comparison
  #aic.est <- modelcmp.glm(response = "survival.summer1.seedling", candidate = models.1, data=d.est, family="gaussian")
  aic.est <- modelcmp.glm(response = "survival.summer1.seedling", candidate = models.1, data=d.est, family="binomial")
  aic.est$response <- "establishment"
  print("establishment")
  #}
  
  #**************************************************
  # 7 recruit size----
  # species level using greenhouse plants
  
  #**************************************************
  # 8 size range: L, U----
  # only compare the biggest individuals
  d.siz <- d %>% 
    filter(focal.species == fc &
             !is.na(size.autumn1) &
             !(!is.na(note.size.autumn0) & note.size.autumn0 == "planting.autumn0" & year == 2020) &   
             background.species %in% sps)
  
  if(nrow(sample.pair(d.siz)) > 5) {
    # AIC comparison
    aic.siz <- modelcmp.lm(response = "size.autumn1", candidate = models.1, data=d.siz)
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
# ******Germination using FuNiche data******------
#****************************************************************************
seeds <- read_excel("/Users/slyu/LVSM/NICH/Data/germination/FunNiche/FunNiche Seed Counts.xlsx")
seeds
colnames(seeds) <- c("rep","Plal", "Poal", "Anal", "Trba", "Seca", "Asal", "Armo", "Plla", "Brer", "Melu", "Trca", "Crbi", "Potr", "Sapr")

funiche <- read_excel("/Users/slyu/LVSM/NICH/Data/germination/FunNiche/FunNiche_data_2019_2020_dec2020.xlsx", na="NA")
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

# seeds per stick------------------------------------------
seeds %>% 
  pivot_longer(cols = 2:15, names_to = "species", values_to = "seeds.per.stick") %>%
  ggplot(aes(x=species, y=seeds.per.stick)) + geom_boxplot() + geom_point()

sn <- seeds %>% 
  pivot_longer(cols = 2:15, names_to = "species", values_to = "seeds.per.stick") %>%
  split(.$species) %>%
  purrr::map_dbl(~mean(.$seeds.per.stick, na.rm=TRUE))
sn

# germination and seedling establishment ---------------------------------------------
funiche$seeds_autumn_18 <- as.numeric(NA)
funiche$germination_spring_19 <- as.numeric(NA)
funiche$establishement_summer_19 <- as.numeric(NA)

# calculate
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
  ggplot(aes(y=germination_spring_19, x=factor(elevation), col=factor(elevation))) + geom_boxplot()+ geom_jitter(height=0.1) + facet_wrap(~species)

# germination against elevation for all sites
funiche %>% 
  ggplot(aes(y=germination_spring_19, x=elevation, col=elevation)) + geom_jitter(height=0.1) + facet_wrap(~species, scales="free")

# plot
funiche %>% 
  filter(site %in% c(4,13,19)) %>%
  ggplot(aes(y=establishement_summer_19, x=factor(elevation), col=factor(elevation))) + geom_boxplot()+ geom_jitter(height=0.1) + facet_wrap(~species)

funiche %>% 
  ggplot(aes(y=establishement_summer_19, x=elevation, col=elevation)) +  geom_jitter(height=0.1) + facet_wrap(~species, scales = "free_y") + stat_summary(col="red")

# plot
funiche %>% 
  filter(site %in% c(13,19)) %>%
  filter(species == "Seca") %>%
  .$germination_spring_19 %>%
  mean()

# species or site level? ---------------------------------------------
# germination
funiche$site <- as.character(funiche$site)
sps <- c("Anal","Armo","Asal","Plal", "Poal", "Seca", "Trba",        # 7 alpine species 
         "Brer", "Crbi", "Trca","Melu" , "Plla", "Potr", "Sapr") # 7 lowland species)
model.site <-  c("1", "site")

aic.ger <- aic.ger.best <- NULL
for(sp in sps) {
  #sp = "Anal"
  aic.germination <- funiche %>%
    filter(species == sp &
             site %in% c("4","13","19")) %>%
    modelcmp.lm(data=., response = "germination_spring_19", candidate = model.site)
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

# write.csv(aic.funiche, "/Users/slyu/LVSM/NICH/Results/R output_ok to delete/aic.funiche_20211117.csv")
# write.csv(aic.best.funiche, "/Users/slyu/LVSM/NICH/Results/R output_ok to delete/aic.funiche.best.csv")
# write.csv(funiche, "/Users/slyu/LVSM/NICH/Data/germination/FunNiche/FunNiche_data_2019_2020_dec2020_20210114.csv")

# recruit size---------------------------------------------
funiche$leaf_summer_19 <- as.numeric(NA)
fn <- subset(funiche, site %in% c("4","13","19") & num_flowerheads_summer_19 == 0)
fn$leaf <- as.numeric(fn$LLL_mean_summer_19) * as.numeric(fn$num_leaves_summer_19)/10
fn$site <- dplyr::recode(fn$site, "4" = "Les Posses", "13" = "Solalex", "19"="Anzeindaz")
fn$leaf
fn$sample <- "FuNiche"

# compare with greenhouse and transplanted plants
# greenhouse
biomass <- read_excel("/Users/slyu/LVSM/NICH/Data/Biomass/focal biomass/biomass data for correlation_2020.xlsx",na="NA", col_names = TRUE)
biomass$greenhouse <- ifelse(biomass$site == "greenhouse" | biomass$site == "chngreenhouse", "yes", "no")
gh <- subset(biomass, greenhouse == "yes")
gh$site = "greenhouse"
gh$sample <- "greenhouse"
gh$leaf <- gh$leaf.length.cm * gh$leaf.number

# transplanted
d <- read_excel("/Users/slyu/LVSM/NICH/Data/Biomass/data_181920_v5.xlsx", na="NA",col_names = TRUE)
d.rec <- d %>%
  filter(
    (!is.na(planting.spring1) & planting.spring1=="yes") &
      !(!is.na(note.size.autumn0) & note.size.autumn0 == "planting.autumn0" & year == 2020) &   
      !is.na(size.autumn1))
d.rec$leaf <- d.rec$leaf.number * d.rec$leaf.length.cm
d.rec$sample <- "transplant"
d.rec$species <- d.rec$focal.species

# function
mean_sd <- function (x, mult = 1) 
{
  x <- stats::na.omit(x)
  se <- mult * sd(x)
  mean <- mean(x)
  re <- data_frame(y = mean, ymin = mean - se, ymax = mean + se)
  as.data.frame(re)
}

rec <- bind_rows(fn, gh, d.rec)
rec %>%
  ggplot(aes(y=leaf, x=sample, col=site)) + geom_jitter() + stat_summary(fun.data=mean_sd, color="black") + facet_wrap(~species, scales = "free_y")

ggplot2::mean_ci

