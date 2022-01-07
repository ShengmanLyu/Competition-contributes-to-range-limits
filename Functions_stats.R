#########################################################################
# Functions for IPMs and coexistence analyses
# Author: Shengman Lyu (shengman.lyu@usys.ethz.ch)
# Date: 07.01.2022
#########################################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# To extrac stats from lmer----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(car)
coef.lmer <- function(model, ...) {
  fit <- Anova(model, ...)
  return(as.data.frame(fit))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# To calculate confidence interval using linear model----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# same as: Hmisc::mean_cl_normal
mean_ci <- function(x) {
  fit <- lm(x~1)
  c(mean(x), confint(fit)[1, ])
  }

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# To calculate confidence interval using SE----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mean_ci_se <- function(x) {
  u = mean(x)
  se = sqrt(var(x)/length(x))
  data.frame(y=u, ymin=u-1.96*se, ymax=u+1.96*se)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# To calculate confidence interval using SD----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mean_ci_sd <- function(x, ...) {
  u = mean(x, ...)
  sd = sqrt(var(x, ...))
  data.frame(y=u, ymin=u-1.96*sd, ymax=u+1.96*sd)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# To calcualte confidence interval using quantitle----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mean_ci_hanna <- function(x, ...) {
  u = median(x, ...)
  x.min = quantile(x, probs = 0.025, ...)
  x.max = quantile(x, probs = 0.975, ...)
  data.frame(y=u, ymin=x.min, ymax=x.max) 
}
