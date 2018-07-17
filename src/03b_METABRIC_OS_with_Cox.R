#######################################################################################################
# CoxPH Modeling in METABRIC data set
# Script authors: David Chen
# Script maintainer: David Chen
# Notes:
#######################################################################################################

rm(list=ls())
library(survival)
library(survminer)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path));
source("helper_functions.R");
sample_clinical <- loadMETABRIC4CoxOS(ADMIN_CENSOR=60, receptorPosOnly=TRUE, HER2negOnly=TRUE);

## Covariate-adjusted CoxPH:
res.cox <- coxph(
  Surv(time=OS_MONTHS, event=event) ~ group + AGE_AT_DIAGNOSIS + Stage,
  data = sample_clinical
);

summary(res.cox)

## OS curve:
fit <- survfit(
  Surv(time=OS_MONTHS, event=event) ~ SVM_BRCA1, 
  sample_clinical
);

ggsurvplot(
  fit,
  palette = c("mediumorchid","darkolivegreen3"), 
  risk.table = FALSE, 
  pval = TRUE
) +
  scale_x_continuous(breaks=seq(0,60,by=10)) +
  scale_y_continuous(limits=c(0.6,1)) +
  ggtitle("METABRIC ER+/PR+ & HER2- tumors") +
  labs(x="Time (months)")
