################################################################################################
# Survival Modeling by SVM BRCA1-like Binary Status in ER+/PR+ & HER2- Tumors
# Script author: David Chen
# Script maintainer: David Chen
# Notes:
################################################################################################

rm(list=ls())
library(survival)
library(survminer)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path));
source("../helper_functions.R");
my_samples <- loadTCGA4CoxOS(ADMIN_CENSOR=60, receptorPosOnly=TRUE, HER2negOnly=TRUE);
  
## CoxPH:
res.cox <- coxph(
  Surv(time=OS_MONTHS, event=event) ~ group + Age + Stage, 
  data = my_samples
);

summary(res.cox)

## Kaplan-Meier:
fit <- survfit(
  Surv(time=OS_MONTHS, event=event) ~ SVM_BRCA1, 
  my_samples
);

ggsurvplot(
  fit,
  palette = c("mediumorchid","darkolivegreen3"),
  risk.table = FALSE,
  pval = TRUE
) + 
  scale_x_continuous(breaks=seq(0,60,by=10)) +
  scale_y_continuous(limits=c(0.6,1)) +
  ggtitle("TCGA ER+/PR+, HER2- tumors") +
  labs(x="Time (months)")

