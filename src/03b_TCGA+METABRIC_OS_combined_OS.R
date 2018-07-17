#######################################################################################################
# CoxPH analysis combining TCGA & METABRIC data sets
# Script authors: David Chen
# Script maintainer: David Chen
# Notes:
#######################################################################################################

rm(list=ls())
library(survival)
library(survminer)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path));
source("helper_functions.R");

## Data loading & exploration:
ADMIN_CENSOR <- 12*5; #in months
AGE_THRESH <- 50; 

sampsCombined <- loadCombinedTCGAandMETABRIC(receptorPosOnly=TRUE, HER2negOnly=TRUE, ADMIN_CENSOR=ADMIN_CENSOR);
checkAgeStrata(sampsCombined, thresh=AGE_THRESH, runFisher=FALSE)

mySubset <- subset(sampsCombined, Age <= AGE_THRESH); 
# mySubset <- subset(sampsCombined, Age > AGE_THRESH); 
# mySubset <- sampsCombined; #no age stratification

calcDescrpStatsByStrata(mySubset)

## CoxPH Regression:
res.cox <- coxph(
  Surv(time=OS_MONTHS, event=event) ~ group + Age + Stage, 
  data = mySubset
);

summary(res.cox)

## Kaplan-Meier:
kmFit <- survfit(
  Surv(time=OS_MONTHS, event=event) ~ SVM_BRCA1, 
  data = mySubset
);

ggsurvplot(
  kmFit,
  palette = c("mediumorchid","darkolivegreen3"),
  risk.table = FALSE,
  pval = TRUE
) + 
  scale_x_continuous(breaks=seq(0,60,by=10)) +
  scale_y_continuous(limits=c(0.6,1)) +
  ggtitle("ER+/PR+, HER2- tumors in TCGA+METABRIC, Age 50 or younger") +
  # ggtitle("ER+/PR+, HER2- tumors in TCGA+METABRIC, Age over 50") +
  # ggtitle("ER+/PR+, HER2- tumors in TCGA and METABRIC combined") +
  labs(x="Time (months)")


