#######################################################################################################
# CoxPH combining TCGA & METABRIC
# Script authors: David Chen
# Script maintainer: David Chen
# Notes:
#######################################################################################################

rm(list=ls())
library(survival)
library(survminer)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path));
source("../helper_functions.R");
ADMIN_CENSOR <- 12*5; #in months

## Combine data sets:
sampsTCGA <- loadTCGA4CoxOS(ADMIN_CENSOR=ADMIN_CENSOR, receptorPosOnly=TRUE, HER2negOnly=TRUE);
sampsTCGA <- sampsTCGA[ , c("patients","SVM_BRCA1","group","Age","Stage","ER","PR","event","OS_MONTHS")]; 
colnames(sampsTCGA) 

sampsMETABRIC <- loadMETABRIC4CoxOS(ADMIN_CENSOR=ADMIN_CENSOR, receptorPosOnly=TRUE, HER2negOnly=TRUE);
sampsMETABRIC <- sampsMETABRIC[ , c("SAMPLE_ID","SVM_BRCA1","group","AGE_AT_DIAGNOSIS","Stage","ER_STATUS","PR_STATUS","event","OS_MONTHS")];
sampsMETABRIC[sampsMETABRIC=="+"] <- "Positive";
sampsMETABRIC[sampsMETABRIC=="-"] <- "Negative";
colnames(sampsMETABRIC) <- colnames(sampsTCGA);

sampsTCGA$dataset <- "TCGA";
sampsMETABRIC$dataset <- "METABRIC";
sampsCombined <- rbind(sampsTCGA, sampsMETABRIC);

table(sampsCombined$SVM_BRCA1)
table(sampsCombined$Stage)
table(sampsCombined$ER)
table(sampsCombined$PR)
table(sampsCombined$event)
table(sampsCombined$dataset)
hist(sampsCombined$OS_MONTHS)
hist(sampsCombined$Age)

## Cox Model:
res.cox <- coxph(
  Surv(time=OS_MONTHS, event=event) ~ group + Age + Stage, 
  data = sampsCombined
);

summary(res.cox)

## Kaplan-Meier Curve:
kmFit <- survfit(
  Surv(time=OS_MONTHS, event=event) ~ SVM_BRCA1, 
  sampsCombined
);

ggsurvplot(
  kmFit,
  palette = c("mediumorchid","darkolivegreen3"),
  risk.table = FALSE,
  pval = TRUE
) + 
  scale_x_continuous(breaks=seq(0,60,by=10)) +
  scale_y_continuous(limits=c(0.6,1)) +
  ggtitle("ER+/PR+, HER2- tumors in TCGA and METABRIC combined") +
  labs(x="Time (months)")


