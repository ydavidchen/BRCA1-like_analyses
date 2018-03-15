#######################################################################################################
# CoxPH Modeling, METABRIC data set
# Script authors: David Chen
# Notes:
#######################################################################################################

rm(list = ls())

library(survival)
library(survminer)

## Load Curtis et al. clinical data
## Survival data already included:
sample_clinical <- read.csv("~/repos/BRCA1ness_by_SVM/annotations_and_backups/030718_METABRIC_study_population.txt",sep="\t",stringsAsFactors=F);

## Define main variable of interest:
sample_clinical$group[sample_clinical$SVM_BRCA1=="BRCA1-like"] <- 1;
sample_clinical$group[sample_clinical$SVM_BRCA1=="non-BRCA1-like"] <- 0;
sample_clinical$group <- as.factor(sample_clinical$group);

## Define clinical subtypes based on O'Sullivan & Johnson et al.:
sample_clinical$Subtype[sample_clinical$ER_STATUS=="+" & sample_clinical$HER2_STATUS=="-"] <- "ER+";
sample_clinical$Subtype[sample_clinical$HER2_STATUS=="+"] <- "HER2+";
sample_clinical$Subtype[sample_clinical$ER_STATUS=="-" & sample_clinical$PR_STATUS=="-" & sample_clinical$HER2_STATUS=="-"] <- "TNBC";
sample_clinical$Subtype[is.na(sample_clinical$Subtype)] <- "Other";
table(sample_clinical$Subtype)

## Define main outcome of interest:
sample_clinical$event <- sample_clinical$OS_STATUS == "DECEASED";

## (Adjust based on hypothesis) Impose administrative censoring:
adminCensor <- 12*5; #in months
sample_clinical$event[sample_clinical$OS_MONTHS >= adminCensor] <- FALSE; #not dead by definition
sample_clinical$OS_MONTHS[sample_clinical$OS_MONTHS >= adminCensor] <- adminCensor;

## CoxPH:
res.cox.other <- coxph(
  Surv(time=OS_MONTHS, event=event) ~ group + AGE_AT_DIAGNOSIS + Stage + Subtype,
  data = subset(sample_clinical, TNBC=="Non-TNBC")
);
summary(res.cox.other)

## OS curves:
fit.other <- survfit(
  Surv(time=OS_MONTHS, event=event) ~ SVM_BRCA1, 
  subset(sample_clinical, TNBC=="Non-TNBC")
);
ggsurvplot(
  fit.other, 
  data = subset(sample_clinical, TNBC=="Non-TNBC"), 
  palette = c("mediumorchid","darkolivegreen3"), 
  risk.table = TRUE, 
  pval = "CoxPH P = 0.21"
) +
  ggtitle("METABRIC receptor positive tumors") +
  labs(x="Time (months)")
