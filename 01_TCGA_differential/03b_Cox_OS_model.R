################################################################################################
# Cox Models by Binary SVM Status, Stratified by TNBC Status, in TCGA
# Script author: David Chen
# Date: 02/28/2018; 03/04/2018
# Notes:
################################################################################################

rm(list = ls())
library(survival)
library(survminer)
source("~/repos/Repos_for_Manuscript_Code/BRCA1-like_analyses/helper_functions.R");
my_samples <- loadReceptorPositiveTumors(receptorPosOnly=TRUE); 

## Merge in survival time (as data-availability check):
clin.breast <- read.csv("~/repos/BRCA1ness_by_SVM/annotations_and_backups/TCGA-BRCA_clinical.csv", stringsAsFactors=F);
colnames(clin.breast)[1] <- "patients";
my_samples <- merge(
  my_samples[ , c("patients","SVM_BRCA1","Age","Stage","ER","PR","HER2","TNBC","PAM50","PAM50lite")],
  clin.breast[ , c("patients","days_to_last_follow_up","days_to_death","vital_status")]
);

## For patients who died, fill in days to last follow-up:
condition <- is.na(my_samples$days_to_last_follow_up) & my_samples$vital_status=="dead";
sum(condition)
my_samples$days_to_last_follow_up[condition] <- my_samples$days_to_death[condition];

## Define main variable of interest:
my_samples$group[my_samples$SVM_BRCA1=="BRCA1-like"] <- 1;
my_samples$group[my_samples$SVM_BRCA1=="non-BRCA1-like"] <- 0;
my_samples$group <- as.factor(my_samples$group);
table(my_samples$group)
table(my_samples$SVM_BRCA1)

## Code in Subtype based on O'Sullivan & Johnson et al.:
my_samples$Subtype[my_samples$TNBC=="TNBC"] <- "TNBC";
my_samples$Subtype[my_samples$HER2=="Positive"] <- "HER2+";
my_samples$Subtype[my_samples$ER=="Positive" & my_samples$HER2=="Negative"] <- "ER+";
my_samples$Subtype[is.na(my_samples$Subtype)] <- "Other"; #will be dropped
table(my_samples$Subtype)

## Define main outcome of interest:
my_samples$event <- my_samples$vital_status == "dead";

## Important: Convert days_to_follow-up to months:
my_samples$OS_MONTHS <- my_samples$days_to_last_follow_up / 30;
my_samples$days_to_last_follow_up <- my_samples$days_to_death <- NULL;

## (Adjust as necessary) Impose administrative censoring:
adminCensor <- 12*5; #in months
my_samples$event[my_samples$OS_MONTHS >= adminCensor] <- FALSE; #not dead by definition
my_samples$OS_MONTHS[my_samples$OS_MONTHS >= adminCensor] <- adminCensor;

## Exclude tumors without stage:
my_samples <- subset(my_samples, ! is.na(Stage) );

## CoxPH:
res.cox.other <- coxph(
  Surv(time=OS_MONTHS, event=event) ~ group + Age + Stage + Subtype, 
  data = my_samples
);
summary(res.cox.other) 

## Kaplan-Meier:
fit.other <- survfit(
  Surv(time=OS_MONTHS, event=event) ~ SVM_BRCA1, 
  my_samples
);
ggsurvplot(
  fit.other,
  palette = c("mediumorchid","darkolivegreen3"),
  risk.table = TRUE,
  pval = TRUE
) + 
  ggtitle("TCGA receptor positive tumors") +
  labs(x="Time (months)")

