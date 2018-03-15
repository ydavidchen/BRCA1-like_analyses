##########################################################################################################
# Relation of Mutational Rates with SVM BRCA1-like Status
# Script author: David Chen
# Notes:
##########################################################################################################

rm(list=ls())

library(gdata)
library(lmtest)
library(MASS)
library(matrixStats)
library(sandwich)

## Clinical annotation for the study population:
my_samples <- read.csv("~/repos/BRCA1ness_by_SVM/annotations_and_backups/030418_TCGA_study_population.txt", header=T, sep="\t", stringsAsFactors=F);

## Load mutational count & rates data (Kandoth et al. 2013):
mutMeasure <- read.xls("~/Dropbox (Christensen Lab)/Pan-cancer-analyses/Kandoth2013Nature/Supplementary_Table_3a.xls", stringsAsFactors=F);
mutMeasure <- subset(mutMeasure, Cancer.Type == "brca");
mutMeasure$patients <- substr(mutMeasure$TCGA.ID, 1, 12);
stopifnot( anyDuplicated(mutMeasure$patients) == 0 )

sum(my_samples$patients %in% mutMeasure$patients)

## Join with sample annotation:
comp.mat <- merge(
  my_samples[ , c("patients","SVM_BRCA1","Age","Stage","ER", "PR", "HER2","TNBC")], 
  mutMeasure, 
  by = "patients"
);
table(comp.mat$SVM_BRCA1[comp.mat$TNBC=="TNBC"])
table(comp.mat$SVM_BRCA1[comp.mat$TNBC=="Non-TNBC"])
table(comp.mat$SVM_BRCA1)

## Statistical test:
comp.mat$group[comp.mat$SVM_BRCA1=="BRCA1-like"] <- 1;
comp.mat$group[comp.mat$SVM_BRCA1=="non-BRCA1-like"] <- 0;
comp.mat$group <- as.factor(comp.mat$group);

fit.rates <- rlm(
  `Mutation.Rate...Mbp.` ~ group + Age + Stage + ER + PR + HER2,
  data = subset(comp.mat, TNBC=="Non-TNBC")
);
summary(fit.rates)
svar.rlm <- sandwich(fit.rates);
coeftest(fit.rates, vcov.=svar.rlm);
