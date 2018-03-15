##########################################################################################################
# Relation of HRD Index and LST Score with SVM BRCA1-like Probability Scores
# Script author: David Chen
# Notes:
##########################################################################################################

rm(list=ls())

## Clinical annotation for the study population:
my_samples <- read.csv("~/repos/BRCA1ness_by_SVM/annotations_and_backups/030418_TCGA_study_population.txt", header=T, sep="\t", stringsAsFactors=F);
table(my_samples$TNBC, useNA="ifany")

## Load published genomic alteration measures:
MarquardHRD <- read.table("~/Dropbox (Christensen Lab)/Pan-cancer-analyses/Marquard AM pan-cancer HRD index/40364_2015_33_MOESM8_ESM.txt", header=T);
MarquardHRD <- subset(MarquardHRD, Tumor %in% my_samples$patients);
colnames(MarquardHRD)[colnames(MarquardHRD)=="Tumor"] <- "patients"; #for merging

## Merge:
comp.mat <- merge(
  my_samples[ , c("patients","BRCA1_prob","SVM_BRCA1","Age","Stage","ER", "PR", "HER2","TNBC")], 
  MarquardHRD, 
  by = "patients"
);
table(comp.mat$SVM_BRCA1[comp.mat$TNBC=="TNBC"])
table(comp.mat$SVM_BRCA1[comp.mat$TNBC=="Non-TNBC"])
comp.mat$Subtype[comp.mat$TNBC=="TNBC"] <- "TNBC (64 BRCA1-like, 7 non-BRCA1-like)";
comp.mat$Subtype[comp.mat$TNBC=="Non-TNBC"] <- "Non-TNBC (143 BRCA1-like, 503 non-BRCA1-like)";

#-------------------------------------------Correlation analysis: HRD index-------------------------------------------
stopifnot( sum(is.na(comp.mat$HRD.LOH)) == 0 ) #no missing value in outcome
fitHRD <- lm(
  HRD.LOH ~ BRCA1_prob,
  data = comp.mat
);
summary(fitHRD)

#-------------------------------------------Correlation analysis: LST Score-------------------------------------------
stopifnot( sum(is.na(comp.mat$LST)) == 0 ) #no missing value in outcome
fitLST <- lm(
  LST ~ BRCA1_prob,
  data = comp.mat
);
summary(fitLST)
