####################################################################################################
# Correlation of WTSI mutational signature 3 with BRCA1-like status in TCGA
# Script author: David Chen
# Notes:
####################################################################################################

rm(list=ls())

## Clinical annotation for study population:
my_samples <- read.csv("~/repos/BRCA1ness_by_SVM/annotations_and_backups/030418_TCGA_study_population.txt", header=T, sep="\t", stringsAsFactors=F);
my_samples <- subset(my_samples, TNBC=="Non-TNBC");

## Load mutational signatures:
rosenthalTCGA <- read.table("~/Dropbox (Christensen Lab)/Pan-cancer-analyses/TCGA_pan_cancer_signature_Rosenthal.txt", header=T);
rosenthalTCGA$id <- gsub(".", "-", rosenthalTCGA$id, fixed=TRUE);
rosenthalTCGA <- subset(rosenthalTCGA, cancer == "BRCA");
colnames(rosenthalTCGA)[1] <- "patients"; #for merging

## Merge:
comp.mat <- merge(my_samples, rosenthalTCGA, by="patients");
comp.mat <- subset(comp.mat, method == "WTSI");

## Statistical test:
comp.mat$group[comp.mat$SVM_BRCA1=="BRCA1-like"] <- 1;
comp.mat$group[comp.mat$SVM_BRCA1=="non-BRCA1-like"] <- 0;
comp.mat$group <- as.factor(comp.mat$group);

fitSig3 <- lm(
  Signature.3 ~ group + Age + Stage + ER + PR + HER2, 
  data = comp.mat
);
summary(fitSig3)
