######################################################################################################
# Comparison of Horvath DNA methylation age
# Script author: David Chen
# Notes:
######################################################################################################

rm(list=ls())
library(gdata)
library(wateRmelon)

## Load data & results:
load("~/Dropbox (Christensen Lab)/Christensen Lab - 2017/BRCA1ness_TCGA_all_types/betas.TCGAC450k.RData");
load("~/Dropbox (Christensen Lab)/Christensen Lab - 2017/BRCA1ness_TCGA_all_types/030418_DMRcate_Non-TNBCs.RData");

## Check if all Non-TNBCs:
table(my_samples$TNBC)

## Mutual subsetting
betas_TCGA450k <- betas_TCGA450k[ , colnames(betas_TCGA450k) %in% my_samples$patients];
dim(betas_TCGA450k)

## Execute methylation age calculation:
DNAmAge <- agep(betas_TCGA450k);
DNAmAge <- as.data.frame(DNAmAge);
colnames(DNAmAge) <- "Horvath";
DNAmAge$patients <- rownames(DNAmAge);

## Merge into sample annotation:
my_samples <- merge(my_samples, DNAmAge, by="patients"); 

## Linear models controlling for covariates:
my_samples$group[my_samples$SVM_BRCA1=="BRCA1-like"] <- 1;
my_samples$group[my_samples$SVM_BRCA1=="non-BRCA1-like"] <- 0;
my_samples$group <- as.factor(my_samples$group);
table(my_samples$group)
table(my_samples$SVM_BRCA1)

horvath.comp <- lm(Horvath ~ group + Stage + ER + PR + HER2, data=my_samples);
summary(horvath.comp)

chrono.comp <- lm(Age ~ group + Stage + ER + PR + HER2, data=my_samples);
summary(chrono.comp)
