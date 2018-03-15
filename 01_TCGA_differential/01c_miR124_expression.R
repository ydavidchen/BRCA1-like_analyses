######################################################################################################
# Verifying expression of miR with differential methylation
# Script author: David Chen
# Notes:
# 1. miR124-2 / miR-124a2: HYPERmethylated in BRCA1-like; modulates proliferation (PMID:25731732)
######################################################################################################

rm(list=ls())

library(data.table)
library(matrixStats)

#-----------------------------------------------Load data for the study population-----------------------------------------------
## Load clinical data & define study population:
my_samples <- read.csv("~/repos/BRCA1ness_by_SVM/annotations_and_backups/030418_TCGA_study_population.txt", header=T, sep="\t", stringsAsFactors=F);
my_samples <- subset(my_samples, TNBC=="Non-TNBC");

## miRNA data:
miRExpr <- fread("~/Dropbox (Christensen Lab)/Breast_Cancer_Data_Sets/All_TCGA_breast_RNAseq/bcgsc.ca_BRCA_IlluminaGA_miRNASeq.miRNAExp.tsv", stringsAsFactors=F)
class(miRExpr) <- "data.frame";
rownames(miRExpr) <- miRExpr$V1; 
miRExpr$V1 <- NULL;
miRExpr <- as.matrix(miRExpr);
stopifnot( mode(miRExpr) == "numeric" )

## Select primary tumors: 
table(substr(colnames(miRExpr), 14, 16))
miRExpr <- miRExpr[ , substr(colnames(miRExpr), 14, 15)=="01"];
colnames(miRExpr) <- substr(colnames(miRExpr), 1, 12);

## Convert to Z-score:
miRExpr <- scale(miRExpr, scale=TRUE, center=TRUE);

## Subset:
my_miRs <- subset(miRExpr, rownames(miRExpr) %in% c("hsa-mir-124-2"));

#-----------------------------------------------Statistical tests-----------------------------------------------
my_miRs <- as.data.frame(t(my_miRs));
my_miRs$patients <- rownames(my_miRs);
my_miRs <- merge(
  my_miRs, 
  my_samples[ , c("patients","SVM_BRCA1","TNBC","Age","Stage","ER","PR","HER2")], 
  by="patients"
);

my_miRs$group[my_miRs$SVM_BRCA1 == "BRCA1-like"] <- 1;
my_miRs$group[my_miRs$SVM_BRCA1 == "non-BRCA1-like"] <- 0;
table(my_miRs$group)
table(my_miRs$SVM_BRCA1)

fit_miR124a2_other <- lm(
  `hsa-mir-124-2` ~ group + Age + Stage + ER + PR + HER2,
  data = my_miRs
);
summary(fit_miR124a2_other)
