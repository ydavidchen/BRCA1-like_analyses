#####################################################################################################
# Comparison of Ki-67 Microarray Gene Expression
# Script author: David Chen
# Notes:
#####################################################################################################

rm(list=ls())

library(ggplot2)
library(data.table)
library(RColorBrewer); gradient_cols <- brewer.pal(12, "Paired")
library(reshape2)

setwd("~/Dropbox (Christensen Lab)/Breast_Cancer_Data_Sets/METABRIC_data_set_cBio/")

## Load clinical data:
sample_clinical <- read.csv("~/repos/BRCA1ness_by_SVM/annotations_and_backups/030718_METABRIC_study_population.txt",sep="\t",stringsAsFactors=F);
sample_clinical <- subset(sample_clinical, TNBC=="Non-TNBC");

## Gene expression data loading:
arrayExpr <- fread("data_expression.txt");
class(arrayExpr) <- "data.frame";
anyDuplicated(arrayExpr$Hugo_Symbol)
rownames(arrayExpr) <- arrayExpr$Hugo_Symbol;
arrayExpr$Hugo_Symbol <- arrayExpr$Entrez_Gene_Id <- NULL;
arrayExpr <- as.matrix(arrayExpr);
stopifnot( mode(arrayExpr) == "numeric" )

## Subset to Ki67:
plt.Ki67 <- arrayExpr[rownames(arrayExpr) == "MKI67", , drop=F];
plt.Ki67 <- as.data.frame( t(plt.Ki67) );
plt.Ki67$SAMPLE_ID <- rownames(plt.Ki67);

## Merge in covariates to adjust:
plt.Ki67 <- merge(
  plt.Ki67,
  sample_clinical[ , c("SAMPLE_ID","SVM_BRCA1","AGE_AT_DIAGNOSIS","Stage","ER_STATUS","PR_STATUS","HER2_STATUS","TNBC")],
  by = "SAMPLE_ID"
);

## Statistical test:
plt.Ki67$group[plt.Ki67$SVM_BRCA1=="BRCA1-like"] <- 1;
plt.Ki67$group[plt.Ki67$SVM_BRCA1=="non-BRCA1-like"] <- 0;
plt.Ki67$group <- as.factor(plt.Ki67$group);
table(plt.Ki67$group)

fit_other <- lm(
  MKI67 ~ group + AGE_AT_DIAGNOSIS + Stage + ER_STATUS + PR_STATUS + HER2_STATUS,
  data = plt.Ki67
);
summary(fit_other)
