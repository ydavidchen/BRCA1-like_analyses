#######################################################################################################
# Explore clinical & survival data of METABRIC for Table 1/S1
# Script authors: David Chen
# Date: 03/13/2018
# Notes:
#######################################################################################################

rm(list=ls())

library(ggplot2)
library(matrixStats)
library(reshape2)
library(tableone)
library(WriteXLS)
library(doParallel); registerDoParallel(detectCores() - 1)

## Load clinical data:
sample_clinical <- read.csv("~/repos/BRCA1ness_by_SVM/annotations_and_backups/030718_METABRIC_study_population.txt",sep="\t",stringsAsFactors=F);
sampNoExclu <- read.csv("~/repos/BRCA1ness_by_SVM/annotations_and_backups/031318_METABRIC_NoExclusion.txt",sep="\t",stringsAsFactors=F);
sampNoExclu$isIncludedInBRCA1analysis <- sampNoExclu$SAMPLE_ID %in% sample_clinical$SAMPLE_ID;
sampNoExclu$anyRecepPos <- sampNoExclu$ER_STATUS=="+" | sampNoExclu$PR_STATUS=="+" | sampNoExclu$HER2_STATUS=="+";
table(sampNoExclu$TNBC, useNA="ifany")

#---------------------------------------Distribution of BRCA1-like probability---------------------------------------
png("~/Downloads/BRCA1ness_figures/Figure1D.png", res=300, units="in", height=8.27, width=5.84);
ggplot(sampNoExclu, aes(x=reorder(SAMPLE_ID, BRCA1_Prob), y=BRCA1_Prob, fill=TNBC)) + #fill=DNAm_avail
  geom_bar(stat="identity") +
  geom_hline(yintercept=0.50, linetype="dashed") + 
  labs(x="Sample", y="SVM BRCA1-like probability", title="METABRIC") +
  scale_x_discrete(expand=c(0, 0)) +
  scale_y_continuous(breaks=seq(0,1,0.2), expand=c(0,0), limits=c(0,1.05)) +
  theme_classic() +
  theme(axis.text.x=element_blank(), axis.text.y=element_text(size=20,color="black"),
        axis.ticks.x=element_blank(),
        axis.title=element_text(size=25, color="black"), title=element_text(size=25, color="black"),
        legend.position="top",legend.title=element_blank(),legend.text=element_text(size=15,color="black") )
dev.off();

#---------------------------------------------Table One---------------------------------------------
## Table One on overall: 
myFactorVars <- c("Stage","ER_STATUS","PR_STATUS","HER2_STATUS","anyRecepPos","isIncludedInBRCA1analysis");
myVars <- c("AGE_AT_DIAGNOSIS", myFactorVars);
myTableOne <- CreateTableOne(
  data = sampNoExclu, 
  vars = myVars, 
  strata = "SVM_BRCA1",
  factorVars = myFactorVars,
  smd = TRUE,
  test = TRUE,
  testExact = fisher.test,
  argsExact = list(workspace = 2*10^5), #default
  testNormal = oneway.test,
  argsNormal = list(var.equal=TRUE),
  includeNA = TRUE
);
myTableOne

