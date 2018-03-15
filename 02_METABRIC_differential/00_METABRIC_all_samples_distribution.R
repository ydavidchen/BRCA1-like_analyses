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
table(sampNoExclu$TNBC, useNA="ifany")

#---------------------------------------Distribution of BRCA1-like probability---------------------------------------
ggplot(sampNoExclu, aes(x=reorder(SAMPLE_ID, BRCA1_Prob), y=BRCA1_Prob, fill=TNBC)) + #fill=DNAm_avail
  geom_bar(stat="identity") +
  geom_hline(yintercept=0.50, linetype="dashed") + 
  labs(x="Sample", y="SVM BRCA1-like probability", title="Curtis et al. / METABRIC data set") +
  scale_x_discrete(expand=c(0, 0)) +
  scale_y_continuous(breaks=seq(0,1,0.2), expand=c(0,0), limits=c(0,1.05)) +
  theme_classic() +
  theme(axis.text.x=element_blank(), axis.text.y=element_text(size=15,color="black"),
        axis.ticks.x=element_blank(),
        axis.title=element_text(size=15, color="black"),
        legend.position="top",legend.title=element_blank(),legend.text=element_text(size=12,color="black") )

#---------------------------------------------Table One---------------------------------------------
## Assign treatment categories:
sampNoExclu$Treatment[sampNoExclu$CHEMOTHERAPY=="YES" & sampNoExclu$HORMONE_THERAPY=="YES"] <- "Both Hormone & Chemo";
sampNoExclu$Treatment[sampNoExclu$CHEMOTHERAPY=="YES" & sampNoExclu$HORMONE_THERAPY=="NO"] <- "Chemotherapy";
sampNoExclu$Treatment[sampNoExclu$CHEMOTHERAPY=="NO" & sampNoExclu$HORMONE_THERAPY=="YES"] <- "Hormone";
sampNoExclu$Treatment[sampNoExclu$CHEMOTHERAPY=="NO" & sampNoExclu$HORMONE_THERAPY=="NO"] <- "Neither";

table(sampNoExclu$CHEMOTHERAPY, useNA="ifany")
table(sampNoExclu$HORMONE_THERAPY, useNA="ifany")
table(sampNoExclu$Treatment, useNA="ifany")

## Assign clinical subtype (order of assignment matters):
sampNoExclu$Clinical_Subtype[sampNoExclu$ER_STATUS=="+" & sampNoExclu$HER2_STATUS=="-"] <- "ER Positive";
sampNoExclu$Clinical_Subtype[sampNoExclu$HER2_STATUS=="+"] <- "HER2 Positive";
sampNoExclu$Clinical_Subtype[sampNoExclu$TNBC=="TNBC"] <- "TNBC";
sampNoExclu$Clinical_Subtype[is.na(sampNoExclu$Clinical_Subtype)] <- "Other";
sampNoExclu$Clinical_Subtype[is.na(sampNoExclu$TNBC)] <- NA;
table(sampNoExclu$Clinical_Subtype)

## Table One: 
myFactorVars <- c("Stage","Treatment","ER_STATUS","PR_STATUS","HER2_STATUS","Clinical_Subtype","CLAUDIN_SUBTYPE","isIncludedInBRCA1analysis");
myVars <- c("AGE_AT_DIAGNOSIS", myFactorVars);
myTableOne <- CreateTableOne(
  data = sampNoExclu, 
  vars = myVars, 
  strata = "SVM_BRCA1",
  factorVars = myFactorVars,
  smd = FALSE, 
  test = FALSE,
  includeNA = TRUE
); 
y <- print(myTableOne, showAllLevels=TRUE); 

