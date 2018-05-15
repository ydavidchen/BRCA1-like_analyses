#######################################################################################################
# Explore TCGA breast cancer data characteristics
# Script authors: David Chen
# Notes:
#######################################################################################################

rm(list=ls())

library(ggplot2)
library(matrixStats)
library(reshape2)
library(tableone)
library(WriteXLS)
library(doParallel); registerDoParallel(detectCores() - 1)

## Clinical annotation for study population:
my_samples <- read.csv("~/repos/BRCA1ness_by_SVM/annotations_and_backups/030418_TCGA_study_population.txt", header=T, sep="\t", stringsAsFactors=F);
mySampNoExclu <- read.csv("~/repos/BRCA1ness_by_SVM/annotations_and_backups/031318_Complete_Sample_Sheet_No_Exclusion.txt", header=T, sep="\t", stringsAsFactors=F);
mySampNoExclu$includedInBRCA1Analysis <- mySampNoExclu$patients %in% my_samples$patients;

## Any hormone receptor positive:
mySampNoExclu$anyRecepPos <- mySampNoExclu$ER=="Positive" |  mySampNoExclu$PR=="Positive" | mySampNoExclu$HER2=="Positive";

#---------------------------------------Distribution of BRCA1-like probability scores---------------------------------------
ggplot(mySampNoExclu, aes(x=reorder(patients, BRCA1_prob), y=BRCA1_prob, fill=TNBC)) + #fill=DNAm_avail
  geom_bar(stat="identity") +
  geom_hline(yintercept=0.50, linetype="dashed") + 
  labs(x="Sample", y="SVM BRCA1-like probability", title="TCGA breast tumors") +
  scale_x_discrete(expand=c(0, 0)) +
  scale_y_continuous(breaks=seq(0,1,0.2), expand=c(0,0), limits=c(0,1.05)) +
  theme_classic() +
  theme(axis.text.x=element_blank(), axis.text.y=element_text(size=15,color="black"),
        axis.ticks.x=element_blank(),
        axis.title=element_text(size=15, color="black"),
        legend.position="top",legend.title=element_blank(),legend.text=element_text(size=12,color="black") )

#---------------------------------------------Table One (part of)---------------------------------------------
myFactorVars <- c("Stage","ER","PR","HER2","TNBC","PAM50", "race", "anyRecepPos", "includedInBRCA1Analysis");
myVars <- c("Age", myFactorVars);
myTableOne <- CreateTableOne(
  data = mySampNoExclu, 
  vars = myVars, 
  strata = "SVM_BRCA1",
  factorVars = myFactorVars,
  smd = TRUE,
  test = TRUE,
  testExact = fisher.test,
  argsExact = list(workspace = 2*10^5),
  testNormal = oneway.test,
  argsNormal = list(var.equal=TRUE),
  includeNA = TRUE
);
myTableOne
