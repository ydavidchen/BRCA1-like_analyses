#######################################################################################################
# Explore TCGA breast tumor characteristics by SVM-BRCA1-like status
# Script author: Y. David Chen
# Script maintainer: Y. David Chen
# Notes:
#######################################################################################################

rm(list=ls())
library(ggplot2)
library(matrixStats)
library(reshape2)
library(tableone)
library(WriteXLS)
source("~/repos/Repos_for_Manuscript_Code/BRCA1-like_analyses/helper_functions.R"); 
source("~/repos/Repos_for_Manuscript_Code/BRCA1-like_analyses/plot_themes.R");

## Clinical annotation for study population:
my_samples <- read.csv("~/repos/BRCA1ness_by_SVM/annotations_and_backups/030418_TCGA_study_population.txt", header=T, sep="\t", stringsAsFactors=F);
mySampNoExclu <- read.csv("~/repos/BRCA1ness_by_SVM/annotations_and_backups/031318_Complete_Sample_Sheet_No_Exclusion.txt", header=T, sep="\t", stringsAsFactors=F);
mySampNoExclu$includedInBRCA1Analysis <- mySampNoExclu$patients %in% my_samples$patients;
mySampNoExclu$anyRecepPos <- mySampNoExclu$ER=="Positive" |  mySampNoExclu$PR=="Positive" | mySampNoExclu$HER2=="Positive";

#---------------------------------------Distribution of BRCA1-like probability---------------------------------------
mySampNoExclu$recepPos[mySampNoExclu$anyRecepPos] <- "Receptor Positive";
mySampNoExclu$recepPos[! mySampNoExclu$anyRecepPos] <- "TNBC";
mySampNoExclu$recepPos[is.na(mySampNoExclu$anyRecepPos)] <- "Unknown"; 

png("~/Downloads/BRCA1ness_figures/Figure1C.png", res=300, units="in", height=8.27, width=5.84);
ggplot(mySampNoExclu, aes(x=reorder(patients, BRCA1_prob), y=BRCA1_prob, fill=recepPos)) +
  geom_bar(stat="identity") +
  geom_hline(yintercept=0.50, linetype="dashed") + 
  labs(x="Sample", y="SVM BRCA1-like probability", title="TCGA") +
  scale_x_discrete(expand=c(0, 0)) +
  scale_y_continuous(breaks=seq(0,1,0.2), expand=c(0,0), limits=c(0,1.05)) +
  scale_fill_manual(values=c("black", "blue", "lightgray")) +
  theme_classic() +
  theme(axis.text.x=element_blank(), axis.text.y=element_text(size=20,color="black"),
        axis.ticks.x=element_blank(),
        axis.title=element_text(size=25, color="black"), title=element_text(size=25, color="black"),
        legend.position="top",legend.title=element_blank(),legend.text=element_text(size=15,color="black") )
dev.off();

#---------------------------------------------Table One: Overall---------------------------------------------
## Table One
myFactorVars <- c("Stage","ER","PR","HER2","TNBC","PAM50", "race", "anyRecepPos", "includedInBRCA1Analysis");
# myFactorVars <- c(myFactorVars, "Clinical_Subtype");
myVars <- c("Age", myFactorVars);
myTableOne <- CreateTableOne(
  data = mySampNoExclu, 
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

## Visualize & export:
y <- print(myTableOne, showAllLevels=TRUE);
# write.csv(y, file="~/Downloads/051118_TCGA-BRCA_No-exclusion_TableOne.csv",quote=F)
# WriteXLS(mySampNoExclu, ExcelFileName="~/Downloads/031318_Sample_Sheet_TCGA_NoExclu.xls");

#---------------------------------------------Table One: Subset---------------------------------------------
## Table One
myFactorVars <- c("Stage","ER","PR","HER2");
myVars <- c("Age", myFactorVars);
myTableOne <- CreateTableOne(
  data = subset(mySampNoExclu, TNBC=="Non-TNBC" & includedInBRCA1Analysis),
  vars = myVars, 
  strata = "SVM_BRCA1",
  factorVars = myFactorVars,
  smd = FALSE, 
  test = FALSE,
  includeNA = TRUE
);

## Visualize & export:
y <- print(myTableOne, showAllLevels=TRUE); 
write.csv(y, file="~/Downloads/BRCA1ness_figures/041118_TCGA_Subset_Table_one.csv",quote=F)
