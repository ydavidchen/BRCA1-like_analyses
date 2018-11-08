# TCGA breast tumor characteristics by SVM-BRCA1-like status
# Script author: David Chen
# Script maintainer: David Chen
# Notes:

rm(list=ls());
library(ggplot2);
library(matrixStats);
library(reshape2);
library(tableone);
library(WriteXLS);
setwd(dirname(rstudioapi::getActiveDocumentContext()$path));
source("helper_functions.R");
source("plot_themes.R");

## Clinical annotation for study population:
my_samples <- read.csv("~/repos/BRCA1ness_by_SVM/annotations_and_backups/030418_TCGA_study_population.txt", header=T, sep="\t", stringsAsFactors=F);
mySampNoExclu <- read.csv("~/repos/BRCA1ness_by_SVM/annotations_and_backups/031318_Complete_Sample_Sheet_No_Exclusion.txt", header=T, sep="\t", stringsAsFactors=F);
mySampNoExclu$includedInBRCA1Analysis <- mySampNoExclu$patients %in% my_samples$patients;
mySampNoExclu$anyRecepPos <- mySampNoExclu$ER=="Positive" |  mySampNoExclu$PR=="Positive" | mySampNoExclu$HER2=="Positive";

#---------------------------------------Distribution of BRCA1-like probability---------------------------------------
mySampNoExclu$recepPos[mySampNoExclu$anyRecepPos] <- "Receptor Positive";
mySampNoExclu$recepPos[! mySampNoExclu$anyRecepPos] <- "TNBC";
mySampNoExclu$recepPos[is.na(mySampNoExclu$anyRecepPos)] <- "Unknown"; 

png("~/Downloads/BRCA1ness_figures/Figure2B.png", res=300, units="in", height=6, width=11.69);
ggplot(subset(mySampNoExclu, recepPos != "Unknown"), aes(x=reorder(patients, BRCA1_prob), y=BRCA1_prob, fill=recepPos)) +
  geom_bar(stat="identity") +
  geom_hline(yintercept=0.50, linetype="dashed") + 
  labs(x="TCGA Sample", y="SVM BRCA1-like probability", title="TCGA") +
  scale_x_discrete(expand=c(0, 0)) +
  scale_y_continuous(breaks=seq(0,1,0.2), expand=c(0,0), limits=c(0,1.05)) +
  scale_fill_manual(values=c("black","blue","lightgray")) +
  myWaterfallTheme  +
  facet_wrap(~ recepPos, scale="free")
dev.off();

#---------------------------------------------Table One: Overall---------------------------------------------
## Table One
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
  argsExact = list(workspace = 2*10^5), #default
  testNormal = oneway.test,
  argsNormal = list(var.equal=TRUE),
  includeNA = TRUE
);

## Visualize & export:
y <- print(myTableOne, showAllLevels=TRUE);
# write.csv(y, file="~/Downloads/051118_TCGA-BRCA_No-exclusion_TableOne.csv",quote=F);
# WriteXLS(mySampNoExclu, ExcelFileName="~/Downloads/031318_Sample_Sheet_TCGA_NoExclu.xls");

## Response to Reviewer #2 comment 6:
myTableOneERpos <- CreateTableOne(
  data = subset(mySampNoExclu, ER=="Positive"), 
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
yERpos <- print(myTableOneERpos, showAllLevels=TRUE);
# write.csv(yERpos, file="~/Downloads/101818_TCGA-BRCA_ER+BRCA.csv",quote=FALSE);

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
# write.csv(y, file="~/Downloads/BRCA1ness_figures/041118_TCGA_Subset_Table_one.csv",quote=F);

#---------------------------------------------Response to Reviewer #1, comment g---------------------------------------------
contTabRace <- matrix(c(60,61, (201+20+1), (488+36+0)), byrow=TRUE, ncol=2);
contTabRace
sum(contTabRace)
fisher.test(contTabRace)

#---------------------------------------------Response to Reviewer #2, comment 7---------------------------------------------
contTabStageByDataset <- matrix(c((63+165),(26+102),(238+473),(220+1106)), byrow=TRUE, ncol=2);
contTabStageByDataset
sum(contTabStageByDataset)
fisher.test(contTabStageByDataset)

contTabStageTCGA <- matrix(c(63,165,238,473), byrow=TRUE, ncol=2);
contTabStageTCGA
fisher.test(contTabStageTCGA)

contTabStageMETABRIC <- matrix(c(26,102,220,1106), byrow=TRUE, ncol=2);
contTabStageMETABRIC
fisher.test(contTabStageMETABRIC)

png("~/Downloads/BRCA1ness_figures/Age_distribution_TCGA.png", res=300, units="in", height=6, width=11.69);
ggplot(mySampNoExclu, aes(Age, color=SVM_BRCA1, fill=SVM_BRCA1)) +
  geom_density(alpha=0.25) +
  scale_fill_manual(values=c("mediumorchid3","darkolivegreen3")) +
  scale_color_manual(values=c("mediumorchid3","darkolivegreen3")) +
  scale_y_continuous(limits=c(0,0.04)) +
  ggtitle("TCGA") +
  facet_wrap(~ SVM_BRCA1) +
  myScatterTheme
dev.off();
